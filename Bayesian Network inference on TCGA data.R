# Required packages
if (!requireNamespace("bnlearn", quietly=TRUE)) install.packages("bnlearn")
if (!requireNamespace("matrixStats", quietly=TRUE)) install.packages("matrixStats")
if (!requireNamespace("Rgraphviz", quietly=TRUE)) BiocManager::install("Rgraphviz")
if (!requireNamespace("arules", quietly=TRUE)) install.packages("arules")  # used as fallback discretizer

library(bnlearn)
library(matrixStats)
library(Rgraphviz)
library(arules)

# ---------------------------
# Helper functions
# ---------------------------

# Safe discretization function (handles ties / constant columns)
discretize_vector_safe <- function(x, bins = 3, method = c("quantile", "frequency")) {
  method <- match.arg(method)
  # If all NA or length 0
  if (all(is.na(x)) || length(x) == 0) return(rep(NA_integer_, length(x)))
  # If nearly constant => single level
  if (length(unique(na.omit(x))) <= 1) return(rep(1L, length(x)))
  
  # Try unique quantile breaks
  if (method == "quantile") {
    br <- unique(quantile(x, probs = seq(0, 1, length.out = bins + 1), na.rm = TRUE))
    if (length(br) <= 2) {
      # fallback to frequency-based discretize from arules
      return(as.integer(arules::discretize(x, method = "frequency", categories = bins)))
    } else {
      # ensure strictly monotonic breaks for cut
      # expand identical end-points slightly
      if (any(duplicated(br))) {
        br[duplicated(br)] <- br[duplicated(br)] + .Machine$double.eps
      }
      return(as.integer(cut(x, breaks = br, include.lowest = TRUE, labels = FALSE)))
    }
  } else {
    # frequency method (arules)
    return(as.integer(arules::discretize(x, method = "frequency", categories = bins)))
  }
}

# Convert discretized matrix -> factor dataframe required by bnlearn
to_factor_df <- function(df_discrete) {
  df_f <- as.data.frame(df_discrete, stringsAsFactors = FALSE)
  df_f[] <- lapply(df_f, function(col) {
    if (is.numeric(col)) {
      # replace NA with a new level if present
      if (any(is.na(col))) {
        col[is.na(col)] <- max(col, na.rm = TRUE) + 1
      }
      factor(col)
    } else {
      factor(col)
    }
  })
  return(df_f)
}

# Remove columns with single level (no variability)
remove_single_level <- function(df) {
  keep <- sapply(df, function(col) length(levels(col)) > 1)
  removed <- names(df)[!keep]
  if (length(removed) > 0) {
    message("Removing single-level variables: ", paste(removed, collapse = ", "))
  }
  df[, keep, drop = FALSE]
}

# Filter top variable genes (samples x genes numeric matrix)
select_top_variable_genes <- function(expr_numeric, top_n = 200) {
  if (!is.matrix(expr_numeric)) expr_numeric <- as.matrix(expr_numeric)
  if (!is.numeric(expr_numeric)) stop("expr_numeric must be numeric matrix (samples x genes).")
  vars <- colVars(expr_numeric, na.rm = TRUE)
  names(vars) <- colnames(expr_numeric)
  if (is.null(names(vars))) stop("expr_numeric must have column names (gene IDs).")
  top_n <- min(top_n, length(vars))
  top_genes <- names(sort(vars, decreasing = TRUE))[1:top_n]
  expr_numeric[, top_genes, drop = FALSE]
}

# Plot only strong arcs and save to file
plot_strong_arcs <- function(boot_res, threshold = 0.7, out_pdf = "bn_strong_arcs.pdf", width = 10, height = 10) {
  avg_net <- averaged.network(boot_res, threshold = threshold)
  if (length(arcs(avg_net)) == 0) {
    message("No arcs exceed threshold = ", threshold)
    return(NULL)
  }
  pdf(out_pdf, width = width, height = height)
  graphviz.plot(avg_net, main = paste0("Averaged BN (strength >= ", threshold, ")"))
  dev.off()
  message("Saved averaged network (strength >= ", threshold, ") to: ", out_pdf)
  return(avg_net)
}

# ---------------------------
# Main pipeline function
# ---------------------------

bn_pipeline_tcga <- function(expr_matrix,
                             sample_is_rows = TRUE,
                             top_n_genes = 200,
                             discretize_bins = 3,
                             discretize_method = c("quantile", "frequency"),
                             seed = 42,
                             hc_score = "bde",    # use BDe as requested
                             bootstrap_R = 200,
                             bootstrap_threshold = 0.7,
                             out_pdf = "bn_strong_arcs.pdf") {
  set.seed(seed)
  discretize_method <- match.arg(discretize_method)
  
  # 0. Ensure matrix orientation: rows = samples, cols = genes
  if (!sample_is_rows) {
    expr_matrix <- t(expr_matrix)
  }
  if (is.null(colnames(expr_matrix))) {
    stop("expr_matrix must have column names (gene IDs).")
  }
  
  message("1) Selecting top variable genes (n = ", top_n_genes, ") ...")
  expr_top <- select_top_variable_genes(expr_matrix, top_n = top_n_genes)
  
  message("2) Discretizing each gene into ", discretize_bins, " bins (method = ", discretize_method, ") ...")
  disc_list <- apply(expr_top, 2, discretize_vector_safe, bins = discretize_bins, method = discretize_method)
  disc_df <- as.data.frame(disc_list, stringsAsFactors = FALSE)
  colnames(disc_df) <- colnames(expr_top)
  rownames(disc_df) <- rownames(expr_top)
  
  message("3) Converting to factors for bnlearn ...")
  disc_f <- to_factor_df(disc_df)
  
  message("4) Removing single-level factors (non-informative variables) ...")
  disc_f_filtered <- remove_single_level(disc_f)
  if (ncol(disc_f_filtered) < 2) stop("Too few variables after filtering single-level factors.")
  message("Variables kept for BN learning: ", ncol(disc_f_filtered))
  
  message("5) Running hc() with score = '", hc_score, "' ...")
  # For discrete BDe we can do: hc(data, score = "bde")
  bn_model <- hc(disc_f_filtered, score = hc_score)
  message("hc() finished. Number of arcs learned: ", nrow(arcs(bn_model)))
  
  message("6) Bootstrapping to estimate arc strength (R = ", bootstrap_R, ") ...")
  # We pass algorithm.args to use the same score inside bootstrapping if supported
  boot_res <- boot.strength(disc_f_filtered, R = bootstrap_R, algorithm = "hc", algorithm.args = list(score = hc_score))
  
  message("7) Summarizing top bootstrap strengths (head):")
  print(head(boot_res[order(-boot_res$strength), ], 10))
  
  message("8) Building averaged network and plotting edges with strength >= ", bootstrap_threshold, " ...")
  avg_net <- plot_strong_arcs(boot_res, threshold = bootstrap_threshold, out_pdf = out_pdf)
  
  # Return a list of outputs
  return(list(
    expr_top = expr_top,
    disc_df = disc_df,
    disc_f_filtered = disc_f_filtered,
    bn_model = bn_model,
    boot_res = boot_res,
    averaged_network = avg_net
  ))
}

# ---------------------------
# Example usage (uncomment and replace `my_expr` with your numeric matrix)
# ---------------------------
# my_expr: numeric matrix with rows = samples, cols = genes (e.g., TPM/log2TPM normalized)

my_expr <- t(expr_scaled)
#or
my_expr <- expr_mat
result <- bn_pipeline_tcga(my_expr,
                            sample_is_rows = TRUE,
                             top_n_genes = 50,
                            discretize_bins = 3,
                            discretize_method = "quantile",
                            seed = 123,
                            hc_score = "bde",
                            bootstrap_R = 30,
                            bootstrap_threshold = 0.75,
                            out_pdf = "tcga_bn_top150_bde_075.pdf")

# # inspect strong arcs table
 head(result$boot_res[order(-result$boot_res$strength), ], 20)
 
 # Simple bnlearn plotting
 graphviz.plot(avg_net, main = "Model-Averaged Network")
 
 

 library(igraph)
 
 # Example: take the top 20 strongest edges
 edges <- head(result$boot_res[order(-result$boot_res$strength), ], 20)
 
 # Keep only one direction for each edge based on higher direction score
 edges_unique <- edges[edges$direction >= 0.5, ]  # Keep arcs with >50% directional confidence
 
 # Build igraph object (directed = TRUE)
 g <- graph_from_data_frame(edges_unique, directed = TRUE)
 
 # Edge aesthetics
 E(g)$width <- edges_unique$strength * 5              # Line width ~ strength
 E(g)$color <- rgb(1 - edges_unique$direction, 0, edges_unique$direction)  # Red→Blue by direction
 
 # Node aesthetics
 V(g)$size <- 20
 V(g)$color <- "lightblue"
 V(g)$label.color <- "black"
 V(g)$label.cex <- 0.7
 
 # Plot DAG
 plot(
   g,
   layout = layout_with_sugiyama(g)$layout, # DAG-friendly layered layout
   edge.arrow.size = 0.4,
   main = "Bootstrap-Inferred DAG"
 )
 
