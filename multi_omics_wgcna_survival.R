# ----------------------------
# Multi-Omics WGCNA + Survival Analysis Pipeline
# TCGA-BRCA Example
# ----------------------------

# Required packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GO.db")

BiocManager::install("impute")

BiocManager::install("clusterProfiler")

BiocManager::install("org.Hs.eg.db")

library(org.Hs.eg.db)


library(WGCNA)

packages <- c("TCGAbiolinks", "SummarizedExperiment", "WGCNA", "survival", "survminer", 
              "clusterProfiler", "org.Hs.eg.db", "ggplot2")

lapply(packages, function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
})

options(stringsAsFactors = FALSE)
enableWGCNAThreads()

# ----------------------------
# Step 1: Download RNA-seq and Clinical Data
# ----------------------------
# Step 1: Define a small sample size
n_samples <- 300

# Select first 300 barcodes
GDCquery_project_samples <- getResults(query_exp, cols = "cases.submitter_id")
selected_barcodes <- unique(GDCquery_project_samples)[1:n_samples]

# Re-run query using selected barcodes
query_exp <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = "Primary Tumor",
  barcode = selected_barcodes
)


GDCdownload(query_exp)
data.exp <- GDCprepare(query = query_exp)


# Extract expression matrix
expr <- assay(data.exp)
expr <- log2(expr + 1)

# ----------------------------
# Step 2: Clinical Data
# ----------------------------
clinical <- colData(data.exp)
clinical_df <- as.data.frame(clinical)

# Extract relevant clinical info
clinical_traits <- data.frame(
  SampleID = rownames(clinical_df),
  Age = as.numeric(clinical_df$age_at_index),
  Stage = as.character(clinical_df$ajcc_pathologic_stage),
  OS.time = as.numeric(clinical_df$days_to_death),
  OS.status = ifelse(is.na(clinical_df$days_to_death), 0, 1)
)

# With this block (handles censored samples properly):
clinical_traits <- data.frame(
  SampleID = rownames(clinical_df),
  Age = as.numeric(clinical_df$age_at_index),
  Stage = as.character(clinical_df$ajcc_pathologic_stage),
  OS.time = ifelse(
    is.na(clinical_df$days_to_death),
    as.numeric(clinical_df$days_to_last_follow_up),
    as.numeric(clinical_df$days_to_death)
  ),
  OS.status = ifelse(is.na(clinical_df$days_to_death), 0, 1)
)

# Now remove only samples with completely missing survival info
clinical_traits <- clinical_traits[!is.na(clinical_traits$OS.time), ]

# Clean missing OS data
clinical_traits <- clinical_traits[!is.na(clinical_traits$OS.time), ]
samples <- intersect(colnames(expr), clinical_traits$SampleID)
expr <- expr[, samples]
clinical_traits <- clinical_traits[match(samples, clinical_traits$SampleID), ]

# Transpose for WGCNA
# Filter out lowly expressed genes
expr_filtered <- expr[rowMeans(expr > 1) > 0.2, ]
datExpr <- t(expr_filtered)

# ----------------------------
# Step 3: WGCNA Preprocessing
# ----------------------------
dim(datExpr)
# Should return: [number of samples, number of genes]
gsg <- goodSamplesGenes(datExpr, verbose = 3)
datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]

# Sample clustering to detect outliers
sampleTree <- hclust(dist(datExpr), method = "average")
plot(sampleTree, main = "Sample Clustering", sub="", xlab="")


# Choose soft thresholding power
powers <- c(1:20)
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

par(mfrow = c(1,2))
cex1 <- 0.9

# Scale-free topology fit
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R^2",
     type = "n", main = "Scale Independence")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
     labels = powers, cex = cex1, col = "red")
abline(h = 0.9, col = "blue")

# Mean connectivity
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity",
     type = "n", main = "Mean Connectivity")
text(sft$fitIndices[,1], sft$fitIndices[,5], labels = powers, cex = cex1, col = "red")

softPower <- sft$powerEstimate

# ----------------------------
# Step 4: Network Construction
# ----------------------------
net <- blockwiseModules(
  datExpr, power = softPower,
  TOMType = "unsigned", minModuleSize = 30,
  reassignThreshold = 0, mergeCutHeight = 0.25,
  numericLabels = TRUE, pamRespectsDendro = FALSE,
  saveTOMs = FALSE, verbose = 3
)

moduleColors <- labels2colors(net$colors)
MEs <- net$MEs

plotDendroAndColors(
  net$dendrograms[[1]],       # Dendrogram of the first block
  moduleColors[net$blockGenes[[1]]],  # Module colors for that block
  "Module colors",            # Title
  dendroLabels = FALSE,       # Hide gene labels for clarity
  hang = 0.03,                 # Branch height adjustment
  addGuide = TRUE,             # Horizontal guide lines
  guideHang = 0.05
)

dim(MEs)
pairs(MEs[, 1:5],  # first 5 modules for example
      main = "Pairs plot of module eigengenes",
      pch = 19, col = "steelblue")





# ----------------------------
# Step 5: Module-Trait Correlation
# ----------------------------
moduleTraitCor <- cor(MEs, clinical_traits[, c("Age", "OS.time")], use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples = nrow(datExpr))

# Plot heatmap
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(clinical_traits[, c("Age", "OS.time")]),
               yLabels = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               main = "Module-trait relationships")

# ----------------------------
# Step 6: Survival Analysis on Module Eigengene
# ----------------------------
# Pick top correlated module
selectedModule <- names(which.max(abs(moduleTraitCor[, "OS.time"])))
ME_selected <- MEs[[selectedModule]]
clinical_traits$ME <- ME_selected

# Kaplan-Meier Survival
fit <- survfit(Surv(OS.time, OS.status) ~ ME > median(ME), data = clinical_traits)
ggsurvplot(fit, data = clinical_traits, risk.table = TRUE,
           pval = TRUE, title = paste("Survival Curve - Module:", selectedModule))

# ----------------------------
# Step 7: Identify Hub Genes
# ----------------------------
# Recalculate module eigengenes to get proper names (if needed)
MEs_colored <- moduleEigengenes(datExpr, colors = moduleColors)$eigengenes

# Get the column names (should now be something like "MEblue", "MEbrown", ...)
colnames(MEs_colored)

# Suppose "ME22" is column index 3 in your original MEs matrix
column_index <- which(colnames(MEs) == "ME22")

# Get the corresponding color from the corrected eigengenes
correct_module_name <- colnames(MEs_colored)[column_index]  # e.g. "MEgreenyellow"

# Strip off the "ME" prefix to get the module color
selectedModuleColor <- substring(correct_module_name, 3)

moduleGenes <- moduleColors == selectedModuleColor
hubGenes <- rownames(t(datExpr))[moduleGenes]

# ----------------------------
# Step 8: Enrichment of Hub Genes
# ----------------------------
# Remove the dot and everything after it
hubGenes_clean <- sub("\\.\\d+$", "", hubGenes)

entrez_ids <- bitr(hubGenes_clean, fromType = "ENSEMBL",
                   toType = "ENTREZID", OrgDb = org.Hs.eg.db)


enrich <- enrichGO(gene = entrez_ids$ENTREZID,
                   OrgDb = org.Hs.eg.db,
                   keyType = "ENTREZID",
                   ont = "BP", pAdjustMethod = "BH",
                   qvalueCutoff = 0.05, readable = TRUE)

dotplot(enrich, showCategory = 10, title = "GO Enrichment of Hub Genes")

# ----------------------------
# Save Results
# ----------------------------
write.csv(hubGenes, "hub_genes.csv", row.names = FALSE)
saveRDS(net, "wgcna_network.rds")
saveRDS(clinical_traits, "clinical_traits.rds")

cat("Pipeline complete. Top module:", selectedModule, "\n")
