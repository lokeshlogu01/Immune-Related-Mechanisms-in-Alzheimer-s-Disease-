# Load required libraries
library(DESeq2)     # For differential gene expression analysis
library(ggplot2)    # For visualization
library(pheatmap)   # For heatmap visualization
library(limma)      # For additional analysis and normalization
library(edgeR)      # For alternative differential analysis
library(caret)      # For machine learning models
library(glmnet)     # For Lasso regression
library(randomForest) # For Random Forest models
library(pROC)       # For ROC analysis

# File paths (update these to match your dataset locations)
file_path <- "C:/Users/LONOVO/Downloads/GSE126003_raw_counts_GRCh38.p13_NCBI (1).tsv"
sample_description_path <- "C:/Users/LONOVO/Downloads/sample_descrption.csv"

# Read raw counts data (RNA-seq)
raw_counts <- as.matrix(read.table(file_path, header=TRUE, row.names=1, sep="\t"))

# Read sample description file
sample_description <- read.csv(sample_description_path, header=FALSE, stringsAsFactors=FALSE)
colnames(sample_description) <- c("SampleID", "Condition")

# Extract condition (assuming conditions like 'AD', 'NT', 'IL6' etc., adjust based on your dataset)
sample_description$Condition <- factor(sapply(strsplit(sample_description$Condition, "_"), "[", 2),
                                       levels = c("NT", "IL6", "veh", "PBS", "AD", "Control"))

# Check if sample names in raw counts match the sample description
if (!all(colnames(raw_counts) %in% sample_description$SampleID)) {
  stop("Mismatch between sample names in raw counts and sample description file")
}

# Create DESeq2 object for differential expression analysis
dds <- DESeqDataSetFromMatrix(countData = raw_counts,
                              colData = sample_description,
                              design = ~ Condition)

# Perform differential expression analysis using DESeq2
dds <- DESeq(dds)

# Contrast of interest (e.g., AD vs Control)
contrast <- c("Condition", "AD", "Control")

# Get differential expression results
results_deseq2 <- results(dds, contrast = contrast)

# Filter significant DE genes (adjusted p-value < 0.05)
sig_results_deseq2 <- results_deseq2[which(results_deseq2$padj < 0.05), ]

# Save significant DE genes to a CSV file for future analysis
write.csv(sig_results_deseq2, file = "DESeq2_Significant_DE_genes_AD_vs_Control.csv")

# Print summary and top DE genes
summary(results_deseq2)
head(sig_results_deseq2)

# Load data again for limma analysis (this is necessary if you've reloaded or altered data)
raw_counts <- read.table(file_path, header=TRUE, row.names=1, sep="\t")

# Create DGEList object for limma analysis
dge <- DGEList(counts = raw_counts)
dge <- calcNormFactors(dge)

# Create design matrix (customize this based on your conditions like AD, Control, etc.)
design <- model.matrix(~ 0 + sample_description$Condition)
colnames(design) <- levels(sample_description$Condition)

# Estimate dispersion for limma
dge <- estimateDisp(dge, design)

# Fit the model using limma
fit <- glmFit(dge, design)

# Perform likelihood ratio test for differential expression (e.g., AD vs Control)
lrt <- glmLRT(fit, contrast = c(1, 0, -1))

# Get differential expression results
lrt_results <- topTags(lrt, n = nrow(dge))

# Filter significant DE genes based on FDR < 0.05
sig_lrt_results <- lrt_results$table[lrt_results$table$FDR < 0.05, ]

# Save significant DE genes to a CSV file
write.csv(sig_lrt_results, file = "limma_Significant_DE_genes_AD_vs_Control.csv")

# Print summary and top DE genes from limma results
head(sig_lrt_results)

# **Optional: Machine learning analysis**
# Prepare the data for a machine learning model, e.g., random forest
rf_data <- sig_results_deseq2[, -c(1:2)] # Remove non-numeric columns like gene_id

# Random Forest model
set.seed(42)
rf_model <- randomForest(x = rf_data, y = factor(sample_description$Condition), ntree = 500)

# Importance of features (genes)
importance(rf_model)

# **Optional: ROC Curve for model evaluation**
# Predicted probabilities
pred_probs <- predict(rf_model, rf_data, type = "prob")

# Create ROC curve for model evaluation
roc_curve <- roc(sample_description$Condition, pred_probs[, 2])
plot(roc_curve)

# Save ROC plot
dev.copy(png, "ROC_curve.png")
dev.off()
