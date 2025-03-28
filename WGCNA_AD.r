# Install required libraries
install.packages(c("WGCNA", "clusterProfiler", "DESeq2", "limma", "randomForest"))

# Load necessary libraries
library(WGCNA)
library(clusterProfiler)
library(DESeq2)
library(limma)
library(randomForest)

# Load your scRNA-seq or bulk RNA-seq data
data <- read.csv("your_data.csv", header=TRUE, row.names=1)

# Pre-process data (remove lowly expressed genes)
data_clean <- data[rowMeans(data) > 1, ]  # Filter out lowly expressed genes

# Perform WGCNA
# Step 1: Choose soft threshold power (beta)
power = pickSoftThreshold(data_clean, powerVector = c(1:20), networkType = "unsigned")$powerEstimate

# Step 2: Construct the gene co-expression network
net = blockwiseModules(data_clean, power = power, TOMType = "unsigned", minModuleSize = 30, 
                        reassignThreshold = 0, mergeCutHeight = 0.25, numericLabels = TRUE, 
                        pamRespectsDendro = FALSE, verbose = 3)

# Step 3: Module-trait relationships (clinical traits)
module_trait_cor = cor(MEs, traits, use = "p")  # traits = Alzheimer’s Disease status or other clinical trait
module_trait_pvalue = corPvalueStudent(module_trait_cor, nSamples)

# Step 4: Identify hub genes in significant modules
module = net$colors
hubGenes <- names(data_clean)[module == "your_module_of_interest"] # Replace with module name of interest
