# breast-cancer-dta-normalization-in-GSE-from-data-normalization-to-volcano-ma-PLOT-
Data Normalization ,Preprocessing,Heatmap,Hierarichal Dendogram,PCA ,DEA-DEG,DEG significant ,volcano,MA plot
# breast-cancer-dta-normalization-in-GSE-from-data-normalization-to-volcano-ma-PLOT-
Data Normalization ,Preprocessing,Heatmap,Hierarichal Dendogram,PCA ,DEA-DEG,DEG significant ,volcano,MA plot
library(GEOquery)
library(limma)
library(dplyr)
library(ggplot2)

# Load dataset
gse <- getGEO("GSE42568", GSEMatrix = TRUE)
metadata <- pData(gse[[1]])
expression_data <- exprs(gse[[1]])

png(("______/GSe42568/boxplot_before_normalization.png"))
boxplot(expression_data, main = "Before Normalization", las = 2, outline = FALSE)
dev.off()

expression_data <- normalizeBetweenArrays(expression_data, method = "quantile")

png(("_______/GSE42568/boxplot_after_normalization.png"))
boxplot(expression_data, main = "After Normalization", las = 2, outline = FALSE)
dev.off()
library(GEOquery)
library(limma)
library(pheatmap)
🔹 2. Download and Process the Dataset
r
gse <- getGEO("GSE42568", GSEMatrix = TRUE)[[1]]
expression_data <- exprs(gse)
pheno_data <- pData(gse)
🔹 3. Filter Top Variable Genes
r
top_var_genes <- order(apply(expression_data, 1, var), decreasing = TRUE)[1:50]
heatmap_data <- expression_data[top_var_genes, ]
🔹 4. Annotate Samples (Optional)
r
annotation_col <- data.frame(Group = pheno_data$characteristics_ch1.1)
rownames(annotation_col) <- colnames(expression_data)
🔹 5. Generate the Heatmap
r
pheatmap(heatmap_data,
         annotation_col = annotation_col,
         show_rownames = FALSE,
         show_colnames = FALSE,
         main = "Top 50 Most Variable Genes - GSE42568")

# Capture the plot
heatmap_plot <- pheatmap(heatmap_data,
                         annotation_col = annotation_col,
                         show_rownames = FALSE,
                         show_colnames = FALSE,
                         main = "Top 50 Most Variable Genes - GSE42568")

ggsave(filename = "c:/Users/___________/heatmap_top50_genes.png",
       plot = heatmap_plot,
       width = 8,
       height = 10)
dev.off()
dist_matrix <- dist(t(expression_data), method = "euclidean")

# Perform hierarchical clustering
hclust_result <- hclust(dist_matrix, method = "complete") 

# Plot dendrogram
plot(hclust_result,
     labels = pheno_data$characteristics_ch1.1,
     main = "Hierarchical Clustering Dendrogram",
     xlab = "",
     sub = "",
     cex = 0.8)
pdf("c:/Users/_____________/hierarchical_dendrogram.pdf", width = 10, height = 6)
plot(hclust_result,
     labels = pheno_data$characteristics_ch1.1,
     main = "Hierarchical Clustering Dendrogram")
dev.off()
png("c:/Users/______________/hierarchical_dendrogram.png", width = 800, height = 600)
plot(hclust_result,
     labels = pheno_data$characteristics_ch1.1,  # or any sample labels you prefer
     main = "Hierarchical Clustering Dendrogram")
dev.off()

png("c:/Users/___________/hierarchical_dendrogram.png", width = 1600, height = 1200, res = 200)
plot(hclust_result,
     labels = pheno_data$characteristics_ch1.1,  # or any sample labels you prefer
     main = "Hierarchical Clustering Dendrogram")
dev.off()



# Install if needed
if (!requireNamespace("GEOquery", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("GEOquery")
}

library(GEOquery)

# Download the series matrix
gse <- getGEO("GSE42568", GSEMatrix = TRUE)
expr_set <- gse[[1]]  

# Extract expression data
expr_matrix <- exprs(expr_set)

# Check the top of the matrix
head(expr_matrix)
write.csv(expr_matrix, "GSE42568_expr_matrix.csv")
# Create the directory if it doesn't exist
dir.create("C:/Users/_________", showWarnings = FALSE)

# Save the expression matrix as a CSV
write.csv(expr_matrix, "C:/Users/____________/GSE42568_expr_matrix.csv")
# 1. Load the expression matrix (skip this if already loaded)
expr_matrix <- read.csv("C:/Users/__________/GSE42568_expr_matrix.csv", row.names = 1)

# 2. Log-transform if not already normalized
expr_matrix <- log2(expr_matrix + 1)

# 3. Transpose: samples as rows
expr_matrix_t <- t(expr_matrix)

# 4. Run PCA
pca_result <- prcomp(expr_matrix_t, scale. = TRUE)

# 5. Create a mock grouping variable (replace with actual phenodata if available)
Group <- rep("Unknown", nrow(expr_matrix_t))  # Placeholder

pca_df <- data.frame(PC1 = pca_result$x[, 1],
                     PC2 = pca_result$x[, 2],
                     Group = Group)

# 6. Plot using ggplot2
library(ggplot2)

ggplot(pca_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3) +
  labs(title = "PCA Plot - GSE42568",
       x = paste0("PC1 (", round(summary(pca_result)$importance[2, 1] * 100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca_result)$importance[2, 2] * 100, 1), "%)")) +
  theme_minimal()
# Remove genes (columns) with zero variance across samples
expr_matrix_t_clean <- expr_matrix_t[, apply(expr_matrix_t, 2, var) != 0]
pca_result <- prcomp(expr_matrix_t_clean, scale. = TRUE)
Group <- rep("Unknown", nrow(expr_matrix_t_clean))  # Still placeholder
pca_df <- data.frame(PC1 = pca_result$x[, 1],
                     PC2 = pca_result$x[, 2],
                     Group = Group)
library(ggplot2)

ggplot(pca_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3) +
  labs(title = "PCA Plot - GSE42568",
       x = paste0("PC1 (", round(summary(pca_result)$importance[2, 1] * 100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca_result)$importance[2, 2] * 100, 1), "%)")) +
  theme_minimal()
# Open PNG device
png(file=file.path("____/PCA_plot.png", width = 1000, height = 800)
# Open PNG device
png(file = file.path("_____/PCA_plot.png"), width = 1000, height = 800)

# Generate the plot
plot(1:10, 1:10, main = "Sample Plot")  # Replace with your actual plotting code

# Close the device
dev.off()

# Plot the PCA
ggplot(pca_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3) +
  labs(title = "PCA Plot - GSE42568",
       x = paste0("PC1 (", round(summary(pca_result)$importance[2, 1] * 100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca_result)$importance[2, 2] * 100, 1), "%)")) +
  theme_minimal()
png("C:/Users/______________/GSE42568_PCA_plot.png", width = 1000, height = 800)

ggplot(pca_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3) +
  labs(title = "PCA Plot - GSE42568",
       x = paste0("PC1 (", round(summary(pca_result)$importance[2, 1] * 100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca_result)$importance[2, 2] * 100, 1), "%)")) +
  theme_minimal()

dev.off()
pdf("C:/Users/____________/GSE42568_PCA_plot.pdf", width = 8, height = 6)  # dimensions in inches

ggplot(pca_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3) +
  labs(title = "PCA Plot - GSE42568",
       x = paste0("PC1 (", round(summary(pca_result)$importance[2, 1] * 100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca_result)$importance[2, 2] * 100, 1), "%)")) +
  theme_minimal(base_size = 14) +  # increase base font size
  theme(legend.position = "right")

dev.off()

tiff("C:/_______________/GSE42568_PCA_plot.tiff", units = "in", width = 8, height = 6, res = 300)

ggplot(pca_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3) +
  labs(title = "PCA Plot - GSE42568",
       x = paste0("PC1 (", round(summary(pca_result)$importance[2, 1] * 100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca_result)$importance[2, 2] * 100, 1), "%)")) +
  theme_minimal(base_size = 14)
dev.off()

save(expression_data, metadata, subtype,
     file = file.path("_____________/GSE42568/processed_data_step1.RData"))
# Load GEOquery
library(GEOquery)
gse <- getGEO("GSE42568", GSEMatrix = TRUE)
expr_set <- gse[[1]]

# Extract expression matrix and phenotype data
expr_matrix <- exprs(expr_set)
pheno_data <- pData(expr_set)
# Use tissue type as grouping variable
group <- as.factor(make.names(pheno_data$`tissue:ch1`))  # e.g., "Breast.tissue.cancer", "Breast.tissue.normal"

# Build design matrix
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

library(limma)

# Fit linear model
fit <- lmFit(expr_matrix, design)

# Define contrast: cancer vs normal
contrast_matrix <- makeContrasts(Diff = Breast.tissue.cancer - Breast.tissue.normal, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)



levels(group)
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
colnames(design)
contrast_matrix <- makeContrasts(Diff = Breast.tissue.cancer - Breast.tissue.normal,
                                 levels = design)
contrast_matrix <- makeContrasts(Diff = `Breast.tissue.cancer` - `Breast.tissue.normal`,
                                 levels = design)
contrast_matrix <- makeContrasts(Diff = `Breast.tissue.cancer` - `Breast.tissue.normal`,
                                 levels = design)
contrast_matrix <- makeContrasts(Diff = Cancer - Normal, levels = design)
colnames(design) <- c("Cancer", "Normal")
contrast_matrix <- makeContrasts(Diff = Cancer - Normal, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)
top_genes <- topTable(fit2, adjust.method = "fdr", number = Inf)
contrast_matrix <- makeContrasts(Diff = Cancer - Normal, levels = design)
# Apply the contrast to the fitted model
fit2 <- contrasts.fit(fit, contrast_matrix)

# Compute moderated statistics
fit2 <- eBayes(fit2)
top_genes <- topTable(fit2, adjust.method = "fdr", number = Inf)
write.csv(top_genes, "C:/Users/________/GSE42568_DEGs_full.csv", row.names = TRUE)
degs_sig <- top_genes[top_genes$adj.P.Val < 0.05 & abs(top_genes$logFC) > 1, ]
write.csv(degs_sig, "____________/GSE42568/GSE42568_DEGs_significant.csv")
degs_sig <- top_genes[top_genes$adj.P.Val < 0.05 & abs(top_genes$logFC) > 1, ]
write.csv(degs_sig, "____________/GSE42568/GSE42568_DEGs_significant.csv")
# Prepare data
top_genes$logP <- -log10(top_genes$P.Value)
top_genes$Significance <- ifelse(top_genes$adj.P.Val < 0.05 & abs(top_genes$logFC) > 1,
                                 "Significant", "Not Significant")

# Plot
library(ggplot2)

ggplot(top_genes, aes(x = logFC, y = logP, color = Significance)) +
  geom_point(alpha = 0.8, size = 1.8) +
  scale_color_manual(values = c("gray", "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(title = "Volcano Plot - GSE42568",
       x = "log2 Fold Change",
       y = "-log10(P-value)") +
  theme_minimal(base_size = 14)
png("__________________/GSE42568_volcano_plot.png", width = 1200, height = 1000, res = 150)

ggplot(top_genes, aes(x = logFC, y = -log10(P.Value), color = Significance)) +
  geom_point(alpha = 0.8, size = 1.8) +
  scale_color_manual(values = c("gray", "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(title = "Volcano Plot - GSE42568", x = "log2 Fold Change", y = "-log10(P-value)") +
  theme_minimal(base_size = 14)

dev.off()
# Calculate average expression
top_genes$AveExpr <- rowMeans(expr_matrix[rownames(top_genes), ])

ggplot(top_genes, aes(x = AveExpr, y = logFC, color = Significance)) +
  geom_point(alpha = 0.8, size = 1.8) +
  scale_color_manual(values = c("gray", "blue")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "MA Plot - GSE42568",
       x = "Average Expression",
       y = "log2 Fold Change") +
  theme_minimal(base_size = 14)

# Calculate average expression if not already done
top_genes$AveExpr <- rowMeans(expr_matrix[rownames(top_genes), ])

# Save MA Plot
tiff("C:_________________/GSE42568/GSE42568_MA_plot.tiff", units = "in", width = 8, height = 6, res = 300)
ggplot(top_genes, aes(x = AveExpr, y = logFC, color = Significance)) +
  geom_point(alpha = 0.8, size = 1.8) +
  scale_color_manual(values = c("gray", "blue")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "MA Plot - GSE42568", x = "Average Expression", y = "log2 Fold Change") +
  theme_minimal(base_size = 14)

dev.off()
