# Description:
# This script processes gene expression data from a WT (wild-type) uninjected sample group and performs several analyses to 
# visualize gene relationships. The data is filtered and selected for Gli target and stress related genes (see Stress Gene Definitions.R). 
# It categorizes the genes into subtypes and generates an annotation for each gene based on its category. The correlation matrix 
# between genes is calculated and visualized through a heatmap, showing gene correlations across samples. 
# Additionally, a significance matrix is created by testing correlations between genes, and a binary heatmap is plotted to highlight 
# significant relationships (p < 0.05) between genes. 
# These steps allow for exploration of gene expression patterns and relationships within the WT uninjected sample group.
#
# Author: Sam Mindham
# Date: 2025-12-03
#
# Inputs:
# - counts_noMito_noZeros_Mutant_vs_WT.txt: Raw gene count data for WT and mutant samples.
# - Gen_Wide_Gli.csv: Gene list for Gli-Target-Genes.
# - stress_genes: List of stress-related genes.
# - significant_stress: Dataframe with stress gene subtypes.
#
# Outputs:
# - Correlation heatmap (UWTucHM.png): A visualization of the correlation matrix for WT uninjected genes.
# - Significance matrix heatmap (UWTucHM_significance.png): A binary heatmap showing significant gene correlations (p < 0.05).
#
# Dependencies:
# - R packages: pheatmap, dplyr.
--------------------------------------------------------------------------------
  library(pheatmap)

# Group: WT Uninjected

# Step 1: Read raw count data
WT_vs_mutant_rawcounts_WT_Uninjected <- read.delim('counts_noMito_noZeros_Mutant_vs_WT.txt', sep = '')

# Step 2: Select relevant columns and prepare the data
WT_rawcounts_WT_Uninjected <- WT_vs_mutant_rawcounts_WT_Uninjected[, c(2, 3, 5)]

# Add row names to the data frame
WT_rawcounts_WT_Uninjected <- data.frame(WT_rawcounts_WT_Uninjected, rownames(WT_rawcounts_WT_Uninjected))

# Clean up row names
WT_rawcounts_WT_Uninjected$rownames.WT_rawcounts_WT_Uninjected <- sub('.*_', '', WT_rawcounts_WT_Uninjected$rownames.WT_rawcounts_WT_Uninjected)
WT_rawcounts_WT_Uninjected$rownames.WT_rawcounts_WT_Uninjected <- make.unique(WT_rawcounts_WT_Uninjected$rownames.WT_rawcounts_WT_Uninjected)

# Remove row names column and set new row names
WT_rawcounts_WT_Uninjected <- WT_rawcounts_WT_Uninjected %>% 
  remove_rownames() %>% 
  column_to_rownames(var = 'rownames.WT_rawcounts_WT_Uninjected')

# Step 3: Filter out genes with all zero counts
WT_rawcounts_WT_Uninjected <- WT_rawcounts_WT_Uninjected[rowSums(WT_rawcounts_WT_Uninjected != 0) > 0, ]

# Step 4: Load gene lists
hedgehog_genes_WT_Uninjected <- read.csv('Gen_Wide_Gli.csv')$Gene
stress_genes_WT_Uninjected <- read.csv('Significant_Stress.csv')
stress_genes_WT_Uninjected <- rownames(stress_genes_WT_Uninjected)
# Step 5: Select genes present in the dataset
selected_genes_WT_Uninjected <- unique(c(hedgehog_genes_WT_Uninjected, stress_genes_WT_Uninjected))
selected_genes_WT_Uninjected <- selected_genes_WT_Uninjected[selected_genes_WT_Uninjected %in% rownames(WT_rawcounts_WT_Uninjected)]

# Step 6: Assign categories (Hedgehog, Stress with subtype, or Ciliary)
gene_category_WT_Uninjected <- data.frame(
  Gene = selected_genes_WT_Uninjected,
  Category = ifelse(
    selected_genes_WT_Uninjected %in% hedgehog_genes_WT_Uninjected, "Gli Target", "Stress")
)

# Step 7: Map stress subtypes
stress_genes <- rownames(significant_stress)
stress_type_map_WT_Uninjected <- setNames(significant_stress$Stress_Type, stress_genes)
gene_category_WT_Uninjected$Stress_Type <- stress_type_map_WT_Uninjected[gene_category_WT_Uninjected$Gene]

# Step 8: Define combined annotation column
gene_category_WT_Uninjected$Combined_Annotation <- ifelse(
  gene_category_WT_Uninjected$Category == "Stress",
  paste("Stress -", gene_category_WT_Uninjected$Stress_Type),
  gene_category_WT_Uninjected$Category
)

# Step 9: Prepare annotation row
annotation_row_WT_Uninjected <- data.frame(Combined_Annotation = gene_category_WT_Uninjected$Combined_Annotation)
rownames(annotation_row_WT_Uninjected) <- gene_category_WT_Uninjected$Gene
colnames(annotation_row_WT_Uninjected) <- 'Gene Type'

# Step 10: Sort genes by category for better visualization
sorted_gene_order_WT_Uninjected <- gene_category_WT_Uninjected$Gene[order(gene_category_WT_Uninjected$Category)]
WT_rawcounts_sorted_WT_Uninjected <- WT_rawcounts_WT_Uninjected[sorted_gene_order_WT_Uninjected, ]

# Remove fourth column and filter out genes with zero counts
WT_rawcounts_sorted_WT_Uninjected <- WT_rawcounts_sorted_WT_Uninjected[, -4]
WT_rawcounts_sorted_WT_Uninjected <- WT_rawcounts_sorted_WT_Uninjected[rowSums(WT_rawcounts_sorted_WT_Uninjected != 0) > 0, ]

# Step 11: Compute correlation matrix
# Ensure the data is numeric
WT_rawcounts_sorted_WT_Uninjected[] <- lapply(WT_rawcounts_sorted_WT_Uninjected, as.numeric)

# Compute correlation matrix
cor_matrix_WT_Uninjected <- cor(t(WT_rawcounts_sorted_WT_Uninjected), use = "pairwise.complete.obs")

# Step 12: Plot correlation heatmap without clustering
UWTucHM <- pheatmap(
  cor_matrix_WT_Uninjected, 
  annotation_row = annotation_row_WT_Uninjected, 
  annotation_col = annotation_row_WT_Uninjected, 
  cluster_rows = FALSE, 
  cluster_cols = FALSE, 
  main = 'WT Uninjected Correlation Matrix', 
  fontsize_row = 6, 
  fontsize_col = 6
)

# Step 13: Compute p-value matrix
WT_rawcounts_sorted_WT_Uninjected_t <- t(WT_rawcounts_sorted_WT_Uninjected)
p_value_matrix_WT_Uninjected <- matrix(NA, nrow = ncol(WT_rawcounts_sorted_WT_Uninjected_t), ncol = ncol(WT_rawcounts_sorted_WT_Uninjected_t))

# Loop through columns and perform correlation tests
for (i in 1:ncol(WT_rawcounts_sorted_WT_Uninjected_t)) {
  for (j in 1:ncol(WT_rawcounts_sorted_WT_Uninjected_t)) {
    if (i != j) {
      cor_test_result <- cor.test(WT_rawcounts_sorted_WT_Uninjected_t[, i], WT_rawcounts_sorted_WT_Uninjected_t[, j])
      p_value_matrix_WT_Uninjected[i, j] <- cor_test_result$p.value
    }
  }
}

# Step 14: Generate binary significance matrix based on p-value threshold
significance_matrix_WT_Uninjected <- ifelse(p_value_matrix_WT_Uninjected < 0.05, 1, 0)
rownames(significance_matrix_WT_Uninjected) <- rownames(WT_rawcounts_sorted_WT_Uninjected)
colnames(significance_matrix_WT_Uninjected) <- rownames(WT_rawcounts_sorted_WT_Uninjected)

# Define color palette for the heatmap
color_palette <- colorRampPalette(c("grey", "green"))(2)

# Step 15: Create annotation for the significance heatmap
annotation_row_WT_Uninjected <- data.frame(Combined_Annotation = gene_category_WT_Uninjected$Combined_Annotation)
rownames(annotation_row_WT_Uninjected) <- gene_category_WT_Uninjected$Gene
colnames(annotation_row_WT_Uninjected) <- 'Gene Type'

# Step 16: Ensure row and column names of the annotation match the significance matrix
annotation_row_WT_Uninjected <- annotation_row_WT_Uninjected[rownames(significance_matrix_WT_Uninjected), , drop = FALSE]
annotation_col_WT_Uninjected <- annotation_row_WT_Uninjected[colnames(significance_matrix_WT_Uninjected), , drop = FALSE]

# Step 17: Plot the significance heatmap
UWTucHM_significance <- pheatmap(
  significance_matrix_WT_Uninjected, 
  annotation_row = annotation_row_WT_Uninjected, 
  annotation_col = annotation_col_WT_Uninjected, 
  cluster_rows = FALSE, 
  cluster_cols = FALSE, 
  main = 'Significance Matrix (p < 0.05) WT Uninjected', 
  fontsize_row = 6, 
  fontsize_col = 6, 
  color = color_palette,  # Red for non-significant (0) and Green for significant (1)
  legend = FALSE,         # No legend as the colors are binary
  display_numbers = FALSE,    # No numbers displayed
  number_format = "%.0f",     # Display numbers as 0 or 1
  breaks = c(-Inf, 0.5, Inf)  # Breaks to separate 0 and 1
)
