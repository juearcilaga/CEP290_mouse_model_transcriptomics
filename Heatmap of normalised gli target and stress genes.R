# Description:
# This script processes RNA-seq count data to generate a heatmap of normalized gene expression.
# It reads raw counts, normalizes them using Variance Stabilizing Transformation (VST),
# filters for common genes across datasets, selects genes of interest, and visualizes the results
# as a heatmap with annotations.
#
# Author: Sam Mindham
# Date: 2025-03-12
#
# Inputs:
# - counts_Mutant_Purmorphamine_WT_Purmorphamine_noMito_noZeros.txt: Raw count data for mutant vs. WT under Purmorphamine.
# - counts_noMito_noZeros_Mutant_vs_WT.txt: Raw count data for mutant vs. WT.
# - METADATA_Mutant_Purmorphamine_WT_Purmorphamine.txt: Metadata file for Purmorphamine-treated samples.
# - METADATA_Mutant_vs_WT.txt: Metadata file for general mutant vs. WT comparison.
# - Significant_Stress.csv: DeSeq data of significantly different stress related genes (Mutant vs WT)
# Outputs:
# - Heatmap of Normalized Gene Expression: A visualization of Z-score normalized expression values for selected genes.
#
# Dependencies:
# - R packages: DESeq2, pheatmap, RColorBrewer.
--------------------------------------------------------------------------------
  
# ------------------------------------------------------------------------------
# Load required libraries
# ------------------------------------------------------------------------------
library(DESeq2)
library(pheatmap)
library(RColorBrewer)

# ------------------------------------------------------------------------------
# Set working directory
# ------------------------------------------------------------------------------
setwd('C:/Users/sammi/OneDrive - Newcastle University/Shared_Folder_UG_proj_2025/DEG_Analysis/only_males/')

# ------------------------------------------------------------------------------
# Read in raw count data
# ------------------------------------------------------------------------------
PM_WT_vs_PM_Mutant_rawcounts <- read.delim(
  'outputs_DEG/counts_Mutant_Purmorphamine_WT_Purmorphamine_noMito_noZeros.txt', sep = ""
)
WT_vs_mutant_rawcounts <- read.delim(
  'outputs_DEG/counts_noMito_noZeros_Mutant_vs_WT.txt', sep = "\t"
)

# Read in metadata
metadata <- read.delim('Metadata_Files/METADATA_Mutant_Purmorphamine_WT_Purmorphamine.txt')
metadata2 <- read.delim('Metadata_Files/METADATA_Mutant_vs_WT.txt')

# ------------------------------------------------------------------------------
# Preprocess raw count data
# ------------------------------------------------------------------------------
# Select dataset to process
raw_counts <- PM_WT_vs_PM_Mutant_rawcounts
raw_counts2 <- WT_vs_mutant_rawcounts   # Change if needed

# Format gene names and make them unique
clean_gene_names <- function(df) {
  df <- data.frame(df, rownames(df))
  df$rownames.df. <- sub('.*_', '', df$rownames.df.)
  df$rownames.df. <- make.unique(df$rownames.df.)
  rownames(df) <- df$rownames.df.
  df <- df[, -ncol(df)]  # Remove the last column (gene names)
  return(df)
}

raw_counts <- clean_gene_names(raw_counts)
raw_counts2 <- clean_gene_names(raw_counts2)

# ------------------------------------------------------------------------------
# Normalize counts using Variance Stabilizing Transformation (VST)
# ------------------------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(countData = round(raw_counts),
                              colData = metadata,
                              design = ~ Genotype)
dds <- estimateSizeFactors(dds)
vst_counts <- assay(dds)

dds2 <- DESeqDataSetFromMatrix(countData = round(raw_counts2),
                               colData = metadata2,
                               design = ~ Genotype)
dds2 <- estimateSizeFactors(dds2)
vst_counts2 <- assay(dds2)

# ------------------------------------------------------------------------------
# Identify common genes between datasets
# ------------------------------------------------------------------------------
common_genes <- intersect(rownames(vst_counts), rownames(vst_counts2))

# Subset normalized counts to common genes
vst_counts_filtered <- vst_counts[common_genes, ]
vst_counts2_filtered <- vst_counts2[common_genes, ]

# Combine datasets
combined_vst_counts <- cbind(vst_counts_filtered, vst_counts2_filtered)
combined_vst_counts <- vst(round(combined_vst_counts))

# ------------------------------------------------------------------------------
# Select genes of interest
# ------------------------------------------------------------------------------
gwg <- read.csv('Gen_Wide_Gli.csv')
significant_stress <- read.csv('Significant_Stress.csv')
stress_genes <- rownames(significant_stress)

selected_genes <- unique(c(gwg, stress_genes))
selected_genes <- selected_genes[selected_genes %in% rownames(combined_vst_counts)]

# Subset normalized counts to selected genes
heatmap_data <- combined_vst_counts[selected_genes, ]

# ------------------------------------------------------------------------------
# Scale data for better visualization (Z-score normalization per gene)
# ------------------------------------------------------------------------------
heatmap_data_scaled <- t(scale(t(heatmap_data)))

# ------------------------------------------------------------------------------
# Define annotation for genes
# ------------------------------------------------------------------------------
annotation_row <- data.frame(Combined_Annotation = gene_category_WT_Uninjected$Combined_Annotation)
rownames(annotation_row) <- gene_category_WT_Uninjected$Gene
colnames(annotation_row) <- 'Gene Type'

# ------------------------------------------------------------------------------
# Define heatmap colors
# ------------------------------------------------------------------------------
heatmap_colors <- colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100)

# ------------------------------------------------------------------------------
# Plot the heatmap
# ------------------------------------------------------------------------------
pheatmap(
  heatmap_data_scaled,
  annotation_row = annotation_row,
  cluster_rows = FALSE,   # Do not cluster genes
  cluster_cols = FALSE,   # Do not cluster samples
  color = heatmap_colors,
  show_rownames = TRUE,   # Display gene names
  show_colnames = TRUE,   # Display sample names
  fontsize_row = 8,
  fontsize_col = 10,
  main = "Heatmap of Normalized Gene Expression"
)


