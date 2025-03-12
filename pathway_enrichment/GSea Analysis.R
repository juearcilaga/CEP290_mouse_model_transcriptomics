#-------------------------------------------------------------------------------
# Description:
# This script performs Gene Set Enrichment Analysis (GSEA) on differential 
# expression results from a DESeq2 analysis to identify biological processes 
# associated with the differentially expressed genes between mutant and wild-type 
# samples. It ranks genes by log2 fold change and uses the clusterProfiler package 
# to conduct the analysis. The results are visualized as a dot plot to highlight 
# the top enriched biological processes.
#
# Author: Sam Mindham
# Date: 2025/19/02
#
# Inputs:
# - Deseq_results_significant_all_mutant_vs_wildtype_padj0.5_l2fc0.3219 - Filtered.csv: 
#   CSV file containing DESeq2 differential expression results with significant genes 
#   based on a padj cutoff of 0.5 and log2 fold change cutoff of 0.3219.
#
# Outputs:
# - GSEA Results: Enriched biological processes based on the ranked gene list.
# - Dotplot: A visualization of the top 20 enriched biological processes.
#-------------------------------------------------------------------------------
# Set Working Directory
#-------------------------------------------------------------------------------
setwd('Project_Cep290_GT_AA/')
#-------------------------------------------------------------------------------
# Load Required Libraries
#-------------------------------------------------------------------------------
library(DESeq2)
library(clusterProfiler)
library(org.Mm.eg.db) 
library(tidyverse)

#-------------------------------------------------------------------------------
# Import Differential Expression Results from DESeq2
#-------------------------------------------------------------------------------
# Read in the filtered DESeq2 results, containing significant genes
res <- read.csv('Deseq_results_significant_all_mutant_vs_wildtype_padj0.5_l2fc0.3219 - Filtered.csv')

# Set the 'Gene' column as rownames and remove it from the data frame
res <- res %>% remove_rownames %>% column_to_rownames(var='Gene')
res <- res[, -1]  # Remove any extra column that might exist

#-------------------------------------------------------------------------------
# Prepare Gene List for GSEA
#-------------------------------------------------------------------------------
# Remove any rows with missing values (NA)
res <- na.omit(res)

# Rank genes based on log2 fold change
ranked_genes <- res$log2FoldChange
names(ranked_genes) <- rownames(res)

#-------------------------------------------------------------------------------
# Perform Gene Set Enrichment Analysis (GSEA)
#-------------------------------------------------------------------------------
# Run GSEA using the ranked gene list and the clusterProfiler package
gsea_results <- gseGO(geneList = ranked_genes,
                      OrgDb = org.Mm.eg.db,
                      keyType = "SYMBOL",
                      ont = "BP",  # Biological Process ontology
                      minGSSize = 10,  # Minimum gene set size
                      maxGSSize = 500,  # Maximum gene set size
                      pvalueCutoff = 0.05)  # p-value threshold for significance

#-------------------------------------------------------------------------------
# Visualize GSEA Results
#-------------------------------------------------------------------------------
# Create a dotplot to display the top 20 enriched biological processes
dotplot(gsea_results, showCategory = 20)
