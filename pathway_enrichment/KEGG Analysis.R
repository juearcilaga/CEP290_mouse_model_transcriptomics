#-------------------------------------------------------------------------------
# Description:
# This script performs KEGG pathway enrichment analysis using Entrez gene IDs 
# obtained from gene symbols. The analysis identifies enriched pathways based 
# on a set of genes, and the results are visualized in a dot plot to highlight 
# the top 20 enriched pathways.
#
# Author: Sam Mindham
# Date: 2025/19/02
#
# Inputs:
# - Deseq_results_significant_all_mutant_vs_wildtype_padj0.5_l2fc0.3219 - 
#   Filtered.csv: CSV file containing DESeq2 differential expression results
#   with significant genes based on a padj cutoff of 0.5 and log2 fold change
#   cutoff of 0.3219.
#
# Outputs:
# - KEGG Pathway Enrichment Results: Enriched KEGG pathways based on the provided 
#   gene list.
# - Dotplot: A visualization of the top 20 enriched KEGG pathways.
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Load Required Libraries
#-------------------------------------------------------------------------------
library(clusterProfiler)

# Import Differential Expression Results from DESeq2
#-------------------------------------------------------------------------------
# Read in the filtered DESeq2 results, containing significant genes
res <- read.csv('Inputs/Deseq_results_significant_all_mutant_vs_wildtype_padj0.5_l2fc0.3219 - Filtered.csv')

# Set the 'Gene' column as rownames and remove it from the data frame
res <- res %>% remove_rownames %>% column_to_rownames(var='Gene')
res <- res[, -1]  # Remove any extra column that might exist

#-------------------------------------------------------------------------------
# Prepare Gene List for KEGG
#-------------------------------------------------------------------------------
# Remove any rows with missing values (NA)
res <- na.omit(res)

# Rank genes based on log2 fold change
ranked_genes <- res$log2FoldChange
names(ranked_genes) <- rownames(res)
#-------------------------------------------------------------------------------
# Map Gene Symbols to Entrez IDs
#-------------------------------------------------------------------------------
# Extract gene symbols from the row names of the result data frame
genes_symbols <- rownames(res)

# Map gene symbols to Entrez IDs using org.Mm.eg.db package
entrez_ids <- mapIds(org.Mm.eg.db, 
                     keys = genes_symbols, 
                     column = "ENTREZID",  # Map to Entrez gene IDs
                     keytype = "SYMBOL",   # Gene symbols as the input keys
                     multiVals = "first")  # Select the first Entrez ID in case of multiple matches

#-------------------------------------------------------------------------------
# Remove Unmapped Genes (NAs)
#-------------------------------------------------------------------------------
# Remove any NA values from the Entrez ID list, which correspond to genes that 
# could not be mapped to an Entrez ID
entrez_ids <- na.omit(entrez_ids)

#-------------------------------------------------------------------------------
# Perform KEGG Pathway Enrichment Analysis
#-------------------------------------------------------------------------------
# Run KEGG pathway enrichment analysis using the Entrez gene IDs
kegg_results <- enrichKEGG(gene = entrez_ids,
                           organism = "mmu",  # Mouse organism code
                           keyType = "kegg",  # KEGG pathway ID
                           pvalueCutoff = 0.05)  # p-value threshold for significance

#-------------------------------------------------------------------------------
# Visualize KEGG Enrichment Results
#-------------------------------------------------------------------------------
# Create a dotplot to display the top 20 enriched KEGG pathways
dotplot(kegg_results, showCategory = 20, font.size = 8)
