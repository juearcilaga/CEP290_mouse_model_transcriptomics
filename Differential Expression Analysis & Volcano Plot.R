#------------------------------------------------------------------------------
# Differential Expression Analysis & Volcano Plot
# Author: Sam Mindham
# Date: 2025-02-09
# Description: This script performs differential expression analysis using DESeq2
#              and generates a volcano plot for visualization.
# Inputs: 
#   - Rawcounts file: counts_noMito_no_Zeros_Mutant_Purmorphamine_WT_Purmorphamine.txt
#   - Metadata file: METADATA_Mutant_Purmorphamine_WT_Purmorphamine.txt
# Outputs:
#   - CSV files with differential expression results
#       - All results: Deseq_results_Mutant_Purmorphamine_WT_Purmorphamine.csv
#       - Padj <0.05, L2FC > 0.3219: Deseq_results_significant_Mutant_Purmorphamine_WT_Purmorphamine_padj0.5_l2fc0.3219.csv
#   - Volcano plots (SVG format)
#       - With labels: outputs_DEG/Mutant_Purmorphamine_WT_Purmorphamine_Volcano.svg
#       - Without labels: outputs_DEG/Mutant_Purmorphamine_WT_Purmorphamine_Volcano2.svg
#------------------------------------------------------------------------------


setwd('C:\\Users\\sammi\\OneDrive - Newcastle University\\Shared_Folder_UG_proj_2025\\DEG_Analysis')
library(biomaRt)
library(DESeq2)    # For differential expression analysis
library(edgeR)     # For working with count data
library(dplyr)     # For data manipulation
library(ggplot2)   # For plotting
library(ggtext)    # For advanced annotations
library(ggrepel)   #For Advanced annotations
library(svglite)
#-------------------------------------------------------------------------------
#                           Step 1: Import data
#-------------------------------------------------------------------------------
# Load count data
raw_counts_noMito_noZeros<-read.delim("outputs_DEG/counts_noMito_no_Zeros_Mutant_Purmorphamine_WT_Purmorphamine.txt", sep = " ", head=TRUE)

# Load metadata
metadata <- as.data.frame(read.table("Metadata_Files/METADATA_Mutant_Purmorphamine_WT_Purmorphamine.txt", header = TRUE))

#-------------------------------------------------------------------------------
#                   Step 2: Differential expression analysis
#-------------------------------------------------------------------------------
# Create DESeq2 dataset

dds <- DESeqDataSetFromMatrix(countData = round(raw_counts_noMito_noZeros),
                              colData = metadata,
                              design = ~ Genotype)

# Pre-filter genes with low counts (retain genes with at least 10 counts in ???3 samples)
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep, ]


# Define genotype levels for comparison
dds$Genotype <- factor(dds$Genotype, levels = c("GT/AA", "WT/WT"))

## Normalize and estimate size factors
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized = TRUE)

# Perform differential expression analysis
dds <- DESeq(dds)
res <- results(dds, name = "Genotype_WT.WT_vs_GT.AA") 
res$ensembl = rownames(res)

# Order the results by p-value
resOrdered <- res[order(res$pvalue), ]
dfres <- as.data.frame(resOrdered)



# Save full results
#write.csv(dfres, file = "outputs_DEG/Deseq_results_Mutant_Purmorphamine_WT_Purmorphamine.csv")

# Filter significant results (padj < 0.05) and calculate fold change (FC)
resSig <- subset(resOrdered, padj < 0.5 & abs(log2FoldChange) > 0.32192809488)
dim(resSig) # 196 7
dfressig <- as.data.frame(resSig)
dfressig$FC <- 2^dfressig$log2FoldChange

## Uncomment the following for further filtering
## Apply additional filter based on log fold change standard error
# Further filter based on fold change thresholds (|FC| > 2)
# dfressig <- dfressig[dfressig$FC > 2 | dfressig$FC < 0.5, ]
# dfressig <- dfressig[dfressig$lfcSE < 0.5, ]

# Save significant results
write.csv(dfressig, file = "Deseq_results_significant_Mutant_Purmorphamine_WT_Purmorphamine_padj0.5_l2fc0.3219.csv")



# Identify upregulated and downregulated genes
# down <- unique(dfressig[dfressig$FC < 0.5, ]$ensembl)
# up <- unique(dfressig[dfressig$FC > 2, ]$ensembl)

#Save the upregulated and downregulated genes to text files
# write.table(down, "results_down_XXX_deseq.txt", sep = "\t", quote = FALSE, row.names = FALSE)
# write.table(up, "results_up_XXX_deseq.txt", sep = "\t", quote = FALSE, row.names = FALSE)
-------------------------------------------------------------------------------
  # Step 3: Volcano Plot
  #-------------------------------------------------------------------------------
Deseq_results<-dfres%>%filter(padj!=0)

# Define gene significance as Log??? FC > 0.58 (FC > 1.5) 
Deseq_results$Significance <-
  ifelse(
    Deseq_results$padj < 0.05 & abs(Deseq_results$log2FoldChange) > 0.58496250072,
    ifelse(
      Deseq_results$log2FoldChange > 1,
      "Upregulated",
      "Downregulated"
    ),
    "Not Significant"
  )

# Adjust for y-axis values
Deseq_results <- Deseq_results %>%
  mutate(
    neg_log10_padj = -log10(padj),
    adjusted_y = if_else(
      neg_log10_padj > 10,
      10 + log2(1+ neg_log10_padj),
      neg_log10_padj
    ))

Deseq_results$gene_name<-gsub(".*_","",Deseq_results$ensembl)

# Identify most significant genes
top_genes <- Deseq_results %>%
  filter(Significance != "Not Significant") %>%
  arrange(padj) %>%
  slice_head(n = 25) %>%
  bind_rows(
    Deseq_results %>%
      filter(Significance != "Not Significant") %>%
      arrange(padj) %>%
      filter(log2FoldChange < 0) %>%
      slice_head(n = 25)
  )


# Generate Volcano Plot
volcano_plot <-
  ggplot(Deseq_results, aes(x = log2FoldChange, y = adjusted_y)) +
  geom_point(
    aes(fill = Significance, color = Significance),
    size = 3,
    alpha = 0.4,
    shape = 21,
    stroke = 0.9
  ) +
  scale_fill_manual(
    values = c(
      "Upregulated" = "#C0272C",
      "Downregulated" = "#395A88",
      "Not Significant" = "grey"
    )
  ) +
  scale_color_manual(
    values = c(
      "Upregulated" = "#C0272C",
      "Downregulated" = "#395A88",
      "Not Significant" = "grey"
    )
  ) +
  geom_text_repel(
    data = top_genes,
    aes(label = gene_name),
    size = 3,
    max.overlaps = Inf,
    nudge_y = 1
  ) +
  geom_vline(xintercept = c(-1, 1), linetype = "dotted", color = "black", linewidth = 0.8) +  # FC thresholds
  geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "black", linewidth = 0.8) +   # padj threshold
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10 padj (Compressed above 10)") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5))

# Print Volcano plot
print(volcano_plot)

# Save Volcano Plot into svg format, keep same size for consistency
#ggsave(file="outputs_DEG/Mutant_Purmorphamine_WT_Purmorphamine_Volcano.svg", plot=volcano_plot, width=5, height=4)

# Volcano Plot without labels
volcano_plot2 <-
  ggplot(Deseq_results, aes(x = log2FoldChange, y = adjusted_y)) +
  geom_point(
    aes(fill = Significance, color = Significance),
    size = 3,
    alpha = 0.4,
    shape = 21,
    stroke = 0.9
  ) +
  scale_fill_manual(
    values = c(
      "Upregulated" = "#C0272C",
      "Downregulated" = "#395A88",
      "Not Significant" = "grey"
    )
  ) +
  scale_color_manual(
    values = c(
      "Upregulated" = "#C0272C",
      "Downregulated" = "#395A88",
      "Not Significant" = "grey"
    )
  ) +
  geom_vline(xintercept = c(-1, 1), linetype = "dotted", color = "black" ) +  # FC thresholds
  geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "black") +   # padj threshold
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10 padj (Compressed above 10)") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5))
print(volcano_plot2)

#ggsave(file="outputs_DEG/Mutant_Purmorphamine_WT_Purmorphamine_Volcano2.svg", plot=volcano_plot2, width=5, height=4)


