# Description:
# This script processes RNA-Seq data to perform differential gene expression analysis and 
# data visualization. It imports gene expression counts, normalizes the data, and filters 
# lowly expressed genes. The script then performs principal component analysis (PCA) to 
# visualize sample clustering and variation across experimental conditions (e.g., control 
# and Purmorphamine-treated groups). Gene annotations are also incorporated using Ensembl 
# biomart for further analysis. The final output includes processed and normalized gene 
# expression data and a PCA plot visualizing the relationship between samples.
#
# Author: [Your Name]
# Date: [Date]
#
# Inputs:
# - METADATA_WT_Purmorphamine_WT_untreated.txt: Metadata file containing sample information.
# - expression_quant_gene: Directory containing gene expression count files.
# - Annotations_GRCm39/mart_Genes_113.tsv.txt: Gene annotation file from Ensembl.
# - Annotations_GRCm39/Gene_annotations_curated.rds: Curated gene annotation file.
#
# Outputs:
# - cpm_WT_Purmorphamine_WT_untreated.txt: Counts per million (CPM) normalized expression data.
# - lcpm_WT_Purmorphamine_WT_untreated.txt: Log-transformed CPM data for visualization.
# - PCA plot (PCA_plot.png): Visualization of principal component analysis (PCA) results.
#
# Dependencies:
# - R packages: biomaRt, clusterProfiler, tidyverse, pheatmap, RColorBrewer, dplyr, writexl, edgeR, ggplot2.
--------------------------------------------------------------------------------


setwd('Project_Cep290_GT_AA/')
library(biomaRt)
library(clusterProfiler)
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(dplyr)
library(writexl)
library(edgeR)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(edgeR)


#Step 1. Importing metadata from data from Metadata_Files Folder: 
metadata<-read.delim("only_males/Metadata_Files/METADATA_WT_Purmorphamine_WT_untreated.txt", header=TRUE)


#Step 2. Importing gene counts data from expression_quant_gene Folder: 
files <- dir("expression_quant_gene/",  full.names = TRUE)
names(files) <- gsub('\\.genes.results', '', dir("expression_quant_gene/"))
#Filter only your files of interest
files <- files[grepl(paste(metadata$Sample, collapse="|"), files)]

counts<-readDGE(files, columns=c(1,5)) ##tags samples (6 TPM, 5 Expected counts)
colnames(counts)<-names(files)


#Step 3: Organising sample information 
counts$samples$Sample<-rownames(counts$samples)
counts$samples<-left_join(counts$samples,metadata, by = "Sample")
group<-as.factor(counts$samples$Genotype)


#Step4 Gene annotation
#We used Ensembl biomart
#Ensembl Genes 113
#Reference genes, Mouse genes (GRCm39)
#Filters
#All mice chr[1:19,X,Y,MT] ,no alt versions.....
#Attributes
#Gene stable ID
#Gene start (bp)
#Gene end (bp)
#Strand
#Karyotype band
#GO term name
#GO term definition
#GO term accession
#Phenotype description
#Gene name
#Save results as .TSV file




#Step 5 - Normalise read count data

################################################################################
#-----------------------Follow tutorial from here-------------------------------

#Step 4: Data pre-processing: Transformations from the raw-scale 

cpm <- cpm(counts)
lcpm <- cpm(counts, log=TRUE)
# Export this for visualisation 
#write.table(cpm, "outputs_DEG/cpm_WT_Purmorphamine_WT_untreated.txt", sep="\t", row.names = TRUE)
#write.table(lcpm, "outputs_DEG/lcpm_WT_Purmorphamine_WT_untreated.txt", sep="\t", row.names = TRUE)



L <- mean(counts$samples$lib.size) * 1e-6   #(Juliana- 14.68002)    #Mine- 12.6453757683333
M <- median(counts$samples$lib.size) * 1e-6 #(Juliana- 15.87456)    #Mine- 14.55102834
print("Mean and median library size values")
print(summary(counts$samples$lib.size))
#write(summary(counts$samples$lib.size), "outputs_DEG/summary_libsize_WT_Purmorphamine_WT_untreated.txt", sep="\t")

print(c(L, M))
print(summary(lcpm))
#write(summary(lcpm), "outputs_DEG/sumary_lcpm_WT_Purmorphamine_WT_untreated.txt", sep="\t")

print("Max and min lcpm values")
print(max(lcpm))   #Juliana- 14.44895    #Sam- 14.44895
print(min(lcpm))   #Juliana- -2.875783   #Sam- -2.660538

## Step 5. Removing genes that are lowly expressed  ###if genes not expressed in kidneys then they will be lowly expressed. Or genes that arent similarly expressed between replicates.
nsamples<-nrow(counts$samples)
zeros<-table(rowSums(counts$counts==0)==nsamples)
print("Features with zero counts across all samples:")
print(zeros)

#           FALSE  TRUE 
#Juliana -> 33270 44969
#Sam ->     33832 44407
keep.exprs <- filterByExpr(counts, group=group)
counts <- counts[keep.exprs,, keep.lib.sizes=FALSE]
dim(counts) #Juliana -> 17262     6
#Sam -> 14770     6  


lcpm.cutoff <- log2(10/M + 2/L)
nsamples <- ncol(counts)

####work done on Thursday 6th february 2025
Gene_annotations <- read.delim("Annotations_GRCm39/mart_Genes_113.tsv.txt", header = TRUE, sep ="\t")

colnames(Gene_annotations) <- c("Gene.stable.ID","Gene.start", "Gene.end", "Chromosome","Strand","Karyotype.band","Gene.type","Gene.name","Entrezgene.ID")

Gene_annotations_curated <- readRDS("Annotations_GRCm39/Gene_annotations_curated.rds")
counts_noMito <- counts$counts[gsub("\\..*", "", rownames(counts$counts)) %in% Gene_annotations_curated$Gene.stable.ID, ]


cpm_noMito <- cpm(counts_noMito)
lcpm_noMito <- cpm(counts_noMito, log=TRUE)
# Export this for visualisation 
#write.table(cpm_noMito, "outputs_DEG/cpm_WT_Purmorphamine_WT_untreated_noMito.txt", sep="\t", row.names = TRUE)
#write.table(lcpm_noMito, "outputs_DEG/lcpm_WT_Purmorphamine_WT_untreated_noMito.txt", sep="\t", row.names = TRUE)


#The values are same with mito
L_noMito <- mean(counts$samples$lib.size) * 1e-6   #(Juliana- 14.68002)    #Mine- 12.57562613
M_noMito <- median(counts$samples$lib.size) * 1e-6 #(Juliana- 15.87456)    #Mine- 14.467971115
print("Mean and median library size values")
print(summary(counts$samples$lib.size))
#write(summary(counts$samples$lib.size), "outputs_DEG/summary_libsize_WT_Purmorphamine_WT_untreated_noMito.txt", sep="\t")

print(c(L_noMito, M_noMito))
print(summary(lcpm))
#No need to write lib sum
#write(summary(lcpm), "outputs_DEG/sumary_lcpm_WT_Purmorphamine_WT_untreated.txt", sep="\t")

print("Max and min lcpm values")
print(max(lcpm_noMito))      #Sam -> 13.23661
print(min(lcpm_noMito))      #Ben -> -0.1406331

## Step 5. Removing genes that are lowly expressed  ###if genes not expressed in kidneys then they will be lowly expressed. Or genes that arent similarly expressed between replicates.
nsamples<-nrow(counts$samples)
zeros<-table(rowSums(counts_noMito==0)==nsamples)
print("Features with zero counts across all samples:")
print(zeros)

#           FALSE  TRUE 
#Sam ->     14755
# No mito


#Remove genes with zero counts across all samples
counts_noMito_noZeros <- counts_noMito[(rowSums(counts_noMito == 0) == nsamples) == FALSE,]
dim(counts_noMito_noZeros) #14755     6

cpm_noMito_noZeros <- cpm(counts_noMito_noZeros)
lcpm_noMito_noZeros <- cpm(counts_noMito_noZeros, log=TRUE)
# Export this for visualisation 
#write.table(cpm_noMito_noZeros, "outputs_DEG/cpm_WT_Purmorphamine_WT_untreated_noMito_noZeros.txt", sep="\t", row.names = TRUE)
#write.table(lcpm_noMito_noZeros, "outputs_DEG/lcpm_WT_Purmorphamine_WT_untreated_noMito_noZeros.txt", sep="\t", row.names = TRUE)
#write.table(counts_noMito_noZeros, 'outputs_DEG/counts_noMito_noZeros_WT_Purmorphamine_WT_untreated.txt')

#Remove genes with at least 10 counts across all samples
counts_noMito_ten_three <- counts_noMito[(rowSums(counts_noMito > 10)>3),]
dim(counts_noMito_ten_three) #14754     6

cpm_noMito_ten_three<- cpm(counts_noMito_ten_three)
lcpm_noMito_ten_three <- cpm(counts_noMito_ten_three, log=TRUE)
# Export this for visualisation 
#write.table(cpm_noMito_noZeros, "outputs_DEG/cpm_WT_Purmorphamine_WT_untreated_noMito_ten_three.txt", sep="\t", row.names = TRUE)
#write.table(lcpm_noMito_noZeros, "outputs_DEG/lcpm_WT_Purmorphamine_WT_untreated_untreated_noMito_ten_three.txt", sep="\t", row.names = TRUE)
#write.table(counts_noMito_noZeros, 'outputs_DEG/counts_noMito_no_Zeros_WT_Purmorphamine_WT_untreated.txt')


#---------------------------------------------------------------------------------###


#-------------------------------------------------------------------------------
# PCA Analysis of Count Data
#-------------------------------------------------------------------------------
log_counts <- log2(cpm_noMito_ten_three[, c("U14M", "U3M", "U5M", "P12M", "P26M", "P29M")] +1)


log_counts <- log_counts[apply(log_counts, 1, var) > 0, ]  # Remove zero variance rows

# Perform PCA
pca <- prcomp(t(log_counts), scale. = TRUE)  # Transpose so samples are rows
pca_data <- as.data.frame(pca$x)
pca_data$Sample <- rownames(pca_data)  # Add sample names

# Calculate variance explained by each PC
percentVar <- pca$sdev^2 / sum(pca$sdev^2) * 100
pc1_var <- round(percentVar[1], 2)
pc2_var <- round(percentVar[2], 2)

# Create PCA plot
pca_data$Type <- c("Control", "Control", "Control", "Purmorphamine", "Purmorphamine", "Purmorphamine")

library(ggplot2)
PCA_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, label = Sample)) +
  geom_point(aes(color = Type), size = 2) + 
  geom_text(vjust = -1.5, size = 3) +
  labs(title = "PCA of Count Data", x = paste0("PC1: ", pc1_var, "% variance"), 
       y = paste0("PC2: ", pc2_var, "% variance")) +
  theme_minimal() + scale_color_manual(values = c(Control = "#367156", Purmorphamine = "#b87800"))

# Print PCA plot
print(PCA_plot)
