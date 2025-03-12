#-------------------------------------------------------------------------------
# Description:
# This script filters a set of genes related to various stress responses 
# (oxidative stress, ER stress, hypoxic stress, and inflammatory stress) from 
# differential expression results. The analysis involves identifying genes 
# associated with each stress type and further filtering based on significance 
# criteria (adjusted p-value < 0.5 and log2 fold change > 0.32). The results are 
# categorized by stress type and prepared for downstream analysis.
#
# Author: Sam Mindham
# Date: 2025/19/02
#
# Inputs:
# - Mutant_vs_WT: A data frame containing differential expression results, 
#   with gene symbols as row names. The dataset is derived from Untreated Mutant
#   vs Untreated Wild Type Deseq Results
#
# Outputs:
# - Significant_Stress.csv: Subsets of the differential expression results 
#   containing genes related to oxidative stress, ER stress, hypoxic stress, 
#   and inflammatory stress, with additional information on the stress type.
#-------------------------------------------------------------------------------
# Set working directory
setwd('Project_Cep290_GT_AA/')

# Load required libraries
library(tidyverse)

# Read DESeq2 results
Mutant_vs_WT <- read.csv('Deseq_results_Mutant_vs_WT.csv')

# Clean up gene names
Mutant_vs_WT$X <- sub('.*_', '', Mutant_vs_WT$X)
Mutant_vs_WT$X <- make.unique(Mutant_vs_WT$X)
Mutant_vs_WT <- Mutant_vs_WT %>% remove_rownames() %>% column_to_rownames(var='X')



# Define oxidative stress-related genes
oxidative_stress_genes <- c(
  "ABCB11", "ABCC1", "ABCD1", "ABL1", "ADAM9", "ADCYAP1R1", "ADIPOQ", "ADPRS", "AGAP3", "AIF1", 
  "AIFM1", "AIFM2", "AKT1", "ALAD", "ALDH3B1", "ALOX5", "ALS2", "ANGPTL7", "ANKZF1", "APOA4", 
  "APOD", "APOE", "APP", "AQP1", "AREG", "ARL6IP5", "ARNT", "ATF2", "ATF4", "ATM", "ATOX1", "ATP13A2", 
  "ATP2A2", "ATP7A", "AXL", "BAK1", "BANF1", "BCL2", "BECN1", "BMAL1", "BMP7", "BNIP3", "BRF2", "BTK", 
  "C19orf12", "CAMKK2", "CAPN2", "CASP3", "CAT", "CBX8", "CCS", "CD2AP", "CD36", "CD38", "CDK1", "CHCHD2", 
  "CHCHD4", "CHRNA4", "CHUK", "COA8", "COL1A1", "COL6A1", "COMT", "CRK", "CRYAB", "CRYGD", "CYB5B", "CYGB", 
  "CYP1B1", "DAPK1", "DDR2", "DGKK", "DHFR", "DHFRP1", "DHRS2", "DUOX1", "DUOX2", "ECT2", "EDN1", "EDNRA", 
  "EGFR", "EIF2S1", "EPAS1", "EPX", "ERCC1", "ERCC2", "ERCC3", "ERCC6", "ERCC6L2", "ERCC8", "ERMP1", "ERN1", 
  "ERO1A", "ETFDH", "ETV5", "EZH2", "FABP1", "FADS2", "FANCC", "FANCD2", "FBLN5", "FER", "FGF8", "FKBP1B", 
  "FOS", "FOSL1", "FOXO1", "FOXO3", "FOXO4", "FOXP1", "FUT8", "FXN", "FYN", "G6PD", "GCH1", "GCLC", "GCLM", 
  "GGT7", "GJB2", "GLRX2", "GPR37", "GPR37L1", "GPX1", "GPX2", "GPX3", "GPX4", "GPX5", "GPX6", "GPX7", 
  "GPX8", "GSR", "GSS", "GSTP1", "HAO1", "HBA1", "HBA2", "HBB", "HDAC2", "HDAC6", "HGF", "HIF1A", "HM13", 
  "HMOX1", "HMOX2", "HP", "HSF1", "HSPA1A", "HSPA1B", "HSPA8", "HTRA2", "HYAL1", "HYAL2", "IDH1", "IL18BP", 
  "IL18RAP", "IL1A", "IL6", "IPCEF1", "JAK2", "JUN", "KAT2B", "KCNA5", "KDM6B", "KEAP1", "KLF2", "KLF4", 
  "KRT1", "LIAS", "LONP1", "LPO", "LRRK2", "MACROH2A1", "MAP1LC3A", "MAP2K4", "MAP3K5", "MAPK1", "MAPK13", 
  "MAPK3", "MAPK7", "MAPK8", "MAPK9", "MAPKAP1", "MAPT", "MB", "MBL2", "MCTP1", "MDM2", "MEAK7", "MET", 
  "MGAT3", "MGST1", "MICB", "MIR103A1", "MIR107", "MIR132", "MIR133A1", "MIR135A1", "MIR17", "MIR21", 
  "MIR34A", "MIR92A1", "MIRLET7B", "MMP14", "MMP2", "MMP3", "MMP9", "MPO", "MPV17", "MSRA", "MSRB2", 
  "MSRB3", "MT-CO1", "MT-ND1", "MT-ND3", "MT-ND5", "MT-ND6", "MT3", "MTF1", "MYB", "NAGLU", "NAPRT", "NCF1", 
  "NCOA7", "NDUFA12", "NDUFA6", "NDUFB4", "NDUFS2", "NDUFS8", "NEIL1", "NET1", "NFE2L2", "NME8", "NOS3", 
  "NQO1", "NR4A2", "NUDT1", "NUDT15", "NUDT2", "OGG1", "OSER1", "OXR1", "OXSR1", "PARK7", "PARP1", "PAWR", 
  "PCGF2", "PCNA", "PDCD10", "PDGFD", "PDGFRA", "PDK2", "PDLIM1", "PENK", "PEX10", "PEX12", "PEX13", "PEX14", 
  "PEX2", "PEX5", "PINK1", "PJVK", "PKD2", "PLA2R1", "PLEKHA1", "PLK3", "PNKP", "PNPLA8", "PNPT1", "PPARGC1A", 
  "PPARGC1B", "PPIA", "PPIF", "PPP1R15B", "PPP2CB", "PRDX1", "PRDX2", "PRDX3", "PRDX4", "PRDX5", "PRDX6", 
  "PRKAA1", "PRKAA2", "PRKCD", "PRKD1", "PRKN", "PRKRA", "PRNP", "PRR5L", "PSEN1", "PSIP1", "PSMB5", "PTGS1", 
  "PTGS2", "PTK2B", "PTPRK", "PTPRN", "PXDN", "PXDNL", "PXN", "PYCR1", "PYCR2", "PYROXD1", "RACK1", "RAD52", 
  "RBM11", "RBPMS", "RBX1", "RCAN1", "RELA", "RGS14", "RHOB", "RIPK1", "RIPK3", "RLIG1", "RNF112", "ROMO1", 
  "RPS3", "RRM2B", "RWDD1", "S100A7", "SCARA3", "SCGB1A1", "SDC1", "SELENOK", "SELENON", "SELENOP", "SELENOS", 
  "SESN1", "SESN2", "SESN3", "SETX", "SGK2", "SIN3A", "SIRPA", "SIRT1", "SIRT2", "SLC1A1", "SLC23A2", "SLC25A14", 
  "SLC25A24", "SLC4A11", "SLC7A11", "SMPD3", "SNCA", "SOD1", "SOD2", "SOD3", "SP1", "SPHK1", "SRC", "SRXN1", 
  "STAT1", "STAT6", "STAU1", "STAU2", "STC2", "STK24", "STK25", "STK26", "STOX1", "STX2", "STX4", "SUMO4", 
  "TACR1", "TAT", "TBC1D24", "TET1", "THG1L", "TLDC2", "TMEM161A", "TMIGD1", "TNFAIP3", "TOP2B", "TOR1A", 
  "TP53", "TP53INP1", "TPM1", "TPO", "TRAP1", "TREX1", "TRIM25", "TRIM32", "TRPA1", "TRPM2", "TXN", "TXN2", 
  "TXNIP", "TXNRD2", "UBE3A", "UCN", "UCP1", "UCP2", "UCP3", "USP25", "VKORC1L1", "VNN1", "VRK2", "WNT16", 
  "WRN", "XRCC1", "ZC3H12A", "ZFAND1", "ZNF277", "ZNF580"
)

# Convert gene names to proper case
oxidative_stress_genes <- tools::toTitleCase(tolower(oxidative_stress_genes))

# Filter oxidative stress genes from DESeq2 results
de_oxidative <- Mutant_vs_WT %>% filter(rownames(Mutant_vs_WT) %in% oxidative_stress_genes)

# Further filter for significance and fold change
filtered_de_oxidative <- de_oxidative %>%
  filter(padj < 0.5, log2FoldChange > 0.32) %>%
  mutate(Stress_Type = 'Oxidative')

# Define Endoplasmic Reticulum Stress related genes
ER_stress_genes <- c(
  "ABCA7", "ABL1", "AFF4", "AGR2", "AIFM1", "ALOX15", "ALOX5", "AMFR", "ANKS4B", "ANKZF1",
  "APAF1", "AQP11", "ATF3", "ATF4", "ATF6", "ATF6B", "ATG10", "ATP2A1", "ATP2A2", "ATP2A3",
  "ATXN3", "AUP1", "BAG6", "BAK1", "BAX", "BBC3", "BCAP31", "BCL2", "BCL2L1", "BCL2L11",
  "BFAR", "BHLHA15", "BOK", "BRSK2", "CALR", "CALR3", "CANX", "CASP4", "CAV1", "CCDC47",
  "CCND1", "CDK5RAP3", "CEBPB", "CERT1", "CFTR", "CHAC1", "CLGN", "CLU", "COPS5", "CREB3",
  "CREB3L1", "CREB3L2", "CREB3L3", "CREBZF", "CTH", "CXCL8", "DAB2IP", "DDIT3", "DDRGK1",
  "DDX3X", "DERL1", "DERL2", "DERL3", "DNAJB12", "DNAJB2", "DNAJB9", "DNAJC10", "DNAJC3",
  "ECPAS", "EDEM1", "EDEM2", "EDEM3", "EIF2AK2", "EIF2AK3", "EIF2AK4", "EIF2B5", "EIF2S1",
  "EIF4G1", "ELAVL4", "ERLEC1", "ERLIN1", "ERLIN2", "ERMP1", "ERN1", "ERN2", "ERO1A", "ERP27",
  "ERP29", "ERP44", "FAF1", "FAF2", "FAM8A1", "FBXO17", "FBXO2", "FBXO27", "FBXO44", "FBXO6",
  "FCGR2B", "FGF21", "FICD", "FLOT1", "FOXRED2", "GET4", "GORASP2", "GRINA", "GSK3B", "HERPUD1",
  "HERPUD2", "HM13", "HSP90B1", "HSPA1A", "HSPA5", "HYOU1", "ITPR1", "JKAMP", "JUN", "KCNJ8",
  "LPCAT3", "LRRK2", "MAGEA3", "MAN1A1", "MAN1A2", "MAN1B1", "MAN1C1", "MANF", "MAP3K5", "MARCHF6",
  "MARCKS", "MBTPS1", "MBTPS2", "MIR199A1", "MIR200C", "NCCRP1", "NCK1", "NCK2", "NFE2L2", "NHLRC1",
  "NIBAN1", "NOD1", "NOD2", "NPLOC4", "NR1H2", "NR1H3", "NRBF2", "NUPR1", "OPA1", "OS9", "P4HB",
  "PARK7", "PARP16", "PARP6", "PARP8", "PDIA2", "PDIA3", "PDIA4", "PDIA6", "PDX1", "PIGBOS1",
  "PIK3R1", "PIK3R2", "PMAIP1", "PML", "PPP1R15A", "PPP1R15B", "PPP2CB", "PRKN", "PSMC6", "PTPN1",
  "PTPN2", "QRICH1", "RACK1", "RASGRF1", "RASGRF2", "RCN3", "RHBDD1", "RHBDD2", "RNF103", "RNF121",
  "RNF139", "RNF175", "RNF183", "RNF185", "RNF186", "RNF5", "RNFT1", "RNFT2", "RPAP2", "SCAMP5",
  "SDF2L1", "SEC16A", "SEC61B", "SEL1L", "SEL1L2", "SELENOK", "SELENOS", "SERINC3", "SERP1", "SERP2",
  "SESN2", "SGF29", "SGTA", "SIRT1", "SRPX", "STC2", "STT3B", "STUB1", "SVIP", "SYVN1", "TARDBP",
  "TBL2", "THBS1", "THBS4", "TMBIM6", "TMCO1", "TMED2", "TMEM117", "TMEM129", "TMEM238L", "TMEM258",
  "TMEM259", "TMEM33", "TMEM67", "TMTC3", "TMTC4", "TMUB1", "TMUB2", "TMX1", "TNFRSF10B", "TOR1A",
  "TP53", "TRAF2", "TRIB3", "TRIM13", "TRIM25", "TTC23L", "TXNDC12", "UBA5", "UBAC2", "UBE2G2",
  "UBE2J1", "UBE2J2", "UBE4A", "UBE4B", "UBQLN1", "UBQLN2", "UBXN1", "UBXN10", "UBXN2A", "UBXN4",
  "UBXN6", "UBXN8", "UFC1", "UFD1", "UFL1", "UFM1", "UGGT1", "UGGT2", "UMOD", "USP13", "USP14",
  "USP19", "USP25", "VAPB", "VCP", "WFS1", "XBP1", "YOD1"
)

# Convert gene names to proper case
ER_stress_genes <- tools::toTitleCase(ER_stress_genes)

# Filter ER stress genes from DESeq2 results
de_ER <- Mutant_vs_WT %>% filter(rownames(Mutant_vs_WT) %in% ER_stress_genes)

# Further filter for significance and fold change
filtered_de_ER <- de_ER %>%
  filter(padj < 0.5, log2FoldChange > 0.32) %>%
  mutate(Stress_Type = 'ER')


#Define Hypoxic Stress related genes
hypoxic_stress_genes <- c(
  "Acaa2", "Adam8", "Ado", "Adrb2", "Ajuba", "Ak4", "Akt1", "Aqp1", "Aqp3", "Atf2", 
  "Bad", "Bcl2", "Bmyc", "Bnip3", "Bnip3l", "Cbs", "Chchd2", "Chchd2-ps", "Cited2", 
  "Commd1", "Cpeb1", "Cpeb2", "Cr1l", "Cybb", "Daxx", "Ddah1", "Dnmt3a", "Drd1", 
  "Drd2", "Egln1", "Egln2", "Egln3", "Eif4ebp1", "Eno1", "Eno1b", "Ep300", "Epas1", 
  "Epha4", "Fabp1", "Fam162a", "Fmn2", "Fndc1", "Gata6", "Gnb1", "Gngt1", "Hif1a", 
  "Hif3a", "Higd1a", "Hp1bp3", "Hyou1", "Ireb2", "Kcnd2", "Kcnk2", "Kcnk3", "Kdm6a", 
  "Map2k1", "Mgarp", "Mir199a-2", "Mir214", "Mir668", "Mir874", "Mlst8", "Mtor", "Myc", 
  "Ndnf", "Ndp", "Nfe2l2", "Ngb", "Nkx3-1", "Nol3", "Nop53", "Notch1", "Npepps", "Ogt", 
  "Oprd1", "P4hb", "Pdk1", "Pdk3", "Pgk1", "Pik3cb", "Pink1", "Plk3", "Ppard", "Pparg", 
  "Prkce", "Pten", "Ptgis", "Rgcc", "Rock2", "Rora", "Rptor", "Rtn4", "Rwdd3", "Scn2a", 
  "Sdhd", "Sirt1", "Sirt2", "Slc2a4", "Slc8a3", "Slc9a1", "Stat3", "Stub1", "Suv39h1", 
  "Suv39h2", "Tbl2", "Tert", "Tgfb1", "Tigar", "Tmbim6", "Trem2", "Trp53", "Twist1", 
  "Ubqln1", "Usp19", "Vasn", "Vegfa", "Vhl", "Vldlr", "Zfas1", "Zfp36l1"
)

# Filter Hypoxic stress genes from DESeq2 results
de_hypoxic <- Mutant_vs_WT %>% filter(rownames(Mutant_vs_WT) %in% hypoxic_stress_genes)

# Further filter for significance and fold change
filtered_de_hypoxic <- de_hypoxic %>%
  filter(padj < 0.5, log2FoldChange > 0.32) %>%
  mutate(Stress_Type = 'Hypoxic')

#Define Inflammatory Stress Related Genes
inflammatory_stress_genes <- c(
  "A2M", "ABCC1", "ABCD1", "ABCD2", "ABCF1", "ABHD12", "ABHD17A", "ACE2", "ACER3", "ACKR1",
  "ACKR2", "ACOD1", "ACP5", "ACVR1", "ADA", "ADAM8", "ADAMTS12", "ADCY1", "ADCY7", "ADCY8",
  "ADGRE2", "ADGRE5", "ADIPOQ", "ADM", "ADORA1", "ADORA2A", "ADORA2B", "ADORA3", "AFAP1L2", "AGER",
  "AGR2", "AGT", "AGTR1", "AGTR2", "AHR", "AHSG", "AIF1", "AIM2", "AIMP1", "AKNA", "AKT1", "ALOX15",
  "ALOX5", "ALOX5AP", "ANO6", "ANXA1", "AOAH", "AOC3", "AP3B1", "APCS", "APOA1", "APOD", "APOE",
  "APOL2", "APOL3", "APP", "APPL1", "APPL2", "AREL1", "ARMH4", "ARNT", "ASH1L", "ASS1", "ATAT1", "ATM",
  "ATRN", "AXL", "AZU1", "B4GALT1", "BAP1", "BCL6", "BCR", "BDKRB1", "BDKRB2", "BIRC2", "BIRC3", "BLNK",
  "BMP2", "BMP6", "BMPR1B", "BPGM", "BRCC3", "BRD4", "BST1", "BTK", "C1QA", "C1QTNF12", "C1QTNF3",
  "C2CD4A", "C2CD4B", "C3", "C3AR1", "C4A", "C4B", "C5", "C5AR1", "C5AR2", "CALCA", "CALHM2", "CAMK1D",
  "CAMK2N1", "CAMK4", "CARD16", "CARD18", "CARD8", "CARD9", "CASP1", "CASP12", "CASP4", "CASP5", "CCL1",
  "CCL11", "CCL13", "CCL14", "CCL15", "CCL16", "CCL17", "CCL18", "CCL19", "CCL2", "CCL20", "CCL21",
  "CCL22", "CCL23", "CCL24", "CCL25", "CCL26", "CCL3", "CCL3L3", "CCL4", "CCL4L2", "CCL5", "CCL7",
  "CCL8", "CCN3", "CCN4", "CCR1", "CCR2", "CCR3", "CCR4", "CCR5", "CCR6", "CCR7", "CCRL2", "CD14",
  "CD163", "CD180", "CD200", "CD200R1", "CD200R1L", "CD28", "CD2AP", "CD36", "CD40", "CD40LG", "CD44",
  "CD47", "CD5L", "CD6", "CD68", "CD81", "CD96", "CDH5", "CDO1", "CEBPA", "CEBPB", "CELA1", "CELF1",
  "CERS6", "CHI3L1", "CHIA", "CHID1", "CHST1", "CHST2", "CHST4", "CHUK", "CLEC7A", "CLOCK", "CLU",
  "CMA1", "CMKLR1", "CNR1", "CNR2", "CNTF", "COL6A1", "CPTP", "CREB3L3", "CRH", "CRHBP", "CRP", "CSF1",
  "CSF1R", "CSNK1A1", "CSRP3", "CST7", "CTNNBIP1", "CTSC", "CUEDC2", "CUL3", "CX3CL1", "CX3CR1", "CXCL1",
  "CXCL10", "CXCL11", "CXCL13", "CXCL17", "CXCL2", "CXCL3", "CXCL5", "CXCL6", "CXCL8", "CXCL9", "CXCR2",
  "CXCR3", "CXCR4", "CXCR6", "CYBA", "CYBB", "CYLD", "CYP19A1", "CYP26B1", "CYSLTR1", "DAB2IP", "DAGLA",
  "DAGLB", "DDT", "DDX3X", "DEFB114", "DHRS7B", "DHX33", "DHX9", "DNASE1", "DNASE1L3", "DPEP1", "DROSHA",
  "DUOXA1", "DUOXA2", "DUSP10", "ECM1", "EDNRB", "EGFR", "EIF2AK1", "EIF2AK2", "ELANE", "ELF3", "ELF4",
  "ENPP3", "EPHA2", "EPHB2", "EPO", "ESR1", "ETS1", "EXT1", "EXTL3", "EZH2", "F11R", "F12", "F2", "F2R",
  "F2RL1", "F3", "F8", "FABP4", "FADS2", "FANCA", "FANCD2", "FASN", "FBXL2", "FCER1A", "FCGR1A", "FCGR1BP",
  "FCGR2A", "FCGR2B", "FCGR2C", "FCGR3A", "FCGR3B", "FEM1A", "FFAR2", "FFAR3", "FFAR4", "FGR", "FKRP",
  "FN1", "FNDC4", "FOLR2", "FOS", "FOSL1", "FOSL2", "FOXF1", "FOXP1", "FOXP3", "FPR1", "FPR2", "FPR3",
  "FURIN", "FUT4", "FUT7", "FXR1", "FYN", "GATA3", "GBA1", "GBP2", "GBP5", "GGT1", "GGT2P", "GGT3P", "GGT5",
  "GHRL", "GHSR", "GIT1", "GNAT2", "GPER1", "GPR17", "GPR31", "GPR32", "GPR32P1", "GPR33", "GPR4", "GPR68",
  "GPRC5B", "GPS2", "GPSM3", "GPX1", "GRN", "GSDMD", "GSTP1", "H2BC1", "HAMP", "HAVCR2", "HCK", "HDAC4",
  "HDAC5", "HDAC9", "HGF", "HIF1A", "HK1", "HLA-DRB1", "HLA-E", "HMGB1", "HMGB2", "HMOX1", "HNRNPA0", "HP",
  "HRH1", "HRH4", "HSPA8", "HSPG2", "HTR2A", "HYAL1", "HYAL2", "HYAL3", "IDO1", "IER3", "IFI16", "IFI35",
  "IFNA2", "IFNG", "IFNGR1", "IFNGR2", "IGF1", "IGHE", "IGHG1", "IKBKB", "IKBKG", "IL10", "IL10RA", "IL10RB",
  "IL12B", "IL13", "IL15", "IL16", "IL17A", "IL17B", "IL17C", "IL17D", "IL17F", "IL17RA", "IL17RC", "IL17RE",
  "IL18", "IL18R1", "IL18RAP", "IL1A", "IL1B", "IL1F10", "IL1R1", "IL1R2", "IL1RAP", "IL1RL1", "IL1RL2",
  "IL1RN", "IL2", "IL20", "IL20RB", "IL21", "IL22", "IL22RA1", "IL22RA2", "IL23A", "IL23R", "IL25", "IL27",
  "IL2RA", "IL31RA", "IL33", "IL34", "IL36A", "IL36B", "IL36G", "IL36RN", "IL37", "IL4", "IL4R", "IL5",
  "IL5RA", "IL6", "IL6R", "IL6ST", "IL9", "INS", "IRAK2", "IRF3", "IRF5", "IRGM", "ISL1", "ITCH", "ITGAL",
  "ITGAM", "ITGB1", "ITGB2", "ITGB6", "ITIH4", "JAK2", "JAM3", "KARS1", "KCNJ8", "KDM4D", "KDM6B", "KIT",
  "KLF4", "KLKB1", "KLRG1", "KNG1", "KPNA6", "KRT1", "KRT16", "LACC1", "LARGE1", "LAT", "LCK", "LDLR",
  "LEF1", "LGR5", "LIF", "LILRA1", "LITAF", "LPS", "LTA", "LTA4H", "LTF", "LXRα", "LXRβ", "MAL", "MAP3K7",
  "MAPK1", "MAPK3", "MAPK8", "MAPK9", "MAPK10", "MAPK14", "MARCH9", "MALT1", "MAP2K4", "MCP-1", "MDA5",
  "MELK", "MIF", "MMP2", "MMP9", "MMP12", "MMP13", "MMP14", "MMP15", "MMP16", "MMP19", "MMP3", "MMP7",
  "MMP8", "MMP9", "MSC", "MS4A1", "MSL1", "MT1A", "MT2A", "MTA1", "MTHFR", "NAB2", "NLRP1", "NLRP3",
  "NLRP6", "NR4A1", "NR4A2", "NRF2", "NUPR1", "OAS1", "OAS2", "OAS3", "ODC1", "OLFM1", "OMT2", "PADI4",
  "PARP1", "PAWR", "PBX1", "PCDHA10", "PCDHA5", "PCDHB4", "PCLAF", "PDGFB", "PDK1", "PGC1A", "PGD",
  "PFN1", "PFKFB3", "PHF6", "PLXDC2", "PMAIP1", "PPARG", "PRKCA", "PRKCD", "PRKCZ", "PTGS2", "PTK2",
  "PTH", "PTPN11", "PTPRC", "PXDN", "PYCARD", "QKI", "RAC1", "RASA1", "RASGRF1", "RAP1B", "RAPGEF2",
  "RELB", "REL", "RELB", "RGS4", "RIPK1", "RIPK2", "RLBP1", "RHOA", "RPS3A", "RPS6KA1", "RTKN", "SAA1",
  "SAA2", "SAA3", "SAMD9", "SASH1", "SELE", "SEMA3A", "SEMA3B", "SEMA3C", "SEMA3D", "SERPINE1", "SERPINB2",
  "SH2D2A", "SH2D3C", "SH3BP5", "SIRPB1", "SIRT1", "SIRT2", "SIRPA", "SKIL", "SLAMF1", "SLAMF7", "SLC12A3",
  "SLC16A1", "SLC22A1", "SLC22A5", "SLC26A1", "SLC29A1", "SLC34A1", "SLC4A4", "SLC40A1", "SLC41A1",
  "SLC51A", "SLC51B", "SLC6A20", "SLC7A1", "SLC9A1", "SLC9A3R1", "SMAD3", "SMAD4", "SMURF1", "SOCS1",
  "SOCS3", "SOD2", "SOD3", "SST", "SSX2", "STAT1", "STAT3", "STAT5A", "STAT5B", "STAT6", "STING1", "SYNJ1",
  "TGFBR1", "TGFBR2", "TLR1", "TLR2", "TLR3", "TLR4", "TLR5", "TLR6", "TLR7", "TLR8", "TLR9", "TNF",
  "TNFRA", "TNFRSF11A", "TNFRSF1A", "TNFRSF10A", "TNFRSF10B", "TNFRSF12A", "TNFRSF14", "TNFRSF18",
  "TNFRSF21", "TNFRSF25", "TRAF3", "TRAF6", "TREM1", "TREM2", "TRIM21", "TRIM5", "TRIP10", "TSG101",
  "TYK2", "UBE2L6", "VDR", "VEGFA", "VEGFB", "VIM", "VPS18", "VPS37A", "WAS", "XAF1", "XBP1"
)

# Convert gene names to proper case
inflammatory_stress_genes <- tolower(inflammatory_stress_genes)
inflammatory_stress_genes <- tools::toTitleCase(inflammatory_stress_genes)

# Filter inflammatory stress genes from DESeq2 results
de_inflam <- Mutant_vs_WT %>% filter(rownames(Mutant_vs_WT) %in% inflammatory_stress_genes)

# Further filter for significance and fold change
filtered_de_inflam <- de_inflam %>%
  filter(padj < 0.5, log2FoldChange > 0.32) %>%
  mutate(Stress_Type = 'Inflammatory')


#Define Chemical Stress related genes
chemical_stress_genes <- c(
  "ABCC1", "ABCD1", "ABL1", "ADPRS", "AGAP3", "AIF1", "AIFM1", "AIFM2", "AKR1B1", "AKT1", 
  "ALDH3B1", "ALOX5", "ANKZF1", "APOA4", "AQP1", "AQP5", "ARHGEF2", "ARL6IP5", "ARNT", "ATF2", 
  "ATF4", "ATG5", "ATG7", "ATM", "ATP13A2", "ATP1A1", "ATP2A2", "ATP7A", "AXL", "BAD", 
  "BDKRB2", "BECN1", "BMAL1", "BMP7", "BNIP3", "BRF2", "BTK", "CAB39", "CAMKK2", "CAPN3", 
  "CASP1", "CASP3", "CAT", "CAV1", "CBX8", "CCS", "CD36", "CDK1", "CHCHD2", "CHCHD4", "CHUK", 
  "CLN3", "CRYGD", "CYGB", "CYP1B1", "DAPK1", "DDIT3", "DDR2", "DDX3X", "DHFR", "DHFRP1", 
  "DHRS2", "DNAJA1", "ECT2", "EDN1", "EDNRA", "EFHD1", "EGFR", "EIF2AK3", "EIF2S1", "EPO", 
  "ERCC6L2", "ERMP1", "ERN1", "ETV5", "EZH2", "FABP1", "FADS2", "FANCC", "FANCD2", "FAS", 
  "FBLN5", "FBP1", "FER", "FOS", "FOXO1", "FOXO3", "FOXP1", "FUT8", "FXN", "FYN", "G6PD", 
  "GCH1", "GDF15", "GJB2", "GLRX2", "GPR37", "GPR37L1", "GPX1", "GPX5", "GPX7", "GPX8", "GSR", 
  "HDAC2", "HDAC6", "HGF", "HIF1A", "HM13", "HSF1", "HSPA1A", "HSPA1B", "HSPA8", "HTRA2", 
  "IL18BP", "IL18RAP", "IL6", "JUN", "KAT2B", "KDM6B", "KEAP1", "KLF2", "KLF4", "LETM1", 
  "LONP1", "LRRC8C", "LRRC8D", "LRRC8E", "LRRK2", "MAP1LC3A", "MAP2K4", "MAP3K5", "MAPK1", 
  "MAPK13", "MAPK3", "MAPK7", "MAPK8", "MAPK9", "MAPKAP1", "MAPT", "MB", "MDM2", "MEAK7", 
  "MET", "MGAT3", "MGST1", "MICU1", "MIR103A1", "MIR107", "MIR132", "MIR133A1", "MIR135A1", 
  "MIR17", "MIR21", "MIR34A", "MIR92A1", "MIRLET7B", "MLST8", "MMP2", "MMP3", "MMP9", "MPO", 
  "MPV17", "MSRA", "MT3", "MTOR", "MYB", "MYLK", "NAGLU", "NCF1", "NCOA7", "NET1", "NFAT5", 
  "NFE2L2", "NLK", "NLRP3", "NME8", "NOS3", "NOX1", "NQO1", "NR4A2", "NUDT2", "OGG1", "OSER1", 
  "OXR1", "OXSR1", "PARK7", "PARP1", "PAWR", "PCGF2", "PCK1", "PCNA", "PDCD10", "PDGFD", "PDGFRA", 
  "PDK2", "PENK", "PEX10", "PEX12", "PEX13", "PEX14", "PEX2", "PEX5", "PIK3CA", "PINK1", "PJVK", 
  "PKD2", "PLA2R1", "PLEC", "PLEKHA1", "PNPLA8", "PNPT1", "PPARGC1A", "PPARGC1B", "PPIA", "PPIF", 
  "PRDX1", "PRDX2", "PRDX3", "PRDX5", "PRKAA1", "PRKAA2", "PRKCD", "PRKCI", "PRKD1", "PRKN", 
  "PRKRA", "PRR5L", "PTGS2", "PTPRK", "PXN", "PYCARD", "PYCR1", "PYCR2", "PYROXD1", "RACK1", 
  "RAD52", "RBM11", "RBX1", "RCSD1", "RELA", "RELB", "RHOB", "RIPK1", "RIPK3", "ROMO1", "RPS3", 
  "RPTOR", "RWDD1", "SCN2A", "SCN7A", "SELENON", "SELENOS", "SERPINB6", "SESN2", "SETX", "SIN3A", 
  "SIRPA", "SIRT1", "SIRT2", "SLC12A6", "SLC1A1", "SLC25A14", "SLC25A23", "SLC25A24", "SLC2A1", 
  "SLC2A4", "SLC4A11", "SLC7A11", "SMPD3", "SNCA", "SOD1", "SOD2", "SOD3", "SPHK1", "SRC", 
  "SRXN1", "STAT6", "STAU1", "STAU2", "STK24", "STK25", "STK26", "STK39", "STOX1", "STX2", 
  "STX4", "SUMO4", "TBC1D24", "TET1", "TIFAB", "TMEM161A", "TMIGD1", "TNFAIP3", "TOP2B", "TP53", 
  "TP53INP1", "TPM1", "TRAP1", "TREX1", "TRIM21", "TRPA1", "TRPM2", "TRPV4", "TSPO", "TXN", 
  "UCP1", "USP15", "VKORC1L1", "VRK2", "WNK1", "WNK3", "WNT16", "XRCC5", "XRCC6", "YBX3", 
  "ZC3H12A", "ZFAND1", "ZFP36L1", "ZNF277", "ZNF580"
)

# Convert gene names to proper case
chemical_stress_genes <- tolower(chemical_stress_genes)
chemical_stress_genes <- tools::toTitleCase(chemical_stress_genes)

# Filter chemical stress genes from DESeq2 results
de_chemical <- Mutant_vs_WT %>% filter(rownames(Mutant_vs_WT) %in% chemical_stress_genes)

# Further filter for significance and fold change
filtered_de_chemical <- de_chemical %>%
  filter(padj < 0.5, log2FoldChange > 0.32) %>%
  mutate(Stress_Type = 'Chemical')

# Combine all filtered stress-related genes
significant_stress <- rbind(filtered_de_chemical, filtered_de_ER, filtered_de_hypoxia, filtered_de_inflam, filtered_de_oxidative)
significant_stress <- significant_stress %>%
  distinct(ensembl, .keep_all = TRUE)

#Write results as .csv
write.csv(significant_stress, 'Significant_Stress.csv')
