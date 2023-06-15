# Timberlake Laboratory Incubations - Metagenomic Analysis
# By Cliff Bueno de Mesquita, JGI, Spring 2023
# Using code by Wyatt Hartman, JGI, 2018
# Wyatt's original code here: "~/Documents/GitHub/Timberlake/Timberlake_KO_Deseq2_tests.ipynb"

# DESeq2 testing of Timberlake KOs via phyloseq, Dev
# 2nd major draft, here normalize all genes, deprioritize FC and volcano plots.
# Incorporates more of DESeq2 vignette along with Phyloseq tutorial materials. 
# Here dropping source soil data entirely.



#### 1. Setup ####
# Libraries
library(phyloseq)
library(plyr)
library(tidyverse)
library(RColorBrewer)
library(stringr)
library(scales)
library(cowplot)
library(plotly)
suppressMessages(library(vegan))
suppressMessages(library(DESeq2))
suppressMessages(library(gplots))
library(pheatmap)
library(reshape2)
library(readxl)
setwd("~/Documents/GitHub/Timberlake/")

# Functions
find_hull <- function(df) df[chull(df$Axis01, df$Axis02),]
find_hullj <- function(df) df[chull(df$Axis01j, df$Axis02j),]

# Import and pre-process files
## KO table import, filter out field samples                           
TL_KO_raw <- read.csv("data/Timberlake_MG_KO_data.csv", header=TRUE) %>%
  select(1:19, Fxn)

# Check pmoABC trends (this was an issue in SF Salinity Gradient)
checkPMO <- TL_KO_raw %>%
  select(-Fxn) %>%
  column_to_rownames(var = "KO") %>%
  t() %>%
  as.data.frame() %>%
  select(K10944, K10945, K10946)
ggplot(checkPMO, aes(K10944, K10945)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dotted") +
  geom_point() +
  geom_smooth(method = "lm", alpha = 0.2) +
  labs(x = "pmoA", y = "pmoB") +
  theme_classic()
ggplot(checkPMO, aes(K10944, K10946)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dotted") +
  geom_point() +
  geom_smooth(method = "lm", alpha = 0.2) +
  labs(x = "pmoA", y = "pmoC") +
  theme_classic()
ggplot(checkPMO, aes(K10945, K10946)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dotted") +
  geom_point() +
  geom_smooth(method = "lm", alpha = 0.2) +
  labs(x = "pmoB", y = "pmoC") +
  theme_classic()
hist(checkPMO$K10944)
hist(checkPMO$K10945)
hist(checkPMO$K10946)
checkPMO_long <- melt(checkPMO, measure.vars = c("K10944", "K10945", "K10946"))
ggplot(checkPMO_long, aes(x = value, colour = variable)) +
  geom_density() +
  labs(x = "Count", y = "Density", colour = "KO") +
  theme_classic()
  
# df of KO gene fxns
TL_KO_fxn <- TL_KO_raw %>%
  select(KO, Fxn) %>%
  column_to_rownames(var = "KO")

# KO fxns as Phyloseq object
suppressWarnings(KOTaxTable <- tax_table(TL_KO_fxn)) # Fxn as PHYLOSEQ tax table 
colnames(KOTaxTable) <- colnames(TL_KO_fxn) # Add back colnames
row.names(KOTaxTable) <- row.names(TL_KO_fxn) # Add back rownames

# Gene Ontology and fxn Colors import
BGC_ont <- read.csv("data/Ontology_KO_CNPSch4_Fm_whh.csv", header=TRUE) # Ontology
BGC_colors <- read.csv("data/Ontol_KO_L2_Color_KEY_whh.csv") # Colors
BGC_ont_col <- merge(BGC_ont, BGC_colors, by = "L2", all.x =TRUE) %>% # Merge
  arrange(Index.x) %>% # Sort for readability
  filter(L2 != "Fermentation") # Drop fermentation

# Only CNPSch4 genes from Ontology
bgc_KOu <- unique(data.frame(KO = BGC_ont$KO)) # only KOs from BGC_ont
TL_bgc_KO <- merge(bgc_KOu, TL_KO_raw) # Merge with TL KO

# Sample Mapping import
map_MG_only <- read.table("data/Timberlake_sample_map_both.txt", sep="\t", header = T) %>%
  drop_na(MG_name) %>% # drop NA in MG samps 
  arrange(itag_meta_order2) %>% # sort by index (itag intentional)
  select(MG_name, Treat, Depth) %>% # keep only relevant cols
  column_to_rownames(var = "MG_name") %>% # Sample as row names, # Drop sample col
  filter(Treat != "SourceSoil") # Drop field samples, incubations only

## Phyloseq to DESeq2
# Prep KO count data for OTU table 
No_fxn <- data.frame(unique(TL_KO_raw)) %>% # DF for KOin, unique used to drop redund.
  column_to_rownames(var = "KO") %>% # Make rownames KO numbers, Drop KO column
  select(-Fxn) %>% # Drop Fxn
  as.matrix()

# Make PHYLOSEQ KO -> "OTU Table"
class(No_fxn) <- "numeric" # Make numeric for phyloseq
KO_otu <- otu_table(No_fxn, taxa_are_rows = TRUE, errorIfNULL = TRUE) # Make OTU table phyloseq object

# Make Phyloseq sample data
Map <- sample_data(map_MG_only) # map

### Make combined Phyloseq OBJECT
physeq = phyloseq(KO_otu, Map, KOTaxTable)

# phyloseq patch before DESeq
sample_data(physeq)$Treat <- relevel(as.factor(get_variable(physeq, "Treat")), 
                                     ref="Control")

# Export phyloseq to DESeq
KO_phy2des <- phyloseq_to_deseq2(physeq, ~ Treat)

## DESeq2 for Diff abund
# Wald test
KO_phy2des <- DESeq(KO_phy2des, 
                    test = "Wald", 
                    fitType = "parametric") # Test

# Inspect DESeq fmt data set, includes estimated dispersions
KO_phy2des

# Overall results  -- why nothing sig when contrasts are?
res <- results(KO_phy2des)
resultsNames(KO_phy2des)
summary(res)
plotMA(res)

# Likelihood ratio test, "ANOVA - like" for multiple factors 
# Note here some are sig again, unlike above Wald test
lrt <- DESeq(KO_phy2des, 
             test = "LRT", 
             reduced = ~1)
res_LRT <- results(lrt)
res_LRT
plotMA(res_LRT)



#### 2. DESeq2 Linear contrasts ####
### Function for linear contrasts with p > filter, FC cutoff
# Function for linear contrasts, results cutoff and merge with Ontology
DSq_cntrstF = function(Treat_col, Treat1, Treat2) {                      # Treat_col: column with treatment data
    alpha = 0.05                                                         # Treat1: for increase compared with Ref
    FC = 1                                                               # Treat2: Ref for comparison  

    T1_T2_De <- results(KO_phy2des,
                        contrast = c(Treat_col, Treat1, Treat2),
                        independentFiltering = FALSE) # DESeq2 results
    T1_T2_DeS <- T1_T2_De[which(T1_T2_De$padj < alpha), ]                # filter by alpha 
    T1_T2_DeSQ <- data.frame(T1_T2_DeS)                                   # make df 
    T1_T2_DeSQ$KO <- row.names(T1_T2_DeSQ)                                # add KO as rownames

    # Merge with CNP ontology, filter FC > 1 
    T1_T2_DeSQ <- merge(T1_T2_DeSQ, BGC_ont_col, by ='KO')               # merge w Ontology
    T1_T2_DeSQ_FC1 <- T1_T2_DeSQ[abs(T1_T2_DeSQ$log2FoldChange) > FC,]   # filter FC > 1 
    return(T1_T2_DeSQ_FC1)
}

# Test contrasts fxn: 
ASW_Ctrl_DeSQ_FC <- DSq_cntrstF("Treat", "ASW", "Control")

# Create ALL Contrasts results sets
ASW_Ctrl_FC <- DSq_cntrstF("Treat", "ASW", "Control")
ASW0S_Ctrl_FC <- DSq_cntrstF("Treat", "ASW_noS", "Control")
SO4_Ctrl_FC <- DSq_cntrstF("Treat", "SO4", "Control")
ASW_SO4_FC <- DSq_cntrstF("Treat", "ASW", "SO4")
ASW_ASW0S_FC <- DSq_cntrstF("Treat","ASW", "ASW_noS")
ASW0S_SO4_FC <- DSq_cntrstF("Treat","ASW_noS","SO4")

### Collect significnant KOs, all treats
# gather KO lists
ASW_Ctrl_KO <- data.frame(KO = ASW_Ctrl_FC$KO)
ASW0S_Ctrl_KO <- data.frame(KO = ASW0S_Ctrl_FC$KO)
SO4_Ctrl_KO <- data.frame(KO = SO4_Ctrl_FC$KO)
ASW_SO4_KO <- data.frame(KO = ASW_SO4_FC$KO)
ASW_ASW0S_KO <- data.frame(KO = ASW_ASW0S_FC$KO)
ASW0S_SO4_KO <- data.frame(KO = ASW0S_SO4_FC$KO)

# Combine lists, make unique
Resp_KO <- (rbind(ASW_Ctrl_KO, ASW0S_Ctrl_KO, SO4_Ctrl_KO, ASW_SO4_KO, ASW_ASW0S_KO, ASW0S_SO4_KO))
Resp_KOu <- unique(Resp_KO)
row.names(Resp_KOu) <- Resp_KOu$KO

## Get DESeq2 variance stablized data and Counts
# From updated phyloseq to DESeq2 tutorial at https://github/Joey711/phyloseq/issues/283
# Adding size factors and dispersion estimates before variance stabilization, as instructed
KO_phy2des_vs <- estimateSizeFactors(KO_phy2des)
KO_phy2des_vs <- estimateDispersions(KO_phy2des)
KO_phy2des_VST <- getVarianceStabilizedData(KO_phy2des)
KO_table_VST <- data.frame(KO_phy2des_VST)    # shorten name
KO_table_VST$KO <- row.names(KO_table_VST)    # index rows 

# Check VST df 
head(KO_table_VST)

## Replace VST data with COUNTS
KO_vs_counts <- counts(KO_phy2des_vs, normalized = TRUE)     # get counts
KO_table_VST <- data.frame(KO_vs_counts)                   # make df, shorten name
head(KO_table_VST)
colSums(KO_table_VST[,1:ncol(KO_table_VST)-1])

# Scale OTU VST table counts to CPM (Counts per million)
sampTots <- colSums(KO_table_VST)                          # Get sample totals   #[,1:ncol(OTU_table_VST)-1]) 
KO_VST_CPM <- sweep(KO_table_VST*1000000, 2, sampTots, '/')   # Sweep matrix by:[div. by samp total * 1e+06]    
KO_VST_CPM$KO <- row.names(KO_VST_CPM)                   # Add KO as column for later merges



#### 3. Heatmap ####
# Prep heatmap - Count data w only responsive KOs

# get KO Counts from DE KOs by MERGE
KO_table_VST_respM <- merge(Resp_KOu, KO_VST_CPM, by = "KO") %>%   # merge
  column_to_rownames(var = "KO") # Drop KO and add as rowname

# Correlations of responsive KOs and CH4 flux
nc_mg <- readRDS("data/input_filt_rare_mTAGs.rds")
CH4 <- nc_mg$map_loaded %>%
  select(Treatment, CH4_ug_m2_h) %>%
  filter(Treatment != "Field") %>%
  select(-Treatment) %>%
  rownames_to_column(var = "sampleID")
map_MG_only <- map_MG_only %>%
  mutate(sampleID = rownames(.)) %>%
  left_join(., CH4, by = "sampleID")
KO_table_VST_respM_t <- as.data.frame(t(KO_table_VST_respM))
sum(map_MG_only$sampleID != rownames(KO_table_VST_respM_t))
cor_df <- as.data.frame(matrix(NA, ncol(KO_table_VST_respM_t), 7)) %>%
  set_names(c("KO", "r", "PearsonP", "rho", "SpearmanP", "tau", "KendallP"))
for (i in 1:ncol(KO_table_VST_respM_t)) {
  m <- cor.test(map_MG_only$CH4_ug_m2_h, KO_table_VST_respM_t[,i], method = "pearson", na.rm = T)
  m1 <- cor.test(map_MG_only$CH4_ug_m2_h, KO_table_VST_respM_t[,i], method = "spearman", na.rm = T)
  m2 <- cor.test(map_MG_only$CH4_ug_m2_h, KO_table_VST_respM_t[,i], method = "kendall", na.rm = T)
  cor_df$KO[i] <- names(KO_table_VST_respM_t)[i]
  cor_df$r[i] <- m$estimate
  cor_df$PearsonP[i] <- m$p.value
  cor_df$rho[i] <- m1$estimate
  cor_df$SpearmanP[i] <- m1$p.value
  cor_df$tau[i] <- m2$estimate
  cor_df$KendallP[i] <- m2$p.value
}
cor_df <- cor_df %>%
  mutate(Pearsonfdr = p.adjust(PearsonP, method = "fdr")) %>%
  mutate(Spearmanfdr = p.adjust(SpearmanP, method = "fdr")) %>%
  mutate(Kendallfdr = p.adjust(KendallP, method = "fdr")) %>%
  mutate(PearsonPcut = factor(ifelse(Pearsonfdr < 0.05, 
                                     "Pfdr < 0.05", "Pfdr > 0.05"),
                              levels = c("Pfdr < 0.05","Pfdr > 0.05"))) %>%
  mutate(SpearmanPcut = factor(ifelse(Spearmanfdr < 0.05, 
                                     "Pfdr < 0.05", "Pfdr > 0.05"),
                              levels = c("Pfdr < 0.05","Pfdr > 0.05"))) %>%
  mutate(KendallPcut = factor(ifelse(Kendallfdr < 0.05, 
                                     "Pfdr < 0.05", "Pfdr > 0.05"),
                              levels = c("Pfdr < 0.05","Pfdr > 0.05"))) %>%
  mutate(CH4_dir = ifelse(rho < 0, "Negative", "Positive")) %>%
  mutate(CH4_sig = ifelse(SpearmanPcut == "Pfdr < 0.05", "Significant", "NS")) %>%
  mutate(CH4_cor = paste(CH4_dir, CH4_sig, sep = "")) %>%
  mutate(CH4_cor = as.factor(CH4_cor)) %>%
  mutate(CH4_cor = recode(CH4_cor, "NegativeSignificant" = "Negative", "NegativeNS" = "None",
                           "PositiveSignificant" = "Positive", "PositiveNS" = "None"))
hist(map_MG_only$CH4_ug_m2_h) # Normal
hist(KO_table_VST_respM_t$K00192) # Not Normal
# Use Spearman and KOs aren't normal. Pearson yielded only 1 significant correlation anyway

# Get z-score data, here of log10
KO_table_VST_respM[KO_table_VST_respM == 0] <- 0.5           # Replace 0 values with psuedo counts for LOG transf
KO_table_VST_respZ <- data.frame(t(scale(t(log10(KO_table_VST_respM)), 
                                         center = TRUE, 
                                         scale = TRUE))) # z-scores, log10 
KO_table_VST_respZ$KO <-row.names(KO_table_VST_respZ) # KO as rowname

## Add Ontology hier / colors
# Merge colors with Ontology
BGC_ont_colors <- merge(BGC_colors[,-1], BGC_ont, by = "L2", all.x = TRUE) # merge
BGC_ont_colors <- BGC_ont_colors[order(BGC_ont_colors$KO),]              # Sort by KO

# Merge Ontology with Responsive KO z-scores
KO_table_VSTrespZc <- merge(BGC_ont_colors, KO_table_VST_respZ, by = 'KO') %>%       # merge
  arrange(Index) # sort by index

# Merge Ontology with Spearman correlations
cor_table <- merge(BGC_ont_colors, cor_df, by = 'KO') %>%       # merge
  arrange(Index) # sort by index

# Get only data, add names to rows
KO_Resp_HeatDS0 <- KO_table_VSTrespZc                   # rename
KO_Resp_HeatDS <- KO_Resp_HeatDS0[,13:30]               # data only
KO_Resp_HeatDSt <- t(KO_Resp_HeatDS)                    # transpose (non-unique rows)
colnames(KO_Resp_HeatDSt) <- KO_Resp_HeatDS0$sm_name    # add sm names for FXN 
KO_Resp_HeatDStt <- t(KO_Resp_HeatDSt)                  # retranspose, keeps non-unique rownames

# Get colors for heatmaps
L2_colors <- as.character(KO_Resp_HeatDS0$color)

## pheatmap
# Put ASW on the right
KO_Resp_HeatDStt <- as.data.frame(KO_Resp_HeatDStt) %>%
  select(-ASW_1, ASW_1) %>%
  select(-ASW_2, ASW_2) %>%
  select(-ASW_3, ASW_3) %>%
  select(-ASW_4, ASW_4) %>%
  select(-ASW_5, ASW_5) %>%
  as.matrix()
L2_cols <- KO_Resp_HeatDS0 %>%
  select(L2, color) %>%
  group_by(L2) %>%
  dplyr::slice(1) %>%
  mutate(L2 = factor(L2,
                     levels = L2))
print(L2_cols, n = 22)
ann_cols <- data.frame(row.names = colnames(KO_Resp_HeatDStt),
                       "Treatment" = c(rep("Control", 5),
                                       rep("SO4", 5),
                                       rep("ASW-SO4", 3), 
                                       rep("ASW", 5)))
ann_rows <- data.frame(row.names = rownames(KO_Resp_HeatDStt),
                       "CH4_cor" = cor_table$CH4_cor,
                       "L2" = KO_Resp_HeatDS0$L2,
                       "L1" = KO_Resp_HeatDS0$L1)
ann_colors <- list(Treatment = c(Control = "#3B528BFF", 
                                 SO4 = "#21908CFF",
                                 `ASW-SO4` = "#5DC863FF",
                                 ASW = "#FDE725FF"),
                   CH4_cor = c(Positive = "red", Negative = "blue", None = "grey95"),
                   L1 = c(`Carbon ` = "brown", 
                          Nitrogen = "forestgreen", 
                          Phosphorus = "deepskyblue3", 
                          Sulfur = "magenta", 
                          CH4_cycling = "grey40", 
                          Fermentation = "darkorange"),
                   L2 = c(Sugars = L2_cols$color[22],
                          Polymers = L2_cols$color[18],
                          Aromatic = L2_cols$color[1],
                          NO3_reduction = L2_cols$color[14],
                          NH3_oxidation = L2_cols$color[11],
                          `NO3_A.reduction` = L2_cols$color[13],
                          NH4_assimilation = L2_cols$color[12],
                          N2_fixation = L2_cols$color[10],
                          P_regulation = L2_cols$color[15],
                          `PolyP-ases ` = L2_cols$color[17],
                          Phn_transport = L2_cols$color[16],
                          `CH3-phosphonate` = L2_cols$color[2],
                          S2O3_oxidation = L2_cols$color[19],
                          SO4_A.reduction = L2_cols$color[20],
                          SO4_D.reduction = L2_cols$color[21],
                          CH4_oxidation = L2_cols$color[7],
                          CH4_H2_reduction = L2_cols$color[4],
                          CH4_acetate = L2_cols$color[5],
                          CH4_methylotroph = L2_cols$color[6],
                          CH4_Archaeal = L2_cols$color[3],
                          H2_production = L2_cols$color[9],
                          Fermentation = L2_cols$color[8]))
colnames(KO_Resp_HeatDStt)[6:10] <- c("SO4_2", "SO4_3", "SO4_4", "SO4_5A", "SO4_5B")
pheatmap(KO_Resp_HeatDStt,
         legend = T,
         legend_breaks = c(-2, -1, 0, 1, 2, max(KO_Resp_HeatDStt)),
         legend_labels = c("-2","-1","0","1","2","Abund.\n"),
         main = "",
         border_color = NA,
         scale = "none",
         angle_col = 315,
         fontsize = 8,
         annotation_col = ann_cols,
         annotation_row = ann_rows,
         annotation_colors = ann_colors,
         cluster_rows = F,
         cluster_cols = F,
         gaps_row = c(12, 19, 25, 31, 50),
         gaps_col = c(5, 10, 13),
         #filename = "InitialFigs/KO_heatmap.png",
         filename = "FinalFigs/Figure7.png",
         width = 7,
         height = 7)
dev.off()
dev.set(dev.next())
dev.set(dev.next())



#### 4. NMDS ####
# Define plotting function
biom_plot_cats_nmds <- function(biom, group, env) {
  biomT <- t(biom)    # Transpose so cols are vars
  
  # Test ADONIS models, extract params
  bT_adonis <- adonis(biomT  ~ group, permutations = 999, method = "bray");
  R2 <- round(bT_adonis$aov.tab$R2[1], digits=3)
  P <- bT_adonis$aov.tab$"Pr(>F)"[1]
  
  # NMDS
  bT_mds <- suppressMessages(metaMDS(biomT[,-1], 
                                     distance = "bray", 
                                     k = 3, 
                                     trymax = 10)); # run NMDS
  bT_mds_DF = data.frame(MDS1 = bT_mds$points[,1], MDS2 = bT_mds$points[,2], group)  # Make data frame
  stress = round(bT_mds$stress, digits = 3)                                      # get NMDS stress
  
  # Plot NMDS by assembly group
  pA <- ggplot(bT_mds_DF, aes(x = MDS1, y = MDS2, color = group)) + 
    geom_point() +
    stat_ellipse(level=0.95);
  pB <-pA + annotate("text", 
                     x = (0.9*(max(bT_mds_DF$MDS1))), 
                     y = max(bT_mds_DF$MDS2), 
                     label = paste("italic(R) ^ 2 ==", R2),
                     parse = TRUE);
  pC <-pB + annotate("text", 
                     x = (0.7*(min(bT_mds_DF$MDS1))), 
                     y = min(bT_mds_DF$MDS2), 
                     label = paste("italic(Stress):", stress),
                     parse = TRUE);
  pD <- pC + theme(legend.title = element_blank())
  
  return(pD)
}

###  New Adaptation
## ALL GENES
# Get counts per sample, from KO_VST_CPM  not  KO_table_VST
KO_table_VST_nk <- KO_VST_CPM[, 1:(ncol(KO_VST_CPM)-1)]  # Drop KO column, previously needed for merge  # 
#KO_table_VST_nk <- KO_table_VST[, 1:(ncol(KO_table_VST)-1)]  # Drop KO column, previously needed for merge  # 

# Get treatment vector
Treat <- map_MG_only$Treat[1:ncol(KO_table_VST_nk)] # Treat

## PLOT
#suppressWarnings(biom_plot_cats_nmds(KO_table_VST_nk, Treat))
suppressMessages(biom_plot_cats_nmds(KO_table_VST_nk, Treat)) 



#### 5. PCoA ####
# All genes
# Bray-Curtis and Jaccard
ko_meta <- map_MG_only %>%
  mutate(Treat = factor(Treat,
                        levels = c("Control", "SO4", "ASW_noS", "ASW"))) %>%
  mutate(Treat = recode(Treat, "ASW_noS" = "ASW-SO4"))
biomT <- t(KO_table_VST_nk)
sum(rownames(ko_meta) != rownames(biomT))
bc_ko <- vegdist(biomT, method = "bray")
pcoa_ko <- cmdscale(bc_ko, k = nrow(ko_meta) - 1, eig = T)
pcoaA1 <- round((eigenvals(pcoa_ko)/sum(eigenvals(pcoa_ko)))[1]*100, digits = 1)
pcoaA2 <- round((eigenvals(pcoa_ko)/sum(eigenvals(pcoa_ko)))[2]*100, digits = 1)
ko_meta$Axis01 <- scores(pcoa_ko)[,1]
ko_meta$Axis02 <- scores(pcoa_ko)[,2]
micro.hulls <- ddply(ko_meta, c("Treat"), find_hull)
g1_ko <- ggplot(ko_meta, aes(Axis01, Axis02, colour = Treat)) +
  geom_polygon(data = micro.hulls, aes(colour = Treat, fill = Treat),
               alpha = 0.1, show.legend = F) +
  geom_point(size = 3, alpha = 1, shape = 17) +
  scale_fill_manual(values = viridis_pal()(5)[2:5]) +
  scale_colour_manual(values = viridis_pal()(5)[2:5]) +
  labs(x = paste("PC1: ", pcoaA1, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2, "%", sep = ""),
       colour = "Treatment",
       title = "KO Bray-Curtis") +
  theme_bw() +  
  theme(legend.position = "none",
        axis.title = element_text(face = "bold", size = 14), 
        axis.text = element_text(size = 12),
        plot.title = element_text(vjust = 0))
g1_ko

jac_ko <- vegdist(biomT, method = "jaccard")
pcoa1_ko <- cmdscale(jac_ko, k = nrow(ko_meta) - 1, eig = T)
pcoa1A1 <- round((eigenvals(pcoa1_ko)/sum(eigenvals(pcoa1_ko)))[1]*100, digits = 1)
pcoa1A2 <- round((eigenvals(pcoa1_ko)/sum(eigenvals(pcoa1_ko)))[2]*100, digits = 1)
ko_meta$Axis01j <- scores(pcoa1_ko)[,1]
ko_meta$Axis02j <- scores(pcoa1_ko)[,2]
micro.hullsj <- ddply(ko_meta, c("Treat"), find_hullj)
g2_ko <- ggplot(ko_meta, aes(Axis01j, Axis02j, colour = Treat)) +
  geom_polygon(data = micro.hullsj, aes(colour = Treat, fill = Treat),
               alpha = 0.1, show.legend = F) +
  geom_point(size = 3, alpha = 1, shape = 17) +
  scale_fill_manual(values = viridis_pal()(5)[2:5]) +
  scale_colour_manual(values = viridis_pal()(5)[2:5]) +
  labs(x = paste("PC1: ", pcoa1A1, "%", sep = ""), 
       y = paste("PC2: ", pcoa1A2, "%", sep = ""),
       colour = "Treatment",
       title = "KO Jaccard") +
  theme_bw() +  
  theme(legend.position = "right",
        axis.title = element_text(face = "bold", size = 14), 
        axis.text = element_text(size = 12),
        plot.title = element_text(vjust = 0))
g2_ko

pdf("InitialFigs/KO_PCoA.pdf", width = 7, height = 5)
plot_grid(g1_ko, g2_ko, ncol = 2, rel_widths = c(1,1.515))
dev.off()
ggplotly(g1_ko)



#### 6. FC Plot ####
# Plot fold change by organism, ontology colors
fold_plot = function(sigtab){

    # L2 order
    x = tapply(sigtab$log2FoldChange, sigtab$L2, function(x) max(x))
    x = sort(x, TRUE)
    sigtab$L2 <- factor(as.character(sigtab$L2), levels = names(x))
    # sm_name order
    x = tapply(sigtab$log2FoldChange, sigtab$sm_name, function(x) max(x))
    x = sort(x, TRUE)
    sigtab$sm_name <- factor(as.character(sigtab$sm_name), levels = names(x))
    
    sigtab <- sigtab[order(sigtab$Index.y),]                                     # Sort Table by factor ordering index
    sigtab$L2 <- factor(sigtab$L2, levels = unique(sigtab$L2[order(sigtab$Index.y)]))  # Reorder factor levels

    # Get factor colors as vector
    colmat <-unique(data.frame(L2=sigtab$L2, color=sigtab$color))
    w = c(unlist(as.character(colmat$color)))
    
        
    #pt_size = (sigtab$baseMean/100)
    p <- ggplot(sigtab, aes(x = sm_name, y = log2FoldChange, color = L2)) + 
      geom_point(aes(size = baseMean)) +
      labs(x = "KO", y = "log2 fold change", 
           size = "Mean count", color = "Gene class") +
      scale_color_manual(values = w) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
    p
}

fold_plot(ASW_Ctrl_DeSQ_FC)



#### 7. Volcano FC Plot ####
# Volcano FC plot, colored by ontology
volc_plot = function(sigtab) {
    sigtab <- sigtab[order(sigtab$Index.y),]                                  # Sort Table by factor ordering index
    sigtab$L2 <- factor(sigtab$L2, levels = unique(sigtab$L2[order(sigtab$Index.y)]))  # Reorder factor levels

    # Get factor colors as vector
    colmat <-unique(data.frame(L2=sigtab$L2, color=sigtab$color))
    w = c(unlist(as.character(colmat$color)))
    
    # PLOT
    ggplot(sigtab, aes(x = log2FoldChange, y = baseMean, group = L2, color = L2)) + 
      geom_point(size = 2) + 
      scale_color_manual(values = w)
}

# TEST volc_plot FUNCTION
volc_plot(ASW_Ctrl_DeSQ_FC)

# vignette("DESeq2")


#### 8. COGs ####
# Import 
IMG_meta <- read.delim("data/IMGmetadata.txt") %>%
  separate(`Genome.Name...Sample.Name`, into = c("Name", "sampleID"),
           sep = " - ", remove = F) %>%
  mutate(sampleID = gsub("SO4_5", "SO4_5B", sampleID)) %>%
  mutate(sampleID = gsub("SO4_4", "SO4_5A", sampleID)) %>%
  mutate(sampleID = gsub("SO4_3", "SO4_4", sampleID)) %>%
  mutate(sampleID = gsub("SO4_2", "SO4_3", sampleID)) %>%
  mutate(sampleID = gsub("SO4_1", "SO4_2", sampleID)) %>%
  separate(sampleID, into = c("Treatment", "Replicate"), sep = "_", remove = F) %>%
  filter(Treatment != "SourceSoil") %>%
  mutate(Treatment = factor(Treatment,
                            levels = c("Control", "SO4", "ASW-S", "ASW"))) %>%
  mutate(Treatment = recode_factor(Treatment, "ASW-S" = "ASW-SO4")) %>%
  mutate(Treatment = factor(Treatment,
                            levels = c("Control", "SO4", "ASW-SO4", "ASW"))) %>%
  arrange(Treatment, Replicate)

cog <- read.delim("data/TL_COG_cat/UI_data_output.txt") %>%
  column_to_rownames(var = "FeatureName") %>%
  select(9:31) %>%
  set_names(substring(names(.), first = 1, last = 11)) %>%
  set_names(gsub("X", "", names(.))) %>%
  select(as.character(IMG_meta$taxon_oid)) %>%
  t() %>%
  as.data.frame()

sum(IMG_meta$taxon_oid != rownames(cog)) # 0, good

# Normalize
dds_input_1 <- DESeqDataSetFromMatrix(countData = t(cog),
                                      colData = IMG_meta,
                                      design = ~ 1)
dds_input_SF <- estimateSizeFactors(dds_input_1)
dds_input_D <- estimateDispersions(dds_input_SF)
cog_DESeq <- as.data.frame((counts(dds_input_D, normalized = T))) %>%
  t() %>%
  as.data.frame()
sum(rownames(cog_DESeq) != IMG_meta$taxon_oid)

# Test
dds_input_da <- DESeqDataSetFromMatrix(countData = t(cog),
                                       colData = IMG_meta,
                                       design = ~ Treatment)
wald <- DESeq(object = dds_input_da, test = "Wald", fitType = "parametric")
res_w <- results(wald,
                 contrast = c("Treatment", "ASW", "Control"),
                 independentFiltering = FALSE,
                 pAdjustMethod = "fdr")
summary(res_w)
plotMA(res_w)

lrt <- DESeq(object = dds_input_da, test = "LRT", reduced = ~1)
res_l <- results(lrt,
                 contrast = c("Treatment", "ASW", "Control"),
                 independentFiltering = FALSE,
                 pAdjustMethod = "fdr")
summary(res_l)
plotMA(res_l)

res_p <- data.frame(LRT.p = res_l$padj,
                    Wald.p = res_w$padj) %>%
  mutate(LRT = ifelse(LRT.p < 0.05, "Pfdr < 0.05", "Pfdr > 0.05")) %>%
  mutate(Wald = ifelse(Wald.p < 0.05, "Pfdr < 0.05", "Pfdr > 0.05"))

# Plot
cog_DESeq <- cog_DESeq %>%
  t() %>%
  as.data.frame() %>%
  set_names(IMG_meta$sampleID)
ann_cols <- data.frame(row.names = colnames(cog_DESeq),
                       "Treatment" = IMG_meta$Treatment)
ann_rows <- data.frame(row.names = rownames(cog_DESeq),
                       "Wald" = res_p$Wald,
                       "LRT" = res_p$LRT)
ann_colors <- list(Treatment = c(Control = "#3B528BFF", 
                                 SO4 = "#21908CFF",
                                 `ASW-SO4` = "#5DC863FF",
                                 ASW = "#FDE725FF"),
                   Wald = c(`Pfdr < 0.05` = "black",
                            `Pfdr > 0.05` = "white"),
                   LRT = c(`Pfdr < 0.05` = "black",
                           `Pfdr > 0.05` = "white"))
pheatmap(cog_DESeq,
         legend = T,
         legend_breaks = c(-2, -1, 0, 1, 2, 2.5),
         legend_labels = c("-2","-1","0","1","2","Abund.\n"),
         main = "",
         border_color = NA,
         scale = "row",
         angle_col = 315,
         fontsize = 8,
         fontsize_row = 6,
         fontsize_col = 6,
         annotation_col = ann_cols,
         annotation_row = ann_rows,
         annotation_colors = ann_colors,
         cluster_rows = T,
         cluster_cols = F,
         gaps_col = c(5, 10, 13),
         filename = "InitialFigs/COG_heatmap.png",
         width = 7,
         height = 7)
dev.off()
dev.set(dev.next())
dev.set(dev.next())



#### 9. Pfams ####
pfam <- read.delim("data/TL_PFam_cat/UI_data_output.txt") %>%
  column_to_rownames(var = "FeatureName") %>%
  select(9:31) %>%
  set_names(substring(names(.), first = 1, last = 11)) %>%
  set_names(gsub("X", "", names(.))) %>%
  select(as.character(IMG_meta$taxon_oid)) %>%
  t() %>%
  as.data.frame()

sum(IMG_meta$taxon_oid != rownames(pfam)) # 0, good

# Normalize
dds_input_1 <- DESeqDataSetFromMatrix(countData = t(pfam),
                                      colData = IMG_meta,
                                      design = ~ 1)
dds_input_SF <- estimateSizeFactors(dds_input_1)
dds_input_D <- estimateDispersions(dds_input_SF)
pfam_DESeq <- as.data.frame((counts(dds_input_D, normalized = T))) %>%
  t() %>%
  as.data.frame()
sum(rownames(pfam_DESeq) != IMG_meta$taxon_oid)

# Test
dds_input_da <- DESeqDataSetFromMatrix(countData = t(pfam),
                                       colData = IMG_meta,
                                       design = ~ Treatment)
wald <- DESeq(object = dds_input_da, test = "Wald", fitType = "parametric")
res_w <- results(wald,
                 contrast = c("Treatment", "ASW", "Control"),
                 independentFiltering = FALSE,
                 pAdjustMethod = "fdr")
summary(res_w)
plotMA(res_w)

lrt <- DESeq(object = dds_input_da, test = "LRT", reduced = ~1)
res_l <- results(lrt,
                 contrast = c("Treatment", "ASW", "Control"),
                 independentFiltering = FALSE,
                 pAdjustMethod = "fdr")
summary(res_l)
plotMA(res_l)

res_p <- data.frame(LRT.p = res_l$padj,
                    Wald.p = res_w$padj) %>%
  mutate(LRT = ifelse(LRT.p < 0.05, "Pfdr < 0.05", "Pfdr > 0.05")) %>%
  mutate(Wald = ifelse(Wald.p < 0.05, "Pfdr < 0.05", "Pfdr > 0.05"))

# Plot
pfam_DESeq <- pfam_DESeq %>%
  t() %>%
  as.data.frame() %>%
  set_names(IMG_meta$sampleID)
ann_cols <- data.frame(row.names = colnames(pfam_DESeq),
                       "Treatment" = IMG_meta$Treatment)
ann_rows <- data.frame(row.names = rownames(pfam_DESeq),
                       "Wald" = res_p$Wald,
                       "LRT" = res_p$LRT)
ann_colors <- list(Treatment = c(Control = "#3B528BFF", 
                                 SO4 = "#21908CFF",
                                 `ASW-SO4` = "#5DC863FF",
                                 ASW = "#FDE725FF"),
                   Wald = c(`Pfdr < 0.05` = "black",
                            `Pfdr > 0.05` = "white"),
                   LRT = c(`Pfdr < 0.05` = "black",
                           `Pfdr > 0.05` = "white"))
pheatmap(pfam_DESeq,
         legend = T,
         legend_breaks = c(-3, -2, -1, 0, 1, 2, 3, 3.3),
         legend_labels = c("-3","-2","-1","0","1","2","3","Abund.\n"),
         main = "",
         border_color = NA,
         scale = "row",
         angle_col = 315,
         fontsize = 8,
         fontsize_row = 6,
         fontsize_col = 6,
         annotation_col = ann_cols,
         annotation_row = ann_rows,
         annotation_colors = ann_colors,
         cluster_rows = T,
         cluster_cols = F,
         gaps_col = c(5, 10, 13),
         filename = "InitialFigs/Pfam_heatmap.png",
         width = 7,
         height = 7)
dev.off()
dev.set(dev.next())
dev.set(dev.next())



#### 10. KEGG Module ####
mod <- read.delim("data/TL_KEGG_mod/UI_data_output.txt") %>%
  mutate(FeatureName = paste(Feature, FeatureName, sep = " ")) %>%
  column_to_rownames(var = "FeatureName") %>%
  select(9:31) %>%
  set_names(substring(names(.), first = 1, last = 11)) %>%
  set_names(gsub("X", "", names(.))) %>%
  select(as.character(IMG_meta$taxon_oid)) %>%
  t() %>%
  as.data.frame()
sum(IMG_meta$taxon_oid != rownames(mod)) # 0, good

# Normalize
dds_input_1 <- DESeqDataSetFromMatrix(countData = t(mod),
                                      colData = IMG_meta,
                                      design = ~ 1)
dds_input_SF <- estimateSizeFactors(dds_input_1)
dds_input_D <- estimateDispersions(dds_input_SF)
mod_DESeq <- as.data.frame((counts(dds_input_D, normalized = T))) %>%
  t() %>%
  as.data.frame()
sum(rownames(mod_DESeq) != IMG_meta$taxon_oid)

# Test
dds_input_da <- DESeqDataSetFromMatrix(countData = t(mod),
                                       colData = IMG_meta,
                                       design = ~ Treatment)
wald <- DESeq(object = dds_input_da, test = "Wald", fitType = "parametric")
res_w <- results(wald,
                 contrast = c("Treatment", "ASW", "Control"),
                 independentFiltering = FALSE,
                 pAdjustMethod = "fdr")
summary(res_w)
plotMA(res_w)

lrt <- DESeq(object = dds_input_da, test = "LRT", reduced = ~1)
res_l <- results(lrt,
                 contrast = c("Treatment", "ASW", "Control"),
                 independentFiltering = FALSE,
                 pAdjustMethod = "fdr")
summary(res_l)
plotMA(res_l)

res_p <- data.frame(mod = colnames(mod),
                    LRT.p = res_l$padj,
                    Wald.p = res_w$padj) %>%
  mutate(LRT = ifelse(LRT.p < 0.05, "Pfdr < 0.05", "Pfdr > 0.05")) %>%
  mutate(Wald = ifelse(Wald.p < 0.05, "Pfdr < 0.05", "Pfdr > 0.05")) %>%
  mutate(Sig = ifelse(LRT.p < 0.05 & Wald.p < 0.05, "Sig", "NS")) %>%
  filter(Sig == "Sig") # Subset to significant (because so many)




# Plot
mod_DESeq <- mod_DESeq %>%
  t() %>%
  as.data.frame() %>%
  filter(rownames(.) %in% res_p$mod) %>%
  set_names(IMG_meta$sampleID)
ann_cols <- data.frame(row.names = colnames(mod_DESeq),
                       "Treatment" = IMG_meta$Treatment)
ann_rows <- data.frame(row.names = rownames(mod_DESeq),
                       "Wald" = res_p$Wald,
                       "LRT" = res_p$LRT)
ann_colors <- list(Treatment = c(Control = "#3B528BFF", 
                                 SO4 = "#21908CFF",
                                 `ASW-SO4` = "#5DC863FF",
                                 ASW = "#FDE725FF"),
                   Wald = c(`Pfdr < 0.05` = "black",
                            `Pfdr > 0.05` = "white"),
                   LRT = c(`Pfdr < 0.05` = "black",
                           `Pfdr > 0.05` = "white"))
pheatmap(mod_DESeq,
         legend = T,
         legend_breaks = c(-3, -2, -1, 0, 1, 2, 3, 3.3),
         legend_labels = c("-3","-2","-1","0","1","2","3","Abund.\n"),
         main = "",
         border_color = NA,
         scale = "row",
         angle_col = 315,
         fontsize = 8,
         fontsize_row = 6,
         fontsize_col = 6,
         annotation_col = ann_cols,
         #annotation_row = ann_rows,
         annotation_colors = ann_colors,
         cluster_rows = T,
         cluster_cols = F,
         gaps_col = c(5, 10, 13),
         filename = "InitialFigs/KEGGmod_heatmap.png",
         width = 8,
         height = 7)
dev.off()
dev.set(dev.next())
dev.set(dev.next())



#### 11. Methanogenesis ####
# Make separate heatmaps for each pathway
# Some overlapping genes, that's okay. Can repeat them if present in multiple pathways
# Transformed KO table DESeq + CPM + Log10 + z-score!
KO_VST_CPM[KO_VST_CPM == 0] <- 0.5 # Replace 0 values with psuedo counts for LOG transform
KO_VST_CPM$KO <- NULL
KO_VST_CPM_Log_Z <- data.frame(t(scale(t(log10(KO_VST_CPM)), 
                                         center = TRUE, 
                                         scale = TRUE))) # z-scores, log10 
ko_tab <- KO_VST_CPM_Log_Z %>%
  select(-ASW_1, ASW_1) %>%
  select(-ASW_2, ASW_2) %>%
  select(-ASW_3, ASW_3) %>%
  select(-ASW_4, ASW_4) %>%
  select(-ASW_5, ASW_5) %>%
  rownames_to_column(var = "KO")
names(ko_tab)[7:11] <- c("SO4_2", "SO4_3", "SO4_4", "SO4_5A", "SO4_5B")

# Read in pathway KOs
mko_ace <- read_excel("data/MethanogenesisKOs.xlsx", sheet = 1) %>%
  mutate(Name = ifelse(is.na(Name), "", Name)) %>%
  mutate(KO_def = paste(KEGG_KO, Name, sep = " ")) %>%
  mutate_if(is.character, as.factor) %>%
  arrange(Pathway_Order, Reaction_Order) %>%
  select(-ModelSeed_ID)
mko_co2 <- read_excel("data/MethanogenesisKOs.xlsx", sheet = 2) %>%
  mutate(Name = ifelse(is.na(Name), "", Name)) %>%
  mutate(KO_def = paste(KEGG_KO, Name, sep = " ")) %>%
  mutate_if(is.character, as.factor) %>%
  arrange(Pathway_Order, Reaction_Order) %>%
  select(-ModelSeed_ID)
mko_metD <- read_excel("data/MethanogenesisKOs.xlsx", sheet = 3) %>%
  mutate(Name = ifelse(is.na(Name), "", Name)) %>%
  mutate(KO_def = paste(KEGG_KO, Name, sep = " ")) %>%
  mutate_if(is.character, as.factor) %>%
  arrange(Pathway_Order, Reaction_Order) %>%
  select(-ModelSeed_ID)
mko_metR <- read_excel("data/MethanogenesisKOs.xlsx", sheet = 4) %>%
  mutate(Name = ifelse(is.na(Name), "", Name)) %>%
  mutate(KO_def = paste(KEGG_KO, Name, sep = " ")) %>%
  mutate_if(is.character, as.factor) %>%
  arrange(Pathway_Order, Reaction_Order) %>%
  select(-ModelSeed_ID)

# Heatmap data and metadata for each
mko_ace_mat <- left_join(mko_ace, ko_tab, by = c("KEGG_KO" = "KO")) %>%
  select(7:25) %>%
  column_to_rownames(var = "KO_def") %>%
  as.matrix()
mko_ace_meta <- left_join(mko_ace, ko_tab, by = c("KEGG_KO" = "KO")) %>%
  select(1:7)

mko_co2_mat <- left_join(mko_co2, ko_tab, by = c("KEGG_KO" = "KO")) %>%
  select(7:25) %>%
  column_to_rownames(var = "KO_def") %>%
  as.matrix()
mko_co2_meta <- left_join(mko_co2, ko_tab, by = c("KEGG_KO" = "KO")) %>%
  select(1:7)

mko_metD_mat <- left_join(mko_metD, ko_tab, by = c("KEGG_KO" = "KO")) %>%
  select(7:25) %>%
  column_to_rownames(var = "KO_def") %>%
  as.matrix()
mko_metD_meta <- left_join(mko_metD, ko_tab, by = c("KEGG_KO" = "KO")) %>%
  select(1:7)

mko_metR_mat <- left_join(mko_metR, ko_tab, by = c("KEGG_KO" = "KO")) %>%
  select(11:29) %>%
  column_to_rownames(var = "KO_def") %>%
  as.matrix()
mko_metR_meta <- left_join(mko_metR, ko_tab, by = c("KEGG_KO" = "KO")) %>%
  select(1:11)

## Heatmaps
# Acetate
ann_cols <- data.frame(row.names = colnames(mko_ace_mat),
                       "Treatment" = c(rep("Control", 5),
                                       rep("SO4", 5),
                                       rep("ASW-SO4", 3), 
                                       rep("ASW", 5)))
ann_rows <- data.frame(row.names = rownames(mko_ace_mat), 
                       Pathway = mko_ace_meta$Pathway_Specific,
                       Process = mko_ace_meta$Pathway_General)
ann_colors <- list(Pathway = c(`Acetate cleavage` = "forestgreen",
                               `Wood-Ljungdahl - CO2 out` = "orange",
                               `Na out` = "purple",
                               "Methyl-CoM reduction" = "#21908CFF",
                               Methanophenazine = brewer.pal(9, "YlOrRd")[1],
                               Ferredoxin = brewer.pal(9, "YlOrRd")[2],
                               Hydrogenase = brewer.pal(9, "YlOrRd")[3]),
                   Process = c("Methyl-CoM formation" = "#440154FF",
                               "Methanogenesis" = "#21908CFF",
                               "CoB - CoM regeneration" = "#FDE725FF"),
                   Treatment = c(Control = "#3B528BFF", 
                                 SO4 = "#21908CFF",
                                 `ASW-SO4` = "#5DC863FF",
                                 ASW = "#FDE725FF"))
pheatmap(mko_ace_mat,
         legend = T,
         legend_breaks = c(-3, -2, -1, 0, 1, 2, 3, max(mko_ace_mat, na.rm = T)),
         legend_labels = c("-3","-2","-1","0","1","2","3","Abund.\n"),
         main = "Acetoclastic Methanogenesis",
         border_color = NA,
         scale = "none",
         angle_col = 315,
         fontsize = 8,
         annotation_col = ann_cols,
         annotation_row = ann_rows,
         annotation_colors = ann_colors,
         cluster_rows = F,
         cluster_cols = F,
         gaps_col = c(5, 10, 13),
         filename = "InitialFigs/KO_heatmap_acetate.png",
         width = 7,
         height = 7)
dev.off()
dev.set(dev.next())
dev.set(dev.next())

# H2/CO2
ann_cols <- data.frame(row.names = colnames(mko_co2_mat),
                       "Treatment" = c(rep("Control", 5),
                                       rep("SO4", 5),
                                       rep("ASW-SO4", 3), 
                                       rep("ASW", 5)))
ann_rows <- data.frame(row.names = rownames(mko_co2_mat), 
                       Pathway = mko_co2_meta$Pathway_Specific,
                       Process = mko_co2_meta$Pathway_General)
ann_colors <- list(Pathway = c(`Wood-Ljungdahl - CO2 in` = "orange",
                               `Na out` = "purple",
                               "Methyl-CoM reduction" = "#21908CFF",
                               F420 = brewer.pal(9, "YlOrRd")[1],
                               Ferredoxin = brewer.pal(9, "YlOrRd")[2],
                               "Ferredoxin-F420-H2" = brewer.pal(9, "YlOrRd")[3],
                               H2 = brewer.pal(9, "YlOrRd")[4],
                               `Formate-F420` = brewer.pal(9, "YlOrRd")[5]),
                   Process = c("Methyl-CoM formation" = "#440154FF",
                               "Methanogenesis" = "#21908CFF",
                               "CoB - CoM regeneration" = "#FDE725FF"),
                   Treatment = c(Control = "#3B528BFF", 
                                 SO4 = "#21908CFF",
                                 `ASW-SO4` = "#5DC863FF",
                                 ASW = "#FDE725FF"))
pheatmap(mko_co2_mat,
         legend = T,
         legend_breaks = c(-3, -2, -1, 0, 1, 2, 3, max(mko_co2_mat, na.rm = T)),
         legend_labels = c("-3","-2","-1","0","1","2","3","Abund.\n"),
         main = "Hydrogenotrophic Methanogenesis",
         border_color = NA,
         scale = "none",
         angle_col = 315,
         fontsize = 8,
         annotation_col = ann_cols,
         annotation_row = ann_rows,
         annotation_colors = ann_colors,
         cluster_rows = F,
         cluster_cols = F,
         gaps_col = c(5, 10, 13),
         filename = "InitialFigs/KO_heatmap_h2_co2.png",
         width = 7,
         height = 7)
dev.off()
dev.set(dev.next())
dev.set(dev.next())
  
# Methyl dismutation
ann_cols <- data.frame(row.names = colnames(mko_metD_mat),
                       "Treatment" = c(rep("Control", 5),
                                       rep("SO4", 5),
                                       rep("ASW-SO4", 3), 
                                       rep("ASW", 5)))
ann_rows <- data.frame(row.names = rownames(mko_metD_mat), 
                       Pathway = mko_metD_meta$Pathway_Specific,
                       Process = mko_metD_meta$Pathway_General)
ann_colors <- list(Pathway = c(Trimethylamine = brewer.pal(7, "Purples")[1], 
                               Dimethylamine = brewer.pal(7, "Purples")[2], 
                               Methylamine = brewer.pal(7, "Purples")[3],
                               Methylamines = brewer.pal(7, "Purples")[4],
                               Methanol = brewer.pal(7, "Purples")[5],
                               Betaine = brewer.pal(7, "Purples")[6],
                               `DMS/MeSH/MMPA` = brewer.pal(7, "Purples")[7],
                               `Wood-Ljungdahl - CO2 out` = "orange",
                               `Na in` = "purple",
                               "Methyl-CoM reduction" = "#21908CFF",
                               Methanophenazine = brewer.pal(9, "YlOrRd")[1],
                               Ferredoxin = brewer.pal(9, "YlOrRd")[2],
                               Hydrogenase = brewer.pal(9, "YlOrRd")[3]),
                   Process = c("Methyl-CoM formation" = "#440154FF",
                               "Methanogenesis" = "#21908CFF",
                               "CoB - CoM regeneration" = "#FDE725FF"),
                   Treatment = c(Control = "#3B528BFF", 
                                 SO4 = "#21908CFF",
                                 `ASW-SO4` = "#5DC863FF",
                                 ASW = "#FDE725FF"))
pheatmap(mko_metD_mat,
         legend = T,
         legend_breaks = c(-3, -2, -1, 0, 1, 2, 3, max(mko_metD_mat, na.rm = T)),
         legend_labels = c("-3","-2","-1","0","1","2","3","Abund.\n"),
         main = "Methyl-dismutating Methanogenesis",
         border_color = NA,
         scale = "none",
         angle_col = 315,
         fontsize = 8,
         annotation_col = ann_cols,
         annotation_row = ann_rows,
         annotation_colors = ann_colors,
         cluster_rows = F,
         cluster_cols = F,
         gaps_col = c(5, 10, 13),
         filename = "InitialFigs/KO_heatmap_methylDis.png",
         width = 7,
         height = 7)
dev.off()
dev.set(dev.next())
dev.set(dev.next())

# Methyl reduction
ann_cols <- data.frame(row.names = colnames(mko_metR_mat),
                       "Treatment" = c(rep("Control", 5),
                                       rep("SO4", 5),
                                       rep("ASW-SO4", 3), 
                                       rep("ASW", 5)))
ann_rows <- data.frame(row.names = rownames(mko_metR_mat), 
                       Pathway = mko_metR_meta$Pathway_Specific,
                       Process = mko_metR_meta$Pathway_General,
                       `M_blatticola` = mko_metR_meta$`M. blatticola`,
                       `M_stadtmanae` = mko_metR_meta$`M. stadtmanae`,
                       `M_luminyensis` = mko_metR_meta$`M. luminyensis`,
                       `M_thermophilum` = mko_metR_meta$`M. thermophilum`)
ann_colors <- list(Pathway = c(Trimethylamine = brewer.pal(7, "Purples")[1], 
                               Dimethylamine = brewer.pal(7, "Purples")[2], 
                               Methylamine = brewer.pal(7, "Purples")[3],
                               Methylamines = brewer.pal(7, "Purples")[4],
                               Methanol = brewer.pal(7, "Purples")[5],
                               Betaine = brewer.pal(7, "Purples")[6],
                               `DMS/MeSH/MMPA` = brewer.pal(7, "Purples")[7],
                               `Wood-Ljungdahl - CO2 out` = "orange",
                               `Na in` = "purple",
                               "Methyl-CoM reduction" = "#21908CFF",
                               Methanophenazine = brewer.pal(9, "YlOrRd")[1],
                               Ferredoxin = brewer.pal(9, "YlOrRd")[2],
                               Hydrogenase = brewer.pal(9, "YlOrRd")[3],
                               Formate = brewer.pal(9, "YlOrRd")[4]),
                   Process = c("Methyl-CoM formation" = "#440154FF",
                               "Methanogenesis" = "#21908CFF",
                               "CoB - CoM regeneration" = "#FDE725FF"),
                   `M_blatticola` = c("Present" = "black",
                                       "Absent" = "white"),
                   `M_stadtmanae` = c("Present" = "black",
                                       "Absent" = "white"),
                   `M_luminyensis` = c("Present" = "black",
                                        "Absent" = "white"),
                   `M_thermophilum` = c("Present" = "black",
                                         "Absent" = "white"),
                   Treatment = c(Control = "#3B528BFF", 
                                 SO4 = "#21908CFF",
                                 `ASW-SO4` = "#5DC863FF",
                                 ASW = "#FDE725FF"))
pheatmap(mko_metR_mat,
         legend = T,
         legend_breaks = c(-3, -2, -1, 0, 1, 2, 3, max(mko_metR_mat, na.rm = T)),
         legend_labels = c("-3","-2","-1","0","1","2","3","Abund.\n"),
         main = "Methyl-reducing Methanogenesis",
         border_color = NA,
         scale = "none",
         angle_col = 315,
         fontsize = 8,
         annotation_col = ann_cols,
         annotation_row = ann_rows,
         annotation_colors = ann_colors,
         cluster_rows = F,
         cluster_cols = F,
         gaps_col = c(5, 10, 13),
         filename = "InitialFigs/KO_heatmap_methylRed.png",
         width = 7,
         height = 7)
dev.off()
dev.set(dev.next())
dev.set(dev.next())



#### 12. Salt Tolerance ####
salt_orig <- read_excel("~/Desktop/Methanolobus/SaltGenes.xlsx") %>%
  filter(KO != "NA") %>%
  mutate(Type = as.factor(Type),
         Solute = as.factor(Solute)) %>%
  group_by(KO) %>%
  slice_head(n = 1)
salt_orig$Order[28] <- 66.5 # Fix glutamine order
KO_table <- KO_VST_CPM %>% # DESeq + CPM normalized KO table
  select(-KO) %>%
  mutate(sum = rowSums(.)) %>% # Get row sums
  filter(sum > 0) %>% # Remove zeroes
  mutate(KO = rownames(.)) %>%
  select(-sum)
s_meta <- KO_table %>%
  filter(KO %in% salt_orig$KO) %>%
  left_join(., salt_orig, by = "KO") %>%
  arrange(Order) %>%
  droplevels()
salt <- read_excel("~/Desktop/Methanolobus/SaltGenes.xlsx") %>%
  filter(KO != "NA") %>%
  group_by(KO) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  arrange(Order) %>%
  select(KO, Pathway_general, Code, Order)
salt$Order[46] <- 66.5 # fix glutamine order
salt <- salt %>%
  arrange(Order)
s <- KO_table %>%
  filter(KO %in% salt$KO) %>%
  left_join(., salt, by = "KO") %>%
  mutate(KO = paste(KO, Code, sep = " ")) %>%
  mutate(KO = paste(KO, Pathway_general, sep = "; ")) %>%
  arrange(Order)
s_mat <- s %>%
  select(-Pathway_general, -Code, -Order) %>%
  column_to_rownames(var = "KO") %>%
  select(-ASW_1, ASW_1) %>%
  select(-ASW_2, ASW_2) %>%
  select(-ASW_3, ASW_3) %>%
  select(-ASW_4, ASW_4) %>%
  select(-ASW_5, ASW_5) %>%
  as.matrix()
ann_rows <- data.frame(row.names = rownames(s_mat), Type = s_meta$Type, Solute = s_meta$Solute)
ann_colors <- list(Type = c(Biosynthesis = "#440154FF", Transport = "#FDE725FF"),
                   Solute = c(Betaine = hue_pal()(length(levels(s_meta$Solute)))[1],
                              Cation = hue_pal()(length(levels(s_meta$Solute)))[2],
                              Choline = hue_pal()(length(levels(s_meta$Solute)))[3],
                              Ectoine = hue_pal()(length(levels(s_meta$Solute)))[4],
                              Glutamate = hue_pal()(length(levels(s_meta$Solute)))[5],
                              Glutamine = hue_pal()(length(levels(s_meta$Solute)))[6],
                              Hydroxyectoine = hue_pal()(length(levels(s_meta$Solute)))[7],
                              Proline = hue_pal()(length(levels(s_meta$Solute)))[8],
                              Sucrose = hue_pal()(length(levels(s_meta$Solute)))[9],
                              Trehalose = hue_pal()(length(levels(s_meta$Solute)))[10]))
pheatmap(s_mat,
         legend = T,
         legend_breaks = c(-4,-3,-2,-1,0,1,2,3,4),
         legend_labels = c("-4","-3","-2","-1","0","1","2","3","Abund.\n"),
         main = "",
         gaps_col = c(5, 10, 13),
         gaps_row = c(16, 28, 31, 36, 41, 42, 43, 48, 49),
         #color = bluered(100),
         border_color = NA,
         scale = "row",
         angle_col = 315,
         annotation_row = ann_rows,
         annotation_colors = ann_colors,
         fontsize = 8,
         fontsize_row = 6,
         cluster_rows = F,
         cluster_cols = F,
         filename = "InitialFigs/KO_heatmap_salt.png",
         width = 8,
         height = 8
         )
dev.off()
dev.set(dev.next())
dev.set(dev.next())
# Note, this shows 75 KOs present in at least 1 metagenome. The list had 78, so 3 weren't present in any!



#### End Script ####