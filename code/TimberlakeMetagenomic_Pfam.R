# Timberlake Laboratory Incubations - Metagenomic Analysis
# Pfam gene copies and spore forming gene abundances
# By Cliff Bueno de Mesquita, JGI, Fall 2023

# Procedure:
## 1. Setup (libraries, functions, working directory etc.)
## 2. Import (import and manipulate date)
## 3. Normalize
## 4. Test
## 5. Plot Spore Forming



#### 1. Setup ####
# Libraries
library(plyr)
library(tidyverse)
library(RColorBrewer)
library(DESeq2)
library(pheatmap)
library(vegan)
library(readr)
library(viridis)
library(mctoolsr)
library(readxl)

# Directory (GitHub repo)
setwd("~/Documents/GitHub/Timberlake/")

# Functions
find_hull <- function(df) df[chull(df$Axis01, df$Axis02),]
`%notin%` <- Negate(`%in%`)
save_pheatmap_pdf <- function(x, filename, width = 7, height = 7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

options(max.print = 10000)



#### 2. Import ####
# Import metadata and clean
# Fix sampleID
# Merge IMG metadata with 16S metadata and GHG
meta_16S <- read.csv("data/metadata_used.csv")
meta_IMG <- read.delim("data/IMGmetadata.txt") %>%
  separate(`Genome.Name...Sample.Name`, 
           into = c("Genome.Name", "SampleID"),
           sep = " - ") %>%
  separate(SampleID,
           into = c("Treatment", "Replicate"),
           sep = "_",
           remove = F) %>%
  arrange(taxon_oid)
meta_IMG_lab <- meta_IMG %>%
  filter(Treatment != "Sourcelab")

# Import spore forming genes from Josep Ramoneda
sf_genes <- read.csv("data/Sporulation_genes_all_pfam.csv") %>%
  mutate(PF = Pfam) %>%
  mutate(Pfam = gsub("PF", "pfam", Pfam))

# Import IMG downloaded functional profile
# Use UI_data_output file, but skip the first 9 columns with the stats output
# The number of columns will vary based on how many genome sets you select in IMG
# In this case it's 9 since 2 sets were input.
# In any case, let's first check the data frame
# Note: you could also use the stats_input file, but we would import that slightly differently
# Note: DO NOT ignore EOF warning - data loads but not all rows! Use quote = "" to avoid.
pf_names <- read.delim("data/TL_Pfam_Copies/UI_data_output.txt", quote = "") %>%
  mutate(X.Feature. = substr(X.Feature., start = 2, stop = 10))
pf <- read.delim("data/TL_Pfam_Copies/UI_data_output.txt", quote = "") %>%
  select(1, 10:ncol(.)) %>%
  select(-X.Feature.) %>%
  apply(., 2, function(y) as.numeric(gsub('"', "", y))) %>%
  as.data.frame() %>%
  mutate(pf = pf_names$X.Feature.) %>%
  column_to_rownames(var = "pf") %>%
  set_names(gsub("X.", "", names(.))) %>%
  set_names(substring(names(.), first = 1, last = 10)) %>%
  t() %>%
  as.data.frame() %>%
  arrange(rownames(.))
pf_lab <- pf %>%
  filter(rownames(.) %in% meta_IMG_lab$taxon_oid)
cs <- as.data.frame(colSums(pf_lab)) %>%
  set_names("sum") %>%
  filter(sum > 0)
pf_lab <- pf_lab %>%
  select(rownames(cs))
  
# Make sure metadata and functional data match (should be 0!)
sum(meta_IMG$taxon_oid != rownames(pf)) # 0, good
sum(meta_IMG_lab$taxon_oid != rownames(pf_lab)) # 0, good



#### 3. Normalize ####
# DESeq normalization, with no design
dds_input_1 <- DESeqDataSetFromMatrix(countData = t(pf),
                                      colData = meta_IMG,
                                      design = ~ 1)
dds_input_SF <- estimateSizeFactors(dds_input_1)
dds_input_D <- estimateDispersions(dds_input_SF)
pf_DESeq <- as.data.frame((counts(dds_input_D, normalized = T))) %>%
  t() %>%
  as.data.frame()
sum(rownames(pf_DESeq) != meta_IMG$taxon_oid)

# Lab only
pf_DESeq_lab <- DESeqDataSetFromMatrix(countData = t(pf_lab),
                                      colData = meta_IMG_lab,
                                      design = ~ 1)
pf_DESeq_lab <- estimateSizeFactors(pf_DESeq_lab)
pf_DESeq_lab <- estimateDispersions(pf_DESeq_lab)
pf_DESeq_lab <- as.data.frame((counts(pf_DESeq_lab, normalized = T))) %>%
  t() %>%
  as.data.frame()
sum(rownames(pf_DESeq_lab) != meta_IMG_lab$taxon_oid)

# Transpose
pf_DESeq_t <- pf_DESeq %>%
  t() %>%
  as.data.frame() %>%
  set_names(meta_IMG$sampleID)

pf_DESeq_lab_t <- pf_DESeq_lab %>%
  t() %>%
  as.data.frame() %>%
  set_names(meta_IMG_lab$sampleID)



#### 4. Test ####
# For differential abundance analysis, we need to give a design
# In this case, we'll use Treatment
# We'll then specify comparisons. Here do RGO vs. DS3
# We'll then run a Wald test and LRT
# We'll adjust p values with FDR
dds_input_da <- DESeqDataSetFromMatrix(countData = t(pf_lab),
                                       colData = meta_IMG_lab,
                                       design = ~ Treatment)

# Wald test - run all contrasts then concatenate
wald <- DESeq(object = dds_input_da, test = "Wald", fitType = "parametric")
cont1 <- results(wald,
                 contrast = c("Treatment", "Control", "SO4"),
                 independentFiltering = FALSE,
                 pAdjustMethod = "fdr")
cont2 <- results(wald,
                 contrast = c("Treatment", "Control", "ASW-S"),
                 independentFiltering = FALSE,
                 pAdjustMethod = "fdr")
cont3 <- results(wald,
                 contrast = c("Treatment", "Control", "ASW"),
                 independentFiltering = FALSE,
                 pAdjustMethod = "fdr")
p1 <- data.frame(Wald.p = cont1$padj,
                 pf = colnames(pf_DESeq_lab))
p2 <- data.frame(Wald.p = cont2$padj,
                 pf = colnames(pf_DESeq_lab))
p3 <- data.frame(Wald.p = cont3$padj,
                 pf = colnames(pf_DESeq_lab))
res_all <- rbind(p1, p2, p3)
res_all_sig <- res_all %>%
  filter(Wald.p < 0.05) %>%
  filter(!duplicated(pf)) %>%
  mutate(Wald = "Pfdr < 0.05") %>%
  column_to_rownames(var = "pf")
res_all_nonsig <- pf_DESeq_lab_t %>%
  set_names(meta_IMG_lab$taxon_oid) %>%
  filter(rownames(.) %notin% rownames(res_all_sig)) %>%
  mutate(Wald.p = NA,
         Wald = "Pfdr > 0.05") %>%
  select(Wald.p, Wald)
sum(nrow(res_all_sig), nrow(res_all_nonsig)) == ncol(pf_lab)
res_all_pcut <- rbind(res_all_sig, res_all_nonsig)



#### 5. Plot Spore Genes ####
# Plot all spore genes, have column showing significant or not.

# Are spore genes in the sig table? Yes - 25 of them!
sum(sf_genes$Pfam %in% rownames(res_all_sig))

# Get significance column
sf_genes_pcut <- res_all_pcut %>%
  filter(rownames(.) %in% sf_genes$Pfam)

# Are all spore genes present? No -  68 of 86
sf_genes_present <- sf_genes %>%
  filter(Pfam %in% colnames(pf_DESeq_lab)) %>%
  droplevels()

# Filter to spore genes
pf_DESeq_lab_t_sf <- pf_DESeq_lab_t %>%
  set_names(meta_IMG_lab$taxon_oid) %>%
  filter(rownames(.) %in% sf_genes_present$Pfam)

# Check order
sum(rownames(pf_DESeq_lab_t_sf) != sf_genes_present$Pfam) # Good
sum(rownames(sf_genes_pcut) != sf_genes_present$Pfam) # Bad

# Reorder and check order
reorder_idx <- match(sf_genes_present$Pfam, rownames(sf_genes_pcut))
sf_genes_pcut <- sf_genes_pcut[reorder_idx, ]
sum(rownames(sf_genes_pcut) != sf_genes_present$Pfam) # Good

# Add definition
sf_genes_present$Def <- paste(sf_genes_present$Pfam, sf_genes_present$Shortname,
                              sep = ": ")

# Plot
pf_DESeq_lab_t_sf <- pf_DESeq_lab_t_sf %>%
  set_names(meta_IMG_lab$SampleID) %>%
  select(Control_1, Control_2, Control_3, Control_4, Control_5,
         SO4_1, SO4_2, SO4_3, SO4_4, SO4_5,
         `ASW-S_2`, `ASW-S_4`, `ASW-S_5`,
         ASW_1, ASW_2, ASW_3, ASW_4, ASW_5)
rownames(pf_DESeq_lab_t_sf) <- sf_genes_present$Def

ann_cols <- data.frame(row.names = colnames(pf_DESeq_lab_t_sf),
                       "Treatment" = c(rep("Control", 5),
                                       rep("SO4", 5),
                                       rep("ASW-S", 3),
                                       rep("ASW", 5)))
ann_rows <- data.frame(row.names = rownames(pf_DESeq_lab_t_sf),
                       Wald = sf_genes_pcut$Wald)
ann_colors <- list(Treatment = c(Control = "#3B528BFF", 
                                 SO4 = "#21908CFF",
                                 `ASW-S` = "#5DC863FF",
                                 ASW = "#FDE725FF"),
                   Wald = c(`Pfdr < 0.05` = "black",
                            `Pfdr > 0.05` = "white"))

pheatmap(pf_DESeq_lab_t_sf, # data
         legend = T, # show legend
         legend_breaks = c(-4, -2, 0, 2, 3.8, 4), # adjust depending on data
         legend_labels = c("-4","-2","0","2","4","Abund.\n"), # adjust depending on data
         main = "", # gives extra space at the top for the legend label
         border_color = NA, # no border color
         scale = "row", # z-scores, by row
         angle_col = 315, # angle of the bottom column labels
         fontsize = 8, # master font size
         fontsize_row = 6, # row font size
         fontsize_col = 6, # column font size
         annotation_col = ann_cols, # column annotations
         annotation_row = ann_rows, # row annotations
         annotation_colors = ann_colors, # annotation colors
         cluster_rows = F, # cluster by rows
         cluster_cols = F, # cluster by column
         filename = "InitialFigs/SporeGenes.png", # file name
         width = 8, # width
         height = 12 # height
         ) 
dev.off()
dev.set(dev.next())
dev.set(dev.next())



#### 6. Firmicutes ####
# Check if spore gene abundance is related to Firmicutes abundance

#### _mTags ####
# Average the DESeq abundance of the 68 spore genes by sample
sf_mean <- colMeans(pf_DESeq_lab_t_sf) %>%
  as.data.frame() %>%
  rename("." = "mean_abund")

# Get mtags Firmicutes
d <- readRDS("data/input_filt_rare_mTAGs.rds")
phyla <- summarize_taxonomy(d, 2, relative = T, report_higher_tax = F) %>%
  select(-Field_1, -Field_2, -Field_3, -Field_4, -Field_5) %>%
  select(Control_1, Control_2, Control_3, Control_4, Control_5,
         SO4_1, SO4_2, SO4_3, SO4_4, SO4_5,
         `ASW_S_2`, `ASW_S_4`, `ASW_S_5`,
         ASW_1, ASW_2, ASW_3, ASW_4, ASW_5) %>%
  t() %>%
  as.data.frame() %>%
  select(Firmicutes) %>%
  mutate(mean_abund = sf_mean$mean_abund) %>%
  mutate("Treatment" = c(rep("Control", 5),
                         rep("SO4", 5),
                         rep("ASW-S", 3),
                         rep("ASW", 5))) %>%
  mutate(Treatment = factor(Treatment,
                            levels = c("Control", "SO4", "ASW-S", "ASW")))

m <- lm(mean_abund ~ Firmicutes, data = phyla)
summary(m)

png("InitialFigs/SporeGenes_Firmicutes.png", width = 7, height = 5, res = 300, units = "in")
ggplot(phyla, aes(Firmicutes*100, mean_abund, colour = Treatment)) +
  geom_point(size = 3, shape = 16, alpha = 1) +
  geom_smooth(aes(Firmicutes*100, mean_abund),
              method = "lm", inherit.aes = F, alpha = 0.1) +
  geom_text(aes(x = 0, y = Inf, hjust = -0.1, vjust = 2, label = "R^2 == 0.39"), 
            parse = TRUE, size = 4,  check_overlap = TRUE, show.legend = F, color = "black") +
  geom_text(aes(x = 0, y = Inf, hjust = -0.1, vjust = 5, label = "p = 0.005"), 
            size = 4, check_overlap = TRUE, show.legend = F, color = "black") +
  scale_colour_manual(values = c("#3B528BFF", "#21908CFF", "#5DC863FF", "#FDE725FF")) +
  labs(x = "Firmicutes % abundance (mTags)",
       y = "Mean spore gene DESeq2 abundance") +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 10))
dev.off()


#### _iTags ####
# Now check iTags too
nc <- readRDS("data/nc.rds")
nc_cn <- read_excel("~/Documents/GitHub/EastCoast/data/Copy of CHN data.xls", sheet = 2) %>%
  #slice(12:91) %>%
  slice_tail(n = 80) %>%
  dplyr::select(Name, `%N`, `%C`, `C:N`) %>%
  mutate(Name = gsub("bottom", "bot", Name)) %>%
  separate(Name, into = c("ID", "Depth"), sep = " ", remove = F) %>%
  mutate(Treatment = substr(ID, start = 1, stop = 1),
         Hydro = substr(ID, start = 2, stop = 2),
         Replicate = substr(ID, start = 3, stop = 3)) %>%
  filter(Hydro == "F") %>%
  mutate(Treatment = recode_factor(Treatment,
                                   "A" = "Control",
                                   "B" = "ASW",
                                   "C" = "ASW_noS",
                                   "D" = "SO4")) %>%
  mutate(Treatment = factor(Treatment,
                            levels = c("ASW", "ASW_noS", "Control", "SO4"))) %>%
  arrange(desc(Depth), Treatment) %>%
  filter(Name %notin% c("BF1 top", "BF2 top", "AF1 top", "DF2 top", "BF2 bot", "CF1 bot", "AF2 bot"))
rownames(nc$map_loaded)[11:43]
nc$map_loaded$sed_per_C[11:43] <- nc_cn$`%C`
nc$map_loaded$sed_per_N[11:43] <- nc_cn$`%N`
nc$map_loaded$sed_CN[11:43] <- nc_cn$`C:N`
for (i in 1:nrow(nc$map_loaded)) {
  if (nc$map_loaded$sampleID[i] == "TL_inc_d1_SO4_5B") {
    nc$map_loaded$Depth[i] <- "10-15 cm"
    nc$map_loaded$sampleID_clean[i] <- "+SO4_5B_D2"
  }
}
for (i in 1:nrow(nc$map_loaded)) {
  if (nc$map_loaded$sampleID[i] == "TL_inc_d2_SO4_5B") {
    nc$map_loaded$Depth[i] <- "0-5 cm"
    nc$map_loaded$sampleID_clean[i] <- "+SO4_5B_D1"
  }
}
new_dat <- read_excel("data/Soil pH-Jessie.xls", sheet = 3, na = "NA")
nc$map_loaded[25, 10:24] <- new_dat[1, 2:16]
nc$map_loaded[26, 10:24] <- new_dat[2, 2:16]
nc$map_loaded[42, 10:24] <- new_dat[3, 2:16]
nc$map_loaded[43, 10:24] <- new_dat[4, 2:16]
nc$map_loaded$SO4 <- recode_factor(nc$map_loaded$Treatment,
                                   "+SO4" = "Sulfate",
                                   "+ASW" = "SaltSulfate")
nc$map_loaded$SO4 <- grepl(x = nc$map_loaded$SO4, pattern = "Sulfate")
nc$map_loaded$ASW <- recode_factor(nc$map_loaded$Treatment,
                                   "+ASW-SO4" = "Salt",
                                   "+ASW" = "SaltSulfate")
nc$map_loaded$ASW <- grepl(x = nc$map_loaded$ASW, pattern = "Salt")

# Get just depth 2, lab, match MG
# Note: missing Control_2, ASW_2, SO4_1
# Plus SO4 5B can't match. 
nc <- filter_data(nc, filter_cat = "Depth", keep_vals = "10-15 cm")
nc <- filter_data(nc, filter_cat = "Treatment", filter_vals = "Field")
nc <- filter_data(nc, filter_cat = "sampleID", filter_vals = "TL_inc_d2_ASW_noS_3")

sf_mean_iTag <- sf_mean %>%
  t() %>%
  as.data.frame() %>%
  select(Control_1, Control_3, Control_4, Control_5,
         SO4_2, SO4_3, SO4_4, SO4_5,
         `ASW-S_2`, `ASW-S_4`, `ASW-S_5`,
         ASW_1, ASW_3, ASW_4, ASW_5) %>%
  t() %>%
  as.data.frame()
  

phyla_iTag <- summarize_taxonomy(nc, 2, relative = T, report_higher_tax = F) %>%
  select(TL_inc_d2_DI_ctrl_1, TL_inc_d2_DI_ctrl_3, TL_inc_d2_DI_ctrl_4, TL_inc_d2_DI_ctrl_5,
         TL_inc_d2_SO4_2, TL_inc_d2_SO4_3, TL_inc_d2_SO4_4, TL_inc_d2_SO4_5A,
         `TL_inc_d2_ASW_noS_2`, `TL_inc_d2_ASW_noS_4`, `TL_inc_d2_ASW_noS_5`,
         TL_inc_d2_ASW_1, TL_inc_d2_ASW_3, TL_inc_d2_ASW_4, TL_inc_d2_ASW_5) %>%
  t() %>%
  as.data.frame() %>%
  select(Firmicutes) %>%
  mutate(mean_abund = sf_mean_iTag$mean_abund) %>%
  mutate("Treatment" = c(rep("Control", 4),
                         rep("SO4", 4),
                         rep("ASW-S", 3),
                         rep("ASW", 4))) %>%
  mutate(Treatment = factor(Treatment,
                            levels = c("Control", "SO4", "ASW-S", "ASW")))

m <- lm(mean_abund ~ Firmicutes, data = phyla_iTag)
summary(m)

png("InitialFigs/SporeGenes_Firmicutes.png", width = 7, height = 5, res = 300, units = "in")
ggplot(phyla_iTag, aes(Firmicutes*100, mean_abund, colour = Treatment)) +
  geom_point(size = 3, shape = 16, alpha = 1) +
  geom_smooth(aes(Firmicutes*100, mean_abund),
              method = "lm", inherit.aes = F, alpha = 0.1) +
  geom_text(aes(x = 0, y = Inf, hjust = -0.1, vjust = 2, label = "R^2 == 0.37"), 
            parse = TRUE, size = 4,  check_overlap = TRUE, show.legend = F, color = "black") +
  geom_text(aes(x = 0, y = Inf, hjust = -0.1, vjust = 5, label = "p = 0.017"), 
            size = 4, check_overlap = TRUE, show.legend = F, color = "black") +
  scale_colour_manual(values = c("#3B528BFF", "#21908CFF", "#5DC863FF", "#FDE725FF")) +
  labs(x = "Firmicutes % abundance (iTags)",
       y = "Mean spore gene DESeq2 abundance") +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 10))
dev.off()

# mTag vs iTag Firmicutes
phyla_15 <- phyla %>%
  t() %>%
  as.data.frame() %>%
  select(Control_1, Control_3, Control_4, Control_5,
         SO4_2, SO4_3, SO4_4, SO4_5,
         `ASW_S_2`, `ASW_S_4`, `ASW_S_5`,
         ASW_1, ASW_3, ASW_4, ASW_5) %>%
  t() %>%
  as.data.frame()

phyla_iTag$Firmicutes_mTag <- as.numeric(phyla_15$Firmicutes)

m <- lm(Firmicutes_mTag ~ Firmicutes, data = phyla_iTag)
summary(m)

ggplot(phyla_iTag, aes(Firmicutes*100, Firmicutes_mTag*100, colour = Treatment)) +
  geom_point(size = 3, shape = 16, alpha = 1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_smooth(aes(Firmicutes*100, Firmicutes_mTag*100),
              method = "lm", inherit.aes = F, alpha = 0.1) +
  geom_text(aes(x = 0, y = Inf, hjust = -0.1, vjust = 2, label = "R^2 == 0.95"), 
            parse = TRUE, size = 4,  check_overlap = TRUE, show.legend = F, color = "black") +
  geom_text(aes(x = 0, y = Inf, hjust = -0.1, vjust = 5, label = "p < 0.001"), 
            size = 4, check_overlap = TRUE, show.legend = F, color = "black") +
  scale_colour_manual(values = c("#3B528BFF", "#21908CFF", "#5DC863FF", "#FDE725FF")) +
  labs(x = "Firmicutes % abundance (iTags)",
       y = "Firmicutes % abundance (mTags") +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 10))


#### End Script ####