# Timberlake MAGs Analysis
# by Cliff Bueno de Mesquita, JGI, Spring 2023
# MAGs binned with 3 methods (concoct, metabat, maxbin2), then dereplicated with dRep
# 19 high quality MAGs
# See KBase narrative
# Mean coverage calculated with coverM



#### 1. Setup ####
library(plyr)
library(tidyverse)
library(ape)
library(picante)
library(dendextend)
library(readxl)
library(scales)
library(pheatmap)
library(RColorBrewer)
library(ggtree)
library(reshape2)
library(cowplot)
setwd("~/Documents/GitHub/Timberlake/")



#### 2. Tree ####
# Tree from KBase made with 49 COGs
TL_tree <- read.tree("~/Documents/GitHub/Timberlake/dRep_Tree_20-labels.newick")

# Prune the one non-MAG
TL_tree <- drop.tip(TL_tree, tip = "Frateuriasp.Soil773")

# Trim the tip labels
TL_tree$tip.label <- substr(TL_tree$tip.label, 
                            start = 1, 
                            stop = nchar(TL_tree$tip.label)-11)

# Trim the node labels
TL_tree$node.label <- as.character(round(as.numeric(TL_tree$node.label), 
                                         digits = 2))
TL_tree$node.label <- replace_na(TL_tree$node.label, replace = "")

# Plot with ape
png("~/Desktop/Timberlake/Tree.png", width = 6.5, height = 6, units = "in", res = 300)
par(oma = c(1,0,0,1))
plot.phylo(TL_tree,
           align.tip.label = T,
           no.margin = T,
           font = 1,
           cex = 0.6,
           edge.width = 2,
           show.node.label = T,
           node.pos = 1,
           label.offset = 0.005,
           adj = 0)
add.scale.bar(x = 0.5, y = 0.5)
title("MAG Phylogeny\nfrom 49 COGs", adj = 0.1, line = -3)
text(x = 0.15, y = 1.5, label = "Archaea")
text(x = 0.15, y = 6.6, label = "Bacteria")
dev.off()

# Plot with no title for multipanel
png("Figures/MAG_tree.png", width = 6.5, height = 6, units = "in", res = 300)
par(oma = c(1,0,0,1))
plot.phylo(TL_tree,
           align.tip.label = T,
           no.margin = T,
           font = 1,
           cex = 0.6,
           edge.width = 2,
           show.node.label = T,
           node.pos = 1,
           label.offset = 0.005,
           adj = 0)
add.scale.bar(x = 0.5, y = 0.5)
text(x = 0.15, y = 1.5, label = "Archaea")
text(x = 0.15, y = 6.6, label = "Bacteria")
dev.off()



#### 3. Abundance ####
sample_info <- read_excel("MetagenomeFiles.xlsx") %>%
  set_names(c("Sample", "Order", "Assembly", "Metagenome")) %>%
  mutate(Metagenome = gsub("-METAGENOME.fastq.gz", "", Metagenome)) %>%
  mutate(Metagenome = gsub("-", ".", Metagenome))
mag_info <- read.delim("CheckM_summary_table.tsv") %>%
  select(Bin.Name, Rename, Order, Completeness, Contamination) %>%
  mutate(Bin.Name = gsub(".fasta", "", Bin.Name)) %>%
  arrange(Order)
mag_abund <- read.delim("coverm_output/bin1.tsv") %>%
  add_case(read.delim("coverm_output/bin2.tsv")) %>%
  add_case(read.delim("coverm_output/bin3.tsv")) %>%
  add_case(read.delim("coverm_output/bin4.tsv")) %>%
  add_case(read.delim("coverm_output/bin5.tsv")) %>%
  add_case(read.delim("coverm_output/bin6.tsv")) %>%
  add_case(read.delim("coverm_output/bin7.tsv")) %>%
  add_case(read.delim("coverm_output/bin8.tsv")) %>%
  add_case(read.delim("coverm_output/bin9.tsv")) %>%
  add_case(read.delim("coverm_output/bin10.tsv")) %>%
  add_case(read.delim("coverm_output/bin11.tsv")) %>%
  add_case(read.delim("coverm_output/bin12.tsv")) %>%
  add_case(read.delim("coverm_output/bin13.tsv")) %>%
  add_case(read.delim("coverm_output/bin14.tsv")) %>%
  add_case(read.delim("coverm_output/bin15.tsv")) %>%
  add_case(read.delim("coverm_output/bin16.tsv")) %>%
  add_case(read.delim("coverm_output/bin17.tsv")) %>%
  add_case(read.delim("coverm_output/bin18.tsv")) %>%
  add_case(read.delim("coverm_output/bin19.tsv")) %>%
  set_names(gsub(".METAGENOME.fastq.gz.Mean", "", names(.))) %>%
  set_names(gsub("X", "", names(.))) %>%
  left_join(., mag_info, by = c("Genome" = "Rename")) %>%
  arrange(Order) %>%
  select(-Genome, -Order, -Completeness, -Contamination) %>%
  column_to_rownames(var = "Bin.Name") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "Metagenome") %>%
  left_join(., sample_info, by = "Metagenome") %>%
  column_to_rownames(var = "Sample") %>%
  arrange(Order) %>%
  select(-Order, -Metagenome, - Assembly) %>%
  t() %>%
  as.data.frame()
mag_abund_mat <- as.matrix(mag_abund)
mag_abund_mat_log <- log2(mag_abund_mat + 0.0000000001)

# Pretty heatmap
ann_cols <- data.frame(row.names = colnames(mag_abund),
                       "Treatment" = c(rep("Field", 5),
                                       rep("Control", 5),
                                       rep("SO4", 5),
                                       rep("ASW-SO4", 3), 
                                       rep("ASW", 5)))
ann_colors <- list(Treatment = c(Field = "#440154FF", 
                                 Control = "#3B528BFF", 
                                 SO4 = "#21908CFF",
                                 `ASW-SO4` = "#5DC863FF",
                                 ASW = "#FDE725FF"))
pheatmap(mag_abund_mat_log,
         legend = T,
         legend_breaks = c(-3, 0, 3, 4.2),
         legend_labels = c("-3", "0", "3", "Abund.\n"),
         main = "",
         #color = bluered(100),
         border_color = NA,
         scale = "row",
         angle_col = 315,
         fontsize = 12,
         fontsize_row = 9,
         fontsize_col = 8,
         labels_row = c(rep("", 19)),
         annotation_col = ann_cols,
         annotation_colors = ann_colors,
         cluster_rows = F,
         cluster_cols = F,
         cellwidth = 10,
         cellheight = 10,
         filename = "Figures/MAG_abund.png",
         width = 6,
         height = 6)
dev.off()
dev.set(dev.next())
dev.set(dev.next())



#### 4. Comp./Cont. ####
mag_info_mat <- read.delim("CheckM_summary_table.tsv") %>%
  select(Bin.Name, Order, Completeness, Contamination) %>%
  arrange(Order) %>%
  select(-Order) %>%
  column_to_rownames(var = "Bin.Name") %>%
  as.matrix()
pheatmap(mag_info_mat,
         legend = T,
         #color = mycolors,
         border_color = NA,
         scale = "none",
         angle_col = 315,
         fontsize = 12,
         fontsize_row = 9,
         fontsize_col = 8,
         labels_row = c(rep("", 19)),
         cluster_rows = F,
         cluster_cols = F,
         display_numbers = T,
         number_color = "black", 
         fontsize_number = 6,
         cellwidth = 30,
         cellheight = 10,
         filename = "Figures/MAG_CompCont.png",
         width = 3,
         height = 6)
dev.off()
dev.set(dev.next())
dev.set(dev.next())



#### 5. Correlations ####
# Correlations between abundance and CH4 flux
mag_abund_t <- as.data.frame(t(mag_abund_mat_log))
nc <- readRDS("input_filt_rare_mTAGs.rds")
CH4 <- nc$map_loaded %>%
  select(CH4_ug_m2_h) %>%
  t() %>%
  as.data.frame() %>%
  set_names(gsub("ASW_S", "ASW-S", names(.))) %>%
  select(rownames(mag_abund_t)) %>%
  t() %>%
  as.data.frame()
sum(rownames(CH4) != rownames(mag_abund_t))
cor_df <- as.data.frame(matrix(NA, ncol(mag_abund_t), 3)) %>%
  set_names(c("MAG", "r", "p"))
for (i in 1:ncol(mag_abund_t)) {
  m <- cor.test(CH4$CH4_ug_m2_h, mag_abund_t[,i], method = "pearson", na.rm = F)
  cor_df$MAG[i] <- names(mag_abund_t)[i]
  cor_df$r[i] <- m$estimate
  cor_df$p[i] <- m$p.value
}
cor_df <- cor_df %>%
  mutate(Pfdr = p.adjust(p, method = "fdr")) %>%
  mutate(PearsonPcut = factor(ifelse(Pfdr < 0.05, 
                                     "Pfdr < 0.05", "Pfdr > 0.05"),
                              levels = c("Pfdr < 0.05","Pfdr > 0.05")))
cor_df_mat <- cor_df %>%
  column_to_rownames(var = "MAG") %>%
  select(r) %>%
  as.matrix()
grey_red <- c("slategray","white","red4")
grey_red <- colorRampPalette(grey_red)(100)
ann_rows <- data.frame(row.names = rownames(cor_df_mat), 
                       "CH4_sig" = cor_df$PearsonPcut)
ann_colors <- list("CH4_sig" = c(`Pfdr < 0.05` = "black", 
                                 `Pfdr > 0.05` = "white"))
pheatmap(cor_df_mat,
         legend = T,
         color = grey_red,
         legend_breaks = c(-0.2, 0, 0.2, 0.4, 0.6, max(cor_df$r)),
         legend_labels = c("-0.2", "0", "0.2", "0.4", "0.6", "CH4_cor\n"),
         main = "",
         #border_color = NA,
         scale = "none",
         angle_col = 315,
         fontsize = 8,
         fontsize_row = 9,
         fontsize_col = 8,
         annotation_row = ann_rows,
         annotation_colors = ann_colors,
         labels_row = c(rep("", 19)),
         labels_col = "CH4_cor",
         cluster_rows = F,
         cluster_cols = F,
         display_numbers = T,
         number_color = "black", 
         fontsize_number = 6,
         cellwidth = 30,
         cellheight = 10,
         filename = "Figures/MAG_Cor.png",
         width = 3,
         height = 6)
dev.off()
dev.set(dev.next())
dev.set(dev.next())

  

#### 6. KO ####
# Import and merge. For some reason 3 files needed to be converted to csv 
bin1_ko <- read.delim2("MAGs_KO/Bin1_ko.txt", header = F) %>%
    select(-V1) %>%
    filter(V2 != "") %>%
    group_by(V2) %>%
    summarise(bin1 = n()) %>%
    set_names(c("KO", "Bin1"))
bin2_ko <- read.delim2("MAGs_KO/Bin2_ko.txt", header = F) %>%
  select(-V1) %>%
  filter(V2 != "") %>%
  group_by(V2) %>%
  summarise(bin2 = n()) %>%
  set_names(c("KO", "Bin2"))
bin3_ko <- read.csv("MAGs_KO/Bin3_ko.csv", header = F) %>%
  select(-V1) %>%
  filter(V2 != "") %>%
  group_by(V2) %>%
  summarise(bin3 = n()) %>%
  set_names(c("KO", "Bin3"))
bin4_ko <- read.csv("MAGs_KO/Bin4_ko.csv", header = F) %>%
  select(-V1) %>%
  filter(V2 != "") %>%
  group_by(V2) %>%
  summarise(bin4 = n()) %>%
  set_names(c("KO", "Bin4"))
bin5_ko <- read.delim2("MAGs_KO/Bin5_ko.txt", header = F) %>%
  select(-V1) %>%
  filter(V2 != "") %>%
  group_by(V2) %>%
  summarise(bin5 = n()) %>%
  set_names(c("KO", "Bin5"))
bin6_ko <- read.delim2("MAGs_KO/Bin6_ko.txt", header = F) %>%
  select(-V1) %>%
  filter(V2 != "") %>%
  group_by(V2) %>%
  summarise(bin6 = n()) %>%
  set_names(c("KO", "Bin6"))
bin7_ko <- read.delim2("MAGs_KO/Bin7_ko.txt", header = F) %>%
  select(-V1) %>%
  filter(V2 != "") %>%
  group_by(V2) %>%
  summarise(bin7 = n()) %>%
  set_names(c("KO", "Bin7"))
bin8_ko <- read.delim2("MAGs_KO/Bin8_ko.txt", header = F) %>%
  select(-V1) %>%
  filter(V2 != "") %>%
  group_by(V2) %>%
  summarise(bin8 = n()) %>%
  set_names(c("KO", "Bin8"))
bin9_ko <- read.delim2("MAGs_KO/Bin9_ko.txt", header = F) %>%
  select(-V1) %>%
  filter(V2 != "") %>%
  group_by(V2) %>%
  summarise(bin9 = n()) %>%
  set_names(c("KO", "Bin9"))
bin10_ko <- read.delim2("MAGs_KO/Bin10_ko.txt", header = F) %>%
  select(-V1) %>%
  filter(V2 != "") %>%
  group_by(V2) %>%
  summarise(bin10 = n()) %>%
  set_names(c("KO", "Bin10"))
bin11_ko <- read.delim2("MAGs_KO/Bin11_ko.txt", header = F) %>%
  select(-V1) %>%
  filter(V2 != "") %>%
  group_by(V2) %>%
  summarise(bin11 = n()) %>%
  set_names(c("KO", "Bin11"))
bin12_ko <- read.delim2("MAGs_KO/Bin12_ko.txt", header = F) %>%
  select(-V1) %>%
  filter(V2 != "") %>%
  group_by(V2) %>%
  summarise(bin12 = n()) %>%
  set_names(c("KO", "Bin12"))
bin13_ko <- read.delim2("MAGs_KO/Bin13_ko.txt", header = F) %>%
  select(-V1) %>%
  filter(V2 != "") %>%
  group_by(V2) %>%
  summarise(bin13 = n()) %>%
  set_names(c("KO", "Bin13"))
bin14_ko <- read.csv("MAGs_KO/Bin14_ko.csv", header = F) %>%
  select(-V1) %>%
  filter(V2 != "") %>%
  group_by(V2) %>%
  summarise(bin14 = n()) %>%
  set_names(c("KO", "Bin14"))
bin15_ko <- read.delim2("MAGs_KO/Bin15_ko.txt", header = F) %>%
  select(-V1) %>%
  filter(V2 != "") %>%
  group_by(V2) %>%
  summarise(bin15 = n()) %>%
  set_names(c("KO", "Bin15"))
bin16_ko <- read.delim2("MAGs_KO/Bin16_ko.txt", header = F) %>%
  select(-V1) %>%
  filter(V2 != "") %>%
  group_by(V2) %>%
  summarise(bin16 = n()) %>%
  set_names(c("KO", "Bin16"))
bin17_ko <- read.delim2("MAGs_KO/Bin17_ko.txt", header = F) %>%
  select(-V1) %>%
  filter(V2 != "") %>%
  group_by(V2) %>%
  summarise(bin17 = n()) %>%
  set_names(c("KO", "Bin17"))
bin18_ko <- read.delim2("MAGs_KO/Bin18_ko.txt", header = F) %>%
  select(-V1) %>%
  filter(V2 != "") %>%
  group_by(V2) %>%
  summarise(bin18 = n()) %>%
  set_names(c("KO", "Bin18"))
bin19_ko <- read.delim2("MAGs_KO/Bin19_ko.txt", header = F) %>%
  select(-V1) %>%
  filter(V2 != "") %>%
  group_by(V2) %>%
  summarise(bin19 = n()) %>%
  set_names(c("KO", "Bin19"))
KO_table <- full_join(bin1_ko, bin2_ko, by = "KO") %>%
  full_join(., bin3_ko, by = "KO") %>%
  full_join(., bin4_ko, by = "KO") %>%
  full_join(., bin5_ko, by = "KO") %>%
  full_join(., bin6_ko, by = "KO") %>%
  full_join(., bin7_ko, by = "KO") %>%
  full_join(., bin8_ko, by = "KO") %>%
  full_join(., bin9_ko, by = "KO") %>%
  full_join(., bin10_ko, by = "KO") %>%
  full_join(., bin11_ko, by = "KO") %>%
  full_join(., bin12_ko, by = "KO") %>%
  full_join(., bin13_ko, by = "KO") %>%
  full_join(., bin14_ko, by = "KO") %>%
  full_join(., bin15_ko, by = "KO") %>%
  full_join(., bin16_ko, by = "KO") %>%
  full_join(., bin17_ko, by = "KO") %>%
  full_join(., bin18_ko, by = "KO") %>%
  full_join(., bin19_ko, by = "KO") %>%
  column_to_rownames(var = "KO") %>%
  mutate_if(is.character, as.integer) %>%
  select(mag_info$Rename) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "Rename") %>%
  left_join(., mag_info, by = "Rename") %>%
  arrange(Order) %>%
  select(-Order, -Rename, -Completeness, -Contamination) %>%
  column_to_rownames(var = "Bin.Name")
KO_table[is.na(KO_table)] <- 0
KO_table[KO_table > 0] <- 1

# Get list of KOs
# Similar to Jinglie Figure 4?
# Or methanogenesis and phn from Methanolobus and Marivita papers
phn <- read_excel("~/Desktop/Marivita/MethylphosphonateKOs.xlsx", sheet = 1) %>%
  mutate(Name = ifelse(is.na(Name), "", Name)) %>%
  mutate(KO_def = paste(KEGG_KO, Name, sep = " ")) %>%
  mutate_if(is.character, as.factor) %>%
  arrange(Pathway_Order, Reaction_Order)
mg <- read_excel("~/Desktop/MetabolicModels/MethanogenesisRxn_for_Seed.xlsx", 
                 sheet = 3) %>%
  mutate(Name = ifelse(is.na(Name), "", Name)) %>%
  mutate(KO_def = paste(KEGG_KO, Name, sep = " ")) %>%
  mutate_if(is.character, as.factor) %>%
  arrange(Pathway_Order, Reaction_Order) %>%
  mutate(Pathway_General = "")

# Or get Wyatt significant responders - use this!
ontol <- read.csv("Ontology_KO_CNPSch4_Fm_whh.csv")
ko_list <- read.csv("SigKOs.csv") %>%
  left_join(., ontol, by = "KO") %>%
  select(KO, sm_name) %>%
  rownames_to_column(var = "rn") %>%
  filter(rn != "2") %>%
  filter(rn != "8") %>%
  filter(rn != "17") %>%
  filter(rn != "19") %>%
  filter(rn != "21") %>%
  filter(rn != "36") %>%
  filter(rn != "38") %>%
  filter(rn != "40") %>%
  filter(rn != "59") %>%
  filter(rn != "61") %>%
  select(-rn)
ko_meta <- read.csv("SigKOs.csv") %>%
  left_join(., ontol, by = "KO") %>%
  rownames_to_column(var = "rn") %>%
  filter(rn != "2") %>%
  filter(rn != "8") %>%
  filter(rn != "17") %>%
  filter(rn != "19") %>%
  filter(rn != "21") %>%
  filter(rn != "36") %>%
  filter(rn != "38") %>%
  filter(rn != "40") %>%
  filter(rn != "59") %>%
  filter(rn != "61") %>%
  select(-rn) %>%
  mutate(L1 = (gsub("Carbon ", "Carbon", L1)))
KO_table_t <- as.data.frame(t(KO_table)) %>%
  rownames_to_column(var = "KO")
sig_kos <- ko_list %>%
  left_join(., KO_table_t, by = "KO") %>%
  select(-KO) %>%
  column_to_rownames(var = "sm_name") %>%
  t() %>%
  as.data.frame()
sig_kos[is.na(sig_kos)] <- 0

# Gene Ontology and fxn Colors import
BGC_ont <- read.csv("Ontology_KO_CNPSch4_Fm_whh.csv", header = TRUE)
BGC_colors <- read.csv("Ontol_KO_L2_Color_KEY_whh.csv")
BGC_ont_col <- merge(BGC_ont, BGC_colors, by = "L2", all.x =TRUE)
BGC_ont_col <-BGC_ont_col[order(BGC_ont_col$Index.x),]

# Pretty heatmap
ann_cols <- data.frame(row.names = colnames(sig_kos),
                       "Type" = ko_meta$L1)
ann_colors <- list(Type = c(Carbon = "#b8302a", 
                            Nitrogen = "#93ce81", 
                            Phosphorus = "#4e96c8",
                            Sulfur = "#e975a8",
                            CH4_cycling = "#d9d9d9",
                            Fermentation = "#fe9829"))
pheatmap(sig_kos,
         legend = T,
         legend_breaks = c(0, 0.9, 1),
         legend_labels = c("0", "1", "KO\n"),
         main = "",
         color = c("white", "black"),
         border_color = "grey90",
         #scale = "row",
         angle_col = 90,
         fontsize = 12,
         fontsize_row = 9,
         fontsize_col = 8,
         labels_row = c(rep("", 19)),
         annotation_col = ann_cols,
         annotation_colors = ann_colors,
         cluster_rows = F,
         cluster_cols = F,
         filename = "Figures/MAG_KO.png",
         width = 8,
         height = 6)
dev.off()
dev.set(dev.next())
dev.set(dev.next())



#### 7. Combined ####
# Best way to do this would be to make a bunch of ggplots
# All of the data frames should be ready, but need to make them long form
# Tree, KO, Abundance, checkM, Cor, 

#### _KO ####
sig_kos_formelt <- sig_kos %>%
  rownames_to_column(var = "MAG")
ko_long <- melt(sig_kos_formelt, id.vars = "MAG") %>%
  mutate(MAG = gsub(".fasta", "", MAG)) %>%
  mutate(MAG = factor(MAG,
                      levels = mag_info$Bin.Name),
         value = as.factor(value))

ann_colors <- data.frame(L1 = c("Carbon", "Nitrogen", "Phosphorus", "Sulfur", 
                                "CH4_cycling", "Fermentation"),
                         colors = c("#b8302a", "#93ce81", "#4e96c8", "#e975a8",
                                             "#d9d9d9", "#fe9829"))
ko_meta_gg <- ko_meta %>%
  left_join(., ann_colors, by = "L1") %>%
  mutate(sm_name = factor(sm_name,
                          levels = ko_list$sm_name),
         L1 = factor(L1,
                     levels = c("Carbon", "Nitrogen", "Phosphorus", "Sulfur", 
                                "CH4_cycling", "Fermentation")))

p1 <- ggplot(ko_long, aes(variable, MAG, fill = value)) +
  geom_tile(color = "grey90") +
  scale_fill_manual(values = c("white", "black"),
                    labels = c("Absent", "Present")) +
  scale_y_discrete(limits = rev(levels(ko_long$MAG))) +
  labs(fill = "Gene\nP/A") +
  guides(fill = guide_legend(ncol = 1,
                             title.position = "top")) +
  theme_classic() +
  theme(legend.position = "top",
        legend.justification = c(0, 1),
        legend.direction = "vertical",
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_text(size = 6, angle = 90, hjust = 1, vjust = 0.5, 
                                   margin = margin(-2, 0, 0, 0)),
        axis.text.y = element_blank(),
        plot.margin = margin(0,0,0,0))
l1 <- get_legend(p1)
p1 <- p1 + theme(legend.position = "none")
p2 <- ggplot(ko_meta_gg, aes(x = sm_name, y = "Type", fill = L1)) +
  geom_bar(stat='identity',  width = 1) + 
  scale_fill_manual(values = c("#b8302a", "#93ce81", "#4e96c8", "#e975a8",
                                                 "#d9d9d9", "#fe9829")) +
  scale_x_discrete(expand = c(0, 0)) +
  labs(fill = "Gene class") +
  guides(fill = guide_legend(ncol = 1,
                             title.position = "top")) +
  theme(legend.position = "top",
        legend.justification = c(0.3, 1),
        legend.direction = "vertical",
        axis.title = element_blank(),                         
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(-2, 0, -2, 0), "cm"))
l2 <- get_legend(p2)
p2 <- p2 + theme(legend.position = "none")
l12 <- plot_grid(l1, l2)
panel1 <- plot_grid(p2, p1, l12, ncol = 1, align = "v", 
                    rel_heights = c(0.1, 0.75, 0.15))
panel1



#### _Abundance ####
abund_log_scaled <- as.data.frame(t(mag_abund_mat_log)) %>%
  scale() %>%
  t() %>%
  as.data.frame()
abund_formelt <- abund_log_scaled %>%
  rownames_to_column(var = "MAG")
abund_long <- melt(abund_formelt, id.vars = "MAG") %>%
  mutate(MAG = factor(MAG,
                      levels = mag_info$Bin.Name))
abund_meta <- data.frame("sampleID" = colnames(mag_abund),
                         "Treatment" = c(rep("Field", 5),
                                         rep("Control", 5),
                                         rep("SO4", 5),
                                         rep("ASW-SO4", 3), 
                                         rep("ASW", 5))) %>%
  mutate(sampleID = factor(sampleID,
                           levels = colnames(mag_abund))) %>%
  mutate(Treatment = factor(Treatment,
                            levels = c("Field", "Control", "SO4", 
                                       "ASW-SO4", "ASW")))
p3 <- ggplot(abund_long, aes(variable, MAG, fill = value)) +
  geom_tile(color = "grey90") +
  scale_fill_distiller(palette = "RdBu",
                       breaks = c(-2, 0, 2, 4),
                       labels = c("-2", "0", "2", "4")) +
  scale_y_discrete(limits = rev(levels(ko_long$MAG))) +
  labs(fill = "Abund.\nz-score") +
  theme_classic() +
  theme(legend.position = "top",
        legend.justification = c(0, 1),
        legend.direction = "vertical",
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_text(size = 6, angle = 90, hjust = 1, 
                                   vjust = 0.5, margin = margin(-2, 0, 0, 0)),
        axis.text.y = element_blank(),
        plot.margin = margin(0,0,0,0))
l3 <- get_legend(p3)
p3 <- p3 + theme(legend.position = "none")
p4 <- ggplot(abund_meta, aes(x = sampleID, y = "Type", fill = Treatment)) +
  geom_bar(stat='identity',  width = 1) + 
  scale_fill_manual(values = viridis_pal()(5)) +
  scale_x_discrete(expand = c(0, 0)) +
  labs(fill = "Treatment") +
  guides(fill = guide_legend(ncol = 1,
                             title.position = "top")) +
  theme(legend.position = "top",
        legend.justification = c(0.3, 1),
        legend.direction = "vertical",
        axis.title = element_blank(),                         
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(-2, 0, -2, 0), "cm"))
l4 <- get_legend(p4)
p4 <- p4 + theme(legend.position = "none")
l34 <- plot_grid(l3, l4)
panel2 <- plot_grid(p4, p3, l34, ncol = 1, align = "v", 
                    rel_heights = c(0.1, 0.75, 0.15))
panel2



#### _checkM ####
mag_info_mat_formelt <- as.data.frame(mag_info_mat) %>%
  rownames_to_column(var = "MAG")
cm_long <- melt(mag_info_mat_formelt, id.vars = "MAG") %>%
  mutate(MAG = gsub(".fasta", "", MAG)) %>%
  mutate(MAG = factor(MAG,
                      levels = mag_info$Bin.Name))
p5 <- ggplot(cm_long, aes(variable, MAG, fill = value)) +
  geom_tile(color = "grey90") +
  geom_text(aes(label = round(value, 2)), size = 3) +
  scale_fill_distiller(palette = "RdYlBu",
                       limits = c(0, 100),
                       breaks = c(0, 25, 50, 75, 95),
                       labels = c("0", "25", "50", "75", "100")) +
  scale_y_discrete(limits = rev(levels(ko_long$MAG))) +
  labs(fill = "%") +
  theme_classic() +
  theme(legend.position = "top",
        legend.justification = c(0, 1),
        legend.direction = "vertical",
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_text(size = 6, angle = 90, hjust = 1, 
                                   vjust = 0.5, margin = margin(-2, 0, 0, 0)),
        axis.text.y = element_blank(),
        plot.margin = margin(0,0,0,0))
l5 <- get_legend(p5)
p5 <- p5 + theme(legend.position = "none")
panel3 <- plot_grid(p5, l5, ncol = 1, 
                    rel_heights = c(0.85, 0.15))
panel3



#### _cor ####
cor_df_mat_formelt <- as.data.frame(cor_df_mat) %>%
  rownames_to_column(var = "MAG")
cor_long <- melt(cor_df_mat_formelt, id.vars = "MAG") %>%
  mutate(MAG = gsub(".fasta", "", MAG)) %>%
  mutate(MAG = factor(MAG,
                      levels = mag_info$Bin.Name))
ann_rows <- data.frame(row.names = rownames(cor_df_mat), 
                       "CH4_sig" = cor_df$PearsonPcut) %>%
  rownames_to_column(var = "MAG") %>%
  mutate(MAG = gsub(".fasta", "", MAG)) %>%
  mutate(MAG = factor(MAG,
                      levels = mag_info$Bin.Name))
p6 <- ggplot(cor_long, aes(variable, MAG, fill = value)) +
  geom_tile(color = "grey90") +
  geom_text(aes(label = round(value, 2)), size = 3) +
  scale_fill_distiller(palette = "RdGy",
                       breaks = c(-0.3, 0, 0.3, 0.6),
                       labels = c("-0.3", "0", "0.3", "0.6")) +
  scale_y_discrete(limits = rev(levels(ko_long$MAG))) +
  scale_x_discrete(labels = "CH4 cor.") +
  labs(fill = "CH4\ncor.") +
  theme_classic() +
  theme(legend.position = "top",
        legend.justification = c(1, 1),
        legend.direction = "vertical",
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_text(size = 8, angle = 90, hjust = 1, 
                                   vjust = 0.5, margin = margin(-2, 0, 0, 0)),
        axis.text.y = element_blank(),
        plot.margin = margin(0,0,0,0))
l6 <- get_legend(p6)
p6 <- p6 + theme(legend.position = "none")
p7 <- ggplot(ann_rows, aes(x = "Blah", y = MAG, fill = CH4_sig)) +
  geom_tile(stat='identity',  width = 1, color = "grey90") + 
  scale_fill_manual(values = c("black", "white"),
                    labels = c(bquote(""~P[FDR]*"<0.05"), 
                               bquote(""~P[FDR]*">0.05"))) +
  scale_y_discrete(limits = rev(levels(ko_long$MAG)),
                   expand = c(0, 0)) +
  scale_x_discrete(labels = "CH4 sig.") +
  labs(fill = "CH4 sig.") +
  guides(fill = guide_legend(ncol = 1,
                             title.position = "top")) +
  theme(legend.position = "top",
        legend.justification = c(1, 1),
        legend.direction = "vertical",
        axis.title = element_blank(),                         
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 8, angle = 90, hjust = 1, 
                                   vjust = 0.5, margin = margin(-2, 0, 0, 0)),
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
l7 <- get_legend(p7)
p7 <- p7 + theme(legend.position = "none")
l67 <- plot_grid(l7, l6)
panel4 <- plot_grid(p7, p6, l7, l6, ncol = 2, align = "h", axis = "b",
                    rel_heights = c(0.8, 0.2))
panel4


#### _Combine ####
# Tricky to get right!
# Try doing in 3 steps by row, to make sure plot rows line up

# Bars
b <- plot_grid(p2, p4, NULL, NULL, NULL, 
          ncol = 5,
          rel_widths = c(0.45, 0.3, 0.15, 0.05, 0.05))

# Plots
p <- plot_grid(p1, p3, p5, p7, p6, 
          ncol = 5, align = "h", 
          rel_widths = c(0.45, 0.3, 0.15, 0.05, 0.05))

# Legends
l <- plot_grid(l12, l34, l5, l7, l6,
          ncol = 5,
          rel_widths = c(0.46, 0.28, 0.09, 0.10, 0.07))

# BPL
bpl <- plot_grid(b, p, l,
                 ncol = 1,
                 rel_heights = c(0.04, 0.71, 0.25))
bpl

# Add tree from base graphics?
brewer_pal(palette = "Paired")(8)
par(mfrow = c(1,1))
par(omi = c(0,0,0,0),
    mai = c(0,0,0,0),
    xpd = NA)
plot.phylo(TL_tree,
           align.tip.label = T,
           no.margin = T,
           font = 1,
           cex = 0.4,
           edge.width = 2,
           show.node.label = T,
           node.pos = 1,
           label.offset = 0.08,
           adj = 1,
           x.lim = 1.8,
           edge.color = c("#A6CEE3", "black", "#1F78B4",
                          "black", "black", "#B2DF8A",
                          "#B2DF8A", "#B2DF8A", "#B2DF8A",
                          "#B2DF8A", "#B2DF8A", "#B2DF8A",
                          "#B2DF8A", "#B2DF8A", "#B2DF8A",
                          "#B2DF8A", "#B2DF8A", "#B2DF8A",
                          "#B2DF8A", "#B2DF8A","#B2DF8A",
                          "black", "black", "#33A02C",
                          "#FB9A99", "#FB9A99", "#FB9A99",
                          "#FB9A99", "#FB9A99", "#FB9A99",
                          "#FB9A99", "black", "#E31A1C",
                          "black", "#FDBF6F", "#FF7F00"))
add.scale.bar(x = 0.5, y = 0.5, cex = 0.4)
text(x = 0.16, y = 1.5, label = "Archaea", cex = 0.5)
text(x = 0.16, y = 6.6, label = "Bacteria", cex = 0.5)
t <- recordPlot()

# Use ggplot to make tree legend
t_df <- data.frame(Phylum = c("Bacteroidota", "Acidobacteriota", "Fimicutes", 
                              "Actinobacteriota", "Gemmatimonadota", "Proteobacteria", 
                              "Verrucomicrobiota", "Halobacteriota"),
                   y = c(1,2,3,4,5,6,7,8)) %>%
  mutate(Phylum = factor(Phylum, levels = Phylum))
t_leg <- get_legend(ggplot(t_df, aes(Phylum, y, fill = Phylum)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = rev(brewer_pal(palette = "Paired")(8))) +
  theme(legend.position = "top",
        legend.justification = c(0, 1),
        legend.direction = "vertical"))


t_multi <- plot_grid(NULL, t, t_leg,
                     ncol = 1,
                     rel_heights = c(0.05, 0.54, 0.42))

png("Figures/MAG_multipanel.png", width = 12, height = 7, units = "in", res = 300)
plot_grid(t_multi, bpl, rel_widths = c(0.2, 0.8))
dev.off()

# With not legend rel_heights = c(0.05, 0.54, 0.42)
tleg <- get_legend(t)
