# Timberlake 16S data analysis
# by Cliff Bueno de Mesquita, Tringe Lab, JGI, Summer/Fall 2022



#### 1. Overview ####
# Samples were sequenced and processed with iTagger at JGI
# Cliff reassigned taxonomy with SILVA v 138.1 using DADA2 in R
# Chloroplast, Mitochondia, Domain level - NAs filtered out
# Input data includes other samples from SF and East Coast, filter them out

# Key Goals:
# Assess microbial taxonomic and functional
# Relate microbial taxa and genes to GHG emissions

# Key papers/experimental information:
# Ardón et al. 2013 Global Change Biology
# Ardón et al. 2016 Biogeochemistry
# Ardón et al. 2018 Biogeochemistry
## 50 samples: 10 initial, 40 experimental (2 x 2 x 2 factorial)
## Control = DI water
## +ASW = artificial salt-water with 5 ppt salinity
## +ASW-SO4 = artificial salt-water with 5 ppt salinity without sulfate
## +SO4 = DI water with same amount of sulfate as ASW, but only sulfate
## Flooded samples (water level maintained at surface)
## There was also a drought treatment but it wasn't sequenced
## 30˚C, 12 weeks

# Analysis Outline
# The analysis contains alpha and beta diversity analyses as well as taxonomic analyses
# Taxonomic analyses include simper and multipatt indicator analyses
# Stacked bar plots are made for each taxonomic level, including functional guilds
# Follow the Document Outline in RStudio (on the right) to navigate among different sections
# Sections are:
## 1. This overview with background information
## 2. Setup (libraries, metadata, guild calling, mctoolsr object, subsetting, rarefaction)
## 3. Alpha
## 4. Beta
## 5. Taxa
## 6. Indicators
## 7. Correlations



#### 2. Setup ####
library(plyr) # Data manipulation
library(tidyverse) # Data manipulation
library(mctoolsr) # Microbial analyses
library(RColorBrewer) # Colors
library(vegan) # Multivariate analyses
library(indicspecies) # Indicator species
library(car) # Stats
library(FSA) # SE
library(magrittr) # Set names
library(PMCMRplus) # Stats
library(readxl) # Excel
library(writexl) # Excel
library(plotly) # Interactive plots
library(ggmap) # Maps
library(ggsn) # Maps
library(multcomp) # Tukey HSD and significance letters
library(emmeans) # Tukey HSD and significance letters
library(scales) # View colors
library(cowplot) # Multipanels
library(qvalue) # q values for indicator species
library(reshape2) # melt
library(gridExtra) # graphs
library(grid) # graphs
library(cowplot) # graphs
library(ggpubr) # graphs
library(ggExtra) # graphs
library(ggh4x) # graphs
library(dendextend) # graphs
library(corrplot) # correlation plots

# Functions
find_hull <- function(df) df[chull(df$Axis01, df$Axis02),]
`%notin%` <- Negate(`%in%`)

# Guild subsetting module from other repo
source("~/Documents/GitHub/SF_microbe_methane/modules/3_OTU_subsetting_modules_v.0.4_strip.r")

# Correlation functions from other repo
source("~/Documents/GitHub/EastCoast/meth_corr_by_taxonomy.R")
source("~/Documents/GitHub/EastCoast/meth_corr_by_bgc.R")

# Repository path
setwd("~/Documents/GitHub/Timberlake/")

# Wyatt Hartman's guild color palette
# Note that MeOB don't exist in this dataset, so removed
# But ANME do exist in this dataset, so add
# Extra methanogen guilds added so colors added too
Guild_cols <- read.table("~/Documents/GitHub/SF_microbe_methane/data/colors/Guild_color_palette.txt",
                         sep='\t') %>%
  dplyr::select(Guild, G_index, color) %>%
  set_names(c("Guild", "Index", "color")) %>%
  mutate(Index = rev(Index)) %>%
  add_row(Guild = "ANME", Index = 10, color = "#836FFF") %>%
  add_row(Guild = "CH4_me", Index = 16, color = "#FDC086") %>%
  add_row(Guild = "CH4_mix", Index = 17, color = "#FFFF99") %>%
  filter(Guild != "MeOB") %>%
  arrange(Index)

# Prepare data. Only need to do once. Then skip to start here.
# Read in combined dataset (with other samples)
input_filt <- readRDS("input_filt_comb_wBGC.rds")
# Get only Timberlake samples
input_filt_nc <- filter_data(input_filt,
                             filter_cat = "Estuary",
                             keep_vals = "Alligator")
# Remove extra samples
input_filt_nc <- filter_data(input_filt_nc,
                             filter_cat = "sampleID",
                             filter_vals = c("TL_nw_d1_DI_ctrl_AF1",
                                             "TL_nw_d1_DI_ctrl_AF3",
                                             "TL_nw_d1_DI_ctrl_AF4",
                                             "TL_nw_d1_ASW_noS_BF3",
                                             "TL_nw_d1_ASW_noS_BF4",
                                             "TL_nw_d1_ASW_noS_BF5",
                                             "TL_inc_d1_SO4_5A",
                                             "TL_inc_d1_SO4_5B",
                                             "TL_inc_d2_SO4_5A",
                                             "TL_inc_d2_SO4_5B"))
# Remove unneeded columns from metadata and format factors
input_filt_nc$map_loaded <- input_filt_nc$map_loaded %>%
  dplyr::select(-Site, -Estuary, -Info, -Conductivity_uS_cm, -CH4_pw_air_ppmv, -NEE_mgC_m2_m,
                -GEP_mgC_m2_m, -PAR_uE_m2_s, -CH4_pot_umol_gdw_h, -CO2_pot_umol_gdw_h,
                -Fe_mgL, -Acetate_mgL, -TotalVFA_uM, -SR_umol_cm3_d, -AMG_umol_cm3_d,
                -c(36:67)) %>%
  rename(Treatment = Detail,
         Salinity = Salinity_ppt_all) %>%
  mutate(Treatment = as.factor(Treatment),
         Depth = as.factor(Depth)) %>%
  mutate(Treatment = recode_factor(Treatment, "Field Reference" = "Field", "DI_ctrl" = "Control", 
                                   "5ppt ASW" = "+ASW", "SW_noSO4" = "+ASW-SO4", "SO4 amended" = "+SO4"),
         Depth = recode_factor(Depth, "0.025" = "0-5 cm", "0.125" = "10-15 cm")) %>%
  mutate(Treatment = factor(Treatment,
                            levels = c("Field", "Control", "+SO4", "+ASW-SO4", "+ASW")),
         TrtDepth = paste(Treatment, Depth, sep = ""))
# Rarefy at 82350
sort(colSums(input_filt_nc$data_loaded))
mean(colSums(input_filt_nc$data_loaded)) # 115713.5
se(colSums(input_filt_nc$data_loaded)) # 2427.631
set.seed(530)
nc <- single_rarefy(input_filt_nc, 82350) # Now n = 343
sort(colSums(nc$data_loaded))
#saveRDS(nc, "nc.rds")

# Start here
nc <- readRDS("nc.rds")



#### 3. Alpha ####
leveneTest(nc$map_loaded$rich ~ nc$map_loaded$Treatment) # Homogeneous
m <- aov(rich ~ Treatment + Depth, data = nc$map_loaded)
Anova(m, type = "II") # Treatment and Depth
m <- aov(rich ~ Treatment, data = nc$map_loaded)
shapiro.test(m$residuals) # Normal
summary(m)
t <- emmeans(object = m, specs = "Treatment") %>%
  cld(object = ., adjust = "Tukey", Letters = letters, alpha = 0.05) %>%
  mutate(name = "rich",
         y = max(nc$map_loaded$rich)+(max(nc$map_loaded$rich)-min(nc$map_loaded$rich))/20)
leveneTest(nc$map_loaded$shannon ~ nc$map_loaded$Treatment) # Homogeneous
m1 <- aov(shannon ~ Treatment + Depth, data = nc$map_loaded)
Anova(m1, type = "II") # Treatment and Depth
m1 <- aov(shannon ~ Treatment, data = nc$map_loaded)
shapiro.test(m1$residuals) # Normal
summary(m1)
t1 <- emmeans(object = m1, specs = "Treatment") %>%
  cld(object = ., adjust = "Tukey", Letters = letters, alpha = 0.05) %>%
  mutate(name = "shannon",
         y = max(nc$map_loaded$shannon)+(max(nc$map_loaded$shannon)-min(nc$map_loaded$shannon))/20)
label_df <- rbind(t, t1)
facet_df <- c("rich" = "(a) Richness",
              "shannon" = "(b) Shannon")
alpha_long <- nc$map_loaded %>%
  pivot_longer(cols = c("rich", "shannon"))
pdf("Figures/Alpha.pdf", width = 6, height = 3)
ggplot(alpha_long, aes(reorder(Treatment, value, mean), value, 
                       colour = Treatment)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 2, alpha = 0.75, width = 0.2, aes(shape = Depth)) +
  geom_text(data = label_df, aes(Treatment, y, label = str_trim(.group)), 
            size = 4, color = "black") +
  labs(x = "Site", y = "Number of OTUs", shape = "Depth (cm)") +
  scale_x_discrete(labels = c("+ASW", "+ASW (no SO4)", "+SO4", "Control", "Field")) +
  scale_colour_viridis_d() +
  guides(colour = "none") +
  facet_wrap(~ name, ncol = 2, scales = "free_y", labeller = as_labeller(facet_df)) +
  theme_bw() +
  theme(legend.position = c(1,0),
        legend.justification = c(1,0),
        legend.background = element_blank(),
        legend.title.align = 0.5,
        axis.title = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 6),
        strip.text = element_text(size = 10))
dev.off()

#### 4. Beta ####
nc_bc <- calc_dm(nc$data_loaded)
set.seed(1150)
adonis2(nc_bc ~ nc$map_loaded$Treatment + nc$map_loaded$Depth) # Both sig
anova(betadisper(nc_bc, nc$map_loaded$Treatment)) # Dispersion not homogeneous
anova(betadisper(nc_bc, nc$map_loaded$Depth)) # Dispersion homogeneous
nc_pcoa <- cmdscale(nc_bc, k = nrow(nc$map_loaded) - 1, eig = T)
pcoaA1 <- round((eigenvals(nc_pcoa)/sum(eigenvals(nc_pcoa)))[1]*100, digits = 1)
pcoaA2 <- round((eigenvals(nc_pcoa)/sum(eigenvals(nc_pcoa)))[2]*100, digits = 1)
nc$map_loaded$Axis01 <- scores(nc_pcoa)[,1]
nc$map_loaded$Axis02 <- scores(nc_pcoa)[,2]
micro.hulls <- ddply(nc$map_loaded, c("Treatment"), find_hull)
pdf("InitialFigs/NC_PCoA.pdf", width = 7, height = 5)
ggplot(nc$map_loaded, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = Treatment, fill = Treatment),
               alpha = 0.1, show.legend = F) +
  geom_point(size = 3, alpha = 1, aes(colour = Treatment, shape = Depth)) +
  labs(x = paste("PC1: ", pcoaA1, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2, "%", sep = ""),
       shape = "Depth (m)") +
  scale_colour_viridis_d(labels = c("Field", "Control", "+SO4", "+ASW (no SO4)", "+ASW"),
                         breaks = c("Field Reference", "DI_ctrl", "SO4 amended",
                                    "SW_noSO4", "5ppt ASW")) +
  scale_fill_viridis_d() +
  guides(colour = guide_legend(override.aes = list(shape = 15),
                               order = 1)) +
  theme_bw() +  
  theme(legend.position = c(0,0),
        legend.justification = c(0,0),
        legend.background = element_blank(),
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        plot.margin = margin(5, 5, 5, 5, "pt"))
dev.off()

#### 5. Taxa ####
nc_phyla <- summarize_taxonomy(nc, level = 2, report_higher_tax = F)
plot_ts_heatmap(nc_phyla, nc$map_loaded, 0.01, 'Treatment', rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
nc_barsP <- plot_taxa_bars(nc_phyla, nc$map_loaded, "Treatment", num_taxa = 12, data_only = T) %>%
  mutate(group_by = factor(group_by, levels = c("Field Reference","DI_ctrl",
                                                "SO4 amended", "SW_noSO4",
                                                "5ppt ASW"))) %>%
  mutate(taxon = fct_rev(taxon))
pdf("InitialFigs/NC_Phyla.pdf", width = 7, height = 5)
ggplot(nc_barsP, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Sample", y = "Relative abundance", fill = "Phylum") +
  scale_fill_manual(values = c("grey90", brewer.pal(12, "Paired")[12:1])) +
  scale_y_continuous(expand = c(0.01, 0.01)) + 
  scale_x_discrete(labels = c("Field", "Control", "+SO4", "+ASW (no SO4)", "+ASW")) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 8))
dev.off()
taxa_summary_by_sample_type(nc_phyla, nc$map_loaded, 'Treatment', 0.01, 'KW')

nc_class <- summarize_taxonomy(nc, level = 3, report_higher_tax = F)
plot_ts_heatmap(nc_class, nc$map_loaded, 0.01, 'Treatment', rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
nc_barsC <- plot_taxa_bars(nc_class, nc$map_loaded, "Treatment", num_taxa = 12, data_only = T) %>%
  mutate(group_by = factor(group_by, levels = c("Field Reference","DI_ctrl",
                                                "SO4 amended", "SW_noSO4",
                                                "5ppt ASW"))) %>%
  mutate(taxon = gsub("NA", "Unclassified", taxon)) %>%
  mutate(taxon = fct_relevel(taxon, "Other", after = Inf)) %>%
  mutate(taxon = fct_relevel(taxon, "Unclassified", after = Inf)) %>%
  mutate(taxon = fct_rev(taxon))
ggplot(nc_barsC, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Sample", y = "Relative abundance", fill = "Class") +
  scale_fill_manual(values = c("grey75", "grey90", brewer.pal(12, "Paired")[1:11])) +
  scale_y_continuous(expand = c(0.01, 0.01)) + 
  scale_x_discrete(labels = c("Field", "Control", "+SO4", "+ASW (no SO4)", "+ASW")) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 8))
taxa_summary_by_sample_type(nc_class, nc$map_loaded, 'Treatment', 0.01, 'KW')

nc_order <- summarize_taxonomy(nc, level = 4, report_higher_tax = F)
plot_ts_heatmap(nc_order, nc$map_loaded, 0.01, 'Treatment', rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
nc_barsO <- plot_taxa_bars(nc_order, nc$map_loaded, "Treatment", num_taxa = 12, data_only = T) %>%
  mutate(group_by = factor(group_by, levels = c("Field Reference","DI_ctrl",
                                                "SO4 amended", "SW_noSO4",
                                                "5ppt ASW"))) %>%
  mutate(taxon = gsub("NA", "Unclassified", taxon)) %>%
  mutate(taxon = fct_relevel(taxon, "Other", after = Inf)) %>%
  mutate(taxon = fct_relevel(taxon, "Unclassified", after = Inf)) %>%
  mutate(taxon = fct_rev(taxon))
ggplot(nc_barsO, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Sample", y = "Relative abundance", fill = "Order") +
  scale_fill_manual(values = c("grey75", "grey90", brewer.pal(12, "Paired")[1:11])) +
  scale_y_continuous(expand = c(0.01, 0.01)) + 
  scale_x_discrete(labels = c("Field", "Control", "+SO4", "+ASW (no SO4)", "+ASW")) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 8))
taxa_summary_by_sample_type(nc_order, nc$map_loaded, 'Treatment', 0.01, 'KW')

nc_family <- summarize_taxonomy(nc, level = 5, report_higher_tax = F)
plot_ts_heatmap(nc_family, nc$map_loaded, 0.01, 'Treatment', rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
nc_barsF <- plot_taxa_bars(nc_family, nc$map_loaded, "Treatment", num_taxa = 12, data_only = T) %>%
  mutate(group_by = factor(group_by, levels = c("Field Reference","DI_ctrl",
                                                "SO4 amended", "SW_noSO4",
                                                "5ppt ASW"))) %>%
  mutate(taxon = gsub("NA", "Unclassified", taxon)) %>%
  mutate(taxon = fct_relevel(taxon, "Other", after = Inf)) %>%
  mutate(taxon = fct_relevel(taxon, "Unclassified", after = Inf)) %>%
  mutate(taxon = fct_rev(taxon))
ggplot(nc_barsF, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Sample", y = "Relative abundance", fill = "Family") +
  scale_fill_manual(values = c("grey75", "grey90", brewer.pal(12, "Paired")[1:11])) +
  scale_y_continuous(expand = c(0.01, 0.01)) + 
  scale_x_discrete(labels = c("Field", "Control", "+SO4", "+ASW (no SO4)", "+ASW")) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 8))
taxa_summary_by_sample_type(nc_family, nc$map_loaded, 'Treatment', 0.01, 'KW')

nc_genus<- summarize_taxonomy(nc, level = 6, report_higher_tax = F)
plot_ts_heatmap(nc_genus, nc$map_loaded, 0.01, 'Treatment', rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
nc_barsG <- plot_taxa_bars(nc_genus, nc$map_loaded, "Treatment", num_taxa = 12, data_only = T) %>%
  mutate(group_by = factor(group_by, levels = c("Field Reference","DI_ctrl",
                                                "SO4 amended", "SW_noSO4",
                                                "5ppt ASW"))) %>%
  mutate(taxon = gsub("NA", "Unclassified", taxon)) %>%
  mutate(taxon = fct_relevel(taxon, "Other", after = Inf)) %>%
  mutate(taxon = fct_relevel(taxon, "Unclassified", after = Inf)) %>%
  mutate(taxon = fct_rev(taxon))
ggplot(nc_barsG, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Sample", y = "Relative abundance", fill = "Genus") +
  scale_fill_manual(values = c("grey75", "grey90", brewer.pal(12, "Paired")[1:11])) +
  scale_y_continuous(expand = c(0.01, 0.01)) + 
  scale_x_discrete(labels = c("Field", "Control", "+SO4", "+ASW (no SO4)", "+ASW")) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 8))
taxa_summary_by_sample_type(nc_genus, nc$map_loaded, 'Treatment', 0.01, 'KW')

nc_guilds <- summarize_taxonomy(nc, level = 9, report_higher_tax = F)
plot_ts_heatmap(nc_guilds, nc$map_loaded, 0.01, 'Treatment', rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
nc_barsGu <- plot_taxa_bars(nc_guilds, nc$map_loaded, "Treatment", 
                            num_taxa = 20, data_only = T) %>%
  filter(taxon != "NA") %>%
  droplevels() %>%
  mutate(taxon = factor(taxon,
                        levels = Guild_cols$Guild)) %>%
  mutate(group_by = factor(group_by, levels = c("Field Reference","DI_ctrl",
                                                "SO4 amended", "SW_noSO4",
                                                "5ppt ASW")))
tallest_bar <- nc_barsGu %>%
  group_by(group_by) %>%
  summarise(sum = sum(mean_value))
pdf("InitialFigs/NC_Guilds.pdf", width = 7, height = 5)
ggplot(nc_barsGu, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Treatment", y = "Relative abundance", fill = "Guild") +
  scale_fill_manual(values = Guild_cols$color) +
  scale_y_continuous(expand = c(max(tallest_bar$sum)/100, max(tallest_bar$sum)/100)) + 
  scale_x_discrete(labels = c("Field", "Control", "+SO4", "+ASW (no SO4)", "+ASW")) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 8))
dev.off()
taxa_summary_by_sample_type(nc_guilds, 
                            nc$map_loaded, 
                            type_header = 'Treatment', 
                            filter_level = 0.01, 
                            test_type = 'KW')

# Guilds by sample (to show variation)
nc$map_loaded$sampleID <- rownames(nc$map_loaded)
nc_barsGu <- plot_taxa_bars(nc_guilds,
                            nc$map_loaded,
                            "sampleID",
                            num_taxa = 20,
                            data_only = TRUE) %>%
  filter(taxon != "NA") %>%
  droplevels() %>%
  mutate(taxon = factor(taxon,
                        levels = Guild_cols$Guild)) %>%
  left_join(., nc$map_loaded, by = c("group_by" = "sampleID")) %>%
  mutate(Treatment = factor(Treatment,
                            levels = c("Field Reference", "DI_ctrl",
                                       "SO4 amended", "SW_noSO4", "5ppt ASW")))
facet_names <- c("0.025" = "Depth 2.5 cm",
                 "0.125" = "Depth 12.5 cm",
                 "Field Reference" = "Field",
                 "DI_ctrl" = "Control",
                 "SO4 amended" = "+SO4",
                 "SW_noSO4" = "+ASW (no SO4)",
                 "5ppt ASW" = "+ASW")
tallest_bar <- nc_barsGu %>%
  group_by(group_by) %>%
  summarise(sum = sum(mean_value))
pdf("InitialFigs/NC_Guilds_samples.pdf", width = 9, height = 6)
ggplot(nc_barsGu, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Sample", y = "Relative abundance", fill = "Guild") +
  scale_fill_manual(values = Guild_cols$color) +
  scale_y_continuous(expand = c(max(tallest_bar$sum)/100, max(tallest_bar$sum)/100)) + 
  facet_nested(~ Depth + Treatment, space = "free", scales = "free_x",
               labeller = as_labeller(facet_names)) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 90, hjust = 1, vjust = 0.5),
        strip.text = element_text(size = 5),
        strip.background = element_rect(size = 0.2),
        axis.line.y = element_blank())
dev.off()

#### 6. Indicators ####
#### _Simper ####
nc_sim <- simper(t(nc$data_loaded), 
                 nc$map_loaded$TrtDepth)
nc_s <- summary(nc_sim)
nc_df1 <- head(nc_s$`5ppt ASW0.025_DI_ctrl0.025`, n = 20) %>%
  mutate(Comparison = "ASW/Control 2.5 cm",
         ASV = rownames(.))
nc_df2 <- head(nc_s$`5ppt ASW0.125_DI_ctrl0.125`, n = 20) %>%
  mutate(Comparison = "ASW/Control 12.5 cm",
         ASV = rownames(.))
nc_simper_results <- rbind(nc_df1, nc_df2) %>%
  left_join(., nc$taxonomy_loaded, by = c("ASV" = "taxonomy8")) %>%
  mutate(SaltResponse = ifelse(ava > avb, "Postive", "Negative")) %>%
  rename("OTU" = "ASV") %>%
  mutate(OTU = gsub("ASV", "OTU", OTU)) %>%
  rename(c("Domain" = "taxonomy1",
           "Phylum" = "taxonomy2",
           "Class" = "taxonomy3",
           "Order" = "taxonomy4",
           "Family" = "taxonomy5",
           "Genus" = "taxonomy6",
           "Species" = "taxonomy7",
           "MeanSalt" = "ava",
           "MeanControl" = "avb",
           "CumulativeContribution" = "cumsum")) %>%
  dplyr::select(Comparison, SaltResponse, Domain, Phylum, Class, Order, Family, Genus,
         Species, OTU, MeanSalt, MeanControl, CumulativeContribution)
write_xlsx(nc_simper_results, 
           "simper_results_nc.xlsx",
           format_headers = F)


#### _Multipatt ####
set.seed(1202)
nc_mp <- multipatt(t(nc$data_loaded), 
                   nc$map_loaded$Treatment, 
                   func = "r.g", 
                   control = how(nperm=999))
multipatt_results <- nc_mp$sign %>%
  mutate(num_sites = rowSums(cbind(`s.5ppt ASW`, `s.DI_ctrl`, `s.Field Reference`,
                                   `s.SO4 amended`, `s.SW_noSO4`)),
         q.value = qvalue(nc_mp$sign$p.value)$qvalues,
         Group = "NA",
         ASV = rownames(nc_mp$sign)) %>%
  filter(q.value < 0.05) %>%
  filter(stat >= 0.65) %>%
  filter(num_sites <= 2) %>%
  left_join(., nc$taxonomy_loaded, by = c("ASV" = "taxonomy8"))
for (i in 1:nrow(multipatt_results)) {
  if (multipatt_results$num_sites[i] == 1 & multipatt_results$`s.5ppt ASW`[i] == 1) {
    multipatt_results$Group[i] <- "+ASW"
  }
}
for (i in 1:nrow(multipatt_results)) {
  if (multipatt_results$num_sites[i] == 1 & multipatt_results$s.DI_ctrl[i] == 1) {
    multipatt_results$Group[i] <- "Control"
  }
}
for (i in 1:nrow(multipatt_results)) {
  if (multipatt_results$num_sites[i] == 1 & multipatt_results$`s.Field Reference`[i] == 1) {
    multipatt_results$Group[i] <- "Field"
  }
}
for (i in 1:nrow(multipatt_results)) {
  if (multipatt_results$num_sites[i] == 1 & multipatt_results$`s.SO4 amended`[i] == 1) {
    multipatt_results$Group[i] <- "+SO4"
  }
}
for (i in 1:nrow(multipatt_results)) {
  if (multipatt_results$num_sites[i] == 1 & multipatt_results$s.SW_noSO4[i] == 1) {
    multipatt_results$Group[i] <- "+ASW (no SO4)"
  }
}
for (i in 1:nrow(multipatt_results)) {
  if (multipatt_results$num_sites[i] == 2 & multipatt_results$`s.5ppt ASW`[i] == 1 &
      multipatt_results$s.SW_noSO4[i] == 1) {
    multipatt_results$Group[i] <- "+ASW, +ASW (no SO4)"
  }
}
table(multipatt_results$Group)
nc_asv <- summarize_taxonomy(nc, level = 8, report_higher_tax = F)
nc_asv_all <- data.frame("RelAbundance" = round(rowMeans(nc_asv) * 100, digits = 4)) %>%
  mutate(ASV = rownames(.))
nc_mp_corrs <- as.data.frame(nc_mp$str) %>%
  dplyr::select(1:5, 9) %>%
  mutate(ASV = rownames(.)) %>%
  filter(ASV %in% multipatt_results$ASV) %>%
  set_names(c("+ASW", "Control", "Field", "+SO4", "+ASW (no SO4)", "+ASW or +ASW (no SO4)", "ASV"))
# Add corrs and taxonomy
multipatt_results <- multipatt_results %>%
  filter(Group == "+ASW (no SO4)" | Group == "+ASW, +ASW (no SO4)" | Group == "+ASW" |
           Group == "+SO4") %>%
  left_join(., nc_asv_all, by = "ASV") %>%
  left_join(., nc_mp_corrs, by = "ASV") %>%
  dplyr::select(taxonomy2, taxonomy6, ASV, Group, "+ASW", "Control", "Field", "+SO4", "+ASW (no SO4)",
                "+ASW or +ASW (no SO4)", "RelAbundance") %>%
  mutate(ASV = gsub("ASV", "OTU", ASV)) %>%
  arrange(Group, taxonomy2, taxonomy6) %>%
  mutate(ASV = factor(ASV, levels = ASV)) %>%
  set_names(c("Phylum", "Genus", "ASV", "Group", "+ASW", "Control", "Field", "+SO4", "+ASW (no SO4)",
              "+ASW/+ASW (no SO4)", "RelAbundance")) %>%
  mutate(taxon = paste(Phylum, Genus, ASV, sep = "; ")) %>%
  mutate(taxon = factor(taxon, levels = taxon))
# Plot heatmap with abundance barplot
hm.melted <- multipatt_results %>%
  dplyr::select(taxon, Field, Control, "+SO4", "+ASW (no SO4)", "+ASW", "+ASW/+ASW (no SO4)") %>%
  melt(., id.vars = c("taxon"))
hm <- ggplot(data = hm.melted, 
             aes(x = factor(taxon), y = variable, fill = value)) + 
  geom_tile() + 
  scale_fill_distiller(name = "Indicator correlation index", palette = "RdBu", 
                       direction = -1, na.value = "transparent", type = "div", limits = c(-0.8, 0.8)) +
  scale_x_discrete(breaks = unique(hm.melted$taxon), labels = unique(hm.melted$taxon),
                   limits = rev(levels(hm.melted$taxon))) +
  scale_y_discrete(position = "right") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        legend.title = element_text(size = 10), legend.key.height = unit(0.5,"cm"),
        legend.key.width = unit(1, "cm"), legend.text = element_text(size = 7)) +
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5)) +
  coord_flip()
l <- get_legend(hm)
hm.clean <- hm +
  theme(axis.title.y = element_blank(), axis.text.y = element_text(margin = margin(c(0,-5,0,0)),
                                                                   size = 7),
        axis.ticks.y = element_blank(), axis.ticks.x = element_blank(), 
        panel.grid.major = element_blank(), legend.position="none",
        axis.text.x.top = element_text(size = 5, margin = margin(c(0,0,-3,0))), 
        plot.margin = margin(c(0,-2,0,5))) +
  labs(y = "Group")
bp.y <- ggplot(data = multipatt_results, aes(x = taxon, y = RelAbundance)) + 
  geom_bar(stat = "identity", fill = "grey") + 
  scale_y_continuous(position = "right") +
  scale_x_discrete(limits = rev(levels(multipatt_results$taxon))) +
  coord_flip() + 
  theme_minimal() +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position="none", axis.text.x.top = element_text(size = 6, margin = margin(c(0,0,-3,0))), 
        plot.margin = margin(c(0,-2,0,-5))) +
  labs(y = "Rel. abund. (%)")
pdf("InitialFigs/NC_Multipatt.pdf", width = 8, height = 5)
plot_grid(hm.clean, bp.y, l, NULL, nrow = 2, ncol = 2, 
          rel_widths = c(10,2), rel_heights = c(12, 2), align = "hv", axis = "b")
dev.off()

# Do at genus level
nc_tax <- nc$taxonomy_loaded %>%
  group_by(taxonomy6) %>%
  slice_head(n = 1) %>%
  dplyr::select(-taxonomy7, -taxonomy8)
nc_genus <- summarize_taxonomy(nc, level = 6, relative = F, report_higher_tax = F)
set.seed(1202)
nc_mp <- multipatt(t(nc_genus), 
                   nc$map_loaded$Treatment, 
                   func = "r.g", 
                   control = how(nperm=999))
nc_mp_results <- nc_mp$sign %>%
  mutate(num_sites = rowSums(cbind(`s.5ppt ASW`, `s.DI_ctrl`, `s.Field Reference`,
                                   `s.SO4 amended`, `s.SW_noSO4`)),
         q.value = qvalue(nc_mp$sign$p.value)$qvalues,
         Group = "NA",
         Genus = rownames(nc_mp$sign)) %>%
  filter(q.value < 0.05) %>%
  filter(stat >= 0.5) %>%
  filter(num_sites <= 2) %>%
  left_join(., nc_tax, by = c("Genus" = "taxonomy6"))
for (i in 1:nrow(nc_mp_results)) {
  if (nc_mp_results$num_sites[i] == 1 & nc_mp_results$`s.5ppt ASW`[i] == 1) {
    nc_mp_results$Group[i] <- "+ASW"
  }
}
for (i in 1:nrow(nc_mp_results)) {
  if (nc_mp_results$num_sites[i] == 1 & nc_mp_results$s.DI_ctrl[i] == 1) {
    nc_mp_results$Group[i] <- "Control"
  }
}
for (i in 1:nrow(nc_mp_results)) {
  if (nc_mp_results$num_sites[i] == 1 & nc_mp_results$`s.Field Reference`[i] == 1) {
    nc_mp_results$Group[i] <- "Field"
  }
}
for (i in 1:nrow(nc_mp_results)) {
  if (nc_mp_results$num_sites[i] == 1 & nc_mp_results$`s.SO4 amended`[i] == 1) {
    nc_mp_results$Group[i] <- "+SO4"
  }
}
for (i in 1:nrow(nc_mp_results)) {
  if (nc_mp_results$num_sites[i] == 1 & nc_mp_results$s.SW_noSO4[i] == 1) {
    nc_mp_results$Group[i] <- "+ASW (no SO4)"
  }
}
for (i in 1:nrow(nc_mp_results)) {
  if (nc_mp_results$num_sites[i] == 2 & nc_mp_results$`s.5ppt ASW`[i] == 1 &
      nc_mp_results$s.SW_noSO4[i] == 1) {
    nc_mp_results$Group[i] <- "+ASW, +ASW (no SO4)"
  }
}
table(nc_mp_results$Group)
nc_genus_all <- data.frame("RelAbundance" = round(rowMeans(nc_genus)/min(colSums(nc$data_loaded)) * 100, digits = 4)) %>%
  mutate(Genus = rownames(.))
nc_mp_corrs <- as.data.frame(nc_mp$str) %>%
  dplyr::select(1:5, 9) %>%
  mutate(Genus = rownames(.)) %>%
  filter(Genus %in% nc_mp_results$Genus) %>%
  set_names(c("+ASW", "Control", "Field", "+SO4", "+ASW (no SO4)", "+ASW or +ASW (no SO4)", "Genus"))
# Add corrs and taxonomy
nc_mp_results <- nc_mp_results %>%
  filter(Group == "+ASW (no SO4)" | Group == "+ASW, +ASW (no SO4)" | Group == "+ASW" |
           Group == "+SO4") %>%
  left_join(., nc_genus_all, by = "Genus") %>%
  left_join(., nc_mp_corrs, by = "Genus") %>%
  dplyr::select(taxonomy2, Genus, Group, "+ASW", "Control", "Field", "+SO4", "+ASW (no SO4)",
                "+ASW or +ASW (no SO4)", "RelAbundance") %>%
  arrange(Group, taxonomy2) %>%
  mutate(Genus = factor(Genus, levels = Genus)) %>%
  set_names(c("Phylum", "Genus", "Group", "+ASW", "Control", "Field", "+SO4", "+ASW (no SO4)",
              "+ASW/+ASW (no SO4)", "RelAbundance")) %>%
  mutate(taxon = paste(Phylum, Genus, sep = "; ")) %>%
  mutate(taxon = factor(taxon, levels = taxon))
# Plot heatmap with abundance barplot
nc.hm.melted <- nc_mp_results %>%
  dplyr::select(taxon, Field, Control, "+SO4", "+ASW (no SO4)", "+ASW", "+ASW/+ASW (no SO4)") %>%
  melt(., id.vars = c("taxon"))
nc.hm <- ggplot(data = nc.hm.melted, 
                aes(x = factor(taxon), y = variable, fill = value)) + 
  geom_tile() + 
  scale_fill_distiller(name = "Indicator correlation index", palette = "RdBu", 
                       direction = -1, na.value = "transparent", type = "div", limits = c(-0.8, 0.8)) +
  scale_x_discrete(breaks = unique(nc.hm.melted$taxon), labels = unique(nc.hm.melted$taxon),
                   limits = rev(levels(nc.hm.melted$taxon))) +
  scale_y_discrete(position = "right") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        legend.title = element_text(size = 10), legend.key.height = unit(0.5,"cm"),
        legend.key.width = unit(1, "cm"), legend.text = element_text(size = 7)) +
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5)) +
  coord_flip()
nc.l <- get_legend(nc.hm)
nc.hm.clean <- nc.hm +
  theme(axis.title.y = element_blank(), axis.text.y = element_text(margin = margin(c(0,-5,0,0)),
                                                                   size = 6),
        axis.ticks.y = element_blank(), axis.ticks.x = element_blank(), 
        panel.grid.major = element_blank(), legend.position="none",
        axis.text.x.top = element_text(size = 5, margin = margin(c(0,0,-2,0))), 
        plot.margin = margin(c(0,-2,0,5))) +
  labs(y = "Group")
nc.bp.y <- ggplot(data = nc_mp_results, aes(x = taxon, y = RelAbundance)) + 
  geom_bar(stat = "identity", fill = "grey") + 
  scale_y_continuous(position = "right") +
  scale_x_discrete(limits = rev(levels(nc_mp_results$taxon))) +
  coord_flip() + 
  theme_minimal() +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position="none", plot.margin = margin(c(0,5,0,-5)),
        axis.text.x.top = element_text(size = 6, margin = margin(c(0,0,-2,0)))) +
  labs(y = "Rel. abund. (%)")
pdf("InitialFigs/nc_Multipatt_genus.pdf", width = 8, height = 5)
plot_grid(nc.hm.clean, nc.bp.y, nc.l, NULL, nrow = 2, ncol = 2, 
          rel_widths = c(10,2), rel_heights = c(12, 2), align = "hv", axis = "b")
dev.off()

#### 7. BGC ####
# Add biogeochem. analysis
# Get variables, make corrplot, envfit, ordination with vectors, ordistep, taxa correlations
# Note only BGC for D2 samples in SC

# Get variables
env <- nc$map_loaded %>%
  dplyr::select(CH4_ug_m2_h, CO2_ug_m2_h, N2O_ug_m2_h, Salinity_ppt_all, 
                TOC_mgL, TN_mgL, NH4_mgL, PO4_mgL, Cl_mgL, SO4_mgL, Br_mgL, NO3_mgL,
                DIN_mgL, DON_mgL, pH)
env_nona <- na.omit(env)

# Corrplot
C <- cor(env_nona)
corrplot(C, method = "number", type = "lower", 
         tl.cex = 0.5, tl.col = "black", number.cex = 0.5)
pdf("InitialFigs/NC_BGC_CH4.pdf", width = 7, height = 5)
meth_corr_by_bgc(env_nona = env_nona)
dev.off()

# Envfit
pcoa <- cmdscale(nc_bc, k = nrow(nc$map_loaded) - 1, eig = T)
set.seed(100)
ef <- envfit(pcoa, env, permutations = 999, na.rm = TRUE)
ef
arrow_factor <- ordiArrowMul(ef)
manual_factor <- 0.4
vec.df <- as.data.frame(ef$vectors$arrows*sqrt(ef$vectors$r)) %>%
  mutate(Dim1 = Dim1 * manual_factor,
         Dim2 = Dim2 * manual_factor) %>%
  mutate(variables = rownames(.)) %>%
  filter(ef$vectors$pvals < 0.05) %>%
  filter(variables != "Cl_mgL") %>%
  mutate(shortnames = c("CH4", "CO2", "N2O", "Salinity", "NH4", "Br"))

# Plot with significant vectors
ordiplot(pcoa)
plot(ef, p.max = 0.05, cex = 0.5)
pcoaA1 <- round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[1]*100, digits = 1)
pcoaA2 <- round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[2]*100, digits = 1)
nc$map_loaded$Axis01 <- scores(pcoa)[,1]
nc$map_loaded$Axis02 <- scores(pcoa)[,2]
micro.hulls <- ddply(nc$map_loaded, c("Treatment"), find_hull)
g <- ggplot(nc$map_loaded, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = Treatment, fill = Treatment),
               alpha = 0.1, show.legend = F) +
  geom_point(size = 3, alpha = 1, aes(colour = Treatment, shape = Depth)) +
  labs(x = paste("PC1: ", pcoaA1, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2, "%", sep = ""),
       shape = "Depth (m)") +
  scale_colour_viridis_d(labels = c("Field", "Control", "+SO4", "+ASW (no SO4)", "+ASW"),
                         breaks = c("Field Reference", "DI_ctrl", "SO4 amended",
                                    "SW_noSO4", "5ppt ASW")) +
  scale_fill_viridis_d() +
  guides(colour = guide_legend(override.aes = list(shape = 15),
                               order = 1)) +
  coord_fixed() +
  theme_bw() +  
  theme(legend.position = c(0,0),
        legend.justification = c(0,0),
        legend.background = element_blank(),
        legend.key.size = unit(0.5, "cm"),
        legend.margin = margin(2,2,2,2, "pt"),
        legend.spacing.y = unit(0, "cm"),
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        plot.margin = margin(5, 5, 5, 5, "pt")) +
  geom_segment(data = vec.df,
               aes(x = 0, xend = Dim1, y = 0, yend = Dim2),
               arrow = arrow(length = unit(0.5, "cm")),
               colour = "gray", alpha = 0.5,
               inherit.aes = FALSE) + 
  geom_text(data = subset(vec.df, shortnames != "Salinity"),
            aes(x = Dim1, y = Dim2, label = shortnames),
            size = 3, color = "black") +
  geom_text(data = subset(vec.df, shortnames == "Salinity"),
            aes(x = Dim1 - 0.03, y = Dim2, label = shortnames),
            size = 3, color = "black")
g
pdf("InitialFigs/NC_PCoA_wBGC.pdf", width = 6, height = 4)
g
dev.off()

# Ordistep
comm_nona <- as.data.frame(t(nc$data_loaded)) %>%
  filter(rownames(.) %in% rownames(env_nona))
mod0 <- rda(comm_nona ~ 1, env_nona)  # Model with intercept only
mod1 <- rda(comm_nona ~ ., env_nona)  # Model with all explanatory variables
mod <- ordistep(mod0, scope = formula(mod1))
mod
mod$anova # TN, N2O

# Plot correlations for different taxonomic levels at a given % relative abundance threshold
pdf("InitialFigs/NC_Phyla_CH4.pdf", width = 7, height = 5)
meth_corr_by_taxonomy(input = nc, level = 2, threshold = 0.5)
dev.off()
meth_corr_by_taxonomy(input = nc, level = 3, threshold = 0.5)
meth_corr_by_taxonomy(input = nc, level = 4, threshold = 0.5)
meth_corr_by_taxonomy(input = nc, level = 5, threshold = 0.5)
meth_corr_by_taxonomy(input = nc, level = 6, threshold = 0.5)
meth_corr_by_taxonomy(input = nc, level = 8, threshold = 0.5)
pdf("InitialFigs/NC_Guilds_CH4.pdf", width = 7, height = 5)
meth_corr_by_taxonomy(input = nc, level = 9, threshold = 0)
dev.off()



#### End Script ####