# Timberlake 16S rRNA marker gene iTag data analysis
# PCR amplicons. For metagenomic derived 16S, see Timberlake16S_mTAG.R
# by Cliff Bueno de Mesquita, Tringe Lab, JGI, Spring 2023



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
## 7. Correlations with BGC
## 8. NCBI Submission



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
library(pheatmap) # heatmaps
library(zCompositions) # CLR
library(compositions) # Aitchison
library(mobr) # rarefaction curves
library(plotly) # interactive graphs
library(pairwiseAdonis) # pairwise permanova
library(patchwork) # insets
library(janitor) # data cleaning
library(ggrepel) # text

# Functions
find_hull <- function(df) df[chull(df$Axis01, df$Axis02),]
`%notin%` <- Negate(`%in%`)

# Guild subsetting module from other repo
source("~/Documents/GitHub/SF_microbe_methane/modules/3_OTU_subsetting_modules_v.0.4_strip.r")

# Correlation functions from other repo
source("~/Documents/GitHub/EastCoast/code/meth_corr_by_taxonomy.R")
source("~/Documents/GitHub/EastCoast/code/meth_corr_by_bgc.R")

# Plotting functions from other repo
source("~/Documents/GitHub/EastCoast/code/cliffplot_taxa_bars.R")
source("~/Documents/GitHub/Extremophilic_Fungi/code/plot_multipatt.R")

# Effect size from Jack Darcy
source("~/Documents/GitHub/Timberlake/code/effectSize.R")

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
  filter(Guild != "Anamx") %>%
  arrange(Index)



#### _Prepare data ####
# Only need to do once. Then skip to start here.
# Read in combined dataset (with other samples)
input_filt <- readRDS("data/input_filt_comb_wBGC.rds")
# Get only Timberlake samples
input_filt_nc <- filter_data(input_filt,
                             filter_cat = "Estuary",
                             keep_vals = "Alligator")
# Remove extra samples
# Keep 5A and 5B - careful about assigning env data though, or don't assign env data
input_filt_nc <- filter_data(input_filt_nc,
                             filter_cat = "sampleID",
                             filter_vals = c("TL_nw_d1_DI_ctrl_AF1",
                                             "TL_nw_d1_DI_ctrl_AF3",
                                             "TL_nw_d1_DI_ctrl_AF4",
                                             "TL_nw_d1_ASW_noS_BF3",
                                             "TL_nw_d1_ASW_noS_BF4",
                                             "TL_nw_d1_ASW_noS_BF5"#,
                                             #"TL_inc_d1_SO4_5A",
                                             #"TL_inc_d1_SO4_5B",
                                             #"TL_inc_d2_SO4_5A",
                                             #"TL_inc_d2_SO4_5B"
                                             ))
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

# Remove singletons and doubletons
singdoub <- data.frame("count" = rowSums(input_filt_nc$data_loaded)) %>%
  filter(count < 3) %>%
  mutate(ASV = rownames(.))

input_filt_nc <- filter_taxa_from_input(input_filt_nc,
                                        taxa_IDs_to_remove = singdoub$ASV)

# Rarefaction curves
rarecurve(t(input_filt_nc$data_loaded), step = 5000)
rarecurve(t(input_filt_nc$data_loaded), step = 5000, sample = 82312)
inv_mob_in <- make_mob_in(t(input_filt_nc$data_loaded), 
                         input_filt_nc$map_loaded)
plot_rarefaction(inv_mob_in, 'Treatment', 'Field', 'SBR',
                 pooled = TRUE, leg_loc = 'bottomright')
plot_rarefaction(inv_mob_in, 'Treatment', 'Field', 'IBR',
                 pooled = TRUE, leg_loc = 'bottomright') # Not enough memory

# Rarefy at 82316
sort(colSums(input_filt_nc$data_loaded))
mean(colSums(input_filt_nc$data_loaded)) # 115744.3
se(colSums(input_filt_nc$data_loaded)) # 2279.646
set.seed(530)
nc <- single_rarefy(input_filt_nc, 82316) # n = 43 still
sort(colSums(nc$data_loaded))

# Add MG:MT ratio and AO:NOB ratio to metadata
nc_guilds <- summarize_taxonomy(input = nc, level = 9, report_higher_tax = F) %>%
  t() %>%
  as.data.frame() %>%
  mutate(sampleID = rownames(.)) %>%
  mutate(AO_NOB = (AOA + AOB) / NOB) %>%
  mutate(MG = CH4_ac + CH4_H2 + CH4_me + CH4_mix) %>%
  mutate(MT = MOB_I + MOB_II + MOB_IIa + ANME) %>%
  mutate(MG_MT = (CH4_ac + CH4_H2 + CH4_me + CH4_mix)/(ANME + MOB_I + MOB_II + MOB_IIa)) %>%
  left_join(., nc$map_loaded, by = "sampleID")
nc$map_loaded$MG_MT <- nc_guilds$MG_MT
nc$map_loaded$AO_NOB <- nc_guilds$AO_NOB

# Clean up sampleID
nc$map_loaded$sampleID_clean <- c("Field_10_D1", "Field_19_D1", "Field_36_D1", "Field_43_D1", "Field_46_D1",
                            "Field_10_D2", "Field_19_D2", "Field_36_D2", "Field_43_D2", "Field_46_D2",
                            "+ASW_3_D1", "+ASW_4_D1", "+ASW_5_D1",
                            "+ASW-SO4_1_D1","+ASW-SO4_2_D1","+ASW-SO4_3_D1","+ASW-SO4_4_D1","+ASW-SO4_5_D1",
                            "Control_2_D1", "Control_3_D1", "Control_4_D1", "Control_5_D1",
                            "+SO4_3_D1", "+SO4_4_D1", "+SO4_5A_D1", "+SO4_5B_D1",
                            "+ASW_1_D2", "+ASW_3_D2", "+ASW_4_D2", "+ASW_5_D2",
                            "+ASW-SO4_2_D2","+ASW-SO4_3_D2","+ASW-SO4_4_D2","+ASW-SO4_5_D2",
                            "Control_1_D2", "Control_3_D2", "Control_4_D2", "Control_5_D2",
                            "+SO4_2_D2", "+SO4_3_D2", "+SO4_4_D2", "+SO4_5A_D2", "+SO4_5B_D2")

# Update richness and Shannon
nc$map_loaded$rich <- specnumber(nc$data_loaded, MARGIN = 2)
nc$map_loaded$shannon <- diversity(nc$data_loaded, index = "shannon", MARGIN = 2)

# Clean up column order
nc$map_loaded <- nc$map_loaded %>%
  dplyr::select(sampleID, sampleID_clean, Treatment, Depth, TrtDepth, rich, 
                shannon, MG_MT, AO_NOB, everything())

# Save
#saveRDS(nc, "data/nc.rds")



#### _Start here ####
nc <- readRDS("data/nc.rds")

# Add C, N, C:N data (was received later)
nc_cn <- read_excel("~/Documents/GitHub/EastCoast/data/Copy of CHN data.xls", sheet = 2) %>%
  slice(12:91) %>%
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

# After looking at microbial alpha and beta diversity it is clear SO4 5B D1 and D2 are switched!
# Fix
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

# Add +SO4 core 5 data to 5A and 5B
new_dat <- read_excel("data/Soil pH-Jessie.xls", sheet = 3, na = "NA")
# rows 25, 26, 42, 43
# columns 10 - 24
nc$map_loaded[25, 10:24] <- new_dat[1, 2:16]
nc$map_loaded[26, 10:24] <- new_dat[2, 2:16]
nc$map_loaded[42, 10:24] <- new_dat[3, 2:16]
nc$map_loaded[43, 10:24] <- new_dat[4, 2:16]
# To recap, C, N, CN, pH have 5A D1 D2, 5B D1 D2. Porewater just has 5 D1 D2. GHG just 5

# Add two new columns for SO4 and ASW (following N + P type experiments)
# There's SO4, ASW, and both (just like N, P, NP)
nc$map_loaded$SO4 <- recode_factor(nc$map_loaded$Treatment,
                                   "+SO4" = "Sulfate",
                                   "+ASW" = "SaltSulfate")
nc$map_loaded$SO4 <- grepl(x = nc$map_loaded$SO4, pattern = "Sulfate")
nc$map_loaded$ASW <- recode_factor(nc$map_loaded$Treatment,
                                   "+ASW-SO4" = "Salt",
                                   "+ASW" = "SaltSulfate")
nc$map_loaded$ASW <- grepl(x = nc$map_loaded$ASW, pattern = "Salt")

# Remove Methylotenera from MOB guild
for (i in 1:nrow(nc$taxonomy_loaded)) {
  if (nc$taxonomy_loaded$taxonomy6[i] == "Methylotenera") {
    nc$taxonomy_loaded$taxonomy9[i] <- "NA"
  }
}



#### 3. Alpha ####
# Get descriptive info
field <- filter_data(nc, "Treatment", keep_vals = "Field")
min(field$map_loaded$rich)
max(field$map_loaded$rich)
lab <- filter_data(nc, "Treatment", filter_vals = "Field")
min(lab$map_loaded$rich)
max(lab$map_loaded$rich)
nc$map_loaded$Treatment2 <- recode_factor(nc$map_loaded$Treatment,
                                          "Control" = "Lab",
                                          "+SO4" = "Lab",
                                          "+ASW-SO4" = "Lab",
                                          "+ASW" = "Lab")
source("~/Documents/GitHub/EastCoast/code/plot_venn_diagram2.R")
png("InitialFigs/Venn.png", width = 7, height = 5, units = "in", res = 300)
plot_venn_diagram2(nc, "Treatment2", pres_thresh = 0.0000000000000000000000000001)
dev.off()

# Include in manuscript as supplemental figure
png("FinalFigs/FigureS2.png", width = 7, height = 5, units = "in", res = 300)
plot_venn_diagram2(nc, "Treatment2", pres_thresh = 0.0000000000000000000000000001)
dev.off()

# Test and plot
leveneTest(nc$map_loaded$rich ~ nc$map_loaded$Treatment) # Homogeneous
m <- aov(rich ~ Treatment * Depth, data = nc$map_loaded)
Anova(m, type = "III") # No interaction
m <- aov(rich ~ Treatment + Depth, data = nc$map_loaded)
Anova(m, type = "II") # Treatment and Depth
m <- aov(rich ~ Treatment, data = nc$map_loaded)
shapiro.test(m$residuals) # Almost normal
hist(m$residuals)
summary(m)
t <- emmeans(object = m, specs = "Treatment") %>%
  cld(object = ., adjust = "Tukey", Letters = letters, alpha = 0.05) %>%
  mutate(name = "rich",
         y = max(nc$map_loaded$rich)+(max(nc$map_loaded$rich)-min(nc$map_loaded$rich))/20)

leveneTest(nc$map_loaded$shannon ~ nc$map_loaded$Treatment) # Homogeneous
m1 <- aov(shannon ~ Treatment * Depth, data = nc$map_loaded)
Anova(m1, type = "III") # Minor interaction
m1 <- aov(shannon ~ Treatment + Depth, data = nc$map_loaded)
Anova(m1, type = "II") # Treatment and Depth
m1 <- aov(shannon ~ Treatment, data = nc$map_loaded)
shapiro.test(m1$residuals) # Normal
hist(m1$residuals)
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
g1 <- ggplot(alpha_long, aes(reorder(Treatment, value, mean), value, 
                             colour = Treatment)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 2, alpha = 0.75, width = 0.2, aes(shape = Depth)) +
  geom_text(data = label_df, aes(Treatment, y, label = str_trim(.group)), 
            size = 4, color = "black") +
  labs(x = "Site", y = "Number of OTUs", shape = "Depth (cm)") +
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
png("InitialFigs/Alpha.png", width = 6, height = 3, units = "in", res = 300)
g1
dev.off()

# Test just lab with 3-way anova
# Make pie chart with effect size
m1 <- aov(rich ~  Depth + SO4 + ASW, data = lab$map_loaded)
summary(m1)
anova(m1)
Anova(m1, type = "II") # Depth and ASW
eta_sq_m1 <- eta_sq(m1) %>%
  as.data.frame() %>%
  rownames_to_column(var = "group") %>%
  set_names(c("group", "value")) %>%
  mutate(group = c("Depth***", "SO4", "ASW***", "Residuals")) %>%
  arrange(desc(group)) %>%
  mutate(ypos = cumsum(value) - 0.5*value)
g2 <- ggplot(eta_sq_m1, aes(x = "", y = value, fill = group)) +
  geom_bar(stat = "identity", width = 1, color = NA) +
  coord_polar("y", start = 0) +
  geom_text(data = subset(eta_sq_m1, group == "SO4"),
            aes(y = ypos, label = group), color = "white", size = 2, check_overlap = F,
            nudge_x = 0.1) +
  geom_text(data = subset(eta_sq_m1, group != "SO4"),
            aes(y = ypos, label = group), color = "white", size = 2, check_overlap = F) +
  scale_fill_manual(values = c("#5DC863FF", "grey30", "grey70", "#21908CFF")) +
  theme_void() + 
  theme(legend.position = "none")

m2 <- aov(shannon ~  Depth + SO4 + ASW, data = lab$map_loaded)
summary(m2)
anova(m2)
Anova(m2, type = "II") # Depth and ASW
eta_sq_m2 <- eta_sq(m2) %>%
  as.data.frame() %>%
  rownames_to_column(var = "group") %>%
  set_names(c("group", "value")) %>%
  mutate(group = c("Depth**", "SO4", "ASW**", "Residuals")) %>%
  arrange(desc(group)) %>%
  mutate(ypos = cumsum(value) - 0.5*value)
g3 <- ggplot(eta_sq_m2, aes(x = "", y = value, fill = group)) +
  geom_bar(stat = "identity", width = 1, color = NA) +
  coord_polar("y", start = 0) +
  geom_text(data = subset(eta_sq_m2, group == "SO4"),
            aes(y = ypos, label = group), color = "white", size = 2, check_overlap = F,
            nudge_x = 0.1) +
  geom_text(data = subset(eta_sq_m2, group != "SO4"),
            aes(y = ypos, label = group), color = "white", size = 2, check_overlap = F) +
  scale_fill_manual(values = c("#5DC863FF", "grey30", "grey70", "#21908CFF")) +
  theme_void() + 
  theme(legend.position = "none")

# Add as inset, save as Figure 2.
plot.with.inset <-
  ggdraw() +
  draw_plot(g1) +
  draw_plot(g2, x = 0, y = 0.56, width = 0.31, height = 0.31) +
  draw_plot(g3, x = 0.47, y = 0.56, width = 0.31, height = 0.31)
plot.with.inset
png("FinalFigs/Figure2.png", width = 7, height = 4, units = "in", res = 300)
plot.with.inset
dev.off()



#### 4. Beta ####
# Get env. vars with no NA
names(nc$map_loaded)
env_nc <- nc$map_loaded %>%
  dplyr::select(10:27)
env_nona_nc <- na.omit(env_nc)
nrow(env_nona_nc) # n = 25



#### _Aitch ####
# Get env. vars with no NA
env_nc2 <- input_filt_nc$map_loaded %>%
  dplyr::select(10:27)
env_nona_nc2 <- na.omit(env_nc2)
nrow(env_nona_nc2) # n = 25

# Use non-rarefied data, do CLR transformation, Aitchison, PCA
sum(rownames(nc$map_loaded) != rownames(input_filt_nc$map_loaded))
input_filt_nc$map_loaded <- nc$map_loaded
dim(input_filt_nc$data_loaded)

# CLT transformation
otu_czm <- cmultRepl(t((input_filt_nc$data_loaded)), label = 0, method = "CZM")
otu_clr <- clr(otu_czm)
aclr <- compositions::dist(otu_clr)

set.seed(1150)
adonis2(aclr ~ input_filt_nc$map_loaded$Treatment*input_filt_nc$map_loaded$Depth) # Both sig
adonis2(aclr ~ nc$map_loaded$Depth*nc$map_loaded$Treatment) # No effect of order
pairwiseAdonis::pairwise.adonis(aclr, input_filt_nc$map_loaded$Treatment)
anova(betadisper(aclr, nc$map_loaded$Treatment)) # Dispersion homogeneous
anova(betadisper(aclr, nc$map_loaded$Depth)) # Dispersion homogeneous

# Stats and pie for just lab
lab_raw <- filter_data(input_filt_nc,
                       filter_cat = "Treatment", 
                       filter_vals = "Field")
otu_czm_lab <- cmultRepl(t((lab_raw$data_loaded)), label = 0, method = "CZM")
otu_clr_lab <- clr(otu_czm_lab)
aclr_lab <- compositions::dist(otu_clr_lab)
set.seed(1150)
m3 <- adonis2(aclr_lab ~ lab_raw$map_loaded$Depth + lab_raw$map_loaded$SO4 + lab_raw$map_loaded$ASW)
eta_sq_m3 <- eta_sq_adonis(m3) %>%
  as.data.frame() %>%
  rownames_to_column(var = "group") %>%
  filter(group != "Total") %>%
  set_names(c("group", "value")) %>%
  mutate(group = c("Depth***", "SO4", "ASW***", "Residuals")) %>%
  arrange(desc(group)) %>%
  mutate(ypos = cumsum(value) - 0.5*value)
g4 <- ggplot(eta_sq_m3, aes(x = "", y = value, fill = group)) +
  geom_bar(stat = "identity", width = 1, color = NA) +
  coord_polar("y", start = 0) +
  geom_text(data = subset(eta_sq_m3, group == "SO4"),
            aes(y = ypos, label = group), color = "white", size = 2.5, check_overlap = F,
            nudge_x = 0.1) +
  geom_text(data = subset(eta_sq_m3, group != "SO4"),
            aes(y = ypos, label = group), color = "white", size = 2.5, check_overlap = F) +
  scale_fill_manual(values = c("#5DC863FF", "grey30", "grey70", "#21908CFF")) +
  theme_void() + 
  theme(legend.position = "none")



# PCA with vectors
d.pcx <- prcomp(aclr)
set.seed(100)
ef_nc <- envfit(d.pcx, env_nc, permutations = 999, na.rm = TRUE)
ef_nc
biplot(d.pcx)
ordiplot(d.pcx)
plot(ef_nc, p.max = 0.075, cex = 0.5)
manual_factor_nc <- 0.28
vec.df_nc <- as.data.frame(ef_nc$vectors$arrows*sqrt(ef_nc$vectors$r)) %>%
  mutate(PC1 = PC1 * manual_factor_nc,
         PC2 = PC2 * manual_factor_nc) %>%
  mutate(variables = rownames(.)) %>%
  filter(ef_nc$vectors$pvals < 0.075) %>%
  filter(variables != "Cl_mgL") %>%
  mutate(shortnames = c("Salinity", "N2O", "CO2", "NH4", "pH", "Br", "C:N"))
d.mvar <- sum(d.pcx$sdev^2)
PC1 <- paste("PC1: ", round((sum(d.pcx$sdev[1]^2)/d.mvar)*100, 1), "%")
PC2 <- paste("PC2: ", round((sum(d.pcx$sdev[2]^2)/d.mvar)*100, 1), "%")
nc$map_loaded$Axis01 <- d.pcx$rotation[,1]
nc$map_loaded$Axis02 <- d.pcx$rotation[,2]
micro.hulls <- ddply(nc$map_loaded, c("Treatment"), find_hull)
g5 <- ggplot(nc$map_loaded, aes(-Axis01, -Axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = Treatment, fill = Treatment),
               alpha = 0.1, show.legend = F) +
  geom_point(size = 3, alpha = 1, aes(colour = Treatment, shape = Depth)) +
  geom_segment(data = vec.df_nc,
               aes(x = 0, xend = PC1, y = 0, yend = PC2),
               arrow = arrow(length = unit(0.35, "cm")),
               colour = "gray", alpha = 0.6,
               inherit.aes = FALSE) + 
  geom_text(data = vec.df_nc,
            aes(x = PC1, y = PC2, label = shortnames),
            size = 3, color = "black") +
  labs(x = PC1, 
       y = PC2,
       shape = "Depth") +
  scale_colour_viridis_d() +
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
g5
png("InitialFigs/BetaAitch.png", width = 7, height = 5, units = "in", res = 300)
g5
dev.off()

# Add as inset, save as Figure 3.
plot.with.inset2 <-
  ggdraw() +
  draw_plot(g5) +
  draw_plot(g4, x = 0.55, y = 0.7, width = 0.31, height = 0.31)
plot.with.inset2
png("FinalFigs/Figure3.png", width = 7, height = 5, units = "in", res = 300)
plot.with.inset2
dev.off()

# RDA
comm_nona_nc <- as.data.frame(t(input_filt_nc$data_loaded)) %>%
  filter(rownames(.) %in% rownames(env_nona_nc2))
mod0 <- rda(comm_nona_nc ~ 1, env_nona_nc2)  # Model with intercept only
mod1 <- rda(comm_nona_nc ~ ., env_nona_nc2)  # Model with all explanatory variables
mod <- ordistep(mod0, scope = formula(mod1))
mod
mod$anova # DIN, Br, C:N



#### _Bray ####
nc_bc <- calc_dm(nc$data_loaded)
set.seed(1150)
adonis2(nc_bc ~ nc$map_loaded$Treatment + nc$map_loaded$Depth) # Both sig
adonis2(nc_bc ~ nc$map_loaded$Depth + nc$map_loaded$Treatment) # No effect of order
anova(betadisper(nc_bc, nc$map_loaded$Treatment)) # Dispersion not homogeneous
anova(betadisper(nc_bc, nc$map_loaded$Depth)) # Dispersion homogeneous

# PCoA with vectors
nc_pcoa <- cmdscale(nc_bc, k = nrow(nc$map_loaded) - 1, eig = T)
set.seed(100)
ef_nc <- envfit(nc_pcoa, env_nc, permutations = 999, na.rm = TRUE)
ef_nc
ordiplot(nc_pcoa)
plot(ef_nc, p.max = 0.05, cex = 0.5)
manual_factor_nc <- 0.45
vec.df_nc <- as.data.frame(ef_nc$vectors$arrows*sqrt(ef_nc$vectors$r)) %>%
  mutate(Dim1 = Dim1 * manual_factor_nc,
         Dim2 = Dim2 * manual_factor_nc) %>%
  mutate(variables = rownames(.)) %>%
  filter(ef_nc$vectors$pvals < 0.05) %>%
  filter(variables != "Cl_mgL") %>%
  mutate(shortnames = c("Salinity", "CH4", "N2O", "CO2", "NH4", "pH", "Br", "C:N"))
pcoaA1 <- round((eigenvals(nc_pcoa)/sum(eigenvals(nc_pcoa)))[1]*100, digits = 1)
pcoaA2 <- round((eigenvals(nc_pcoa)/sum(eigenvals(nc_pcoa)))[2]*100, digits = 1)
nc$map_loaded$Axis01 <- scores(nc_pcoa)[,1]
nc$map_loaded$Axis02 <- scores(nc_pcoa)[,2]
micro.hulls <- ddply(nc$map_loaded, c("Treatment"), find_hull)
png("InitialFigs/Beta.png", width = 7, height = 5, units = "in", res = 300)
ggplot(nc$map_loaded, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = Treatment, fill = Treatment),
               alpha = 0.1, show.legend = F) +
  geom_point(size = 3, alpha = 1, aes(colour = Treatment, shape = Depth)) +
  geom_segment(data = vec.df_nc,
               aes(x = 0, xend = Dim1, y = 0, yend = Dim2),
               arrow = arrow(length = unit(0.35, "cm")),
               colour = "gray", alpha = 0.6,
               inherit.aes = FALSE) + 
  geom_text(data = vec.df_nc,
            aes(x = Dim1, y = Dim2, label = shortnames),
            size = 3, color = "black") +
  labs(x = paste("PC1: ", pcoaA1, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2, "%", sep = ""),
       shape = "Depth") +
  scale_colour_viridis_d() +
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

# RDA
comm_nona_nc <- as.data.frame(t(nc$data_loaded)) %>%
  filter(rownames(.) %in% rownames(env_nona_nc))
mod0 <- rda(comm_nona_nc ~ 1, env_nona_nc)  # Model with intercept only
mod1 <- rda(comm_nona_nc ~ ., env_nona_nc)  # Model with all explanatory variables
mod <- ordistep(mod0, scope = formula(mod1))
mod
mod$anova # DIN, Br

# Check PCoA at different levels
tax_sum_phyla <- summarize_taxonomy(input = nc, level = 2, report_higher_tax = T)
phy_bc <- calc_dm(tax_sum_phyla)
nc_pcoa <- cmdscale(phy_bc, k = nrow(nc$map_loaded) - 1, eig = T)
pcoaA1 <- round((eigenvals(nc_pcoa)/sum(eigenvals(nc_pcoa)))[1]*100, digits = 1)
pcoaA2 <- round((eigenvals(nc_pcoa)/sum(eigenvals(nc_pcoa)))[2]*100, digits = 1)
nc$map_loaded$Axis01 <- scores(nc_pcoa)[,1]
nc$map_loaded$Axis02 <- scores(nc_pcoa)[,2]
micro.hulls <- ddply(nc$map_loaded, c("Treatment"), find_hull)
g1 <- ggplot(nc$map_loaded, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = Treatment, fill = Treatment),
               alpha = 0.1, show.legend = F) +
  geom_point(size = 3, alpha = 1, aes(colour = Treatment, shape = Depth)) +
  labs(x = paste("PC1: ", pcoaA1, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2, "%", sep = ""),
       shape = "Depth") +
  scale_colour_viridis_d() +
  scale_fill_viridis_d() +
  guides(colour = guide_legend(override.aes = list(shape = 15),
                               order = 1)) +
  theme_bw() +  
  ggtitle("Phylum") +
  theme(legend.position = "right",
        plot.title = element_text(size = 12, vjust = -1),
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        plot.margin = margin(5, 5, 5, 5, "pt"))
leg <- get_legend(g1)
g1 <- g1 + theme(legend.position = "none")

tax_sum_classes <- summarize_taxonomy(input = nc, level = 3, report_higher_tax = T)
cla_bc <- calc_dm(tax_sum_classes)
nc_pcoa <- cmdscale(cla_bc, k = nrow(nc$map_loaded) - 1, eig = T)
pcoaA1 <- round((eigenvals(nc_pcoa)/sum(eigenvals(nc_pcoa)))[1]*100, digits = 1)
pcoaA2 <- round((eigenvals(nc_pcoa)/sum(eigenvals(nc_pcoa)))[2]*100, digits = 1)
nc$map_loaded$Axis01 <- scores(nc_pcoa)[,1]
nc$map_loaded$Axis02 <- scores(nc_pcoa)[,2]
micro.hulls <- ddply(nc$map_loaded, c("Treatment"), find_hull)
g2 <- ggplot(nc$map_loaded, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = Treatment, fill = Treatment),
               alpha = 0.1, show.legend = F) +
  geom_point(size = 3, alpha = 1, aes(colour = Treatment, shape = Depth)) +
  labs(x = paste("PC1: ", pcoaA1, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2, "%", sep = ""),
       shape = "Depth") +
  scale_colour_viridis_d() +
  scale_fill_viridis_d() +
  guides(colour = guide_legend(override.aes = list(shape = 15),
                               order = 1)) +
  theme_bw() + 
  ggtitle("Class") +
  theme(legend.position = "none",
        plot.title = element_text(size = 12, vjust = -1),
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        plot.margin = margin(5, 5, 5, 5, "pt"))

tax_sum_orders <- summarize_taxonomy(input = nc, level = 4, report_higher_tax = T)
ord_bc <- calc_dm(tax_sum_orders)
nc_pcoa <- cmdscale(ord_bc, k = nrow(nc$map_loaded) - 1, eig = T)
pcoaA1 <- round((eigenvals(nc_pcoa)/sum(eigenvals(nc_pcoa)))[1]*100, digits = 1)
pcoaA2 <- round((eigenvals(nc_pcoa)/sum(eigenvals(nc_pcoa)))[2]*100, digits = 1)
nc$map_loaded$Axis01 <- scores(nc_pcoa)[,1]
nc$map_loaded$Axis02 <- scores(nc_pcoa)[,2]
micro.hulls <- ddply(nc$map_loaded, c("Treatment"), find_hull)
g3 <- ggplot(nc$map_loaded, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = Treatment, fill = Treatment),
               alpha = 0.1, show.legend = F) +
  geom_point(size = 3, alpha = 1, aes(colour = Treatment, shape = Depth)) +
  labs(x = paste("PC1: ", pcoaA1, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2, "%", sep = ""),
       shape = "Depth") +
  scale_colour_viridis_d() +
  scale_fill_viridis_d() +
  guides(colour = guide_legend(override.aes = list(shape = 15),
                               order = 1)) +
  theme_bw() +  
  ggtitle("Order") +
  theme(legend.position = "none",
        plot.title = element_text(size = 12, vjust = -1),
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        plot.margin = margin(5, 5, 5, 5, "pt"))

tax_sum_families <- summarize_taxonomy(input = nc, level = 5, report_higher_tax = T)
fam_bc <- calc_dm(tax_sum_families)
nc_pcoa <- cmdscale(fam_bc, k = nrow(nc$map_loaded) - 1, eig = T)
pcoaA1 <- round((eigenvals(nc_pcoa)/sum(eigenvals(nc_pcoa)))[1]*100, digits = 1)
pcoaA2 <- round((eigenvals(nc_pcoa)/sum(eigenvals(nc_pcoa)))[2]*100, digits = 1)
nc$map_loaded$Axis01 <- scores(nc_pcoa)[,1]
nc$map_loaded$Axis02 <- scores(nc_pcoa)[,2]
micro.hulls <- ddply(nc$map_loaded, c("Treatment"), find_hull)
g4 <- ggplot(nc$map_loaded, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = Treatment, fill = Treatment),
               alpha = 0.1, show.legend = F) +
  geom_point(size = 3, alpha = 1, aes(colour = Treatment, shape = Depth)) +
  labs(x = paste("PC1: ", pcoaA1, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2, "%", sep = ""),
       shape = "Depth") +
  scale_colour_viridis_d() +
  scale_fill_viridis_d() +
  guides(colour = guide_legend(override.aes = list(shape = 15),
                               order = 1)) +
  theme_bw() +  
  ggtitle("Family") +
  theme(legend.position = "none",
        plot.title = element_text(size = 12, vjust = -1),
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        plot.margin = margin(5, 5, 5, 5, "pt"))

tax_sum_genera <- summarize_taxonomy(input = nc, level = 6, report_higher_tax = T)
gen_bc <- calc_dm(tax_sum_genera)
nc_pcoa <- cmdscale(gen_bc, k = nrow(nc$map_loaded) - 1, eig = T)
pcoaA1 <- round((eigenvals(nc_pcoa)/sum(eigenvals(nc_pcoa)))[1]*100, digits = 1)
pcoaA2 <- round((eigenvals(nc_pcoa)/sum(eigenvals(nc_pcoa)))[2]*100, digits = 1)
nc$map_loaded$Axis01 <- scores(nc_pcoa)[,1]
nc$map_loaded$Axis02 <- scores(nc_pcoa)[,2]
micro.hulls <- ddply(nc$map_loaded, c("Treatment"), find_hull)
g5 <- ggplot(nc$map_loaded, aes(Axis01, -Axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = Treatment, fill = Treatment),
               alpha = 0.1, show.legend = F) +
  geom_point(size = 3, alpha = 1, aes(colour = Treatment, shape = Depth)) +
  labs(x = paste("PC1: ", pcoaA1, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2, "%", sep = ""),
       shape = "Depth") +
  scale_colour_viridis_d() +
  scale_fill_viridis_d() +
  guides(colour = guide_legend(override.aes = list(shape = 15),
                               order = 1)) +
  theme_bw() +  
  ggtitle("Genus") +
  theme(legend.position = "none",
        plot.title = element_text(size = 12, vjust = -1),
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        plot.margin = margin(5, 5, 5, 5, "pt"))

nc_pcoa <- cmdscale(nc_bc, k = nrow(nc$map_loaded) - 1, eig = T)
pcoaA1 <- round((eigenvals(nc_pcoa)/sum(eigenvals(nc_pcoa)))[1]*100, digits = 1)
pcoaA2 <- round((eigenvals(nc_pcoa)/sum(eigenvals(nc_pcoa)))[2]*100, digits = 1)
nc$map_loaded$Axis01 <- scores(nc_pcoa)[,1]
nc$map_loaded$Axis02 <- scores(nc_pcoa)[,2]
micro.hulls <- ddply(nc$map_loaded, c("Treatment"), find_hull)
g6 <- ggplot(nc$map_loaded, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = Treatment, fill = Treatment),
               alpha = 0.1, show.legend = F) +
  geom_point(size = 3, alpha = 1, aes(colour = Treatment, shape = Depth)) +
  labs(x = paste("PC1: ", pcoaA1, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2, "%", sep = ""),
       shape = "Depth") +
  scale_colour_viridis_d() +
  scale_fill_viridis_d() +
  guides(colour = guide_legend(override.aes = list(shape = 15),
                               order = 1)) +
  theme_bw() + 
  ggtitle("OTU") +
  theme(legend.position = "none",
        plot.title = element_text(size = 12, vjust = -1),
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        plot.margin = margin(5, 5, 5, 5, "pt"))

p1 <- plot_grid(g1, g2, g3, g4, g5, g6, ncol = 3)

png("InitialFigs/Beta_allLevels.png", width = 8, height = 6, units = "in", res = 300)
plot_grid(p1, leg, rel_widths = c(0.85, 0.15))
dev.off()



#### _Jaccard ####
nc_ja <- calc_dm(nc$data_loaded, method = "jaccard")
set.seed(1150)
adonis2(nc_ja ~ nc$map_loaded$Treatment + nc$map_loaded$Depth) # Both sig
adonis2(nc_ja ~ nc$map_loaded$Depth + nc$map_loaded$Treatment) # No effect of order
anova(betadisper(nc_ja, nc$map_loaded$Treatment)) # Dispersion not homogeneous
anova(betadisper(nc_ja, nc$map_loaded$Depth)) # Dispersion homogeneous

# PCoA with vectors
nc_pcoa <- cmdscale(nc_ja, k = nrow(nc$map_loaded) - 1, eig = T)
set.seed(100)
ef_nc <- envfit(nc_pcoa, env_nc, permutations = 999, na.rm = TRUE)
ef_nc
ordiplot(nc_pcoa)
plot(ef_nc, p.max = 0.05, cex = 0.5)
manual_factor_nc <- 0.45
vec.df_nc <- as.data.frame(ef_nc$vectors$arrows*sqrt(ef_nc$vectors$r)) %>%
  mutate(Dim1 = Dim1 * manual_factor_nc,
         Dim2 = Dim2 * manual_factor_nc) %>%
  mutate(variables = rownames(.)) %>%
  filter(ef_nc$vectors$pvals < 0.05) %>%
  filter(variables != "Cl_mgL") %>%
  mutate(shortnames = c("Salinity", "CH4", "N2O", "CO2", "NH4", "pH", "Br"))
pcoaA1 <- round((eigenvals(nc_pcoa)/sum(eigenvals(nc_pcoa)))[1]*100, digits = 1)
pcoaA2 <- round((eigenvals(nc_pcoa)/sum(eigenvals(nc_pcoa)))[2]*100, digits = 1)
nc$map_loaded$Axis01 <- scores(nc_pcoa)[,1]
nc$map_loaded$Axis02 <- scores(nc_pcoa)[,2]
micro.hulls <- ddply(nc$map_loaded, c("Treatment"), find_hull)
#png("InitialFigs/Beta_Jac.png", width = 7, height = 5, units = "in", res = 300)
png("FinalFigs/FigureS4.png", width = 7, height = 5, units = "in", res = 300)
ggplot(nc$map_loaded, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = Treatment, fill = Treatment),
               alpha = 0.1, show.legend = F) +
  geom_point(size = 3, alpha = 1, aes(colour = Treatment, shape = Depth)) +
  geom_segment(data = vec.df_nc,
               aes(x = 0, xend = Dim1, y = 0, yend = Dim2),
               arrow = arrow(length = unit(0.35, "cm")),
               colour = "gray", alpha = 0.6,
               inherit.aes = FALSE) + 
  geom_text(data = vec.df_nc,
            aes(x = Dim1, y = Dim2, label = shortnames),
            size = 3, color = "black") +
  labs(x = paste("PC1: ", pcoaA1, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2, "%", sep = ""),
       shape = "Depth") +
  scale_colour_viridis_d() +
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

# Check PCoA at different levels
tax_sum_phyla <- summarize_taxonomy(input = nc, level = 2, report_higher_tax = T)
phy_ja <- calc_dm(tax_sum_phyla, method = "jaccard")
nc_pcoa <- cmdscale(phy_ja, k = nrow(nc$map_loaded) - 1, eig = T)
pcoaA1 <- round((eigenvals(nc_pcoa)/sum(eigenvals(nc_pcoa)))[1]*100, digits = 1)
pcoaA2 <- round((eigenvals(nc_pcoa)/sum(eigenvals(nc_pcoa)))[2]*100, digits = 1)
nc$map_loaded$Axis01 <- scores(nc_pcoa)[,1]
nc$map_loaded$Axis02 <- scores(nc_pcoa)[,2]
micro.hulls <- ddply(nc$map_loaded, c("Treatment"), find_hull)
g1 <- ggplot(nc$map_loaded, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = Treatment, fill = Treatment),
               alpha = 0.1, show.legend = F) +
  geom_point(size = 3, alpha = 1, aes(colour = Treatment, shape = Depth)) +
  labs(x = paste("PC1: ", pcoaA1, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2, "%", sep = ""),
       shape = "Depth") +
  scale_colour_viridis_d() +
  scale_fill_viridis_d() +
  guides(colour = guide_legend(override.aes = list(shape = 15),
                               order = 1)) +
  theme_bw() +  
  ggtitle("Phylum") +
  theme(legend.position = "right",
        plot.title = element_text(size = 12, vjust = -1),
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        plot.margin = margin(5, 5, 5, 5, "pt"))
leg <- get_legend(g1)
g1 <- g1 + theme(legend.position = "none")

tax_sum_classes <- summarize_taxonomy(input = nc, level = 3, report_higher_tax = T)
cla_ja <- calc_dm(tax_sum_classes, method = "jaccard")
nc_pcoa <- cmdscale(cla_ja, k = nrow(nc$map_loaded) - 1, eig = T)
pcoaA1 <- round((eigenvals(nc_pcoa)/sum(eigenvals(nc_pcoa)))[1]*100, digits = 1)
pcoaA2 <- round((eigenvals(nc_pcoa)/sum(eigenvals(nc_pcoa)))[2]*100, digits = 1)
nc$map_loaded$Axis01 <- scores(nc_pcoa)[,1]
nc$map_loaded$Axis02 <- scores(nc_pcoa)[,2]
micro.hulls <- ddply(nc$map_loaded, c("Treatment"), find_hull)
g2 <- ggplot(nc$map_loaded, aes(Axis01, -Axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = Treatment, fill = Treatment),
               alpha = 0.1, show.legend = F) +
  geom_point(size = 3, alpha = 1, aes(colour = Treatment, shape = Depth)) +
  labs(x = paste("PC1: ", pcoaA1, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2, "%", sep = ""),
       shape = "Depth") +
  scale_colour_viridis_d() +
  scale_fill_viridis_d() +
  guides(colour = guide_legend(override.aes = list(shape = 15),
                               order = 1)) +
  theme_bw() + 
  ggtitle("Class") +
  theme(legend.position = "none",
        plot.title = element_text(size = 12, vjust = -1),
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        plot.margin = margin(5, 5, 5, 5, "pt"))

tax_sum_orders <- summarize_taxonomy(input = nc, level = 4, report_higher_tax = T)
ord_ja <- calc_dm(tax_sum_orders, method = "jaccard")
nc_pcoa <- cmdscale(ord_ja, k = nrow(nc$map_loaded) - 1, eig = T)
pcoaA1 <- round((eigenvals(nc_pcoa)/sum(eigenvals(nc_pcoa)))[1]*100, digits = 1)
pcoaA2 <- round((eigenvals(nc_pcoa)/sum(eigenvals(nc_pcoa)))[2]*100, digits = 1)
nc$map_loaded$Axis01 <- scores(nc_pcoa)[,1]
nc$map_loaded$Axis02 <- scores(nc_pcoa)[,2]
micro.hulls <- ddply(nc$map_loaded, c("Treatment"), find_hull)
g3 <- ggplot(nc$map_loaded, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = Treatment, fill = Treatment),
               alpha = 0.1, show.legend = F) +
  geom_point(size = 3, alpha = 1, aes(colour = Treatment, shape = Depth)) +
  labs(x = paste("PC1: ", pcoaA1, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2, "%", sep = ""),
       shape = "Depth") +
  scale_colour_viridis_d() +
  scale_fill_viridis_d() +
  guides(colour = guide_legend(override.aes = list(shape = 15),
                               order = 1)) +
  theme_bw() +  
  ggtitle("Order") +
  theme(legend.position = "none",
        plot.title = element_text(size = 12, vjust = -1),
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        plot.margin = margin(5, 5, 5, 5, "pt"))

tax_sum_families <- summarize_taxonomy(input = nc, level = 5, report_higher_tax = T)
fam_ja <- calc_dm(tax_sum_families, method = "jaccard")
nc_pcoa <- cmdscale(fam_ja, k = nrow(nc$map_loaded) - 1, eig = T)
pcoaA1 <- round((eigenvals(nc_pcoa)/sum(eigenvals(nc_pcoa)))[1]*100, digits = 1)
pcoaA2 <- round((eigenvals(nc_pcoa)/sum(eigenvals(nc_pcoa)))[2]*100, digits = 1)
nc$map_loaded$Axis01 <- scores(nc_pcoa)[,1]
nc$map_loaded$Axis02 <- scores(nc_pcoa)[,2]
micro.hulls <- ddply(nc$map_loaded, c("Treatment"), find_hull)
g4 <- ggplot(nc$map_loaded, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = Treatment, fill = Treatment),
               alpha = 0.1, show.legend = F) +
  geom_point(size = 3, alpha = 1, aes(colour = Treatment, shape = Depth)) +
  labs(x = paste("PC1: ", pcoaA1, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2, "%", sep = ""),
       shape = "Depth") +
  scale_colour_viridis_d() +
  scale_fill_viridis_d() +
  guides(colour = guide_legend(override.aes = list(shape = 15),
                               order = 1)) +
  theme_bw() +  
  ggtitle("Family") +
  theme(legend.position = "none",
        plot.title = element_text(size = 12, vjust = -1),
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        plot.margin = margin(5, 5, 5, 5, "pt"))

tax_sum_genera <- summarize_taxonomy(input = nc, level = 6, report_higher_tax = T)
gen_ja <- calc_dm(tax_sum_genera)
nc_pcoa <- cmdscale(gen_ja, k = nrow(nc$map_loaded) - 1, eig = T)
pcoaA1 <- round((eigenvals(nc_pcoa)/sum(eigenvals(nc_pcoa)))[1]*100, digits = 1)
pcoaA2 <- round((eigenvals(nc_pcoa)/sum(eigenvals(nc_pcoa)))[2]*100, digits = 1)
nc$map_loaded$Axis01 <- scores(nc_pcoa)[,1]
nc$map_loaded$Axis02 <- scores(nc_pcoa)[,2]
micro.hulls <- ddply(nc$map_loaded, c("Treatment"), find_hull)
g5 <- ggplot(nc$map_loaded, aes(Axis01, -Axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = Treatment, fill = Treatment),
               alpha = 0.1, show.legend = F) +
  geom_point(size = 3, alpha = 1, aes(colour = Treatment, shape = Depth)) +
  labs(x = paste("PC1: ", pcoaA1, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2, "%", sep = ""),
       shape = "Depth") +
  scale_colour_viridis_d() +
  scale_fill_viridis_d() +
  guides(colour = guide_legend(override.aes = list(shape = 15),
                               order = 1)) +
  theme_bw() +  
  ggtitle("Genus") +
  theme(legend.position = "none",
        plot.title = element_text(size = 12, vjust = -1),
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        plot.margin = margin(5, 5, 5, 5, "pt"))

nc_pcoa <- cmdscale(nc_ja, k = nrow(nc$map_loaded) - 1, eig = T)
pcoaA1 <- round((eigenvals(nc_pcoa)/sum(eigenvals(nc_pcoa)))[1]*100, digits = 1)
pcoaA2 <- round((eigenvals(nc_pcoa)/sum(eigenvals(nc_pcoa)))[2]*100, digits = 1)
nc$map_loaded$Axis01 <- scores(nc_pcoa)[,1]
nc$map_loaded$Axis02 <- scores(nc_pcoa)[,2]
micro.hulls <- ddply(nc$map_loaded, c("Treatment"), find_hull)
g6 <- ggplot(nc$map_loaded, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = Treatment, fill = Treatment),
               alpha = 0.1, show.legend = F) +
  geom_point(size = 3, alpha = 1, aes(colour = Treatment, shape = Depth)) +
  labs(x = paste("PC1: ", pcoaA1, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2, "%", sep = ""),
       shape = "Depth") +
  scale_colour_viridis_d() +
  scale_fill_viridis_d() +
  guides(colour = guide_legend(override.aes = list(shape = 15),
                               order = 1)) +
  theme_bw() + 
  ggtitle("OTU") +
  theme(legend.position = "none",
        plot.title = element_text(size = 12, vjust = -1),
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        plot.margin = margin(5, 5, 5, 5, "pt"))

p1 <- plot_grid(g1, g2, g3, g4, g5, g6, ncol = 3)

png("InitialFigs/Beta_allLevels_Jac.png", width = 8, height = 6, units = "in", res = 300)
plot_grid(p1, leg, rel_widths = c(0.85, 0.15))
dev.off()



#### 5. Taxa ####
cliffplot_taxa_bars(input = nc, level = 1, variable = "Treatment")
cliffplot_taxa_bars(input = nc, level = 2, variable = "Treatment")
cliffplot_taxa_bars(input = nc, level = 3, variable = "Treatment")
cliffplot_taxa_bars(input = nc, level = 4, variable = "Treatment")
cliffplot_taxa_bars(input = nc, level = 5, variable = "Treatment")
cliffplot_taxa_bars(input = nc, level = 6, variable = "Treatment")
cliffplot_taxa_bars(input = nc, level = 9, variable = "Treatment")

# Field SRB
tax_sum_guilds <- summarize_taxonomy(input = field, level = 9, report_higher_tax = F) %>%
  t() %>%
  as.data.frame() %>%
  mutate(SRB_tot = SRB + SRB_syn)
mean(tax_sum_guilds$SRB_tot)
se(tax_sum_guilds$SRB_tot)

# Phyla, all samples
tax_sum_phyla <- summarize_taxonomy(input = nc, level = 2, report_higher_tax = F)
nc$map_loaded$sampleID <- rownames(nc$map_loaded)
barsP <- plot_taxa_bars(tax_sum_phyla,
                        nc$map_loaded,
                        "sampleID_clean",
                        num_taxa = 12,
                        data_only = TRUE) %>%
  mutate(taxon = fct_rev(taxon)) %>%
  left_join(., nc$map_loaded, by = c("group_by" = "sampleID_clean"))
topphy <- barsP %>%
  group_by(taxon) %>%
  summarise(mean = mean(mean_value)) %>%
  filter(taxon != "Other") %>%
  arrange(-mean)
phy <- ggplot(barsP, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Sample", y = "Relative abundance", fill = "Phylum") +
  scale_fill_manual(values = c("grey90", brewer.pal(12, "Paired")[12:1]),
                    labels = c("Other", phy_stats_lab$StarLab)) +
  scale_y_continuous(expand = c(0.01, 0.01)) +  
  facet_nested(~ Treatment + Depth, space = "free", scales = "free_x") +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        strip.text = element_text(size = 7),
        strip.background = element_rect(linewidth = 0.2),
        axis.line.y = element_blank(),
        legend.margin = margin(0, 0, 0, 5, unit = "pt"),
        legend.box.margin = margin(0, 0, 0, 5, unit = "pt"),
        legend.key.size = unit(0.3, "cm"),
        panel.spacing.x = unit(c(0.25, 0.5, 0.25, 0.5, 0.25, 0.5, 0.25, 0.5, 0.25), "lines"))

# Guilds all samples
tax_sum_guilds <- summarize_taxonomy(input = nc, level = 9, report_higher_tax = F)
nc$map_loaded$sampleID <- rownames(nc$map_loaded)
barsG <- plot_taxa_bars(tax_sum_guilds,
                        nc$map_loaded,
                        "sampleID_clean",
                        num_taxa = 20,
                        data_only = TRUE) %>%
  filter(taxon != "NA") %>%
  droplevels() %>%
  mutate(taxon = factor(taxon,
                        levels = Guild_cols$Guild)) %>%
  left_join(., nc$map_loaded, by = c("group_by" = "sampleID_clean"))
tallest_bar <- barsG %>%
  group_by(group_by) %>%
  summarise(sum = sum(mean_value))
topgui <- barsG %>%
  group_by(taxon) %>%
  summarise(mean = mean(mean_value)) %>%
  arrange(-mean)
gui <- ggplot(barsG, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Sample", y = "Relative abundance", fill = "Guild") +
  scale_fill_manual(values = Guild_cols$color,
                    labels = gui_stats_lab$StarLab) +
  scale_y_continuous(expand = c(max(tallest_bar$sum)/100, max(tallest_bar$sum)/100)) + 
  facet_nested(~ Treatment + Depth, space = "free", scales = "free_x") +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5, hjust = 1),
        axis.ticks.x = element_blank(),
        strip.text = element_blank(),
        strip.background = element_rect(linewidth = 0.2),
        axis.line.y = element_blank(),
        legend.margin = margin(0, 0, 0, -10, unit = "pt"),
        legend.box.margin = margin(0, 0, 0, -10, unit = "pt"),
        legend.key.size = unit(0.3, "cm"),
        panel.spacing.x = unit(c(0.25, 0.5, 0.25, 0.5, 0.25, 0.5, 0.25, 0.5, 0.25), "lines"))

plot_grid(phy, gui, ncol = 1, rel_heights = c(0.45, 0.55), align = "v")

png("InitialFigs/PhylaGuilds.png", width = 8, height = 6, units = "in", res = 300)
plot_grid(phy, gui, ncol = 1, rel_heights = c(0.45, 0.55), align = "v")
dev.off()

png("FinalFigs/Figure4.png", width = 8, height = 6, units = "in", res = 300)
plot_grid(phy, gui, ncol = 1, rel_heights = c(0.45, 0.55), align = "v", labels = "auto")
dev.off()

# Check SRB abundance range
srbs <- as.data.frame(t(tax_sum_guilds)) %>%
  mutate(srbs = SRB + SRB_syn)
min(srbs$srbs)
max(srbs$srbs)

# Stats
# All data
phy_stats <- taxa_summary_by_sample_type(tax_sum_phyla, 
                                         nc$map_loaded, 
                                         type_header = 'Treatment', 
                                         filter_level = 0.01, 
                                         test_type = 'KW') %>%
  filter(rownames(.) %in% barsP$taxon) %>%
  arrange(desc(rownames(.))) %>%
  mutate(Sig = ifelse(pvalsFDR < 0.05, "Pfdr < 0.05", "Pfdr > 0.05"))

gui_stats <- taxa_summary_by_sample_type(tax_sum_guilds, 
                                         nc$map_loaded, 
                                         type_header = 'Treatment', 
                                         filter_level = 0, 
                                         test_type = 'KW') %>%
  filter(rownames(.) %in% barsG$taxon) %>%
  arrange(desc(rownames(.))) %>%
  mutate(Sig = ifelse(pvalsFDR < 0.05, "Pfdr < 0.05", "Pfdr > 0.05"))

# But this should be refined. Let's test field vs. lab control. And just lab treatments
fiecon <- filter_data(nc,
                      "Treatment",
                      keep_vals = c("Field", "Control"))
sum_phy_fiecon <- summarize_taxonomy(input = fiecon, level = 2, report_higher_tax = F)
phy_stats_fiecon <- taxa_summary_by_sample_type(sum_phy_fiecon, 
                                                fiecon$map_loaded, 
                                                type_header = 'Treatment', 
                                                filter_level = 0, 
                                                test_type = 'MW') %>%
  filter(rownames(.) %in% barsP$taxon) %>%
  arrange(desc(rownames(.))) %>%
  mutate(Sig = ifelse(pvalsFDR < 0.05, "Pfdr < 0.05", "Pfdr > 0.05")) %>%
  mutate(Star = ifelse(pvalsFDR < 0.05, "*", "")) %>%
  rownames_to_column(var = "taxon") %>%
  arrange(match(taxon, levels(barsP$taxon))) %>%
  mutate(StarLab = paste(Guild, Star, sep = " "))

sum_gui_fiecon <- summarize_taxonomy(input = fiecon, level = 9, report_higher_tax = F)
gui_stats_fiecon <- taxa_summary_by_sample_type(sum_gui_fiecon, 
                                                fiecon$map_loaded, 
                                                type_header = 'Treatment', 
                                                filter_level = 0, 
                                                test_type = 'MW') %>%
  filter(rownames(.) %in% barsG$taxon) %>%
  arrange(desc(rownames(.))) %>%
  mutate(Sig = ifelse(pvalsFDR < 0.05, "Pfdr < 0.05", "Pfdr > 0.05"))

lab <- filter_data(nc,
                   "Treatment",
                   filter_vals = "Field")
sum_phy_lab <- summarize_taxonomy(input = lab, level = 2, report_higher_tax = F)
phy_stats_lab <- taxa_summary_by_sample_type(sum_phy_lab, 
                                                lab$map_loaded, 
                                                type_header = 'Treatment', 
                                                filter_level = 0, 
                                                test_type = 'KW') %>%
  filter(rownames(.) %in% barsP$taxon) %>%
  arrange(desc(rownames(.))) %>%
  mutate(Sig = ifelse(pvalsFDR < 0.05, "Pfdr < 0.05", "Pfdr > 0.05")) %>%
  mutate(Star = ifelse(pvalsFDR < 0.05, "*", "")) %>%
  rownames_to_column(var = "taxon") %>%
  mutate(StarLab = paste(taxon, Star, sep = " "))

sum_gui_lab <- summarize_taxonomy(input = lab, level = 9, report_higher_tax = F)
gui_stats_lab <- taxa_summary_by_sample_type(sum_gui_lab, 
                                                lab$map_loaded, 
                                                type_header = 'Treatment', 
                                                filter_level = 0, 
                                                test_type = 'KW') %>%
  filter(rownames(.) %in% barsG$taxon) %>%
  arrange(desc(rownames(.))) %>%
  mutate(Sig = ifelse(pvalsFDR < 0.05, "Pfdr < 0.05", "Pfdr > 0.05")) %>%
  mutate(Star = ifelse(pvalsFDR < 0.05, "*", "")) %>%
  rownames_to_column(var = "Guild") %>%
  left_join(., Guild_cols, by = "Guild") %>%
  arrange(Index) %>%
  mutate(StarLab = paste(Guild, Star, sep = " "))
  


#### _Ratios ####
MG_MT <- summarize_taxonomy(nc, 9, report_higher_tax = F) %>%
  t() %>%
  as.data.frame() %>%
  mutate(MT = MOB_I + MOB_II + MOB_IIa + ANME)
sum(rownames(MG_MT) != rownames(nc$map_loaded))
nc$map_loaded$MT <- MG_MT$MT

summary(lm(log(MG_MT) ~ log(AO_NOB), data = nc$map_loaded)) # NSD

a <- ggplot(nc$map_loaded, aes(AO_NOB, MG_MT)) +
  geom_point(size = 2, aes(colour = Treatment, shape = Depth)) +
  geom_smooth(method = "lm", size = 0.5, alpha = 0.1, linetype = "dotted", se = F) +
  labs(x = "Ammonia oxidizers: Nitrite oxidizing bacteria",
       y = "Methanogens: Methanotrophs",
       colour = "Treatment",
       shape = "Depth") +
  scale_colour_viridis_d() +
  scale_x_continuous(trans = 'log10') +
  scale_y_continuous(trans = 'log10',
                     breaks = c(0.01, 0.1, 1, 10, 100),
                     labels = c(0.01, 0.1, 1, 10, 100)) +
  guides(colour = guide_legend(shape = c(17, 16, 18, 18, 15))) +
  theme_bw()
l <- get_legend(a)
a <- a + theme(legend.position = "none",
               axis.title = element_text(size = 10))

summary(lm(MT ~ log(AO_NOB), data = nc$map_loaded))

b <- ggplot(nc$map_loaded, aes(AO_NOB, MT*100)) +
  geom_point(size = 2, aes(colour = Treatment, shape = Depth)) +
  geom_smooth(method = "lm", size = 0.5, alpha = 0.1, linetype = "dotted", se = F) +
  labs(x = "Ammonia oxidizers: Nitrite oxidizing bacteria",
       y = "Methanotroph % abundance",
       colour = "Treatment",
       shape = "Depth") +
  scale_colour_viridis_d() +
  scale_x_continuous(trans = 'log10') +
  theme_bw() +
  theme(legend.position = "none",
        axis.title = element_text(size = 10))

png("InitialFigs/Ratios.png", width = 8, height = 4, units = "in", res = 300)
plot_grid(a, b, l, ncol = 3, rel_widths = c(0.43, 0.43, 0.14), labels = c("a", "b", ""))
dev.off()

png("FinalFigs/FigureS7.png", width = 8, height = 4, units = "in", res = 300)
plot_grid(a, b, l, ncol = 3, rel_widths = c(0.43, 0.43, 0.14), labels = c("a", "b", ""))
dev.off()

# MG:MT
leveneTest(nc$map_loaded$MG_MT ~ nc$map_loaded$Treatment) # Homogeneous
m <- aov(MG_MT ~ Treatment + Depth, data = nc$map_loaded)
Anova(m, type = "II") # Treatment
m <- aov(MG_MT ~ Treatment, data = nc$map_loaded)
shapiro.test(m$residuals) # Not normal
summary(m)
TukeyHSD(m)
t <- emmeans(object = m, specs = "Treatment") %>%
  cld(object = ., adjust = "Tukey", Letters = letters, alpha = 0.05) %>%
  mutate(name = "MG_MT",
         y = max(nc$map_loaded$MG_MT)+(max(nc$map_loaded$MG_MT)-min(nc$map_loaded$MG_MT))/2)

png("InitialFigs/MG_MT.png", width = 7, height = 5, units = "in", res = 300)
ggplot(nc$map_loaded, aes(reorder(Treatment, MG_MT, mean), MG_MT)) +
  geom_hline(yintercept = 1, linetype = "dotted") +
  geom_boxplot(aes(colour = Treatment), outlier.shape = NA) +
  geom_jitter(size = 3, width = 0.2, aes(colour = Treatment, shape = Depth)) + 
  geom_text(data = t, aes(Treatment, y, label = str_trim(.group)), 
            size = 4, color = "black") +
  labs(x = "Treatment",
       y = "Methanogens:Methanotrophs",
       shape = "Depth") +
  scale_y_continuous(trans = 'log10',
                     breaks = c(0.01, 0.1, 1, 10, 100),
                     labels = c(0.01, 0.1, 1, 10, 100)) +
  scale_colour_viridis_d() +
  guides(colour = "none") +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.position = c(1,0),
        legend.justification = c(1,0),
        legend.background = element_blank())
dev.off()



#### _Methanogens ####
# Need to fix SILVA typo: Halobacterota -> Halobacteriota
for (i in 1:nrow(nc$taxonomy_loaded)) {
  if (nc$taxonomy_loaded$taxonomy2[i] == "Halobacterota") {
    nc$taxonomy_loaded$taxonomy2[i] <- "Halobacteriota"
  }
}

nc_mg <- filter_taxa_from_input(nc,
                                taxa_to_keep = c("CH4_ac", "CH4_H2", "CH4_me", "CH4_mix"),
                                at_spec_level = 9)
tax_sum_mg <- summarize_taxonomy(input = nc_mg, 
                                 level = 5, 
                                 report_higher_tax = T, 
                                 relative = F) %>%
  mutate_all(funs((./82312)*100))
tax_sum_mg <- summarize_taxonomy(input = nc_mg, 
                                 level = 5, 
                                 report_higher_tax = F, 
                                 relative = F) %>%
  mutate_all(funs((./82312)*100))
barsMG <- plot_taxa_bars(tax_sum_mg,
                         nc_mg$map_loaded,
                         "sampleID_clean",
                         num_taxa = 20,
                         data_only = TRUE) %>%
  mutate(taxon = gsub("NA", "Unclassified", taxon)) %>%
  mutate(taxon = fct_relevel(taxon, "Unclassified", after = Inf)) %>%
  mutate(taxon = fct_rev(taxon)) %>%
  left_join(., nc_mg$map_loaded, by = c("group_by" = "sampleID_clean"))
topmg <- barsMG %>%
  group_by(taxon) %>%
  summarise(mean = mean(mean_value)) %>%
  arrange(-mean)
barsMG <- barsMG %>%
  mutate(taxon = factor(taxon, levels = topmg$taxon)) %>%
  mutate(taxon = fct_relevel(taxon, "Unclassified", after = Inf)) %>%
  mutate(taxon = fct_rev(taxon))
png("InitialFigs/Methanogens.png", width = 7, height = 5, units = "in", res = 300)
fig5 <- ggplot(barsMG, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Sample", y = "% Abundance", fill = "Family") +
  scale_fill_manual(values = c("grey90", brewer_pal(palette = "Paired")(8))) +
  scale_y_continuous(expand = c(0.01, 0.01)) +  
  facet_nested(~ Treatment + Depth, space = "free", scales = "free_x") +
  guides(fill = guide_legend(ncol = 1)) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5),
        axis.ticks.x = element_blank(),
        strip.text = element_text(size = 6),
        strip.background = element_rect(size = 0.2),
        axis.line.y = element_blank(),
        legend.position = c(0,1),
        legend.justification = c(0,1),
        legend.background = element_blank(),
        legend.key.size = unit(0.3, "cm"),
        panel.spacing.x = unit(c(0.25, 0.5, 0.25, 0.5, 0.25, 0.5, 0.25, 0.5, 0.25), "lines"))
fig5
dev.off()



#### _Methanotrophs ####
nc_mt <- filter_taxa_from_input(nc,
                                taxa_to_keep = c("ANME", "MOB_I", "MOB_II", "MOB_IIa"),
                                at_spec_level = 9)
tax_sum_mt <- summarize_taxonomy(input = nc_mt, 
                                 level = 6, 
                                 report_higher_tax = F, 
                                 relative = F) %>%
  # filter(rownames(.) != "NA") %>%
  mutate_all(funs((./82312)*100))
barsMT <- plot_taxa_bars(tax_sum_mt,
                         nc_mt$map_loaded,
                         "sampleID_clean",
                         num_taxa = 20,
                         data_only = TRUE) %>%
  mutate(taxon = gsub("NA", "Unclassified", taxon)) %>%
  mutate(taxon = fct_relevel(taxon, "Other", after = Inf)) %>%
  mutate(taxon = fct_relevel(taxon, "Unclassified", after = Inf)) %>%
  mutate(taxon = fct_rev(taxon)) %>%
  left_join(., nc_mt$map_loaded, by = c("group_by" = "sampleID_clean"))
topmt <- barsMT %>%
  group_by(taxon) %>%
  summarise(mean = mean(mean_value)) %>%
  arrange(-mean)
barsMT <- barsMT %>%
  mutate(taxon = factor(taxon, levels = topmt$taxon)) %>%
  mutate(taxon = fct_relevel(taxon, "Other", after = Inf)) %>%
  mutate(taxon = fct_relevel(taxon, "Unclassified", after = Inf)) %>%
  mutate(taxon = fct_rev(taxon))
nb.cols <- 19
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)

# Stats
lab <- filter_data(nc,
                   "Treatment",
                   filter_vals = "Field")
sum_gen_lab <- summarize_taxonomy(input = lab, level = 6, report_higher_tax = F)
gen_stats_lab <- taxa_summary_by_sample_type(sum_gen_lab,
                                             lab$map_loaded, 
                                             type_header = 'Treatment', 
                                             filter_level = 0, 
                                             test_type = 'KW') %>%
  filter(rownames(.) %in% barsMT$taxon) %>%
  arrange(desc(rownames(.))) %>%
  mutate(Pfdr = p.adjust(pvals, method = "fdr")) %>%
  mutate(Sig = ifelse(Pfdr < 0.05, "Pfdr < 0.05", "Pfdr > 0.05")) %>%
  mutate(Star = ifelse(Pfdr < 0.05, "*", "")) %>%
  rownames_to_column(var = "taxon") %>%
  arrange(match(taxon, levels(barsMT$taxon))) %>%
  mutate(StarLab = paste(taxon, Star, sep = " "))
statsMT <- barsMT %>%
  left_join(., gen_stats_lab, by = "taxon") %>%
  mutate(StarLab = ifelse(is.na(StarLab), taxon, StarLab)) %>%
  mutate(taxon = fct_relevel(taxon, "Other", after = Inf)) %>%
  mutate(taxon = fct_relevel(taxon, "Unclassified", after = Inf)) %>%
  mutate(taxon = fct_rev(taxon)) %>%
  group_by(taxon) %>%
  slice_head(n = 1)
stats_lab <- data.frame(taxon = levels(barsMT$taxon)) %>%
  left_join(., statsMT, by = "taxon")

png("InitialFigs/Methanotrophs.png", width = 6.5, height = 6, units = "in", res = 300)
fig6 <- ggplot(barsMT, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Sample", y = "% Abundance", fill = "Genus") +
  scale_fill_manual(values = c("grey75", "grey90", mycolors),
                    labels = stats_lab$StarLab) +
  scale_y_continuous(expand = c(0.01, 0.01)) +  
  facet_nested(~ Treatment + Depth, space = "free", scales = "free_x") +
  guides(fill = guide_legend(ncol = 1)) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5, hjust = 1),
        axis.ticks.x = element_blank(),
        strip.text = element_text(size = 7),
        strip.background = element_rect(size = 0.2),
        axis.line.y = element_blank(),
        legend.position = c(1,1),
        legend.justification = c(1,1),
        legend.background = element_blank(),
        legend.key.size = unit(0.3, "cm"),
        panel.spacing.x = unit(c(0.25, 0.5, 0.25, 0.5, 0.25, 0.5, 0.25, 0.5, 0.25), "lines"))
fig6
dev.off()

png("FinalFigs/Figure6.png", width = 6.5, height = 6, units = "in", res = 300)
fig6
dev.off()


#### _SRB ####
# Need to get top SRB taxa
# Then look at cultured isolates to check salinity tolerance
nc_srb <- filter_taxa_from_input(nc,
                                taxa_to_keep = c("SRB", "SRB_syn"),
                                at_spec_level = 9)
tax_sum_srb <- summarize_taxonomy(input = nc_srb, 
                                 level = 6, 
                                 report_higher_tax = F, 
                                 relative = F) %>%
  # filter(rownames(.) != "NA") %>%
  mutate_all(funs((./82312)*100))
barssrb <- plot_taxa_bars(tax_sum_srb,
                         nc_srb$map_loaded,
                         "sampleID_clean",
                         num_taxa = 20,
                         data_only = TRUE) %>%
  mutate(taxon = gsub("NA", "Unclassified", taxon)) %>%
  mutate(taxon = fct_relevel(taxon, "Other", after = Inf)) %>%
  mutate(taxon = fct_relevel(taxon, "Unclassified", after = Inf)) %>%
  mutate(taxon = fct_rev(taxon)) %>%
  left_join(., nc_srb$map_loaded, by = c("group_by" = "sampleID_clean"))
topsrb <- barssrb %>%
  group_by(taxon) %>%
  summarise(mean = mean(mean_value)) %>%
  arrange(-mean)
barssrb <- barssrb %>%
  mutate(taxon = factor(taxon, levels = topsrb$taxon)) %>%
  mutate(taxon = fct_relevel(taxon, "Other", after = Inf)) %>%
  mutate(taxon = fct_relevel(taxon, "Unclassified", after = Inf)) %>%
  mutate(taxon = fct_rev(taxon))
nb.cols <- 19
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)

# Stats
lab <- filter_data(nc,
                   "Treatment",
                   filter_vals = "Field")
sum_gen_lab <- summarize_taxonomy(input = lab, level = 6, report_higher_tax = F)
gen_stats_lab <- taxa_summary_by_sample_type(sum_gen_lab,
                                             lab$map_loaded, 
                                             type_header = 'Treatment', 
                                             filter_level = 0, 
                                             test_type = 'KW') %>%
  filter(rownames(.) %in% barssrb$taxon) %>%
  arrange(desc(rownames(.))) %>%
  mutate(Pfdr = p.adjust(pvals, method = "fdr")) %>%
  mutate(Sig = ifelse(Pfdr < 0.05, "Pfdr < 0.05", "Pfdr > 0.05")) %>%
  mutate(Star = ifelse(Pfdr < 0.05, "*", "")) %>%
  rownames_to_column(var = "taxon") %>%
  arrange(match(taxon, levels(barssrb$taxon))) %>%
  mutate(StarLab = paste(taxon, Star, sep = " "))
statssrb <- barssrb %>%
  left_join(., gen_stats_lab, by = "taxon") %>%
  mutate(StarLab = ifelse(is.na(StarLab), taxon, StarLab)) %>%
  mutate(taxon = fct_relevel(taxon, "Other", after = Inf)) %>%
  mutate(taxon = fct_relevel(taxon, "Unclassified", after = Inf)) %>%
  mutate(taxon = fct_rev(taxon)) %>%
  group_by(taxon) %>%
  slice_head(n = 1)
stats_lab <- data.frame(taxon = levels(barssrb$taxon)) %>%
  left_join(., statssrb, by = "taxon")

ggplot(barssrb, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Sample", y = "% Abundance", fill = "Genus") +
  scale_fill_manual(values = c("grey75", "grey90", mycolors),
                    labels = stats_lab$StarLab) +
  scale_y_continuous(expand = c(0.01, 0.01)) +  
  facet_nested(~ Treatment + Depth, space = "free", scales = "free_x") +
  guides(fill = guide_legend(ncol = 1)) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5, hjust = 1),
        axis.ticks.x = element_blank(),
        strip.text = element_text(size = 7),
        strip.background = element_rect(size = 0.2),
        axis.line.y = element_blank(),
        legend.position = c(0,1),
        legend.justification = c(0,1),
        legend.background = element_blank(),
        legend.key.size = unit(0.3, "cm"),
        panel.spacing.x = unit(c(0.25, 0.5, 0.25, 0.5, 0.25, 0.5, 0.25, 0.5, 0.25), "lines"))



#### _SRA ####
# Check if any sulfate reducing archaea (Archaeoglobus, Caldivirga)
nc_sra <- filter_taxa_from_input(nc,
                                 taxa_to_keep = c("Archaeoglobus", "Caldivirga"),
                                 at_spec_level = 6)
# Not found!

#### 6. Indicators ####
#### _Simper ####
nc$map_loaded$TrtDepth <- paste(nc$map_loaded$Treatment, nc$map_loaded$Depth, sep = "")
nc_sim <- simper(t(nc$data_loaded), 
                 nc$map_loaded$TrtDepth)
nc_s <- summary(nc_sim)
nc_df1 <- head(nc_s$`Control0-5 cm_+SO40-5 cm`, n = 10) %>%
  mutate(Comparison = "Control_+SO4",
         ASV = rownames(.),
         Depth = "0-5 cm") %>%
  rename("MeanControl" = ava,
         "MeanTrt" = avb)
nc_df2 <- head(nc_s$`+ASW-SO40-5 cm_Control0-5 cm`, n = 10) %>%
  mutate(Comparison = "Control_+ASW-SO4",
         ASV = rownames(.),
         Depth = "0-5 cm") %>%
  rename("MeanControl" = avb,
         "MeanTrt" = ava)
nc_df3 <- head(nc_s$`+ASW0-5 cm_Control0-5 cm`, n = 10) %>%
  mutate(Comparison = "Control_+ASW",
         ASV = rownames(.),
         Depth = "0-5 cm") %>%
  rename("MeanControl" = avb,
         "MeanTrt" = ava)
nc_df4 <- head(nc_s$`+SO410-15 cm_Control10-15 cm`, n = 10) %>%
  mutate(Comparison = "Control_+SO4",
         ASV = rownames(.),
         Depth = "10-15 cm") %>%
  rename("MeanControl" = avb,
         "MeanTrt" = ava)
nc_df5 <- head(nc_s$`+ASW-SO410-15 cm_Control10-15 cm`, n = 10) %>%
  mutate(Comparison = "Control_+ASW-SO4",
         ASV = rownames(.),
         Depth = "10-15 cm") %>%
  rename("MeanControl" = avb,
         "MeanTrt" = ava)
nc_df6 <- head(nc_s$`+ASW10-15 cm_Control10-15 cm`, n = 10) %>%
  mutate(Comparison = "Control_+ASW",
         ASV = rownames(.),
         Depth = "10-15 cm") %>%
  rename("MeanControl" = avb,
         "MeanTrt" = ava)
nc_simper_results <- rbind(nc_df1, nc_df2, nc_df3, nc_df4, nc_df5, nc_df6) %>%
  left_join(., nc$taxonomy_loaded, by = c("ASV" = "taxonomy8")) %>%
  mutate(Response = ifelse(MeanTrt > MeanControl, "Positive", "Negative")) %>%
  rename("OTU" = "ASV") %>%
  mutate(OTU = gsub("ASV", "OTU", OTU)) %>%
  rename(c("Domain" = "taxonomy1",
           "Phylum" = "taxonomy2",
           "Class" = "taxonomy3",
           "Order" = "taxonomy4",
           "Family" = "taxonomy5",
           "Genus" = "taxonomy6",
           "Species" = "taxonomy7",
           "Guild" = "taxonomy9",
           "CumulativeContribution" = "cumsum")) %>%
  mutate(MeanTrt = round((MeanTrt/82312)*100, digits = 2),
         MeanControl = round((MeanControl/82312)*100, digits = 2),
         CumulativeContribution = round(CumulativeContribution*100, digits = 2)) %>%
  unite(Taxonomy, Phylum, Class, Order, Family, Genus, Species, OTU,
                          sep = "; ") %>%
  mutate(Taxonomy = make.unique(Taxonomy)) %>%
  dplyr::select(Comparison, Depth, Response, Domain, Guild, Taxonomy, MeanTrt,
                MeanControl, CumulativeContribution)

simper_meta <- nc_simper_results %>%
  column_to_rownames(var = "Taxonomy") %>%
  dplyr::select(Guild, Domain, Response, Depth, Comparison)
simper_mat <- nc_simper_results %>%
  column_to_rownames(var = "Taxonomy") %>%
  dplyr::select(MeanControl, MeanTrt, CumulativeContribution) %>%
  as.matrix()
ann_rows <- data.frame(row.names = rownames(simper_meta),
                       "Guild" = simper_meta$Guild,
                       "Domain" = simper_meta$Domain,
                       "Response" = simper_meta$Response,
                       "Depth" = simper_meta$Depth, 
                       "Comparison" = simper_meta$Comparison)
ann_colors <- list(Guild = c(AOA = "#7CFC00",
                             FeOB = "#CD6600",
                             SRB = "#8B008B",
                             SRB_syn = "#CD2990",
                             `NA` = "white"),
                   Domain = c(Bacteria = "#440154FF",
                              Archaea = "#FCFDBFFF"),
                   Response = c(Positive = "red",
                                Negative = "blue"),
                   Depth = c(`0-5 cm` = "grey",
                             `10-15 cm` = "black"),
                   Comparison = c(`Control_+SO4` = "#21908CFF", 
                                  `Control_+ASW-SO4` = "#5DC863FF", 
                                  `Control_+ASW` = "#FDE725FF"))

pheatmap(simper_mat,
         legend = T,
         legend_breaks = c(0, 5, 10, 15, 20, 25, max(na.omit(simper_mat))),
         legend_labels = c("0   ", "5   ", "10   ", "15   ", "20   ", "25   ", "%\n"),
         main = "",
         border_color = NA,
         scale = "none",
         angle_col = 315,
         fontsize = 6,
         fontsize_row = 6,
         fontsize_col = 8,
         annotation_row = ann_rows,
         annotation_colors = ann_colors,
         cluster_rows = F,
         cluster_cols = F,
         display_numbers = T,
         gaps_row = c(10, 20, 30, 40, 50, 60),
         filename = "InitialFigs/Simper.png",
         width = 8,
         height = 7)
dev.off()
dev.set(dev.next())
dev.set(dev.next())

# Get % rel. abundance for the top 10 contributors of each comparison
tax_sum_ASV <- summarize_taxonomy(input = nc, level = 8, report_higher_tax = F)
nc_df1_abund <- tax_sum_ASV %>%
  filter(rownames(.) %in% nc_df1$ASV) %>%
  colSums(.) %>%
  as.data.frame() %>%
  set_names("RelAbund") %>%
  rownames_to_column(var = "sampleID") %>%
  mutate(Comparison = "Comp. 1")
nc_df2_abund <- tax_sum_ASV %>%
  filter(rownames(.) %in% nc_df2$ASV) %>%
  colSums(.) %>%
  as.data.frame() %>%
  set_names("RelAbund") %>%
  rownames_to_column(var = "sampleID") %>%
  mutate(Comparison = "Comp. 2")
nc_df3_abund <- tax_sum_ASV %>%
  filter(rownames(.) %in% nc_df3$ASV) %>%
  colSums(.) %>%
  as.data.frame() %>%
  set_names("RelAbund") %>%
  rownames_to_column(var = "sampleID") %>%
  mutate(Comparison = "Comp. 3")
nc_df4_abund <- tax_sum_ASV %>%
  filter(rownames(.) %in% nc_df4$ASV) %>%
  colSums(.) %>%
  as.data.frame() %>%
  set_names("RelAbund") %>%
  rownames_to_column(var = "sampleID") %>%
  mutate(Comparison = "Comp. 4")
nc_df5_abund <- tax_sum_ASV %>%
  filter(rownames(.) %in% nc_df5$ASV) %>%
  colSums(.) %>%
  as.data.frame() %>%
  set_names("RelAbund") %>%
  rownames_to_column(var = "sampleID") %>%
  mutate(Comparison = "Comp. 5")
nc_df6_abund <- tax_sum_ASV %>%
  filter(rownames(.) %in% nc_df6$ASV) %>%
  colSums(.) %>%
  as.data.frame() %>%
  set_names("RelAbund") %>%
  rownames_to_column(var = "sampleID") %>%
  mutate(Comparison = "Comp. 6")
check_abund <- rbind(nc_df1_abund, nc_df2_abund, nc_df3_abund,
                     nc_df4_abund, nc_df5_abund, nc_df6_abund)
check_abund_min_max <- check_abund %>%
  group_by(Comparison) %>%
  summarise(min = min(RelAbund),
            max = max(RelAbund))
mean(check_abund$RelAbund)
se(check_abund$RelAbund)




#### _Multipatt ####
tax_sum_OTU <- summarize_taxonomy(nc, level = 8, report_higher_tax = F, relative = T)
set.seed(425)
mp <- multipatt(t(nc$data_loaded), 
                nc$map_loaded$Treatment, 
                func = "r.g", 
                control = how(nperm=999))
summary(mp) # many

png("InitialFigs/Multipatt.png", width = 6, height = 8, units = "in", res = 300)
plot_multipatt(mp_obj = mp, 
               input = nc,
               tax_sum = tax_sum_OTU,
               group = "Treatment",
               filter = TRUE,
               filter_vals = "Field",
               abund = "% Rel. Abund.",
               qcut = 0.05,
               rcut = 0.65)
dev.off()

# Redo with only abundant taxa. Looks like the others are rare
# Remember rarefied to 82316
# 0.05 percent is 82316*0.0005 = 41.158
nrow(nc$taxonomy_loaded)
View(nc$data_loaded)
meancounts <- data.frame(meancount = sort(rowMeans(nc$data_loaded)))
nc_abund <- filter_taxa_from_input(nc, filter_thresh = 41.158)
nrow(nc_abund$taxonomy_loaded)
sort(rowMeans(nc_abund$data_loaded)/82316*100)
# This is the top 318 taxa

tax_sum_OTU <- summarize_taxonomy(nc_abund, level = 8, report_higher_tax = F, relative = T)
sort(rowMeans(tax_sum_OTU)*100)
set.seed(425)
mp <- multipatt(t(nc_abund$data_loaded), 
                nc_abund$map_loaded$Treatment, 
                func = "r.g", 
                control = how(nperm=999))
summary(mp) # Number of species associated to 1 group: 114

png("InitialFigs/Multipatt_abund.png", width = 6, height = 8, units = "in", res = 300)
plot_multipatt(mp_obj = mp, 
               input = nc,
               tax_sum = tax_sum_OTU,
               group = "Treatment",
               filter = TRUE,
               filter_vals = "Field",
               abund = "% Rel. Abund.",
               qcut = 0.05,
               rcut = 0)
dev.off()



#### 7. BGC ####
# Plot BGC variable by treatment, just for this time point

#### _Correlations ####
# Correlations with methane
# Corrplot
env_nona <- env_nona_nc %>%
  dplyr::select(-Cl_mgL)
C <- cor(env_nona)
corrplot(C, method = "number", type = "lower", 
         tl.cex = 0.5, tl.col = "black", number.cex = 0.5)

# Methane correlations
png("InitialFigs/BGC_CH4_cors.png", width = 6, height = 6, units = "in", res = 300)
meth_corr_by_bgc(env_nona = env_nona)
dev.off()



#### _Fluxes ####
# Boxplots by treatment - Fluxes
# Remove Field
nc_lab <- filter_data(nc,
                      "Treatment",
                      filter_vals = "Field")
flux_data <- nc_lab$map_loaded %>%
  group_by(CH4_ug_m2_h) %>%
  slice_head(n = 1)

leveneTest(flux_data$CH4_ug_m2_h ~ flux_data$Treatment) # Homogeneous
m <- aov(CH4_ug_m2_h ~ Treatment, data = flux_data)
shapiro.test(m$residuals) # Normal
hist(m$residuals)
summary(m)
TukeyHSD(m)
t1 <- emmeans(object = m, specs = "Treatment") %>%
  cld(object = ., adjust = "Tukey", Letters = letters, alpha = 0.05) %>%
  mutate(name = "CH4_ug_m2_h",
         y = max(flux_data$CH4_ug_m2_h)+(max(flux_data$CH4_ug_m2_h)-min(flux_data$CH4_ug_m2_h)) + 1000000)

g1 <- ggplot(flux_data, aes(Treatment, CH4_ug_m2_h)) +
  geom_boxplot(aes(colour = Treatment), outlier.shape = NA) +
  geom_jitter(size = 3, width = 0.1, aes(colour = Treatment)) + 
  geom_text(data = t1, aes(Treatment, y, label = str_trim(.group)), 
            size = 4, color = "black") +
  labs(x = "Treatment",
       y = expression(""*CH[4]*" flux (µg/"*m^2*"/h)")) +
  scale_y_continuous(trans = 'log10',
                     breaks = c(0.01, 0.1, 1, 10, 100, 1000, 10000, 100000),
                     labels = c(0.01, 0.1, 1, 10, 100, 1000, 10000, "100000")) +
  scale_colour_manual(values = viridis_pal()(5)[2:5]) +
  guides(colour = "none") +
  theme_bw() +
  theme(axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_blank())
g1

leveneTest(flux_data$CO2_ug_m2_h ~ flux_data$Treatment) # Not Homogeneous
m <- aov(CO2_ug_m2_h ~ Treatment, data = flux_data)
shapiro.test(m$residuals) # Almost Normal
hist(m$residuals)
summary(m)
TukeyHSD(m)
t2 <- emmeans(object = m, specs = "Treatment") %>%
  cld(object = ., adjust = "Tukey", Letters = letters, alpha = 0.05) %>%
  mutate(name = "CO2_ug_m2_h",
         y = max(flux_data$CO2_ug_m2_h)+(max(flux_data$CO2_ug_m2_h)-min(flux_data$CO2_ug_m2_h))/2)

g2 <- ggplot(flux_data, aes(Treatment, CO2_ug_m2_h)) +
  geom_boxplot(aes(colour = Treatment), outlier.shape = NA) +
  geom_jitter(size = 3, width = 0.1, aes(colour = Treatment)) + 
  geom_text(data = t2, aes(Treatment, y, label = str_trim(.group)), 
            size = 4, color = "black") +
  labs(x = "Treatment",
       y = expression(""*CO[2]*" flux (µg/"*m^2*"/h)")) +
  scale_y_continuous(trans = 'log10',
                     breaks = c(30000, 100000, 300000, 1000000),
                     labels = c("30000", "100000", "300000", "1000000")) +
  scale_colour_manual(values = viridis_pal()(5)[2:5]) +
  guides(colour = "none") +
  theme_bw() +
  theme(axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_blank())
g2

leveneTest(flux_data$N2O_ug_m2_h ~ flux_data$Treatment) # Almost Homogeneous
m <- aov(N2O_ug_m2_h ~ Treatment, data = flux_data)
shapiro.test(m$residuals) # Not Normal
hist(m$residuals)
summary(m)
TukeyHSD(m)
t3 <- emmeans(object = m, specs = "Treatment") %>%
  cld(object = ., adjust = "Tukey", Letters = letters, alpha = 0.05) %>%
  mutate(name = "N2O_ug_m2_h",
         y = max(flux_data$N2O_ug_m2_h)+(max(flux_data$N2O_ug_m2_h)-min(flux_data$N2O_ug_m2_h))/2)

g3 <- ggplot(flux_data, aes(Treatment, N2O_ug_m2_h)) +
  geom_boxplot(aes(colour = Treatment), outlier.shape = NA) +
  geom_jitter(size = 3, width = 0.1, aes(colour = Treatment)) + 
  geom_text(data = t3, aes(Treatment, y, label = str_trim(.group)), 
            size = 4, color = "black") +
  labs(x = "Treatment",
       y = expression(""*N[2]*"O flux (µg/"*m^2*"/h)")) +
  scale_y_continuous(trans = 'log10') +
  scale_colour_manual(values = viridis_pal()(5)[2:5]) +
  guides(colour = "none") +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 10))
g3

png("InitialFigs/Fluxes.png", width = 6, height = 6, units = "in", res = 300)
plot_grid(g1, g2, g3, ncol = 1, align = "v", rel_heights = c(0.32, 0.32, 0.36))
dev.off()

# Facet wrap
flux_data_long <- melt(flux_data, measure.vars = c("CH4_ug_m2_h", "CO2_ug_m2_h", "N2O_ug_m2_h"))
t_comb <- rbind(t1, t2, t3) %>%
  rename(variable = name)
facet_labels <- c("CH4_ug_m2_h" = "CH[4]",
                  "CO2_ug_m2_h" = "CO[2]",
                  "N2O_ug_m2_h" = "N[2]*O")
fw <- ggplot(flux_data_long, aes(Treatment, value)) +
  geom_boxplot(aes(colour = Treatment), outlier.shape = NA) +
  geom_jitter(size = 3, width = 0.1, aes(colour = Treatment)) + 
  geom_text(data = t_comb, aes(Treatment, y, label = str_trim(.group)), 
            size = 4, color = "black") +
  labs(x = "Treatment",
       y = expression("Gas flux (µg/"*m^2*"/h)")) +
  scale_y_continuous(trans = 'log10') +
  scale_colour_manual(values = viridis_pal()(5)[2:5]) +
  facet_wrap(~ variable, ncol = 1, scales = "free_y", labeller = as_labeller(facet_labels, default = label_parsed)) +
  guides(colour = "none") +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 10),
        strip.text = element_text(size = 10),
        strip.background = element_rect(fill = "white"))
fw

# Make pie charts with effect size
m1 <- aov(CH4_ug_m2_h ~ SO4 + ASW, data = flux_data)
summary(m1)
anova(m1)
Anova(m1, type = "II") # ASW
eta_sq_m1 <- eta_sq(m1) %>%
  as.data.frame() %>%
  rownames_to_column(var = "group") %>%
  set_names(c("group", "value")) %>%
  mutate(group = c("SO4", "ASW***", "Residuals")) %>%
  arrange(desc(group)) %>%
  mutate(ypos = cumsum(value) - 0.5*value)
p1 <- ggplot(eta_sq_m1, aes(x = "", y = value, fill = group)) +
  geom_bar(stat = "identity", width = 1, color = NA) +
  coord_polar("y", start = 0) +
  geom_text(aes(y = ypos, label = group), color = "white", size = 3, check_overlap = T) +
  scale_fill_manual(values = c("#5DC863FF", "grey70", "#21908CFF")) +
  theme_void() + 
  theme(legend.position = "none")

m2 <- aov(CO2_ug_m2_h ~ SO4 + ASW, data = flux_data)
summary(m2)
anova(m2)
Anova(m2, type = "II") # ASW
eta_sq_m2 <- eta_sq(m2) %>%
  as.data.frame() %>%
  rownames_to_column(var = "group") %>%
  set_names(c("group", "value")) %>%
  mutate(group = c("SO4", "ASW*", "Residuals")) %>%
  arrange(desc(group)) %>%
  mutate(ypos = cumsum(value) - 0.5*value)
p2 <- ggplot(eta_sq_m2, aes(x = "", y = value, fill = group)) +
  geom_bar(stat = "identity", width = 1, color = NA) +
  coord_polar("y", start = 0) +
  geom_text(aes(y = ypos, label = group), color = "white", size = 3, check_overlap = T) +
  scale_fill_manual(values = c("#5DC863FF", "grey70", "#21908CFF")) +
  theme_void() + 
  theme(legend.position = "none")

m3 <- aov(N2O_ug_m2_h ~ SO4 + ASW, data = flux_data)
summary(m3)
anova(m3)
Anova(m3, type = "II") # ASW
eta_sq_m3 <- eta_sq(m3) %>%
  as.data.frame() %>%
  rownames_to_column(var = "group") %>%
  set_names(c("group", "value")) %>%
  mutate(group = c("SO4", "ASW*", "Residuals")) %>%
  arrange(desc(group)) %>%
  mutate(ypos = cumsum(value) - 0.5*value)
p3 <- ggplot(eta_sq_m3, aes(x = "", y = value, fill = group)) +
  geom_bar(stat = "identity", width = 1, color = NA) +
  coord_polar("y", start = 0) +
  geom_text(aes(y = ypos, label = group), color = "white", size = 3, check_overlap = T) +
  scale_fill_manual(values = c("#5DC863FF", "grey70", "#21908CFF")) +
  theme_void() + 
  theme(legend.position = "none")

plot_grid(g1, p1, g2, p2, g3, p3, ncol = 2, 
          rel_heights = c(0.32, 0.32, 0.36), rel_widths = c(0.7, 0.3),
          align = "hv")
pies <- plot_grid(p1, p2, p3, NULL, ncol = 1, rel_heights = c(0.32, 0.32, 0.32, 0.04))
pies
graphs <- plot_grid(g1, g2, g3, ncol = 1, align = "v", rel_heights = c(0.32, 0.32, 0.36))
graphs
r1 <- plot_grid(g1, p1, ncol = 2, align = "h", rel_widths = c(0.7, 0.3))
r1
r2 <- plot_grid(g2, p2, ncol = 2, align = "h", rel_widths = c(0.7, 0.3))
r2
r3 <- plot_grid(g3, p3, ncol = 2, align = "h", rel_widths = c(0.7, 0.3))
r3

png("FinalFigs/Figure1.png", width = 6, height = 6, units = "in", res = 300)
plot_grid(fw, pies, ncol = 2, rel_widths = c(0.7, 0.3))
dev.off()



# Linear regressions and scatterplots - Fluxes
summary(lm(log(N2O_ug_m2_h) ~ log(CO2_ug_m2_h), data = flux_data)) # Sig.
summary(lm(log(CH4_ug_m2_h) ~ log(CO2_ug_m2_h), data = flux_data)) # NSD
summary(lm(log(CH4_ug_m2_h) ~ log(N2O_ug_m2_h), data = flux_data)) # Sig.

g4 <- ggplot(flux_data, aes(CO2_ug_m2_h, N2O_ug_m2_h)) +
  geom_point(size = 3, aes(colour = Treatment)) +
  geom_smooth(method = "lm", alpha = 0.1) +
  geom_text(aes(x = 0, y = Inf, hjust = -0.1, vjust = 2, label = "R^2 == 0.57"), 
            parse = TRUE, size = 3, check_overlap = TRUE) +
  geom_text(aes(x = 0, y = Inf, hjust = -0.1, vjust = 5, label = "p < 0.001"), 
            size = 3, check_overlap = TRUE) +
  labs(x = expression(""*CO[2]*" flux (µg/"*m^2*"/h)"),
       y = expression(""*N[2]*"O flux (µg/"*m^2*"/h)")) +
  scale_x_continuous(trans = 'log10') +
  scale_y_continuous(trans = 'log10') +
  scale_colour_manual(values = viridis_pal()(5)[2:5]) +
  guides(colour = "none") +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 10))

g5 <- ggplot(flux_data, aes(CO2_ug_m2_h, CH4_ug_m2_h)) +
  geom_point(size = 3, aes(colour = Treatment)) +
  geom_smooth(method = "lm", alpha = 0.1, linetype = "dashed") +
  geom_text(aes(x = Inf, y = Inf, hjust = 1.5, vjust = 2, label = "R^2 == 0.09"), 
            parse = TRUE, size = 3, check_overlap = TRUE) +
  geom_text(aes(x = Inf, y = Inf, hjust = 1.5, vjust = 5, label = "p = 0.26 "), 
            size = 3, check_overlap = TRUE) +
  labs(x = expression(""*CO[2]*" flux (µg/"*m^2*"/h)"),
       y = expression(""*CH[4]*" flux (µg/"*m^2*"/h)")) +
  scale_x_continuous(trans = 'log10') +
  scale_y_continuous(trans = 'log10') +
  scale_colour_manual(values = viridis_pal()(5)[2:5]) +
  guides(colour = "none") +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 10))

g6 <- ggplot(flux_data, aes(N2O_ug_m2_h, CH4_ug_m2_h)) +
  geom_point(size = 3, aes(colour = Treatment)) +
  geom_smooth(method = "lm", alpha = 0.1) +
  geom_text(aes(x = Inf, y = Inf, hjust = 1.5, vjust = 2, label = "R^2 == 0.67"), 
            parse = TRUE, size = 3, check_overlap = TRUE) +
  geom_text(aes(x = Inf, y = Inf, hjust = 1.5, vjust = 5, label = "p < 0.001"), 
            size = 3, check_overlap = TRUE) +
  labs(x = expression(""*N[2]*"O flux (µg/"*m^2*"/h)"),
       y = expression(""*CH[4]*" flux (µg/"*m^2*"/h)")) +
  scale_x_continuous(trans = 'log10') +
  scale_y_continuous(trans = 'log10') +
  scale_colour_manual(values = viridis_pal()(5)[2:5]) +
  guides(colour = "none") +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 10))

png("InitialFigs/FluxesScatter.png", width = 4, height = 6, units = "in", res = 300)
plot_grid(g4, g5, g6, ncol = 1, align = "v")
dev.off()

#### _Other BGC ###
# Other BGC by treatment, with stats
bgc <- nc_lab$map_loaded %>%
  dplyr::select(Treatment, Depth,
                Salinity, Br_mgL, SO4_mgL, 
                TOC_mgL, NH4_mgL, NO3_mgL, DIN_mgL, DON_mgL, TN_mgL, PO4_mgL, pH,
                sed_per_C, sed_per_N, sed_CN)
bgc_long <- melt(bgc,
                 id.vars = c("Treatment", "Depth"),
                 measure.vars = c(names(bgc)[3:16]))
bgc_only <- bgc %>%
  dplyr::select(-Treatment, -Depth)
m <- list()
t <- list()
for (i in 1:ncol(bgc_only)) {
  m[[i]] <- aov(bgc_only[,i] ~ bgc$Treatment)
  t[[i]] <- emmeans(object = m[[i]], specs = "Treatment") %>%
    cld(object = ., adjust = "Tukey", Letters = letters, alpha = 0.05) %>%
    mutate(name = names(bgc_only)[i],
           y = max(bgc_only[,i], na.rm = T)+
             (max(bgc_only[,i], na.rm = T)-min(bgc_only[,i], na.rm = T))/10) %>%
    rename(variable = name)
}
t_df <- do.call(rbind.data.frame, t) %>%
  mutate(variable = factor(variable, levels = levels(bgc_long$variable)))

facet_names <- c("Salinity" = "Salinity (ppt)",
                 "Br_mgL" = "Br (mg/L)",
                 "SO4_mgL" = "SO4 (mg/L)",
                 "TOC_mgL" = "TOC (mg/L)",
                 "NH4_mgL" = "NH4 (mg/L)",
                 "NO3_mgL" = "NO3 (mg/L)",
                 "DIN_mgL" = "DIN (mg/L)",
                 "DON_mgL" = "DON (mg/L)",
                 "TN_mgL" = "TN (mg/L)",
                 "PO4_mgL" = "PO4 (mg/L)",
                 "pH" = "pH",
                 "sed_per_C" = "% C",
                 "sed_per_N" = "% N",
                 "sed_CN" = "C:N")
# Figure S1
png("InitialFigs/BGC.png", width = 8, height = 6, units = "in", res = 300)
ggplot(bgc_long, aes(Treatment, value)) +
  geom_boxplot(aes(colour = Treatment), outlier.shape = NA) +
  geom_jitter(size = 2, width = 0.1, aes(colour = Treatment, shape = Depth)) + 
  geom_text(data = t_df, aes(Treatment, y, label = str_trim(.group)), 
            size = 3, color = "black") +
  labs(x = "Treatment",
       y = NULL) +
  scale_colour_manual(values = viridis_pal()(5)[2:5]) +
  facet_wrap(~ variable, scales = "free_y", ncol = 3, labeller = as_labeller(facet_names)) +
  guides(colour = "none") +
  theme_bw() +
  theme(axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        legend.position = c(1,0),
        legend.justification = c(1,0.5))
dev.off()

# Plot correlations for different taxonomic levels at a given % relative abundance threshold
png("InitialFigs/CH4_Phyla.png", width = 7, height = 5, units = "in", res = 300)
meth_corr_by_taxonomy(input = nc, level = 2, threshold = 0.5, data = "No")
dev.off()
meth_corr_by_taxonomy(input = nc, level = 3, threshold = 0.5, data = "No")
meth_corr_by_taxonomy(input = nc, level = 4, threshold = 0.5, data = "No")
meth_corr_by_taxonomy(input = nc, level = 5, threshold = 0.5, data = "No")
meth_corr_by_taxonomy(input = nc, level = 6, threshold = 0.5, data = "No")
meth_corr_by_taxonomy(input = nc, level = 8, threshold = 0.5, data = "No")
png("InitialFigs/CH4_Guilds.png", width = 7, height = 5, units = "in", res = 300)
meth_corr_by_taxonomy(input = nc, level = 9, threshold = 0, data = "No")
dev.off()



#### 8. O2/H2S ####
# Oxygen and Hydrogen sulfide data from Sep 8 (Day 14 of the experiment)
o_stdev <- read_xls("data/microelectrode data.xls", sheet = 8) %>%
  clean_names() %>%
  mutate(depth_cm = depth/1000) %>%
  dplyr::select(treatment, depth_cm, avg_o2_mg_l, stdev_o2_mg_l, 
                avg_h2s_mg_l, stdev_h2s_mg_l) %>%
  pivot_longer(cols = c(stdev_o2_mg_l, stdev_h2s_mg_l)) %>%
  mutate(name = gsub("stdev", "avg", name)) %>%
  mutate(treatment = factor(treatment,
                            levels = c("Control", "+SO4", "+ASW-SO4", "+ASW")))

o <- read_xls("data/microelectrode data.xls", sheet = 8) %>%
  clean_names() %>%
  mutate(depth_cm = depth/1000) %>%
  dplyr::select(treatment, depth_cm, avg_o2_mg_l, stdev_o2_mg_l, 
                avg_h2s_mg_l, stdev_h2s_mg_l) %>%
  pivot_longer(cols = c(avg_o2_mg_l, avg_h2s_mg_l)) %>%
  mutate(stdev = o_stdev$value) %>%
  mutate(se = stdev/sqrt(5)) %>%
  mutate(treatment = factor(treatment,
                            levels = c("Control", "+SO4", "+ASW-SO4", "+ASW")))


facet_names <- c("avg_h2s_mg_l" = "(a) Hydrogen sulfide (mg/L)",
                 "avg_o2_mg_l" = "(b) Oxygen (mg/L)")

png("FinalFigs/FigureS10.png", width = 6, height = 3, units = "in", res = 300)
ggplot(o, aes(x = value, y = depth_cm, colour = treatment, group = treatment)) +
  geom_point(position = position_jitter(seed = 1, width = 0, height = 0.03)) +
  geom_path(position = position_jitter(seed = 1, width = 0, height = 0.03)) +
  geom_errorbarh(aes(xmin = value - se,
                     xmax = value + se,
                     y = depth_cm),
                 position = position_jitter(seed = 1, width = 0, height = 0.03),
                 colour = "black", linewidth = 0.1, alpha = 0.5) +
  scale_y_reverse(limits = c(5, -0.1),
                  breaks = c(5, 4, 3, 2, 1, 0),
                  labels = c(5, 4, 3, 2, 1, 0)) +
  scale_colour_manual(values = viridis_pal()(5)[2:5]) +
  labs(y = "Depth (cm)",
       x = NULL,
       colour = "Treatment") +
  facet_wrap(~ name, scales = "free_x", labeller = as_labeller(facet_names)) +
  theme_bw()
dev.off()



#### 9. NCBI ####
#write.csv(nc$map_loaded, "data/metadata_used.csv")
# Repset is in EastCoast folder, then filter to these samples.
# Also make sample-sequence map file
library(microseq)
f <- microseq::readFasta("~/Documents/GitHub/EastCoast/data/repset.comb.fasta") %>%
  separate(Header, into = c("Header", "Junk"), sep = " ") %>%
  dplyr::select(-Junk) %>%
  filter(Header %in% rownames(nc$data_loaded))
microseq::writeFasta(f, "data/repset_used.fasta")

# Map sequences to samples - Info is in data_loaded
# This takes a while!
info <- nc$data_loaded
names(info) <- nc$map_loaded$sampleID_clean
for (i in 1:ncol(info)) {
  for (j in 1:nrow(info)) {
    ifelse(info[j, i] > 0, info[j, i] <- names(info)[i], info[j, i] <- "")
  }
}

# Merge columns
info_cat <- info
info_cat <- info_cat %>%
  mutate_all(na_if, "") %>%
  mutate(unite(., "sample_name", c(names(info)), sep = ", ")) %>%
  mutate(sample_name = gsub("NA, ", "", sample_name)) %>%
  mutate(sample_name = gsub(", NA", "", sample_name)) %>%
  rownames_to_column(var = "Sequence_ID") %>%
  dplyr::select(Sequence_ID, sample_name)

# They don't accept multiple samples, so make another one with just the first sample
info_first <- info_cat %>%
  separate(sample_name, into = c("sample_name", "Junk"), sep = ", ") %>%
  dplyr::select(Sequence_ID, sample_name)

# Save
write_tsv(info_cat, file = "data/biosample_assignment.tsv")
write_tsv(info_first, file = "data/biosample_assignment_first.tsv")

# Filter out NCBI flagged sequences
flagged <- read.csv("data/NCBI_filter.csv")
f <- readFasta("data/repset_used.fasta") %>%
  filter(Header %notin% flagged$ASV_ID)
microseq::writeFasta(f, "data/repset_used_filtered.fasta")
sum(f$Header %in% flagged$ASV_ID)

info_cat <- read_tsv("data/biosample_assignment.tsv") %>%
  filter(Sequence_ID %notin% flagged$ASV_ID)
info_first <- read_tsv("data/biosample_assignment_first.tsv") %>%
  filter(Sequence_ID %notin% flagged$ASV_ID)
write_tsv(info_cat, file = "data/biosample_assignment_filt.tsv")
write_tsv(info_first, file = "data/biosample_assignment_first_filt.tsv")

# Filter a second time (ASV_32675)
f <- readFasta("data/repset_used_filtered.fasta") %>%
  filter(Header != "ASV_32675")
microseq::writeFasta(f, "data/repset_used_filtered2.fasta")
info_cat <- read_tsv("data/biosample_assignment_filt.tsv") %>%
  filter(Sequence_ID != "ASV_32675")
info_first <- read_tsv("data/biosample_assignment_first_filt.tsv") %>%
  filter(Sequence_ID != "ASV_32675")
write_tsv(info_cat, file = "data/biosample_assignment_filt2.tsv")
write_tsv(info_first, file = "data/biosample_assignment_first_filt2.tsv")


#### End Script ####