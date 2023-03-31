# Timberlake metagenomic 16S data analysis
# SSU and LSU extracted and annotated with mTAGs program
# by Cliff Bueno de Mesquita, Tringe Lab, JGI, Spring 2023



#### 1. Overview ####
# Samples were sequenced and processed at JGI
# SSU and LSU extracted and annotated with mTAGs program
# This script will basically replicate the amplicon analysis Timberlake16S_iTAG.R
# Main difference is that these are only D2 (10-15 cm) depths, so depth will not be a factor here
# Other difference is that some replicates are present that were missing in iTAG (ASW_2, Control_2)
# Also has section to compare iTAG and mTAG

# Key Goals:
# Assess microbial taxonomic and functional response to treatments
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
## 8. Compare iTAG and mTAG



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

# Functions
find_hull <- function(df) df[chull(df$Axis01, df$Axis02),]
`%notin%` <- Negate(`%in%`)

# Guild subsetting module from other repo
source("~/Documents/GitHub/SF_microbe_methane/modules/3_OTU_subsetting_modules_v.0.4_strip.r")
paste_ranks = function(sm_taxa){
  k = data.frame(k ="k_", sm_taxa['taxonomy1'])
  k2 <- do.call(paste, c(k, sep = ""))
  
  p = data.frame(k ="p_", sm_taxa['taxonomy2'])
  p2 <- do.call(paste, c(p, sep = ""))
  
  c = data.frame(k ="c_", sm_taxa['taxonomy3'])
  c2 <- do.call(paste, c(c, sep = ""))
  
  o = data.frame(k ="o_", sm_taxa['taxonomy4'])
  o2 <- do.call(paste, c(o, sep = ""))
  
  f = data.frame(k ="f_", sm_taxa['taxonomy5'])
  f2 <- do.call(paste, c(f, sep = ""))
  
  g = data.frame(k ="g_", sm_taxa['taxonomy6'])
  g2 <- do.call(paste, c(g, sep = ""))
  
  # NOT USING SPECIES HERE! OTU preprocessing doesn't!
  s = data.frame(k ="s_", sm_taxa['taxonomy7'])
  s2 <- do.call(paste, c(s, sep = ""))
  
  # combine all
  lineage_df = data.frame(k2, p2, c2, o2, f2, g2)
  lineage = do.call(paste, c(lineage_df, sep = ';'))
  return(lineage)
}

# Correlation functions from other repo
source("~/Documents/GitHub/EastCoast/meth_corr_by_taxonomy.R")
source("~/Documents/GitHub/EastCoast/meth_corr_by_bgc.R")

# Plotting functions from other repo
source("~/Documents/GitHub/EastCoast/cliffplot_taxa_bars.R")
source("~/Documents/GitHub/Extremophilic_Fungi/plot_multipatt.R")

# Repository path
setwd("~/Documents/GitHub/Timberlake/")

# Wyatt Hartman's guild color palette
# Note that MeOB, ANME, Anamx don't exist in this dataset, so removed
# Extra methanogen guilds added so colors added too
Guild_cols <- read.table("~/Documents/GitHub/SF_microbe_methane/data/colors/Guild_color_palette.txt",
                         sep='\t') %>%
  dplyr::select(Guild, G_index, color) %>%
  set_names(c("Guild", "Index", "color")) %>%
  mutate(Index = rev(Index)) %>%
  #add_row(Guild = "ANME", Index = 10, color = "#836FFF") %>%
  add_row(Guild = "CH4_me", Index = 15, color = "#FDC086") %>%
  add_row(Guild = "CH4_mix", Index = 16, color = "#FFFF99") %>%
  filter(Guild != "MeOB") %>%
  filter(Guild != "ANME") %>%
  filter(Guild != "Anamx") %>%
  arrange(Index)

# Prepare data. Only need to do once. Then skip to start here
# Need to make mctoolsr object from the merged mTAGs output
# Input file and reformat everything
# Separate euks and proks
p <- read.delim("data/merged_profile.otu.tsv") %>%
  set_names(substr(names(.), 1, nchar(names(.))-5)) %>%
  rename(taxonomy = X.ta) %>%
  filter(., !grepl("Eukaryota", taxonomy)) %>%
  mutate(taxonomy = gsub("root__Root;", "", taxonomy)) %>%
  mutate(taxonomy = gsub("domain__", "", taxonomy)) %>%
  mutate(taxonomy = gsub("phylum__", "", taxonomy)) %>%
  mutate(taxonomy = gsub("class__", "", taxonomy)) %>%
  mutate(taxonomy = gsub("order__", "", taxonomy)) %>%
  mutate(taxonomy = gsub("family__", "", taxonomy)) %>%
  mutate(taxonomy = gsub("genus__", "", taxonomy)) %>%
  mutate(taxonomy = gsub("otu__", "", taxonomy)) %>%
  mutate(taxonomy = gsub("silva_138_complink_cons_", "", taxonomy)) %>%
  mutate(taxonomy = gsub("unknown", "NA", taxonomy)) %>%
  mutate(taxonomy = gsub("otu_", "NA;otu_", taxonomy)) %>%
  filter(taxonomy != "Unassigned") %>%
  filter(taxonomy != "Unaligned") %>%
  separate(taxonomy, into = c("a","s","d","f","g","h","j","otu"), remove = F, sep = ";") %>%
  dplyr::select(-c(a,s,d,f,g,h,j)) %>%
  dplyr::select(otu, everything()) %>%
  dplyr::select(-taxonomy, taxonomy)
names(p) <- gsub("ASW.", "ASW_", names(p))
names(p) <- gsub("Control", "Control_", names(p))
names(p) <- gsub("Field", "Field_", names(p))

e <- read.delim("data/merged_profile.otu.tsv") %>%
  set_names(substr(names(.), 1, nchar(names(.))-5)) %>%
  rename(taxonomy = X.ta) %>%
  filter(., grepl("Eukaryota", taxonomy)) %>%
  mutate(taxonomy = gsub("root__Root;", "", taxonomy)) %>%
  mutate(taxonomy = gsub("domain__", "", taxonomy)) %>%
  mutate(taxonomy = gsub("phylum__", "", taxonomy)) %>%
  mutate(taxonomy = gsub("class__", "", taxonomy)) %>%
  mutate(taxonomy = gsub("order__", "", taxonomy)) %>%
  mutate(taxonomy = gsub("family__", "", taxonomy)) %>%
  mutate(taxonomy = gsub("genus__", "", taxonomy)) %>%
  mutate(taxonomy = gsub("otu__", "", taxonomy)) %>%
  mutate(taxonomy = gsub("silva_138_complink_cons_", "", taxonomy)) %>%
  mutate(taxonomy = gsub("unknown", "NA", taxonomy)) %>%
  mutate(taxonomy = gsub("otu_", "NA;otu_", taxonomy)) %>%
  filter(taxonomy != "Unassigned") %>%
  filter(taxonomy != "Unaligned") %>%
  separate(taxonomy, into = c("a","s","d","f","g","h","j","otu"), remove = F, sep = ";") %>%
  dplyr::select(-c(a,s,d,f,g,h,j)) %>%
  dplyr::select(otu, everything()) %>%
  dplyr::select(-taxonomy, taxonomy)
names(e) <- gsub("ASW.", "ASW_", names(e))
names(e) <- gsub("Control", "Control_", names(e))
names(e) <- gsub("Field", "Field_", names(e))

# Sequencing info
mtag_output <- read.delim("data/merged_profile.otu.tsv") %>%
  set_names(substr(names(.), 1, nchar(names(.))-5)) %>%
  rename(taxonomy = X.ta)
mtag_output_t <- mtag_output %>%
  column_to_rownames(var = "taxonomy") %>%
  t() %>%
  as.data.frame()
depth_info <- data.frame("Bacteria+Archaea" = colSums(p[,2:24]),
                         "Eukaryota" = colSums(e[,2:24]),
                         #"Classified" = colSums(mtag_output[1:10316, 2:24]),
                         "Unclassified" = mtag_output_t$Unassigned,
                         "Unaligned" = mtag_output_t$Unaligned,
                         "Treatment" = c(rep("ASW-SO4", 3), 
                                         rep("ASW", 5),
                                         rep("Control", 5),
                                         rep("Field", 5),
                                         rep("SO4", 5))) %>%
  mutate(sampleID = rownames(.)) %>%
  mutate(sampleID = gsub("ASW.", "ASW_", sampleID)) %>%
  mutate(sampleID = gsub("Control", "Control_", sampleID)) %>%
  mutate(sampleID = gsub("Field", "Field_", sampleID)) %>%
  mutate(Treatment = factor(Treatment, levels = c("Field", "Control", "SO4", "ASW-SO4", "ASW")))
mean(depth_info$Bacteria.Archaea)
mean(depth_info$Eukaryota)
mean(depth_info$Unclassified)
mean(depth_info$Unaligned)
depth_info_rel <- depth_info %>%
  mutate(Total = Bacteria.Archaea + Eukaryota + Unclassified + Unaligned) %>%
  mutate(BacArcRel = round(Bacteria.Archaea/Total*100, digits = 2),
         EukRel = round(Eukaryota/Total*100, digits = 2))
depth_info_long <- melt(depth_info, id.vars = c("sampleID", "Treatment"))
png("InitialFigs/mTAG_info.png", width = 6, height = 3, units = "in", res = 300)
ggplot(depth_info_long, aes(sampleID, value, fill = variable)) +
  geom_bar(stat = "identity") +
  labs(x = NULL,
       y = "Reads",
       fill = "Type") +
  facet_grid(~ Treatment, space = "free", scales = "free_x") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

# Check and write files
sum(names(p) %in% depth_info$sampleID)
out_fp <- "~/Documents/GitHub/Timberlake/seqtab_wTax_mctoolsr_mtagProk.txt"
names(p)[1] = "#OTU_ID"
write("#Exported for mctoolsr", out_fp)
suppressWarnings(write.table(p, out_fp, sep = "\t", row.names = FALSE, append = TRUE))


sum(names(e) %in% depth_info$sampleID)
out_fp <- "~/Documents/GitHub/Timberlake/seqtab_wTax_mctoolsr_mtagEuk.txt"
names(e)[1] = "#OTU_ID"
write("#Exported for mctoolsr", out_fp)
suppressWarnings(write.table(e, out_fp, sep = "\t", row.names = FALSE, append = TRUE))
  
# iTags (to compare). Also make mapping file for mTAGs. 
# NB: Metagenomes are from D2!
# NB: Use Wyatt's mapping file
nc <- readRDS("data/nc.rds")

wyatt_map <- read.delim("data/Timberlake_sample_map_both.txt") %>%
  dplyr::select(Itag_sample, MG_name) %>%
  filter(., !grepl("d1", Itag_sample))
wyatt_map <- wyatt_map[rowSums(is.na(wyatt_map)) != ncol(wyatt_map),]

mTAGs_map <- nc$map_loaded %>%
  filter(Depth == "10-15 cm") %>%
  dplyr::select(1:3, 10:24) %>%
  full_join(., wyatt_map, by = c("sampleID" = "Itag_sample")) %>%
  arrange(MG_name) %>%
  drop_na(MG_name) %>%
  dplyr::select(MG_name, everything()) %>%
  dplyr::select(-Treatment, -sampleID) %>%
  mutate(MG_name = gsub("SourceSoil", "Field", MG_name)) %>%
  left_join(., depth_info, by = c("MG_name" = "sampleID")) %>%
  dplyr::select(-Classified, -Unclassified, -Unaligned) %>%
  dplyr::select(MG_name, Treatment, everything())
#write.table(mTAGs_map, "~/Documents/GitHub/Timberlake/mtags_metadata.txt", sep = "\t", row.names = F)
# Add in additional metadata in Excel
# NB: Salinity = Cl_mgL * 0.0018066. 
# https://www.ldeo.columbia.edu/edu/k12/snapshotday/activities/2011/Classroom%20HS%20activity/chloride%20conversion/Chloride%20and%20Salinity.pdf

# Import Data (n = 23), Filter, Rarefy, Calc richness
tax_table_fp <- "data/seqtab_wTax_mctoolsr_mtagProk.txt"
map_fp <- "data/mtags_metadata.txt"
input = load_taxa_table(tax_table_fp, map_fp)

# Filter chloroplast, mitochondria, eukaryotes, unassigned at domain
input_filt <- filter_taxa_from_input(input,
                                     taxa_to_remove = "Chloroplast") # 13 removed
input_filt <- filter_taxa_from_input(input_filt,
                                     taxa_to_remove = "Mitochondria") # 23 removed
input_filt <- filter_taxa_from_input(input_filt,
                                     taxa_to_remove = "Eukaryota") # none
# (those were already filtered)
input_filt <- filter_taxa_from_input(input_filt,
                                     taxa_to_remove = "NA",
                                     at_spec_level = 1) # none 
# (those were already filtered by mTAGs as "Unclassified)

# Remove singletons and doubletons
singdoub <- data.frame("count" = rowSums(input_filt$data_loaded)) %>%
  filter(count < 3) %>%
  mutate(ASV = rownames(.))

input_filt <- filter_taxa_from_input(input_filt,
                                     taxa_IDs_to_remove = singdoub$ASV)

# Guilds
# Use updated Wyatt guild calling script
# Version that pulls out 4 different methanogen guilds
# Need to first reorganize OTU table
# Then place guilds as a column in input$taxonomy_loaded to analyze as any other taxonomic group
# Follow Wyatt's pipeline to prepare an mctoolsr object for guild analysis
taxa = input_filt$taxonomy_loaded
OTUs = input_filt$data_loaded
raw_OTU_table = data.frame(OTUs, taxa)
taxa[taxa=='NA'] <- ""
Consensus.lineage = paste_ranks(taxa)
reformed_OTU_table = data.frame(OTUs, Consensus.lineage) %>%
  mutate_if(is.integer, as.numeric) %>%
  mutate(OTU = rownames(.)) %>%
  dplyr::select(OTU, everything()) %>%
  left_join(., input_filt$taxonomy_loaded, by = c("OTU" = "taxonomy8")) %>%
  dplyr::rename(Kingdom = taxonomy1,
                Phylum= taxonomy2,
                Class = taxonomy3,
                Order = taxonomy4,
                Family = taxonomy5,
                Genus = taxonomy6) %>%
  dplyr::select(-taxonomy7) %>%
  mutate(Taxonomy = Phylum)
rownames(reformed_OTU_table) <- reformed_OTU_table$OTU
ANME <- otu_V[grepl("(ANME|Methanoperedenaceae|Syntrophoarchaeaceae)", otu_V$Consensus.lineage),]
anamox_list <-"(Kuenenia|Anammoxoglobus|Scalindua|Brocadia|Jettenia)"
anamox <- reformed_OTU_table[grepl(anamox_list, reformed_OTU_table$Consensus.lineage),]
Guild_OTUs <- Get_16S_Guilds_alt(reformed_OTU_table)
levels(as.factor(Guild_OTUs$Guild))

# Now add as 9th column to input_filt$taxonomy_loaded
input_filt$taxonomy_loaded <- input_filt$taxonomy_loaded %>%
  left_join(., Guild_OTUs, by = c("taxonomy8" = "OTU")) %>%
  rename(taxonomy9 = Guild) %>%
  mutate(taxonomy9 = as.character(taxonomy9)) %>%
  mutate(taxonomy9 = replace_na(taxonomy9, "NA"))
rownames(input_filt$taxonomy_loaded) <- input_filt$taxonomy_loaded$taxonomy8

# Save
#saveRDS(input_filt, "data/input_filt_mTAGs.rds")

# Rarefy at minimum (4271)
sort(colSums(input_filt$data_loaded))
mean(colSums(input_filt$data_loaded))
se(colSums(input_filt$data_loaded))
# Depth 6554.304 ±281.0974
set.seed(530)
input_filt_rare <- single_rarefy(input_filt, 4271)
sort(colSums(input_filt_rare$data_loaded))

# OTU Richness
input_filt_rare$map_loaded$rich <- specnumber(input_filt_rare$data_loaded, 
                                              MARGIN = 2)

# Shannon diversity
input_filt_rare$map_loaded$shannon <- diversity(input_filt_rare$data_loaded, 
                                                index = "shannon", 
                                                MARGIN = 2)

# Save
#saveRDS(input_filt_rare, "data/input_filt_rare_mTAGs.rds")

# Start here
nc <- readRDS("data/input_filt_rare_mTAGs.rds")
nc$map_loaded <- nc$map_loaded %>%
  mutate("Treatment" = c(rep("+ASW-SO4", 3), 
                         rep("+ASW", 5),
                         rep("Control", 5),
                         rep("Field", 5),
                         rep("+SO4", 5))) %>%
  mutate(Treatment = factor(Treatment, levels = c("Field", "Control", "+SO4", "+ASW-SO4", "+ASW")))
nc$map_loaded$sampleID <- rownames(nc$map_loaded)
nc$map_loaded$sampleID[5] <- "TL_inc_d2_ASW_2"
nc$map_loaded$sampleID[10] <- "TL_inc_d2_DI_ctrl_2"



#### 3. Alpha ####
leveneTest(nc$map_loaded$rich ~ nc$map_loaded$Treatment) # Homogeneous
m <- aov(rich ~ Treatment, data = nc$map_loaded)
shapiro.test(m$residuals) # Almost normal
summary(m)
t <- emmeans(object = m, specs = "Treatment") %>%
  cld(object = ., adjust = "Tukey", Letters = letters, alpha = 0.05) %>%
  mutate(name = "rich",
         y = max(nc$map_loaded$rich)+(max(nc$map_loaded$rich)-min(nc$map_loaded$rich))/20)
leveneTest(nc$map_loaded$shannon ~ nc$map_loaded$Treatment) # Homogeneous
m1 <- aov(shannon ~ Treatment, data = nc$map_loaded)
shapiro.test(m1$residuals) # Almost normal
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
png("InitialFigs/Alpha_mTAG.png", width = 6, height = 3, units = "in", res = 300)
ggplot(alpha_long, aes(reorder(Treatment, value, mean), value, 
                       colour = Treatment)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 2, alpha = 0.75, width = 0.2) +
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
dev.off()
  


#### 4. Beta ####
# Get env. vars with no NA
env_nc <- nc$map_loaded %>%
  dplyr::select(3:17)
env_nona_nc <- na.omit(env_nc)
nrow(env_nona_nc) # n = 14



#### _Aitch ####
# Use non-rarefied data, do CLR transformation, Aitchison, PCA
sum(rownames(nc$map_loaded) != rownames(input_filt$map_loaded))
input_filt$map_loaded <- nc$map_loaded
dim(input_filt$data_loaded)

# CLT transformation
otu_czm <- cmultRepl(t((input_filt$data_loaded)), label = 0, method = "CZM")
otu_clr <- clr(otu_czm)
aclr <- compositions::dist(otu_clr)

set.seed(1150)
adonis2(aclr ~ nc$map_loaded$Treatment) # R2 = 0.41, p = 0.001
anova(betadisper(aclr, nc$map_loaded$Treatment)) # Dispersion homogeneous

# PCA with vectors
d.pcx <- prcomp(aclr)
set.seed(100)
ef_nc <- envfit(d.pcx, env_nc, permutations = 999, na.rm = TRUE)
ef_nc
ordiplot(d.pcx)
plot(ef_nc, p.max = 0.075, cex = 0.5)
manual_factor_nc <- 0.3
vec.df_nc <- as.data.frame(ef_nc$vectors$arrows*sqrt(ef_nc$vectors$r)) %>%
  mutate(PC1 = PC1 * manual_factor_nc,
         PC2 = PC2 * manual_factor_nc) %>%
  mutate(variables = rownames(.)) %>%
  filter(ef_nc$vectors$pvals < 0.075) %>%
  filter(variables != "Cl_mgL") %>%
  mutate(shortnames = c("Salinity", "N2O", "CO2", "NH4", "pH", "Br"))
d.mvar <- sum(d.pcx$sdev^2)
PC1 <- paste("PC1: ", round((sum(d.pcx$sdev[1]^2)/d.mvar)*100, 1), "%")
PC2 <- paste("PC2: ", round((sum(d.pcx$sdev[2]^2)/d.mvar)*100, 1), "%")
nc$map_loaded$Axis01 <- d.pcx$rotation[,1]
nc$map_loaded$Axis02 <- d.pcx$rotation[,2]
micro.hulls <- ddply(nc$map_loaded, c("Treatment"), find_hull)
png("InitialFigs/BetaAitch_mTAG.png", width = 7, height = 5, units = "in", res = 300)
ggplot(nc$map_loaded, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = Treatment, fill = Treatment),
               alpha = 0.1, show.legend = F) +
  geom_point(size = 3, alpha = 1, aes(colour = Treatment)) +
  geom_segment(data = vec.df_nc,
               aes(x = 0, xend = -PC1, y = 0, yend = -PC2),
               arrow = arrow(length = unit(0.35, "cm")),
               colour = "gray", alpha = 0.6,
               inherit.aes = FALSE) + 
  geom_text(data = vec.df_nc,
            aes(x = -PC1, y = -PC2, label = shortnames),
            size = 3, color = "black") +
  labs(x = PC1, 
       y = PC2) +
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



#### _Bray ####
nc_bc <- calc_dm(nc$data_loaded)
set.seed(1150)
adonis2(nc_bc ~ nc$map_loaded$Treatment) # R2 = 0.41, p = 0.001
anova(betadisper(nc_bc, nc$map_loaded$Treatment)) # Dispersion homogeneous

# PCoA with vectors
nc_pcoa <- cmdscale(nc_bc, k = nrow(nc$map_loaded) - 1, eig = T)
set.seed(100)
ef_nc <- envfit(nc_pcoa, env_nc, permutations = 999, na.rm = TRUE)
ef_nc
ordiplot(nc_pcoa)
plot(ef_nc, p.max = 0.1, cex = 0.5)
manual_factor_nc <- 0.45
vec.df_nc <- as.data.frame(ef_nc$vectors$arrows*sqrt(ef_nc$vectors$r)) %>%
  mutate(Dim1 = Dim1 * manual_factor_nc,
         Dim2 = Dim2 * manual_factor_nc) %>%
  mutate(variables = rownames(.)) %>%
  filter(ef_nc$vectors$pvals < 0.1) %>%
  filter(variables != "Cl_mgL") %>%
  mutate(shortnames = c("Salinity", "N2O", "pH", "Br"))
pcoaA1 <- round((eigenvals(nc_pcoa)/sum(eigenvals(nc_pcoa)))[1]*100, digits = 1)
pcoaA2 <- round((eigenvals(nc_pcoa)/sum(eigenvals(nc_pcoa)))[2]*100, digits = 1)
nc$map_loaded$Axis01 <- scores(nc_pcoa)[,1]
nc$map_loaded$Axis02 <- scores(nc_pcoa)[,2]
micro.hulls <- ddply(nc$map_loaded, c("Treatment"), find_hull)
png("InitialFigs/Beta_mTAG.png", width = 7, height = 5, units = "in", res = 300)
ggplot(nc$map_loaded, aes(-Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = Treatment, fill = Treatment),
               alpha = 0.1, show.legend = F) +
  geom_point(size = 3, alpha = 1, aes(colour = Treatment)) +
  geom_segment(data = vec.df_nc,
               aes(x = 0, xend = -Dim1, y = 0, yend = Dim2),
               arrow = arrow(length = unit(0.35, "cm")),
               colour = "gray", alpha = 0.6,
               inherit.aes = FALSE) + 
  geom_text(data = vec.df_nc,
            aes(x = -Dim1, y = Dim2, label = shortnames),
            size = 3, color = "black") +
  labs(x = paste("PC1: ", pcoaA1, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2, "%", sep = ""),
       ) +
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
mod$anova # Cl, CH4

# Check PCoA at different levels
tax_sum_phyla <- summarize_taxonomy(input = nc, level = 2, report_higher_tax = T)
phy_bc <- calc_dm(tax_sum_phyla)
nc_pcoa <- cmdscale(phy_bc, k = nrow(nc$map_loaded) - 1, eig = T)
pcoaA1 <- round((eigenvals(nc_pcoa)/sum(eigenvals(nc_pcoa)))[1]*100, digits = 1)
pcoaA2 <- round((eigenvals(nc_pcoa)/sum(eigenvals(nc_pcoa)))[2]*100, digits = 1)
nc$map_loaded$Axis01 <- scores(nc_pcoa)[,1]
nc$map_loaded$Axis02 <- scores(nc_pcoa)[,2]
micro.hulls <- ddply(nc$map_loaded, c("Treatment"), find_hull)
g1 <- ggplot(nc$map_loaded, aes(-Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = Treatment, fill = Treatment),
               alpha = 0.1, show.legend = F) +
  geom_point(size = 3, alpha = 1, aes(colour = Treatment)) +
  labs(x = paste("PC1: ", pcoaA1, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2, "%", sep = ""),
       ) +
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
g2 <- ggplot(nc$map_loaded, aes(-Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = Treatment, fill = Treatment),
               alpha = 0.1, show.legend = F) +
  geom_point(size = 3, alpha = 1, aes(colour = Treatment)) +
  labs(x = paste("PC1: ", pcoaA1, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2, "%", sep = ""),
       ) +
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
g3 <- ggplot(nc$map_loaded, aes(-Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = Treatment, fill = Treatment),
               alpha = 0.1, show.legend = F) +
  geom_point(size = 3, alpha = 1, aes(colour = Treatment)) +
  labs(x = paste("PC1: ", pcoaA1, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2, "%", sep = ""),
       ) +
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
g4 <- ggplot(nc$map_loaded, aes(-Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = Treatment, fill = Treatment),
               alpha = 0.1, show.legend = F) +
  geom_point(size = 3, alpha = 1, aes(colour = Treatment)) +
  labs(x = paste("PC1: ", pcoaA1, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2, "%", sep = ""),
       ) +
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
g5 <- ggplot(nc$map_loaded, aes(-Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = Treatment, fill = Treatment),
               alpha = 0.1, show.legend = F) +
  geom_point(size = 3, alpha = 1, aes(colour = Treatment)) +
  labs(x = paste("PC1: ", pcoaA1, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2, "%", sep = ""),
       ) +
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
g6 <- ggplot(nc$map_loaded, aes(-Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = Treatment, fill = Treatment),
               alpha = 0.1, show.legend = F) +
  geom_point(size = 3, alpha = 1, aes(colour = Treatment)) +
  labs(x = paste("PC1: ", pcoaA1, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2, "%", sep = ""),
       ) +
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

png("InitialFigs/Beta_allLevels_mTAG.png", width = 8, height = 6, units = "in", res = 300)
plot_grid(p1, leg, rel_widths = c(0.85, 0.15))
dev.off()



#### _Jaccard ####
nc_ja <- calc_dm(nc$data_loaded, method = "jaccard")
set.seed(1150)
adonis2(nc_ja ~ nc$map_loaded$Treatment) # R2 = 0.30, p = 0.001
anova(betadisper(nc_ja, nc$map_loaded$Treatment)) # Dispersion homogeneous

# PCoA with vectors
nc_pcoa <- cmdscale(nc_ja, k = nrow(nc$map_loaded) - 1, eig = T)
set.seed(100)
ef_nc <- envfit(nc_pcoa, env_nc, permutations = 999, na.rm = TRUE)
ef_nc
ordiplot(nc_pcoa)
plot(ef_nc, p.max = 0.1, cex = 0.5)
manual_factor_nc <- 0.35
vec.df_nc <- as.data.frame(ef_nc$vectors$arrows*sqrt(ef_nc$vectors$r)) %>%
  mutate(Dim1 = Dim1 * manual_factor_nc,
         Dim2 = Dim2 * manual_factor_nc) %>%
  mutate(variables = rownames(.)) %>%
  filter(ef_nc$vectors$pvals < 0.1) %>%
  filter(variables != "Cl_mgL") %>%
  mutate(shortnames = c("Salinity", "N2O_ug_m2_h", "pH", "Br"))
pcoaA1 <- round((eigenvals(nc_pcoa)/sum(eigenvals(nc_pcoa)))[1]*100, digits = 1)
pcoaA2 <- round((eigenvals(nc_pcoa)/sum(eigenvals(nc_pcoa)))[2]*100, digits = 1)
nc$map_loaded$Axis01 <- scores(nc_pcoa)[,1]
nc$map_loaded$Axis02 <- scores(nc_pcoa)[,2]
micro.hulls <- ddply(nc$map_loaded, c("Treatment"), find_hull)
png("InitialFigs/Beta_Jac_mTAG.png", width = 7, height = 5, units = "in", res = 300)
ggplot(nc$map_loaded, aes(-Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = Treatment, fill = Treatment),
               alpha = 0.1, show.legend = F) +
  geom_point(size = 3, alpha = 1, aes(colour = Treatment)) +
  geom_segment(data = vec.df_nc,
               aes(x = 0, xend = -Dim1, y = 0, yend = Dim2),
               arrow = arrow(length = unit(0.35, "cm")),
               colour = "gray", alpha = 0.6,
               inherit.aes = FALSE) + 
  geom_text(data = vec.df_nc,
            aes(x = -Dim1, y = Dim2, label = shortnames),
            size = 3, color = "black") +
  labs(x = paste("PC1: ", pcoaA1, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2, "%", sep = ""),
       ) +
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
g1 <- ggplot(nc$map_loaded, aes(-Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = Treatment, fill = Treatment),
               alpha = 0.1, show.legend = F) +
  geom_point(size = 3, alpha = 1, aes(colour = Treatment)) +
  labs(x = paste("PC1: ", pcoaA1, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2, "%", sep = ""),
       ) +
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
g2 <- ggplot(nc$map_loaded, aes(-Axis01, -Axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = Treatment, fill = Treatment),
               alpha = 0.1, show.legend = F) +
  geom_point(size = 3, alpha = 1, aes(colour = Treatment)) +
  labs(x = paste("PC1: ", pcoaA1, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2, "%", sep = ""),
       ) +
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
g3 <- ggplot(nc$map_loaded, aes(-Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = Treatment, fill = Treatment),
               alpha = 0.1, show.legend = F) +
  geom_point(size = 3, alpha = 1, aes(colour = Treatment)) +
  labs(x = paste("PC1: ", pcoaA1, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2, "%", sep = ""),
       ) +
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
g4 <- ggplot(nc$map_loaded, aes(-Axis01, -Axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = Treatment, fill = Treatment),
               alpha = 0.1, show.legend = F) +
  geom_point(size = 3, alpha = 1, aes(colour = Treatment)) +
  labs(x = paste("PC1: ", pcoaA1, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2, "%", sep = ""),
       ) +
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
g5 <- ggplot(nc$map_loaded, aes(-Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = Treatment, fill = Treatment),
               alpha = 0.1, show.legend = F) +
  geom_point(size = 3, alpha = 1, aes(colour = Treatment)) +
  labs(x = paste("PC1: ", pcoaA1, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2, "%", sep = ""),
       ) +
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
g6 <- ggplot(nc$map_loaded, aes(-Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = Treatment, fill = Treatment),
               alpha = 0.1, show.legend = F) +
  geom_point(size = 3, alpha = 1, aes(colour = Treatment)) +
  labs(x = paste("PC1: ", pcoaA1, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2, "%", sep = ""),
       ) +
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

png("InitialFigs/Beta_allLevels_Jac_mTAG.png", width = 8, height = 6, units = "in", res = 300)
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

# Phyla, all samples
tax_sum_phyla <- summarize_taxonomy(input = nc, level = 2, report_higher_tax = F)
nc$map_loaded$sampleID <- rownames(nc$map_loaded)
barsP <- plot_taxa_bars(tax_sum_phyla,
                        nc$map_loaded,
                        "sampleID",
                        num_taxa = 12,
                        data_only = TRUE) %>%
  mutate(taxon = fct_rev(taxon)) %>%
  left_join(., nc$map_loaded, by = c("group_by" = "sampleID"))
topphy <- barsP %>%
  group_by(taxon) %>%
  summarise(mean = mean(mean_value)) %>%
  filter(taxon != "Other") %>%
  arrange(-mean)
phy <- ggplot(barsP, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Sample", y = "Relative abundance", fill = "Phylum") +
  scale_fill_manual(values = c("grey90", brewer.pal(12, "Paired")[12:1])) +
  scale_y_continuous(expand = c(0.01, 0.01)) +  
  facet_grid(~ Treatment, space = "free", scales = "free_x") +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        strip.text = element_text(size = 4),
        strip.background = element_rect(size = 0.2),
        axis.line.y = element_blank(),
        legend.margin = margin(0, 0, 0, 5, unit = "pt"),
        legend.box.margin = margin(0, 0, 0, 5, unit = "pt"),
        legend.key.size = unit(0.4, "cm"))

# Guilds all samples
tax_sum_guilds <- summarize_taxonomy(input = nc, level = 9, report_higher_tax = F)
nc$map_loaded$sampleID <- rownames(nc$map_loaded)
barsG <- plot_taxa_bars(tax_sum_guilds,
                        nc$map_loaded,
                        "sampleID",
                        num_taxa = 20,
                        data_only = TRUE) %>%
  filter(taxon != "NA") %>%
  droplevels() %>%
  mutate(taxon = factor(taxon,
                        levels = Guild_cols$Guild)) %>%
  left_join(., nc$map_loaded, by = c("group_by" = "sampleID"))
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
  scale_fill_manual(values = Guild_cols$color) +
  scale_y_continuous(expand = c(max(tallest_bar$sum)/100, max(tallest_bar$sum)/100)) + 
  facet_nested(~ Treatment, space = "free", scales = "free_x") +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5),
        axis.ticks.x = element_blank(),
        strip.text = element_blank(),
        strip.background = element_rect(size = 0.2),
        axis.line.y = element_blank(),
        legend.margin = margin(0, 0, 0, 5, unit = "pt"),
        legend.box.margin = margin(0, 0, 0, 5, unit = "pt"),
        legend.key.size = unit(0.4, "cm"))

plot_grid(phy, gui, ncol = 1, rel_heights = c(0.45, 0.55), align = "v", axis = "trbl")

png("InitialFigs/PhylaGuilds_mTAG.png", width = 8, height = 6, units = "in", res = 300)
plot_grid(phy, gui, ncol = 1, rel_heights = c(0.45, 0.55), align = "v", axis = "trbl")
dev.off()

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
  mutate(Sig = ifelse(pvalsFDR < 0.05, "Pfdr < 0.05", "Pfdr > 0.05"))

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
  mutate(Sig = ifelse(pvalsFDR < 0.05, "Pfdr < 0.05", "Pfdr > 0.05"))

sum_gui_lab <- summarize_taxonomy(input = lab, level = 9, report_higher_tax = F)
gui_stats_lab <- taxa_summary_by_sample_type(sum_gui_lab, 
                                                lab$map_loaded, 
                                                type_header = 'Treatment', 
                                                filter_level = 0, 
                                                test_type = 'KW') %>%
  filter(rownames(.) %in% barsG$taxon) %>%
  arrange(desc(rownames(.))) %>%
  mutate(Sig = ifelse(pvalsFDR < 0.05, "Pfdr < 0.05", "Pfdr > 0.05"))

#### _Ratios ####
# Add MG:MT ratio and AO:NOB ratio to metadata
nc_guilds <- summarize_taxonomy(input = nc, level = 9, report_higher_tax = F) %>%
  t() %>%
  as.data.frame() %>%
  mutate(AO_NOB = (AOA + AOB) / NOB) %>%
  mutate(MG = CH4_ac + CH4_H2 + CH4_me + CH4_mix) %>%
  mutate(MT = MOB_I + MOB_II + MOB_IIa) %>%
  mutate(MG_MT = (CH4_ac + CH4_H2 + CH4_me + CH4_mix)/(MOB_I + MOB_II + MOB_IIa)) %>%
  mutate_all(function(x) ifelse(is.infinite(x), 0.031608523, x))
  
sum(rownames(nc_guilds) != rownames(nc$map_loaded))
nc$map_loaded$MG_MT <- nc_guilds$MG_MT
nc$map_loaded$AO_NOB <- nc_guilds$AO_NOB

summary(lm(log(MG_MT) ~ log(AO_NOB), data = nc$map_loaded)) # Sig +

a <- ggplot(nc$map_loaded, aes(AO_NOB, MG_MT)) +
  geom_point(size = 2, aes(colour = Treatment)) +
  geom_smooth(method = "lm", size = 0.5, alpha = 0.1, linetype = "solid", se = F) +
  labs(x = "Ammonia oxidizers: Nitrite oxidizing bacteria",
       y = "Methanogens: Methanotrophs",
       colour = "Treatment",
       ) +
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

summary(lm(MT ~ log(AO_NOB), data = nc$map_loaded)) # NS

b <- ggplot(nc$map_loaded, aes(AO_NOB, MT*100)) +
  geom_point(size = 2, aes(colour = Treatment)) +
  geom_smooth(method = "lm", size = 0.5, alpha = 0.1, linetype = "dotted", se = F) +
  labs(x = "Ammonia oxidizers: Nitrite oxidizing bacteria",
       y = "Methanotroph % abundance",
       colour = "Treatment",
       ) +
  scale_colour_viridis_d() +
  scale_x_continuous(trans = 'log10') +
  theme_bw() +
  theme(legend.position = "none",
        axis.title = element_text(size = 10))

png("InitialFigs/Ratios_mTAG.png", width = 8, height = 4, units = "in", res = 300)
plot_grid(a, b, l, ncol = 3, rel_widths = c(0.43, 0.43, 0.14), labels = c("a", "b", ""))
dev.off()

# MG:MT
leveneTest(nc$map_loaded$MG_MT ~ nc$map_loaded$Treatment) # Homogeneous
m <- aov(MG_MT ~ Treatment, data = nc$map_loaded)
shapiro.test(m$residuals) # Not normal
hist(m$residuals)
summary(m) # NSD
TukeyHSD(m)
t <- emmeans(object = m, specs = "Treatment") %>%
  cld(object = ., adjust = "Tukey", Letters = letters, alpha = 0.05) %>%
  mutate(name = "MG_MT",
         y = max(nc$map_loaded$MG_MT)+(max(nc$map_loaded$MG_MT)-min(nc$map_loaded$MG_MT))/2)

png("InitialFigs/MG_MT_mTAG.png", width = 7, height = 5, units = "in", res = 300)
ggplot(nc$map_loaded, aes(reorder(Treatment, MG_MT, mean), MG_MT)) +
  geom_hline(yintercept = 1, linetype = "dotted") +
  geom_boxplot(aes(colour = Treatment), outlier.shape = NA) +
  geom_jitter(size = 3, width = 0.2, aes(colour = Treatment)) + 
  geom_text(data = t, aes(Treatment, y, label = str_trim(.group)), 
            size = 4, color = "black") +
  labs(x = "Treatment",
       y = "Methanogens:Methanotrophs") +
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
nc_mg <- filter_taxa_from_input(nc,
                                  taxa_to_keep = c("CH4_ac", "CH4_H2", "CH4_me", "CH4_mix"),
                                  at_spec_level = 9)
tax_sum_mg <- summarize_taxonomy(input = nc_mg, 
                                 level = 5, 
                                 report_higher_tax = T, 
                                 relative = F) %>%
  mutate_all(funs((./4271)*100))
rownames(tax_sum_mg) <- gsub("Archaea; ", "", rownames(tax_sum_mg))
barsMG <- plot_taxa_bars(tax_sum_mg,
                         nc_mg$map_loaded,
                         "sampleID",
                         num_taxa = 20,
                         data_only = TRUE) %>%
  left_join(., nc_mg$map_loaded, by = c("group_by" = "sampleID"))
topmg <- barsMG %>%
  group_by(taxon) %>%
  summarise(mean = mean(mean_value)) %>%
  arrange(-mean)
barsMG <- barsMG %>%
  mutate(taxon = factor(taxon, levels = topmg$taxon)) %>%
  mutate(taxon = fct_rev(taxon))
png("InitialFigs/Methanogens_mTAG.png", width = 7, height = 5, units = "in", res = 300)
ggplot(barsMG, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Sample", y = "% Abundance", fill = "Phylum; Class; Order; Family") +
  scale_fill_manual(values = brewer_pal(palette = "Paired")(10)) +
  scale_y_continuous(expand = c(0.01, 0.01)) +  
  facet_grid(~ Treatment, space = "free", scales = "free_x") +
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
        legend.key.size = unit(0.25, "cm"),
        legend.text = element_text(size = 6))
dev.off()



#### _Methanotrophs ####
nc_mt <- filter_taxa_from_input(nc,
                                taxa_to_keep = c("MOB_I", "MOB_II", "MOB_IIa"),
                                at_spec_level = 9)
tax_sum_mt <- summarize_taxonomy(input = nc_mt, 
                                 level = 6, 
                                 report_higher_tax = F, 
                                 relative = F) %>%
  mutate_all(funs((./4271)*100))
rownames(tax_sum_mt)[21] <- "Methyloligellaceae; uncultured"
barsMT <- plot_taxa_bars(tax_sum_mt,
                         nc_mt$map_loaded,
                         "sampleID",
                         num_taxa = 20,
                         data_only = TRUE) %>%
  mutate(taxon = gsub("NA", "Unclassified", taxon)) %>%
  mutate(taxon = fct_relevel(taxon, "Other", after = Inf)) %>%
  mutate(taxon = fct_relevel(taxon, "Unclassified", after = Inf)) %>%
  mutate(taxon = fct_rev(taxon)) %>%
  left_join(., nc_mt$map_loaded, by = c("group_by" = "sampleID"))
topmt <- barsMT %>%
  group_by(taxon) %>%
  summarise(mean = mean(mean_value)) %>%
  arrange(-mean)
barsMT <- barsMT %>%
  mutate(taxon = factor(taxon, levels = topmt$taxon)) %>%
  mutate(taxon = fct_relevel(taxon, "Other", after = Inf)) %>%
  mutate(taxon = fct_relevel(taxon, "Unclassified", after = Inf)) %>%
  mutate(taxon = fct_rev(taxon))
nb.cols <- 21
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)
png("InitialFigs/Methanotrophs_mTAG.png", width = 6.5, height = 6, units = "in", res = 300)
ggplot(barsMT, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Sample", y = "% Abundance", fill = "Genus") +
  scale_fill_manual(values = c("grey75", "grey90", mycolors)) +
  scale_y_continuous(expand = c(0.01, 0.01)) +  
  facet_grid(~ Treatment, space = "free", scales = "free_x") +
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
        legend.position = "right",
        #legend.justification = c(1,1),
        #legend.background = element_blank(),
        legend.key.size = unit(0.3, "cm"))
dev.off()



#### 6. Indicators ####
#### _Simper ####
nc_sim <- simper(t(nc$data_loaded), nc$map_loaded$Treatment)
nc_s <- summary(nc_sim)
nc_df1 <- head(nc_s$`Control_+SO4`, n = 10) %>%
  mutate(Comparison = "Control_+SO4",
         ASV = rownames(.)) %>%
  rename("MeanControl" = ava,
         "MeanTrt" = avb)
nc_df2 <- head(nc_s$`+ASW-SO4_Control`, n = 10) %>%
  mutate(Comparison = "Control_+ASW-SO4",
         ASV = rownames(.)) %>%
  rename("MeanControl" = avb,
         "MeanTrt" = ava)
nc_df3 <- head(nc_s$`+ASW_Control`, n = 10) %>%
  mutate(Comparison = "Control_+ASW",
         ASV = rownames(.)) %>%
  rename("MeanControl" = avb,
         "MeanTrt" = ava)
nc_simper_results <- rbind(nc_df1, nc_df2, nc_df3) %>%
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
  mutate(MeanTrt = round((MeanTrt/4271)*100, digits = 2),
         MeanControl = round((MeanControl/4271)*100, digits = 2),
         CumulativeContribution = round(CumulativeContribution*100, digits = 2)) %>%
  unite(Taxonomy, Phylum, Class, Order, Family, Genus, Species, OTU,
                          sep = "; ") %>%
  mutate(Taxonomy = make.unique(Taxonomy)) %>%
  dplyr::select(Comparison, Response, Domain, Guild, Taxonomy, MeanTrt,
                MeanControl, CumulativeContribution)

simper_meta <- nc_simper_results %>%
  column_to_rownames(var = "Taxonomy") %>%
  dplyr::select(Guild, Domain, Response, Comparison)
simper_mat <- nc_simper_results %>%
  column_to_rownames(var = "Taxonomy") %>%
  dplyr::select(MeanControl, MeanTrt, CumulativeContribution) %>%
  as.matrix()
ann_rows <- data.frame(row.names = rownames(simper_meta),
                       "Guild" = simper_meta$Guild,
                       "Domain" = simper_meta$Domain,
                       "Response" = simper_meta$Response,
                       "Comparison" = simper_meta$Comparison)
ann_colors <- list(Guild = c(AOA = "#7CFC00",
                             CH4_H2 = "#CD4F39",
                             `NA` = "white"),
                   Domain = c(Bacteria = "#440154FF",
                              Archaea = "#FCFDBFFF"),
                   Response = c(Positive = "red",
                                Negative = "blue"),
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
         gaps_row = c(10, 20),
         filename = "InitialFigs/Simper_mTAG.png",
         width = 8,
         height = 7)
dev.off()
dev.set(dev.next())
dev.set(dev.next())



#### _Multipatt ####
tax_sum_OTU <- summarize_taxonomy(nc, level = 8, report_higher_tax = F, relative = T)
set.seed(425)
mp <- multipatt(t(nc$data_loaded), 
                nc$map_loaded$Treatment, 
                func = "r.g", 
                control = how(nperm=999))
summary(mp) # 697 associated to 1 group

png("InitialFigs/Multipatt_mTAG.png", width = 6, height = 8, units = "in", res = 300)
plot_multipatt(mp_obj = mp, 
               input = nc,
               tax_sum = tax_sum_OTU,
               group = "Treatment",
               filter = TRUE,
               filter_vals = "Field",
               abund = "% Rel. Abund.",
               qcut = 0.1,
               rcut = 0)
dev.off()

# Redo with only abundant taxa. Looks like the others are rare
# Remember rarefied to 4271
# 0.05 percent is 4271*0.0005 = 2.1355
nrow(nc$taxonomy_loaded)
View(nc$data_loaded)
meancounts <- data.frame(meancount = sort(rowMeans(nc$data_loaded)))
nc_abund <- filter_taxa_from_input(nc, filter_thresh = 2.1355)
nrow(nc_abund$taxonomy_loaded)
sort(rowMeans(nc_abund$data_loaded)/4271*100)
# This is the top 331 taxa

tax_sum_OTU <- summarize_taxonomy(nc_abund, level = 8, report_higher_tax = F, relative = T)
sort(rowMeans(tax_sum_OTU)*100)
set.seed(425)
mp <- multipatt(t(nc_abund$data_loaded), 
                nc_abund$map_loaded$Treatment, 
                func = "r.g", 
                control = how(nperm=999))
summary(mp) # Number of species associated to 1 group: 81 

png("InitialFigs/Multipatt_abund_mTAG.png", width = 6, height = 8, units = "in", res = 300)
plot_multipatt(mp_obj = mp, 
               input = nc,
               tax_sum = tax_sum_OTU,
               group = "Treatment",
               filter = TRUE,
               filter_vals = "Field",
               abund = "% Rel. Abund.",
               qcut = 0.1,
               rcut = 0)
dev.off()



#### 7. BGC ####
# Don't need to repeat the BGC-only stuff
# But do run the BGC-microbe stuff
env_nona <- env_nona_nc %>%
  dplyr::select(-Cl_mgL)

# Plot correlations for different taxonomic levels at a given % relative abundance threshold
png("InitialFigs/CH4_Phyla_mTAG.png", width = 7, height = 5, units = "in", res = 300)
meth_corr_by_taxonomy(input = nc, level = 2, threshold = 0.5, data = "No")
dev.off()
meth_corr_by_taxonomy(input = nc, level = 3, threshold = 0.5, data = "No")
meth_corr_by_taxonomy(input = nc, level = 4, threshold = 0.5, data = "No")
meth_corr_by_taxonomy(input = nc, level = 5, threshold = 0.5, data = "No")
meth_corr_by_taxonomy(input = nc, level = 6, threshold = 0.5, data = "No")
meth_corr_by_taxonomy(input = nc, level = 8, threshold = 0.5, data = "No")
png("InitialFigs/CH4_Guilds_mTAG.png", width = 7, height = 5, units = "in", res = 300)
meth_corr_by_taxonomy(input = nc, level = 9, threshold = 0, data = "No")
dev.off()
# Note not as many sig phyla and no sig guilds, perhaps due to lower n because only D2



#### 8. Compare iTAG and mTAG ####
# Compare abundances of guilds and top phyla

# Input
mtag <- readRDS("input_filt_rare_mTAGs.rds")
itag <- readRDS("nc.rds")

# Get same samples
s <- mtag$map_loaded %>%
  filter(sampleID %in% itag$map_loaded$sampleID)
mtag <- filter_data(mtag,
                    filter_cat = "sampleID",
                    keep_vals = s$sampleID)
itag <- filter_data(itag,
                    filter_cat = "sampleID",
                    keep_vals = s$sampleID)

# Get same sample order
mtag$map_loaded <- mtag$map_loaded %>%
  arrange(sampleID)
mtag$data_loaded <- mtag$data_loaded %>%
  dplyr::select(rownames(mtag$map_loaded))
itag$map_loaded <- itag$map_loaded %>%
  arrange(sampleID)
itag$data_loaded <- itag$data_loaded %>%
  dplyr::select(rownames(itag$map_loaded))
sum(mtag$map_loaded$sampleID != itag$map_loaded$sampleID) # Good

# Phyla
Pmtag <- summarize_taxonomy(input = mtag, level = 2, report_higher_tax = F) %>%
  plot_taxa_bars(.,
                 mtag$map_loaded,
                 "sampleID",
                 num_taxa = 12,
                 data_only = TRUE) %>%
  filter(taxon != "Other") %>%
  pivot_wider(names_from = group_by,
              values_from = mean_value) %>%
  column_to_rownames(var = "taxon") %>%
  t() %>%
  as.data.frame()
Pitag <- summarize_taxonomy(input = itag, level = 2, report_higher_tax = F) %>%
  plot_taxa_bars(.,
                 itag$map_loaded,
                 "sampleID",
                 num_taxa = 20,
                 data_only = TRUE) %>%
  filter(taxon != "Other") %>%
  filter(taxon %in% names(Pmtag)) %>%
  pivot_wider(names_from = group_by,
              values_from = mean_value) %>%
  column_to_rownames(var = "taxon") %>%
  t() %>%
  as.data.frame()

sum(rownames(Pmtag) != rownames(Pitag)) # Good
names_p <- rev(names(Pmtag))
names_p[10] <- "RCP2.54"
names(Pmtag) <- paste0(names(Pmtag),"_m")
names(Pitag) <- paste0(names(Pitag),"_i")
Combined_phyla = data.frame(Pmtag, Pitag)
names(Combined_phyla)

# Function for plotting
compare_abund <- function(data, taxon) {
  data <- data %>%
    dplyr::select(paste0(taxon, "_i"), paste0(taxon, "_m"))
  ggplot(data, aes(x = data[,1], y = data[,2])) + 
    geom_point(size = 2) + 
    labs(x = NULL, y = NULL) +
    stat_smooth(method = "lm", se = F, size = 0.5, linetype = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", size = 0.2, color = "grey10") +
    ggtitle(taxon) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 10, vjust = -1),
          axis.text = element_text(size = 6))
}

# Run loop through the taxa
p <- list()
for (i in 1:length(names_p)) {
  p[[i]] <- compare_abund(data = Combined_phyla, taxon = names_p[i])
}

png("InitialFigs/CompareMethodsPhyla.png", width = 8, height = 6, units = "in", res = 300)
grid.arrange(p[[1]], p[[2]], p[[3]], p[[4]],
             p[[5]], p[[6]], p[[7]], p[[8]],
             p[[9]], p[[10]], p[[11]], p[[12]],
             ncol = 4, nrow = 3,
             bottom = "iTAG PCR amplicon", left = "mTAG metagenomic")
dev.off()

# Guilds
Gmtag <- summarize_taxonomy(input = mtag, level = 9, report_higher_tax = F) %>%
  plot_taxa_bars(.,
                 mtag$map_loaded,
                 "sampleID",
                 num_taxa = 12,
                 data_only = TRUE) %>%
  filter(taxon != "Other") %>%
  filter(taxon != "NA") %>%
  pivot_wider(names_from = group_by,
              values_from = mean_value) %>%
  column_to_rownames(var = "taxon") %>%
  t() %>%
  as.data.frame()
Gitag <- summarize_taxonomy(input = itag, level = 9, report_higher_tax = F) %>%
  plot_taxa_bars(.,
                 itag$map_loaded,
                 "sampleID",
                 num_taxa = 20,
                 data_only = TRUE) %>%
  filter(taxon != "Other") %>%
  filter(taxon != "NA") %>%
  filter(taxon %in% names(Gmtag)) %>%
  pivot_wider(names_from = group_by,
              values_from = mean_value) %>%
  column_to_rownames(var = "taxon") %>%
  t() %>%
  as.data.frame()

sum(rownames(Gmtag) != rownames(Gitag)) # Good
names_g <- rev(names(Gmtag))
names(Gmtag) <- paste0(names(Gmtag),"_m")
names(Gitag) <- paste0(names(Gitag),"_i")
Combined_guilds = data.frame(Gmtag, Gitag)
names(Combined_guilds)

# Run loop through the taxa
g <- list()
for (i in 1:length(names_g)) {
  g[[i]] <- compare_abund(data = Combined_guilds, taxon = names_g[i])
}

png("InitialFigs/CompareMethodsGuilds.png", width = 8, height = 6, units = "in", res = 300)
grid.arrange(g[[1]], g[[2]], g[[3]], g[[4]],
             g[[5]], g[[6]], g[[7]], g[[8]],
             g[[9]], g[[10]], g[[11]],
             ncol = 4, nrow = 3,
             bottom = "iTAG PCR amplicon", left = "mTAG metagenomic")
dev.off()



#### End Script ####