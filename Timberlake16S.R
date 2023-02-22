# East Coast/SF 16S data analysis
# by Cliff Bueno de Mesquita, Tringe Lab, JGI, Summer/Fall 2022



#### 1. Overview ####
# Samples from Weston, Neubauer, Ardón/Bernhardt labs
# Delaware River, Alligator River (NC), and Waccamaw River (SC) estuaries
# Two plates were sequenced and processed with iTagger at JGI
# Cliff reassigned taxonomy with SILVA v 138.1 using DADA2 in R

# Key Goals:
# Compare communities to Tringe Lab SF Bay/Delta samples
# Identify main methanogens and guilds
# Assess role of salinity, other environmental factors, and site
# See if microbial communities can explain discrepancies in CH4/salinity results

# Key papers/experimental information:
# Neubauer 2013 Estuaries and Coasts
## - South Carolina, experiment with Control, +salt, +fresh plots

# Ardón et al. 2013 Global Change Biology
## - North Carolina, experiment with hydrology, salt, sulfate microcosms
## 50 samples: 10 initial, 40 experimental (2 x 2 x 2 factorial)
## Control = DI water
## ASW = artificial salt-water with 5 ppt salinity
## ASW-SO4 = artificial salt-water with 5 ppt salinity without sulfate
## DI+SO4 = DI water with same sulfate as ASW
## Flooded samples (water level maintained at surface)
## There was also a drought treatment but it wasn't sequenced
## 30˚C, 12 weeks

# Weston et al. 2014 Biogeochemistry
## - Delaware River, field sampling
## Tidal freshwater, oligohaline, mesohaline

## - Delaware River, experiment
## Freshwater or ASW, 0, 4 wk, 7 wk, 12 wk sampling

## - Delaware River, transplants
## Transplanted TFM to another TFM, Oligohaline marsh, and Mesohaline marsh
## Done at surface and 40 cm below surface

# Hartman et al. in prep
## - SF Bay/Delta Estuary, field sampling
## - Various restored and reference wetlands ranging from freshwater to oligo/meso/poly

# Analysis Outline
# The analysis contains alpha and beta diversity analyses as well as taxonomic analyses
# Taxonomic analyses include simper and multipatt indicator analyses
# Stacked bar plots are made for each taxonomic level, including functional guilds
# Similar analysis repeated on different subsets of data
# Follow the Document Outline in RStudio (on the right) to navigate among different sections
# Sections are:
## 1. This overview with background information
## 2. Setup (libraries, metadata, guild calling, mctoolsr object, rarefaction)
## 3. East Coast Overview (all east coast samples)
## 4. East Coast Experiments (each each coast experiment individually)
## 5. Comparison Overview (all west coast and east coast samples)
## 6. Comparison Field Control (all unmanipulated field samples, west and east coasts)
## 7. Comparison Lab Inc (all laboratory incubations, Delaware and Alligator)
## 8. Comparison Field Exp (all field manipulations, Delaware and Waccamaw)



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

# Guild subsetting module from other repository
source("~/Documents/GitHub/SF_microbe_methane/modules/3_OTU_subsetting_modules_v.0.4_strip.r")

# Correlations
source("~/Documents/GitHub/EastCoast/meth_corr_by_taxonomy.R")
source("~/Documents/GitHub/EastCoast/meth_corr_by_bgc.R")

# Repository path
setwd("~/Documents/GitHub/EastCoast/")

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



#### _Delaware ####
# Make mapping file using the 2 Excel files Wyatt sent Tijana for iTag sequencing
# Note, already fixed "Freswater" to "Freshwater" typo in Excel
p1 <- read_excel("SPITS Wyatt 1520 itags4.xlsx")
p2 <- read_excel("SPITS Wyatt 1520 itags pl2v2.xlsx")
metadata <- rbind(p1, p2) %>%
  dplyr::select(`Sample Name*`, `Collection Year*`, `Collection Month*`, 
                `Collection Day*`, `Sample Isolated From*`, 
                `Collection Site or Growth Conditions`, `Latitude*`, 
                `Longitude*`, `Altitude or Depth*`) %>%
  set_names(c("sampleID", "Year", "Month", "Day", "Experiment", "Treatment", 
              "Latitude", "Longitude", "Depth")) %>%
  mutate(Estuary = "NA")
for (i in 1:nrow(metadata)) {
  if (metadata$Latitude[i] == 33.5250) {
    metadata$Estuary[i] <- "Waccamaw"
  }
  if (metadata$Latitude[i] == 35.9061) {
    metadata$Estuary[i] <- "Alligator"
  }
  if (metadata$Latitude[i] > 36) {
    metadata$Estuary[i] <- "Delaware"
  }
}
write.table(metadata, "delaware_metadata.txt", sep = "\t", row.names = F)

# Metadata has 184 samples
# ASV table has 177 samples
# 7 samples lost/not sequenced - check which ones
otu_table <- read.table("seqtab_wTax_mctoolsr.txt", header = 2)
missing <- metadata %>%
  filter(sampleID %notin% names(otu_table))
write.table(missing, "not_sequenced.txt", sep = "\t", row.names = F)

# Import Data (n = 177), Filter, Rarefy, Calc richness
tax_table_fp <- "seqtab_wTax_mctoolsr.txt"
map_fp <- "delaware_metadata.txt"
input = load_taxa_table(tax_table_fp, map_fp)

# Filter chloroplast, mitochondria, eukaryotes, unassigned at domain
input_filt <- filter_taxa_from_input(input,
                                     taxa_to_remove = "Chloroplast") # 292 removed
input_filt <- filter_taxa_from_input(input_filt,
                                     taxa_to_remove = "Mitochondria") # 787 removed
input_filt <- filter_taxa_from_input(input_filt,
                                     taxa_to_remove = "Eukarya") # none
input_filt <- filter_taxa_from_input(input_filt,
                                     taxa_to_remove = "NA",
                                     at_spec_level = 1) # 34 removed

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
Guild_OTUs <- Get_16S_Guilds_alt(reformed_OTU_table)
levels(as.factor(Guild_OTUs$Guild))

# Now add as 9th column to input_filt$taxonomy_loaded
input_filt$taxonomy_loaded <- input_filt$taxonomy_loaded %>%
  left_join(., Guild_OTUs, by = c("taxonomy8" = "OTU")) %>%
  rename(taxonomy9 = Guild) %>%
  mutate(taxonomy9 = as.character(taxonomy9)) %>%
  mutate(taxonomy9 = replace_na(taxonomy9, "NA"))
rownames(input_filt$taxonomy_loaded) <- input_filt$taxonomy_loaded$taxonomy8

# Check ANME as this was updated
ANME <- reformed_OTU_table[grepl("(ANME|Methanoperedenaceae|Syntrophoarchaeaceae)",
                                 reformed_OTU_table$Consensus.lineage),]     
ANME <- Subs_to_DF(ANME)
ANME["Guild"] <- "ANME"
ANME.2 <- input_filt$taxonomy_loaded[grepl("ANME",
                                           input_filt$taxonomy_loaded$taxonomy9),]
ANME$OTU %in% ANME.2$taxonomy8

# Save
saveRDS(input_filt, "input_filt.rds")

# Rarefy at minimum (26429)
sort(colSums(input_filt$data_loaded))
mean(colSums(input_filt$data_loaded))
se(colSums(input_filt$data_loaded))
# Depth 92024 ± 2919
set.seed(530)
input_filt_rare <- single_rarefy(input_filt, 26429)
sort(colSums(input_filt_rare$data_loaded))

# OTU Richness
input_filt_rare$map_loaded$rich <- specnumber(input_filt_rare$data_loaded, 
                                              MARGIN = 2)

# Shannon diversity
input_filt_rare$map_loaded$shannon <- diversity(input_filt_rare$data_loaded, 
                                                index = "shannon", 
                                                MARGIN = 2)

# Save
saveRDS(input_filt_rare, "input_filt_rare.rds")



#### _Combined ####
# Make combined metadata table
de <- read.delim("delaware_metadata.txt") %>%
  dplyr::select(sampleID, Experiment, Treatment, Depth, Estuary) %>%
  set_names(c("sampleID", "Site", "Detail", "Depth", "Estuary")) %>%
  mutate(Info = Detail)
sf <- read.delim("~/Desktop/Wyatt Manuscript/metadata_mctoolsr.txt") %>%
  dplyr::select(SampleID, Location, Vegetation, Depth) %>%
  set_names(c("sampleID", "Site", "Detail", "Depth")) %>%
  mutate(Estuary = "SF") %>%
  mutate(Info = Site)
comb <- rbind(de, sf)

# Metadata has 352 samples
# ASV table has 345 samples
# 7 lost from Delaware

write.table(comb, "combined_metadata.txt", sep = "\t", row.names = F)

# Import Data (n = 345), Filter, Rarefy, Calc richness
tax_table_fp <- "seqtab_wTax_mctoolsr_comb.txt"
map_fp <- "combined_metadata.txt"
input = load_taxa_table(tax_table_fp, map_fp)

# Filter chloroplast, mitochondria, eukaryotes, unassigned at domain
input_filt <- filter_taxa_from_input(input,
                                     taxa_to_remove = "Chloroplast") # 368 removed
input_filt <- filter_taxa_from_input(input_filt,
                                     taxa_to_remove = "Mitochondria") # 815 removed
input_filt <- filter_taxa_from_input(input_filt,
                                     taxa_to_remove = "Eukarya") # none
input_filt <- filter_taxa_from_input(input_filt,
                                     taxa_to_remove = "NA",
                                     at_spec_level = 1) # 54 removed

# Guilds
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
Guild_OTUs <- Get_16S_Guilds_alt(reformed_OTU_table)
levels(as.factor(Guild_OTUs$Guild))

# Now add as 9th column to input_filt$taxonomy_loaded
input_filt$taxonomy_loaded <- input_filt$taxonomy_loaded %>%
  left_join(., Guild_OTUs, by = c("taxonomy8" = "OTU")) %>%
  rename(taxonomy9 = Guild) %>%
  mutate(taxonomy9 = as.character(taxonomy9)) %>%
  mutate(taxonomy9 = replace_na(taxonomy9, "NA"))
rownames(input_filt$taxonomy_loaded) <- input_filt$taxonomy_loaded$taxonomy8

# Check ANME as this was updated
ANME <- reformed_OTU_table[grepl("(ANME|Methanoperedenaceae|Syntrophoarchaeaceae)",
                                 reformed_OTU_table$Consensus.lineage),]     
ANME <- Subs_to_DF(ANME)
ANME["Guild"] <- "ANME"
ANME.2 <- input_filt$taxonomy_loaded[grepl("ANME",
                                           input_filt$taxonomy_loaded$taxonomy9),]
ANME$OTU %in% ANME.2$taxonomy8

# OTU Richness
input_filt$map_loaded$rich <- specnumber(input_filt$data_loaded, 
                                         MARGIN = 2)

# Shannon diversity
input_filt$map_loaded$shannon <- diversity(input_filt$data_loaded, 
                                           index = "shannon", 
                                           MARGIN = 2)

# Save
saveRDS(input_filt, "input_filt_comb.rds")

# Rarefy at 26429
sort(colSums(input_filt$data_loaded))
mean(colSums(input_filt$data_loaded))
se(colSums(input_filt$data_loaded))
# Original depth 118864 ± 2323
# Drop Sandmound_TuleB_D1 (1296 reads) and Muzzi_PWB_D2 (5287 reads)
set.seed(530)
input_filt_rare <- single_rarefy(input_filt, 26429) # Now n = 343
sort(colSums(input_filt_rare$data_loaded))

# OTU Richness
input_filt_rare$map_loaded$rich <- specnumber(input_filt_rare$data_loaded, 
                                              MARGIN = 2)

# Shannon diversity
input_filt_rare$map_loaded$shannon <- diversity(input_filt_rare$data_loaded, 
                                                index = "shannon", 
                                                MARGIN = 2)

# Save
saveRDS(input_filt_rare, "input_filt_rare_comb.rds")

#### _Biogeochemistry ####
# Received biogeochemical data later
# Wrangled it in Excel and in PrepBiogeochem.R
# Got sampleID to each row so can merge here
# Reinput the old files and save new files with BGC data
# Make sure to set rownames! dplyr loses them...
biogeochem <- read.csv("biogeochem_all_clean.csv") %>%
  dplyr::select(-X, -Estuary, -Salinity_calcd_ppt, -Salinity)

input_filt_rare <- readRDS("input_filt_rare.rds")
nrow(input_filt_rare$map_loaded) # 177
input_filt_rare_wBGC <- input_filt_rare
input_filt_rare_wBGC$map_loaded <- input_filt_rare_wBGC$map_loaded %>%
  mutate(sampleID = rownames(.)) %>%
  left_join(., biogeochem, by = "sampleID")
rownames(input_filt_rare_wBGC$map_loaded) <- input_filt_rare_wBGC$map_loaded$sampleID
nrow(input_filt_rare_wBGC$map_loaded) # 177
saveRDS(input_filt_rare_wBGC, "input_filt_rare_wBGC.rds")

input_filt <- readRDS("input_filt.rds")
nrow(input_filt$map_loaded)
input_filt_wBGC <- input_filt
input_filt_wBGC$map_loaded <- input_filt_wBGC$map_loaded %>%
  mutate(sampleID = rownames(.)) %>%
  left_join(., biogeochem, by = "sampleID")
rownames(input_filt_wBGC$map_loaded) <- input_filt_wBGC$map_loaded$sampleID
nrow(input_filt_wBGC$map_loaded) # 177
saveRDS(input_filt_wBGC, "input_filt_wBGC.rds")

input_filt_comb <- readRDS("input_filt_comb.rds")
for (i in 1:nrow(input_filt_comb$map_loaded)) {
if (rownames(input_filt_comb$map_loaded)[i] == "Sandmound_Tule_C_D1") {
  rownames(input_filt_comb$map_loaded)[i] <- "Sandmound_TuleC_D1"
 }
}
for (i in 1:nrow(input_filt_comb$map_loaded)) {
  if (rownames(input_filt_comb$map_loaded)[i] == "Sandmound_Tule_C_D2") {
    rownames(input_filt_comb$map_loaded)[i] <- "Sandmound_TuleC_D2"
  }
}
for (i in 1:ncol(input_filt_comb$data_loaded)) {
  if (colnames(input_filt_comb$data_loaded)[i] == "Sandmound_Tule_C_D1") {
    colnames(input_filt_comb$data_loaded)[i] <- "Sandmound_TuleC_D1"
  }
}
for (i in 1:ncol(input_filt_comb$data_loaded)) {
  if (colnames(input_filt_comb$data_loaded)[i] == "Sandmound_Tule_C_D2") {
    colnames(input_filt_comb$data_loaded)[i] <- "Sandmound_TuleC_D2"
  }
}
nrow(input_filt_comb$map_loaded)
input_filt_comb_wBGC <- input_filt_comb
input_filt_comb_wBGC$map_loaded <- input_filt_comb_wBGC$map_loaded %>%
  mutate(sampleID = rownames(.)) %>%
  left_join(., biogeochem, by = "sampleID")
rownames(input_filt_comb_wBGC$map_loaded) <- input_filt_comb_wBGC$map_loaded$sampleID
nrow(input_filt_comb_wBGC$map_loaded) # 345
saveRDS(input_filt_comb_wBGC, "input_filt_comb_wBGC.rds")

input_filt_rare_comb <- readRDS("input_filt_rare_comb.rds")
for (i in 1:nrow(input_filt_rare_comb$map_loaded)) {
  if (rownames(input_filt_rare_comb$map_loaded)[i] == "Sandmound_Tule_C_D1") {
    rownames(input_filt_rare_comb$map_loaded)[i] <- "Sandmound_TuleC_D1"
  }
}
for (i in 1:nrow(input_filt_rare_comb$map_loaded)) {
  if (rownames(input_filt_rare_comb$map_loaded)[i] == "Sandmound_Tule_C_D2") {
    rownames(input_filt_rare_comb$map_loaded)[i] <- "Sandmound_TuleC_D2"
  }
}
for (i in 1:ncol(input_filt_rare_comb$data_loaded)) {
  if (colnames(input_filt_rare_comb$data_loaded)[i] == "Sandmound_Tule_C_D1") {
    colnames(input_filt_rare_comb$data_loaded)[i] <- "Sandmound_TuleC_D1"
  }
}
for (i in 1:ncol(input_filt_rare_comb$data_loaded)) {
  if (colnames(input_filt_rare_comb$data_loaded)[i] == "Sandmound_Tule_C_D2") {
    colnames(input_filt_rare_comb$data_loaded)[i] <- "Sandmound_TuleC_D2"
  }
}
nrow(input_filt_rare_comb$map_loaded)
input_filt_rare_comb_wBGC <- input_filt_rare_comb
input_filt_rare_comb_wBGC$map_loaded <- input_filt_rare_comb_wBGC$map_loaded %>%
  mutate(sampleID = rownames(.)) %>%
  left_join(., biogeochem, by = "sampleID")
rownames(input_filt_rare_comb_wBGC$map_loaded) <- input_filt_rare_comb_wBGC$map_loaded$sampleID
nrow(input_filt_rare_comb_wBGC$map_loaded) # 177
saveRDS(input_filt_rare_comb_wBGC, "input_filt_rare_comb_wBGC.rds")



#### ...................................... ####
#### 3. East Coast Overview ####
# Look at just the East Coast Data, but altogether
# Using the default "Treatment" category to group for first look
# Will change in the combined and detailed analyses further down
# So, don't save figures except for the maps
# Showing top 10 taxa, here, later will show top 12 and format Other and Unclassified
input_filt_rare <- readRDS("input_filt_rare.rds")



#### _Map ####
# Summarize data, get unique coordinates
coords <- input_filt_rare$map_loaded %>%
  group_by(Latitude, Longitude) %>%
  slice_head()
min(coords$Latitude)
max(coords$Latitude)
min(coords$Longitude)
max(coords$Longitude)

del <- get_stamenmap(bbox = c(left = min(coords$Longitude) - 0.5, 
                              bottom = min(coords$Latitude) - 0.5, 
                              right = max(coords$Longitude) + 0.5, 
                              top = max(coords$Latitude) + 0.5),
                     zoom = 10, 
                     maptype = "terrain-background")
del_attributes <- attributes(del)
del_transparent <- matrix(adjustcolor(del, alpha.f = 0.4), nrow = nrow(del))
attributes(del_transparent) <- del_attributes

pdf("InitialFigs/Map_EastCoast.pdf", width = 6, height = 6)
ggmap(del_transparent, extent = "device") + # the base map
  geom_point(data = coords,
             aes(x = Longitude, y = Latitude), size = 4) +
  xlab(NULL) + 
  ylab(NULL) +
  theme(legend.position = "none",
        plot.margin = unit(c(0,-1,0,-1), "cm"),
        axis.text = element_text(size = 8, color = "black"))
dev.off()

# There are 5 Delaware River sites, one North Carolina site and one South Carolina site
# Zoom in on Delaware River
coords_del_only <- input_filt_rare$map_loaded %>%
  group_by(Latitude, Longitude) %>%
  slice_head() %>%
  filter(Latitude > 37)
del_only <- get_stamenmap(bbox = c(left = min(coords_del_only$Longitude) - 0.25, 
                                   bottom = min(coords_del_only$Latitude) - 0.2, 
                                   right = max(coords_del_only$Longitude) + 0.1, 
                                   top = max(coords_del_only$Latitude) + 0.1),
                     zoom = 10, 
                     maptype = "terrain-background")
del_only_attributes <- attributes(del_only)
del_only_transparent <- matrix(adjustcolor(del_only, alpha.f = 0.4), 
                               nrow = nrow(del_only))
attributes(del_only_transparent) <- del_only_attributes
pdf("InitialFigs/Map_Delaware.pdf", width = 6, height = 6)
ggmap(del_only_transparent, extent = "device") + # the base map
  geom_point(data = coords_del_only,
             aes(x = Longitude, y = Latitude), size = 4) +
  geom_text(aes(x = -75.2, y = 39.95, label = "Philadelphia"), 
            colour = "grey40", size = 3, check_overlap = T) +
  geom_text(aes(x = -75.56, y = 39.75, label = "Wilmington"), 
            colour = "grey40", size = 3, check_overlap = T) +
  geom_text(aes(x = -75.39, y = 39.33, label = "Delaware River"), 
            colour = "white", angle = 315, size = 4, fontface = "bold", 
            check_overlap = T) +
  geom_segment(aes(x = -74.8, xend = -74.8, y = 39.25, yend = 39.28), 
               arrow = arrow(length = unit(0.30, "cm"))) +
  geom_text(aes(x = -74.8, y = 39.3, label = "N"), 
            colour = "black", size = 4, check_overlap = T) +
  xlab(NULL) + 
  ylab(NULL) +
  scalebar(x.min = min(coords_del_only$Longitude), 
           y.min = min(coords_del_only$Latitude)-0.15, 
           x.max = max(coords_del_only$Longitude)-0.05, 
           y.max = max(coords_del_only$Latitude), 
           dist = 10, dist_unit = "km", height = 0.02, st.dist = 0.03, st.size = 4,
           transform = TRUE, model = "WGS84", location = "bottomright") +
  theme(legend.position = "none",
        plot.margin = unit(c(0,-1,0,-1), "cm"),
        axis.text = element_text(size = 8, color = "black"))
dev.off()



#### _Alpha ####
# OTU Richness
leveneTest(input_filt_rare$map_loaded$rich ~ input_filt_rare$map_loaded$Treatment)
# Variance homogeneous (p > 0.05)
m <- aov(input_filt_rare$map_loaded$rich ~ input_filt_rare$map_loaded$Treatment)
shapiro.test(m$residuals)
# Residuals not normally distributed (p < 0.05)
summary(m)

ggplot(input_filt_rare$map_loaded, aes(Treatment, rich, colour = Experiment)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 2, alpha = 0.5, width = 0.25) +
  labs(x = "Site", y = "Number of OTUs", colour = "Experiment") +
  facet_wrap(~ Experiment, ncol = 3, scales = "free") +
  scale_colour_viridis_d() +
  theme_bw() +
  theme(legend.position = c(0.75, 0.2),
        axis.title = element_text(face = "bold", size = 14),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))

# Shannon
leveneTest(input_filt_rare$map_loaded$shannon ~ input_filt_rare$map_loaded$Treatment)
# Variance homogeneous (p > 0.05)
m1 <- aov(input_filt_rare$map_loaded$shannon ~ input_filt_rare$map_loaded$Treatment)
shapiro.test(m1$residuals)
# Residuals not normally distributed (p < 0.05)
summary(m1)

ggplot(input_filt_rare$map_loaded, aes(Treatment, shannon, colour = Experiment)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 2, alpha = 0.5, width = 0.25) +
  labs(x = "Treatment", y = "Shannon diversity", colour = "Experiment") +
  facet_wrap(~ Experiment, ncol = 3, scales = "free") +
  scale_colour_viridis_d() +
  theme_bw() +
  theme(legend.position = c(0.75, 0.2),
        axis.title = element_text(face = "bold", size = 14),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))



#### _Beta  ####
bc <- calc_dm(input_filt_rare$data_loaded)
pcoa <- cmdscale(bc, k = nrow(input_filt_rare$map_loaded) - 1, eig = T)
pcoaA1 <- round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[1]*100, digits = 1)
pcoaA2 <- round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[2]*100, digits = 1)
input_filt_rare$map_loaded$Axis01 <- scores(pcoa)[,1]
input_filt_rare$map_loaded$Axis02 <- scores(pcoa)[,2]
micro.hulls <- ddply(input_filt_rare$map_loaded, c("Experiment"), find_hull)
ggplot(input_filt_rare$map_loaded, aes(Axis01, Axis02, colour = Experiment)) +
  geom_polygon(data = micro.hulls, aes(colour = Experiment, fill = Experiment),
               alpha = 0.1, show.legend = F) +
  geom_point(size = 2, alpha = 0.5, aes(shape = Estuary)) +
  labs(x = paste("PC1: ", pcoaA1, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2, "%", sep = ""),
       colour = "Experiment") +
  scale_colour_viridis_d() +
  scale_fill_viridis_d() +
  theme_bw() +  
  theme(legend.position = "right",
        axis.title = element_text(face = "bold", size = 16), 
        axis.text = element_text(size = 14))

# Stats
set.seed(1150)
input_filt_rare$map_loaded$Estuary <- as.factor(input_filt_rare$map_loaded$Estuary)
adonis2(bc ~ Estuary + Experiment + Depth, data = input_filt_rare$map_loaded)
anova(betadisper(bc, input_filt_rare$map_loaded$Estuary)) # Dispersion not homogeneous
anova(betadisper(bc, input_filt_rare$map_loaded$Experiment)) # Dispersion not homogeneous
anova(betadisper(bc, input_filt_rare$map_loaded$Depth)) # Dispersion not homogeneous



#### _Taxa ####
# Prelim exploration but don't save anything. 
# Will redo with SF Bay data included and save figures

#### _Indicators ####
sim <- simper(t(input_filt_rare$data_loaded), 
              input_filt_rare$map_loaded$Experiment)
s <- summary(sim)
head(s$`Soil incubation_Soil Incubation`, n = 10)
head(s$`Soil_Soil Field plots`, n = 10)

# MULTIPATT (list ASVs associated with each group)
set.seed(1202)
mp <- multipatt(t(input_filt_rare$data_loaded), 
                input_filt_rare$map_loaded$Experiment, 
                func = "IndVal.g", 
                control = how(nperm=999))
summary(mp)

#### _Domain ####
tax_sum_domain <- summarize_taxonomy(input_filt_rare, level = 1, 
                                     report_higher_tax = F)
plot_ts_heatmap(tax_sum_domain, 
                input_filt_rare$map_loaded, 
                0.01, 
                'Experiment',
                rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
bars <- plot_taxa_bars(tax_sum_domain,
                       input_filt_rare$map_loaded,
                       "Experiment",
                       num_taxa = 10,
                       data_only = TRUE)
ggplot(bars, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Experiment", y = "Relative abundance", fill = "Domain") +
  scale_fill_manual(values = c("red", "blue")) +
  scale_y_continuous(expand = c(0.01, 0.01)) + 
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 12), 
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
taxa_summary_by_sample_type(tax_sum_domain, 
                            input_filt_rare$map_loaded, 
                            type_header = 'Experiment', 
                            filter_level = 0.01, 
                            test_type = 'KW')

#### _Phylum ####
tax_sum_phyla <- summarize_taxonomy(input_filt_rare, level = 2, 
                                    report_higher_tax = F)
plot_ts_heatmap(tax_sum_phyla, 
                input_filt_rare$map_loaded, 
                0.01, 
                'Experiment',
                rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
bars <- plot_taxa_bars(tax_sum_phyla,
                       input_filt_rare$map_loaded,
                       "Experiment",
                       num_taxa = 10,
                       data_only = TRUE)
ggplot(bars, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Experiment", y = "Relative abundance", fill = "Phylum") +
  scale_fill_brewer(palette = "Paired") +
  scale_y_continuous(expand = c(0.01, 0.01)) + 
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 12), 
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
taxa_summary_by_sample_type(tax_sum_phyla, 
                            input_filt_rare$map_loaded, 
                            type_header = 'Experiment', 
                            filter_level = 0.01, 
                            test_type = 'KW')

# Look at archaeal phyla
tax_sum_phyla_ar <- summarize_taxonomy(input_filt_rare, level = 2, 
                                       report_higher_tax = T)
tax_sum_phyla_ar <- tax_sum_phyla_ar[grep("Archaea", rownames(tax_sum_phyla_ar)),]
bars_ar <- plot_taxa_bars(tax_sum_phyla_ar,
                          input_filt_rare$map_loaded,
                          "Experiment",
                          num_taxa = 13,
                          data_only = TRUE)
nb.cols <- nrow(tax_sum_phyla_ar)
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)
ggplot(bars_ar, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Experiment", y = "Relative abundance", fill = "Phylum") +
  scale_fill_manual(values = mycolors) +
  scale_y_continuous(expand = c(0.01, 0.01)) + 
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 12), 
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))

#### _Class ####
tax_sum_class <- summarize_taxonomy(input_filt_rare, level = 4, report_higher_tax = F)
plot_ts_heatmap(tax_sum_class, 
                input_filt_rare$map_loaded, 
                0.01, 
                'Estuary',
                rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
bars <- plot_taxa_bars(tax_sum_class,
                       input_filt_rare$map_loaded,
                       "Estuary",
                       num_taxa = 10,
                       data_only = TRUE)
ggplot(bars, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Estuary", y = "Relative abundance", fill = "Class") +
  scale_fill_brewer(palette = "Paired") +
  scale_y_continuous(expand = c(0.01, 0.01)) + 
  theme_classic() +
  theme(axis.title = element_text(face = "bold", size = 12), 
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
taxa_summary_by_sample_type(tax_sum_class, 
                            input_filt_rare$map_loaded, 
                            type_header = 'Estuary', 
                            filter_level = 0.01, 
                            test_type = 'KW')

# Check sulfur reducers (not all, just quick Desulfo check)
tax_sum_class_su <- summarize_taxonomy(input_filt_rare, level = 3, 
                                       report_higher_tax = T)
tax_sum_class_su <- tax_sum_class_su[grep("Desulfo", rownames(tax_sum_class_su)),]

bars_su <- plot_taxa_bars(tax_sum_class_su,
                          input_filt_rare$map_loaded,
                          "Experiment",
                          num_taxa = nrow(tax_sum_class_su),
                          data_only = TRUE) %>%
  mutate(taxon = substring(taxon, 11))
tallest_bar <- bars_su %>%
  group_by(Estuary, Salt) %>%
  summarise(sum = sum(mean_value))
nb.cols <- nrow(tax_sum_class_su)
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)
ggplot(bars_su, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Experiment", y = "Relative abundance", fill = "Class") +
  scale_fill_manual(values = mycolors) +
  scale_y_continuous(expand = c(max(tallest_bar$sum)/100, max(tallest_bar$sum)/100)) + 
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 12), 
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))

#### _Order ####
tax_sum_order <- summarize_taxonomy(input_filt_rare, level = 5, 
                                       report_higher_tax = FALSE)
plot_ts_heatmap(tax_sum_order, 
                input_filt_rare$map_loaded, 
                0.01, 
                'Experiment',
                rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
bars <- plot_taxa_bars(tax_sum_order,
                       input_filt_rare$map_loaded,
                       "Experiment",
                       num_taxa = 10,
                       data_only = TRUE)
ggplot(bars, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Experiment", y = "Relative abundance", fill = "Family") +
  scale_fill_brewer(palette = "Paired") +
  scale_y_continuous(expand = c(0.01, 0.01)) + 
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 12), 
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
taxa_summary_by_sample_type(tax_sum_order, 
                            input_filt_rare$map_loaded, 
                            type_header = 'Experiment', 
                            filter_level = 0.01, 
                            test_type = 'KW')

#### _Family ####
tax_sum_families <- summarize_taxonomy(input_filt_rare, level = 5, 
                                       report_higher_tax = FALSE)
plot_ts_heatmap(tax_sum_families, 
                input_filt_rare$map_loaded, 
                0.01, 
                'Experiment',
                rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
bars <- plot_taxa_bars(tax_sum_families,
                       input_filt_rare$map_loaded,
                       "Experiment",
                       num_taxa = 10,
                       data_only = TRUE)
ggplot(bars, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Experiment", y = "Relative abundance", fill = "Family") +
  scale_fill_brewer(palette = "Paired") +
  scale_y_continuous(expand = c(0.01, 0.01)) + 
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 12), 
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
taxa_summary_by_sample_type(tax_sum_families, 
                            input_filt_rare$map_loaded, 
                            type_header = 'Experiment', 
                            filter_level = 0.01, 
                            test_type = 'KW')

# Look at methanogens (careful to do this accurately!)
tax_sum_families_meth <- summarize_taxonomy(input_filt_rare, level = 5, 
                                            report_higher_tax = F)
tax_sum_families_meth <- tax_sum_families_meth[grep("(Methano|Methermicoccaceae|
                                                    Syntrophoarchaeaceae)", 
                                              rownames(tax_sum_families_meth)),]
tax_sum_families_meth <- tax_sum_families_meth[!grepl("Methanoperedenaceae", 
                                              rownames(tax_sum_families_meth)),]
bars_meth <- plot_taxa_bars(tax_sum_families_meth,
                            input_filt_rare$map_loaded,
                            "Experiment",
                            num_taxa = nrow(tax_sum_families_meth),
                            data_only = TRUE)
nb.cols <- nrow(tax_sum_families_meth)
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)
tallest_bar <- bars_meth %>%
  group_by(Estuary, Salt) %>%
  summarise(sum = sum(mean_value))
ggplot(bars_meth, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Experiment", y = "Relative abundance", fill = "Genus") +
  scale_fill_manual(values = mycolors) +
  scale_y_continuous(expand = c(max(tallest_bar$sum)/100, max(tallest_bar$sum)/100)) + 
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 12), 
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))

#### _Genus ####
tax_sum_genera <- summarize_taxonomy(input_filt_rare, level = 6, report_higher_tax = T)
plot_ts_heatmap(tax_sum_genera, 
                input_filt_rare$map_loaded, 
                0.01, 
                'Experiment',
                rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
bars <- plot_taxa_bars(tax_sum_genera,
                       input_filt_rare$map_loaded,
                       "Experiment",
                       num_taxa = 10,
                       data_only = TRUE)
ggplot(bars, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Experiment", y = "Relative abundance", fill = "Genus") +
  scale_fill_brewer(palette = "Paired") +
  scale_y_continuous(expand = c(0.01, 0.01)) + 
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 12), 
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
taxa_summary_by_sample_type(tax_sum_genera, 
                            input_filt_rare$map_loaded, 
                            type_header = 'Experiment', 
                            filter_level = 0.01, 
                            test_type = 'KW')

#### _Guilds ####
tax_sum_guilds <- summarize_taxonomy(input_filt_rare, level = 9, report_higher_tax = F)
plot_ts_heatmap(tax_sum_guilds, 
                input_filt_rare$map_loaded, 
                0, 
                'Experiment',
                rev_taxa = T,
                remove_other = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
bars <- plot_taxa_bars(tax_sum_guilds,
                       input_filt_rare$map_loaded,
                       "Experiment",
                       num_taxa = nrow(tax_sum_guilds),
                       data_only = TRUE) %>%
  filter(taxon != "NA") %>%
  droplevels() %>%
  mutate(taxon = factor(taxon,
                        levels = Guild_cols$Guild))
tallest_bar <- bars %>%
  group_by(Estuary, Salt) %>%
  summarise(sum = sum(mean_value))
ggplot(bars, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Experiment", y = "Relative abundance", fill = "Guild") +
  scale_fill_manual(values = Guild_cols$color) +
  scale_y_continuous(expand = c(max(tallest_bar$sum)/100, max(tallest_bar$sum)/100)) + 
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 12), 
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
taxa_summary_by_sample_type(tax_sum_guilds, 
                            input_filt_rare$map_loaded, 
                            type_header = 'Experiment', 
                            filter_level = 0.01, 
                            test_type = 'KW')

#### _Venn ####
phy <- summarize_taxonomy(input_filt_rare, level = 2, report_higher_tax = F)
cla <- summarize_taxonomy(input_filt_rare, level = 3, report_higher_tax = F)
ord <- summarize_taxonomy(input_filt_rare, level = 4, report_higher_tax = F)
fam <- summarize_taxonomy(input_filt_rare, level = 5, report_higher_tax = F)
gen <- summarize_taxonomy(input_filt_rare, level = 6, report_higher_tax = F)

input_phylum <- input_filt_rare
input_phylum$data_loaded <- phy
input_class <- input_filt_rare
input_class$data_loaded <- cla
input_order <- input_filt_rare
input_order$data_loaded <- ord
input_family <- input_filt_rare
input_family$data_loaded <- fam
input_genus <- input_filt_rare
input_genus$data_loaded <- gen

plot_venn_diagram(input_filt_rare,
                  "Estuary",
                  0.00000000000000001)

plot_grid(plot_venn_diagram(input_phylum, "Estuary", 0.00000000000000001),
          plot_venn_diagram(input_class, "Estuary", 0.00000000000000001),
          plot_venn_diagram(input_order, "Estuary", 0.00000000000000001),
          plot_venn_diagram(input_family, "Estuary", 0.00000000000000001),
          plot_venn_diagram(input_genus, "Estuary", 0.00000000000000001),
          plot_venn_diagram(input_filt_rare, "Estuary", 0.00000000000000001),
          labels = c("(a) Phylum", "(b) Class", "(c) Order", 
                     "(d) Family", "(e) Genus", "(f) OTU"))



#### ...................................... ####
#### 4. East Coast Experiments ####
# Analyze each of the East Coast experiments separately
# Look at microbial responses to the different treatments/time points/depths
# Here plots are annotated with text, later facet_wrap is used



#### _South Carolina ####
# input_filt <- readRDS("input_filt.rds")
input_filt <- readRDS("input_filt_wBGC.rds")
sc <- filter_data(input_filt,
                  filter_cat = "Estuary",
                  keep_vals = "Waccamaw")
set.seed(530)
sc <- single_rarefy(sc, min(colSums(sc$data_loaded))) # 79282
sc$map_loaded <- sc$map_loaded %>%
  mutate(rich = specnumber(sc$data_loaded, MARGIN = 2),
         shannon = diversity(sc$data_loaded, index = "shannon", MARGIN = 2),
         TrtDepth = paste(sc$map_loaded$Treatment, sc$map_loaded$Depth, sep = ""),
         Depth = as.factor(Depth)) %>%
  mutate_if(is.character, as.factor)

#### __Alpha ####
leveneTest(sc$map_loaded$rich ~ sc$map_loaded$Treatment)
m <- aov(rich ~ Treatment + Depth, data = sc$map_loaded)
Anova(m, type = "II") # Treatment, not Depth
m <- aov(rich ~ Treatment, data = sc$map_loaded)
shapiro.test(m$residuals)
summary(m)
t <- emmeans(object = m, specs = "Treatment") %>%
  cld(object = ., adjust = "Tukey", Letters = letters, alpha = 0.05) %>%
  mutate(name = "rich",
         y = max(sc$map_loaded$rich) +
           (max(sc$map_loaded$rich) - min(sc$map_loaded$rich))/20)
leveneTest(sc$map_loaded$shannon ~ sc$map_loaded$Treatment)
m1 <- aov(shannon ~ Treatment + Depth, data = sc$map_loaded)
Anova(m1, type = "II") # Treatment, not Depth
m1 <- aov(shannon ~ Treatment, data = sc$map_loaded)
shapiro.test(m1$residuals)
summary(m1)
t1 <- emmeans(object = m1, specs = "Treatment") %>%
  cld(object = ., adjust = "Tukey", Letters = letters, alpha = 0.05) %>%
  mutate(name = "shannon",
         y = max(sc$map_loaded$shannon) + 
           (max(sc$map_loaded$shannon) - min(sc$map_loaded$shannon))/20)
label_df <- rbind(t, t1)
facet_df <- c("rich" = "(a) Richness",
              "shannon" = "(b) Shannon")
alpha_long <- sc$map_loaded %>%
  pivot_longer(cols = c("rich", "shannon"))
pdf("InitialFigs/SC_Alpha.pdf", width = 6, height = 3)
ggplot(alpha_long, aes(reorder(Treatment, value, mean), value, 
                       colour = Treatment)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 2, alpha = 0.75, width = 0.2, aes(shape = Depth)) +
  geom_text(data = label_df, aes(Treatment, y, label = str_trim(.group)), 
            size = 4, color = "black") +
  labs(x = "Site", y = "Number of OTUs", shape = "Depth (m)") +
  scale_x_discrete(labels = c("+Saltwater", "+Freshwater", "Control")) +
  scale_colour_viridis_d() +
  guides(colour = "none") +
  facet_wrap(~ name, ncol = 2, scales = "free_y", 
             labeller = as_labeller(facet_df)) +
  theme_bw() +
  theme(legend.position = c(1,0),
        legend.justification = c(1,0),
        legend.background = element_blank(),
        legend.title.align = 0.5,
        axis.title = element_blank(),
        axis.text = element_text(size = 10),
        strip.text = element_text(size = 10))
dev.off()

#### __Beta ####
sc_bc <- calc_dm(sc$data_loaded)
set.seed(1150)
adonis2(sc_bc ~ sc$map_loaded$Treatment + sc$map_loaded$Depth) # Both sig
anova(betadisper(sc_bc, sc$map_loaded$Treatment)) # Dispersion not homogeneous
anova(betadisper(sc_bc, sc$map_loaded$Depth)) # Dispersion homogeneous
sc_pcoa <- cmdscale(sc_bc, k = nrow(sc$map_loaded) - 1, eig = T)
pcoaA1 <- round((eigenvals(sc_pcoa)/sum(eigenvals(sc_pcoa)))[1]*100, digits = 1)
pcoaA2 <- round((eigenvals(sc_pcoa)/sum(eigenvals(sc_pcoa)))[2]*100, digits = 1)
sc$map_loaded$Axis01 <- scores(sc_pcoa)[,1]
sc$map_loaded$Axis02 <- scores(sc_pcoa)[,2]
micro.hulls <- ddply(sc$map_loaded, c("TrtDepth"), find_hull)
g <- ggplot(sc$map_loaded, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = TrtDepth, fill = TrtDepth),
               alpha = 0.1, show.legend = F) +
  geom_point(size = 3, alpha = 1, aes(colour = TrtDepth, shape = TrtDepth),
             show.legend = F) +
  labs(x = paste("PC1: ", pcoaA1, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2, "%", sep = "")) +
  scale_colour_manual(values = c("#440154FF", "#440154FF", "#21908CFF", 
                                 "#21908CFF", "#FDE725FF", "#FDE725FF")) +
  scale_fill_manual(values = c("#440154FF", "#440154FF", "#21908CFF", 
                               "#21908CFF", "#FDE725FF", "#FDE725FF")) +
  scale_shape_manual(values = c(16,17,16,17,16,17)) +
  theme_bw() +  
  theme(legend.position = "none",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        plot.margin = margin(5, -5, 5, 5, "pt"))
leg <- get_legend(ggplot(sc$map_loaded, aes(Axis01, Axis02)) +
      geom_point(size = 3, alpha = 1, aes(colour = Treatment, shape = Depth)) +
      theme_bw() +
      scale_colour_viridis_d(labels = c("Control", "+Freshwater", "+Saltwater"),
                             guide = guide_legend(reverse = T,
                                                  override.aes = list(shape = 15))) +
                    labs(shape = "Depth (m)"))
pdf("InitialFigs/SC_PCoA.pdf", width = 6, height = 4)
plot_grid(g, leg, rel_widths = c(3.5, 1))
dev.off()

#### __Taxa ####
bar_text <- data.frame(group_by = c("Freshwater amended0.02", "Freshwater amended0.1"),
                       y = c(1.05, 1.05),
                       label = c("|-----2 cm-----|",
                                 "|-----10 cm-----|"))

sc_phyla <- summarize_taxonomy(sc, level = 2, report_higher_tax = F)
plot_ts_heatmap(sc_phyla, sc$map_loaded, 0.01, 'TrtDepth', rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
sc_barsP <- plot_taxa_bars(sc_phyla, sc$map_loaded, "TrtDepth", num_taxa = 12, data_only = T) %>%
  mutate(group_by = factor(group_by, levels = c("Control0.02","Freshwater amended0.02",
                                                "Saltwater amended0.02", "Control0.1",
                                                "Freshwater amended0.1", "Saltwater amended0.1"))) %>%
  mutate(taxon = fct_rev(taxon))
pdf("InitialFigs/SC_Phyla.pdf", width = 7, height = 5)
ggplot(sc_barsP, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  geom_text(data = bar_text,
            aes(x = group_by, y = y, label = label),
            inherit.aes = F) +
  labs(x = "Sample", y = "Relative abundance", fill = "Phylum") +
  scale_fill_manual(values = c("grey90", brewer.pal(12, "Paired")[12:1])) +
  scale_y_continuous(expand = c(0.01, 0.01)) + 
  scale_x_discrete(labels = c("Control", "+Fresh", "+Salt", "Control", "+Fresh", "+Salt")) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 10))
dev.off()
taxa_summary_by_sample_type(sc_phyla, sc$map_loaded, 'TrtDepth', 0.01, 'KW')

sc_class <- summarize_taxonomy(sc, level = 3, report_higher_tax = F)
plot_ts_heatmap(sc_class, sc$map_loaded, 0.01, 'TrtDepth', rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
sc_barsC <- plot_taxa_bars(sc_class, sc$map_loaded, "TrtDepth", num_taxa = 12, data_only = T) %>%
  mutate(group_by = factor(group_by, 
                           levels = c("Control0.02","Freshwater amended0.02",
                                      "Saltwater amended0.02", "Control0.1",
                                      "Freshwater amended0.1", "Saltwater amended0.1"))) %>%
  mutate(taxon = gsub("NA", "Unclassified", taxon)) %>%
  mutate(taxon = fct_relevel(taxon, "Other", after = Inf)) %>%
  mutate(taxon = fct_relevel(taxon, "Unclassified", after = Inf)) %>%
  mutate(taxon = fct_rev(taxon))
ggplot(sc_barsC, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  geom_text(data = bar_text,
            aes(x = group_by, y = y, label = label),
            inherit.aes = F) +
  labs(x = "Sample", y = "Relative abundance", fill = "Class") +
  scale_fill_manual(values = c("grey75", "grey90", brewer.pal(12, "Paired")[1:11])) +
  scale_y_continuous(expand = c(0.01, 0.01)) + 
  scale_x_discrete(labels = c("Control", "+Fresh", "+Salt", "Control", "+Fresh", "+Salt")) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 10))
taxa_summary_by_sample_type(sc_class, sc$map_loaded, 'TrtDepth', 0.01, 'KW')

sc_order <- summarize_taxonomy(sc, level = 4, report_higher_tax = F)
plot_ts_heatmap(sc_order, sc$map_loaded, 0.01, 'TrtDepth', rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
sc_barsO <- plot_taxa_bars(sc_order, sc$map_loaded, "TrtDepth", num_taxa = 12, data_only = T) %>%
  mutate(group_by = factor(group_by, levels = c("Control0.02","Freshwater amended0.02",
                                                "Saltwater amended0.02", "Control0.1",
                                                "Freshwater amended0.1", "Saltwater amended0.1"))) %>%
  mutate(taxon = gsub("NA", "Unclassified", taxon)) %>%
  mutate(taxon = fct_relevel(taxon, "Other", after = Inf)) %>%
  mutate(taxon = fct_relevel(taxon, "Unclassified", after = Inf)) %>%
  mutate(taxon = fct_rev(taxon))
ggplot(sc_barsO, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  geom_text(data = bar_text,
            aes(x = group_by, y = y, label = label),
            inherit.aes = F) +
  labs(x = "Sample", y = "Relative abundance", fill = "Order") +
  scale_fill_manual(values = c("grey75", "grey90", brewer.pal(12, "Paired")[1:11])) +
  scale_y_continuous(expand = c(0.01, 0.01)) + 
  scale_x_discrete(labels = c("Control", "+Fresh", "+Salt", "Control", "+Fresh", "+Salt")) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 10))
taxa_summary_by_sample_type(sc_order, sc$map_loaded, 'TrtDepth', 0.01, 'KW')

sc_family <- summarize_taxonomy(sc, level = 5, report_higher_tax = F)
plot_ts_heatmap(sc_family, sc$map_loaded, 0.01, 'TrtDepth', rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
sc_barsF <- plot_taxa_bars(sc_family, sc$map_loaded, "TrtDepth", num_taxa = 12, data_only = T) %>%
  mutate(group_by = factor(group_by, levels = c("Control0.02","Freshwater amended0.02",
                                                "Saltwater amended0.02", "Control0.1",
                                                "Freshwater amended0.1", "Saltwater amended0.1"))) %>%
  mutate(taxon = gsub("NA", "Unclassified", taxon)) %>%
  mutate(taxon = fct_relevel(taxon, "Other", after = Inf)) %>%
  mutate(taxon = fct_relevel(taxon, "Unclassified", after = Inf)) %>%
  mutate(taxon = fct_rev(taxon))
ggplot(sc_barsF, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  geom_text(data = bar_text,
            aes(x = group_by, y = y, label = label),
            inherit.aes = F) +
  labs(x = "Sample", y = "Relative abundance", fill = "Family") +
  scale_fill_manual(values = c("grey75", "grey90", brewer.pal(12, "Paired")[1:11])) +
  scale_y_continuous(expand = c(0.01, 0.01)) + 
  scale_x_discrete(labels = c("Control", "+Fresh", "+Salt", "Control", "+Fresh", "+Salt")) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 10))
taxa_summary_by_sample_type(sc_family, sc$map_loaded, 'TrtDepth', 0.01, 'KW')

sc_genus <- summarize_taxonomy(sc, level = 6, report_higher_tax = F)
plot_ts_heatmap(sc_genus, sc$map_loaded, 0.01, 'TrtDepth', rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
sc_barsG <- plot_taxa_bars(sc_genus, sc$map_loaded, "TrtDepth", num_taxa = 12, data_only = T) %>%
  mutate(group_by = factor(group_by, 
                           levels = c("Control0.02","Freshwater amended0.02",
                                      "Saltwater amended0.02", "Control0.1",
                                      "Freshwater amended0.1", "Saltwater amended0.1"))) %>%
  mutate(taxon = gsub("NA", "Unclassified", taxon)) %>%
  mutate(taxon = fct_relevel(taxon, "Other", after = Inf)) %>%
  mutate(taxon = fct_relevel(taxon, "Unclassified", after = Inf)) %>%
  mutate(taxon = fct_rev(taxon))
ggplot(sc_barsG, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  geom_text(data = bar_text,
            aes(x = group_by, y = y, label = label),
            inherit.aes = F) +
  labs(x = "Sample", y = "Relative abundance", fill = "Genus") +
  scale_fill_manual(values = c("grey75", "grey90", brewer.pal(12, "Paired")[1:11])) +
  scale_y_continuous(expand = c(0.01, 0.01)) + 
  scale_x_discrete(labels = c("Control", "+Fresh", "+Salt", "Control", "+Fresh", "+Salt")) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 10))
taxa_summary_by_sample_type(sc_genus, sc$map_loaded, 'TrtDepth', 0.01, 'KW')

sc_guilds <- summarize_taxonomy(sc, level = 9, report_higher_tax = F)
plot_ts_heatmap(sc_guilds, sc$map_loaded, 0.01, 'TrtDepth', rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
sc_barsGu <- plot_taxa_bars(sc_guilds,
                       sc$map_loaded,
                       "TrtDepth",
                       num_taxa = 20,
                       data_only = TRUE) %>%
  mutate(group_by = factor(group_by, 
                           levels = c("Control0.02","Freshwater amended0.02",
                                      "Saltwater amended0.02", "Control0.1",
                                      "Freshwater amended0.1", "Saltwater amended0.1"))) %>%
  filter(taxon != "NA") %>%
  droplevels() %>%
  mutate(taxon = factor(taxon,
                        levels = Guild_cols$Guild))
tallest_bar <- sc_barsGu %>%
  group_by(group_by) %>%
  summarise(sum = sum(mean_value))
bar_textGu <- data.frame(group_by = c("Freshwater amended0.02", "Freshwater amended0.1"),
                         y = c(max(tallest_bar$sum) + 0.01, max(tallest_bar$sum) + 0.01),
                         label = c("|-----2 cm-----|",
                                   "|-----10 cm-----|"))
pdf("InitialFigs/SC_Guilds.pdf", width = 7, height = 5)
ggplot(sc_barsGu, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  geom_text(data = bar_textGu,
            aes(x = group_by, y = y, label = label),
            inherit.aes = F) +
  labs(x = "Sample", y = "Relative abundance", fill = "Guild") +
  scale_fill_manual(values = Guild_cols$color) +
  scale_y_continuous(expand = c(max(tallest_bar$sum)/100, max(tallest_bar$sum)/100)) + 
  scale_x_discrete(labels = c("Control", "+Fresh", "+Salt", "Control", "+Fresh", "+Salt")) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 10))
dev.off()
taxa_summary_by_sample_type(sc_guilds, 
                            sc$map_loaded, 
                            type_header = 'TrtDepth', 
                            filter_level = 0.01, 
                            test_type = 'KW')

# Guilds by sample (to show variation)
sc$map_loaded$sampleID <- rownames(sc$map_loaded)
sc_barsGu <- plot_taxa_bars(sc_guilds,
                            sc$map_loaded,
                            "sampleID",
                            num_taxa = 20,
                            data_only = TRUE) %>%
  filter(taxon != "NA") %>%
  droplevels() %>%
  mutate(taxon = factor(taxon,
                        levels = Guild_cols$Guild)) %>%
  left_join(., sc$map_loaded, by = c("group_by" = "sampleID"))
facet_names <- c("0.02" = "Depth 2 cm",
                 "0.1" = "Depth 10 cm",
                 "Control" = "Control",
                 "Freshwater amended" = "+Freshwater",
                 "Saltwater amended" = "+Saltwater")
tallest_bar <- sc_barsGu %>%
  group_by(group_by) %>%
  summarise(sum = sum(mean_value))
pdf("InitialFigs/SC_Guilds_samples.pdf", width = 9, height = 5)
ggplot(sc_barsGu, aes(group_by, mean_value, fill = taxon)) +
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
        strip.background = element_rect(size = 0.2),
        axis.line.y = element_blank())
dev.off()

#### __Simper ####
sc_sim <- simper(t(sc$data_loaded), 
                 sc$map_loaded$TrtDepth)
sc_s <- summary(sc_sim)
head(sc_s$`Control0.02_Saltwater amended0.02`)
head(sc_s$`Control0.1_Saltwater amended0.1`)
sc_df1 <- head(sc_s$`Control0.02_Saltwater amended0.02`, n = 20) %>%
  mutate(Comparison = "ASW/Control 2 cm",
         ASV = rownames(.))
sc_df2 <- head(sc_s$`Control0.1_Saltwater amended0.1`, n = 20) %>%
  mutate(Comparison = "ASW/Control 10 cm",
         ASV = rownames(.))
sc_simper_results <- rbind(sc_df1, sc_df2) %>%
  left_join(., sc$taxonomy_loaded, by = c("ASV" = "taxonomy8")) %>%
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
write_xlsx(sc_simper_results, 
           "simper_results_sc.xlsx",
           format_headers = F)

#### __Multipatt ####
set.seed(1202)
sc_mp <- multipatt(t(sc$data_loaded), 
                   sc$map_loaded$TrtDepth, 
                   func = "r.g", 
                   control = how(nperm=999))
sc_mp_results <- sc_mp$sign %>%
  mutate(num_sites = rowSums(cbind(`s.Control0.02`, `s.Control0.1`, `s.Freshwater amended0.02`,
                                   `s.Freshwater amended0.1`, `s.Saltwater amended0.02`,
                                   `s.Saltwater amended0.1`)),
         q.value = qvalue(sc_mp$sign$p.value)$qvalues,
         Group = "NA",
         ASV = rownames(sc_mp$sign)) %>%
  filter(q.value < 0.05) %>%
  filter(stat >= 0.85) %>%
  filter(num_sites <= 2) %>%
  left_join(., sc$taxonomy_loaded, by = c("ASV" = "taxonomy8"))
for (i in 1:nrow(sc_mp_results)) {
  if (sc_mp_results$num_sites[i] == 1 & sc_mp_results$s.Control0.02[i] == 1) {
    sc_mp_results$Group[i] <- "Control_0.02"
  }
}
for (i in 1:nrow(sc_mp_results)) {
  if (sc_mp_results$num_sites[i] == 1 & sc_mp_results$s.Control0.1[i] == 1) {
    sc_mp_results$Group[i] <- "Control_0.1"
  }
}
for (i in 1:nrow(sc_mp_results)) {
  if (sc_mp_results$num_sites[i] == 1 & sc_mp_results$`s.Freshwater amended0.02`[i] == 1) {
    sc_mp_results$Group[i] <- "+Fresh_0.02"
  }
}
for (i in 1:nrow(sc_mp_results)) {
  if (sc_mp_results$num_sites[i] == 1 & sc_mp_results$`s.Freshwater amended0.1`[i] == 1) {
    sc_mp_results$Group[i] <- "+Fresh_0.1"
  }
}
for (i in 1:nrow(sc_mp_results)) {
  if (sc_mp_results$num_sites[i] == 1 & sc_mp_results$`s.Saltwater amended0.02`[i] == 1) {
    sc_mp_results$Group[i] <- "+Salt_0.02"
  }
}
for (i in 1:nrow(sc_mp_results)) {
  if (sc_mp_results$num_sites[i] == 1 & sc_mp_results$`s.Saltwater amended0.1`[i] == 1) {
    sc_mp_results$Group[i] <- "+Salt_0.1"
  }
}
for (i in 1:nrow(sc_mp_results)) {
  if (sc_mp_results$num_sites[i] == 2 & sc_mp_results$`s.Saltwater amended0.02`[i] == 1 &
      sc_mp_results$`s.Saltwater amended0.1`[i] == 1) {
    sc_mp_results$Group[i] <- "+Salt_both"
  }
}
table(sc_mp_results$Group)
sc_asv <- summarize_taxonomy(sc, level = 8, report_higher_tax = F)
sc_asv_all <- data.frame("RelAbundance" = round(rowMeans(sc_asv) * 100, digits = 4)) %>%
  mutate(ASV = rownames(.))
sc_mp_corrs <- as.data.frame(sc_mp$str) %>%
  dplyr::select(1:length(levels(sc$map_loaded$TrtDepth)), 
                `Saltwater amended0.02+Saltwater amended0.1`) %>%
  mutate(ASV = rownames(.)) %>%
  filter(ASV %in% sc_mp_results$ASV) %>%
  set_names(c("Control_0.02", "Control_0.1", "+Fresh_0.02", "+Fresh_0.1", "+Salt_0.02", 
              "+Salt_0.1", "+Salt_both", "ASV"))
# Add corrs and taxonomy
sc_mp_results <- sc_mp_results %>%
  filter(Group == "+Salt_0.02" | Group == "+Salt_0.1" | Group == "+Salt_both") %>%
  left_join(., sc_asv_all, by = "ASV") %>%
  left_join(., sc_mp_corrs, by = "ASV") %>%
  dplyr::select(taxonomy2, taxonomy6, ASV, Group, names(sc_mp_corrs)[1:7], "RelAbundance") %>%
  mutate(ASV = gsub("ASV", "OTU", ASV)) %>%
  arrange(Group, taxonomy2, taxonomy6) %>%
  mutate(ASV = factor(ASV, levels = ASV)) %>%
  set_names(c("Phylum", "Genus", "ASV", "Group", names(sc_mp_corrs)[1:7], "RelAbundance")) %>%
  mutate(taxon = paste(Phylum, Genus, ASV, sep = "; ")) %>%
  mutate(taxon = factor(taxon, levels = taxon))
# Plot heatmap with abundance barplot
sc.hm.melted <- sc_mp_results %>%
  dplyr::select(taxon, names(sc_mp_corrs)[1:7]) %>%
  melt(., id.vars = c("taxon"))
sc.hm <- ggplot(data = sc.hm.melted, 
             aes(x = factor(taxon), y = variable, fill = value)) + 
  geom_tile() + 
  scale_fill_distiller(name = "Indicator correlation index", palette = "RdBu", 
                       direction = -1, na.value = "transparent", type = "div", limits = c(-1, 1)) +
  scale_x_discrete(breaks = unique(sc.hm.melted$taxon), labels = unique(sc.hm.melted$taxon),
                   limits = rev(levels(sc.hm.melted$taxon))) +
  scale_y_discrete(position = "right") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        legend.title = element_text(size = 10), legend.key.height = unit(0.5,"cm"),
        legend.key.width = unit(1, "cm"), legend.text = element_text(size = 7)) +
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5)) +
  coord_flip()
sc.l <- get_legend(sc.hm)
sc.hm.clean <- sc.hm +
  theme(axis.title.y = element_blank(), axis.text.y = element_text(margin = margin(c(0,-5,0,0)),
                                                                   size = 6),
        axis.ticks.y = element_blank(), axis.ticks.x = element_blank(), 
        panel.grid.major = element_blank(), legend.position="none",
        axis.text.x.top = element_text(size = 6, margin = margin(c(0,0,-2,0))), 
        plot.margin = margin(c(0,-2,0,5))) +
  labs(y = "Group")
sc.bp.y <- ggplot(data = sc_mp_results, aes(x = taxon, y = RelAbundance)) + 
  geom_bar(stat = "identity", fill = "grey") + 
  scale_y_continuous(position = "right") +
  scale_x_discrete(limits = rev(levels(sc_mp_results$taxon))) +
  coord_flip() + 
  theme_minimal() +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position="none", plot.margin = margin(c(0,5,0,-5)),
        axis.text.x.top = element_text(size = 6, margin = margin(c(0,0,-2,0)))) +
  labs(y = "Rel. abund. (%)")
pdf("InitialFigs/SC_Multipatt.pdf", width = 8, height = 5)
plot_grid(sc.hm.clean, sc.bp.y, sc.l, NULL, nrow = 2, ncol = 2, 
          rel_widths = c(10,2), rel_heights = c(12, 2), align = "hv", axis = "b")
dev.off()

# Problem is that at OTU level, OTUs are too rare, not even 1%.
# So let's see what indicators there are at the genus level.
# First get aggregate taxonomy table
sc_tax <- sc$taxonomy_loaded %>%
  group_by(taxonomy6) %>%
  slice_head(n = 1) %>%
  dplyr::select(-taxonomy7, -taxonomy8)
sc_genus <- summarize_taxonomy(sc, level = 6, relative = F, report_higher_tax = F)
set.seed(1202)
sc_mp <- multipatt(t(sc_genus), 
                   sc$map_loaded$TrtDepth, 
                   func = "r.g", 
                   control = how(nperm=999))
sc_mp_results <- sc_mp$sign %>%
  mutate(num_sites = rowSums(cbind(`s.Control0.02`, `s.Control0.1`, `s.Freshwater amended0.02`,
                                   `s.Freshwater amended0.1`, `s.Saltwater amended0.02`,
                                   `s.Saltwater amended0.1`)),
         q.value = qvalue(sc_mp$sign$p.value)$qvalues,
         Group = "NA",
         Genus = rownames(sc_mp$sign)) %>%
  filter(q.value < 0.05) %>%
  filter(stat >= 0.7) %>%
  filter(num_sites <= 2) %>%
  left_join(., sc_tax, by = c("Genus" = "taxonomy6"))
for (i in 1:nrow(sc_mp_results)) {
  if (sc_mp_results$num_sites[i] == 1 & sc_mp_results$s.Control0.02[i] == 1) {
    sc_mp_results$Group[i] <- "Control_0.02"
  }
}
for (i in 1:nrow(sc_mp_results)) {
  if (sc_mp_results$num_sites[i] == 1 & sc_mp_results$s.Control0.1[i] == 1) {
    sc_mp_results$Group[i] <- "Control_0.1"
  }
}
for (i in 1:nrow(sc_mp_results)) {
  if (sc_mp_results$num_sites[i] == 1 & sc_mp_results$`s.Freshwater amended0.02`[i] == 1) {
    sc_mp_results$Group[i] <- "+Fresh_0.02"
  }
}
for (i in 1:nrow(sc_mp_results)) {
  if (sc_mp_results$num_sites[i] == 1 & sc_mp_results$`s.Freshwater amended0.1`[i] == 1) {
    sc_mp_results$Group[i] <- "+Fresh_0.1"
  }
}
for (i in 1:nrow(sc_mp_results)) {
  if (sc_mp_results$num_sites[i] == 1 & sc_mp_results$`s.Saltwater amended0.02`[i] == 1) {
    sc_mp_results$Group[i] <- "+Salt_0.02"
  }
}
for (i in 1:nrow(sc_mp_results)) {
  if (sc_mp_results$num_sites[i] == 1 & sc_mp_results$`s.Saltwater amended0.1`[i] == 1) {
    sc_mp_results$Group[i] <- "+Salt_0.1"
  }
}
for (i in 1:nrow(sc_mp_results)) {
  if (sc_mp_results$num_sites[i] == 2 & sc_mp_results$`s.Saltwater amended0.02`[i] == 1 &
      sc_mp_results$`s.Saltwater amended0.1`[i] == 1) {
    sc_mp_results$Group[i] <- "+Salt_both"
  }
}
table(sc_mp_results$Group)
sc_genus_all <- data.frame("RelAbundance" = round(rowMeans(sc_genus)/min(colSums(sc$data_loaded)) * 100,
                                                  digits = 4)) %>%
  mutate(Genus = rownames(.))
sc_mp_corrs <- as.data.frame(sc_mp$str) %>%
  dplyr::select(1:length(levels(sc$map_loaded$TrtDepth)), 
                `Saltwater amended0.02+Saltwater amended0.1`) %>%
  mutate(Genus = rownames(.)) %>%
  filter(Genus %in% sc_mp_results$Genus) %>%
  set_names(c("Control_0.02", "Control_0.1", "+Fresh_0.02", "+Fresh_0.1", "+Salt_0.02", 
              "+Salt_0.1", "+Salt_both", "Genus"))
# Add corrs and taxonomy
sc_mp_results <- sc_mp_results %>%
  filter(Group == "+Salt_0.02" | Group == "+Salt_0.1" | Group == "+Salt_both") %>%
  left_join(., sc_genus_all, by = "Genus") %>%
  left_join(., sc_mp_corrs, by = "Genus") %>%
  dplyr::select(taxonomy2, Genus, Group, names(sc_mp_corrs)[1:7], "RelAbundance") %>%
  arrange(Group, taxonomy2) %>%
  mutate(Genus = factor(Genus, levels = Genus)) %>%
  set_names(c("Phylum", "Genus", "Group", names(sc_mp_corrs)[1:7], "RelAbundance")) %>%
  mutate(taxon = paste(Phylum, Genus, sep = "; ")) %>%
  mutate(taxon = factor(taxon, levels = taxon))
# Plot heatmap with abundance barplot
sc.hm.melted <- sc_mp_results %>%
  dplyr::select(taxon, names(sc_mp_corrs)[1:7]) %>%
  melt(., id.vars = c("taxon"))
sc.hm <- ggplot(data = sc.hm.melted, 
                aes(x = factor(taxon), y = variable, fill = value)) + 
  geom_tile() + 
  scale_fill_distiller(name = "Indicator correlation index", palette = "RdBu", 
                       direction = -1, na.value = "transparent", type = "div", limits = c(-1, 1)) +
  scale_x_discrete(breaks = unique(sc.hm.melted$taxon), labels = unique(sc.hm.melted$taxon),
                   limits = rev(levels(sc.hm.melted$taxon))) +
  scale_y_discrete(position = "right") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        legend.title = element_text(size = 10), legend.key.height = unit(0.5,"cm"),
        legend.key.width = unit(1, "cm"), legend.text = element_text(size = 7)) +
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5)) +
  coord_flip()
sc.l <- get_legend(sc.hm)
sc.hm.clean <- sc.hm +
  theme(axis.title.y = element_blank(), 
        axis.text.y = element_text(margin = margin(c(0,-5,0,0)),
                                   size = 6),
        axis.ticks.y = element_blank(), axis.ticks.x = element_blank(), 
        panel.grid.major = element_blank(), legend.position="none",
        axis.text.x.top = element_text(size = 6, margin = margin(c(0,0,-2,0))), 
        plot.margin = margin(c(0,-2,0,5))) +
  labs(y = "Group")
sc.bp.y <- ggplot(data = sc_mp_results, aes(x = taxon, y = RelAbundance)) + 
  geom_bar(stat = "identity", fill = "grey") + 
  scale_y_continuous(position = "right") +
  scale_x_discrete(limits = rev(levels(sc_mp_results$taxon))) +
  coord_flip() + 
  theme_minimal() +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position="none", plot.margin = margin(c(0,5,0,-5)),
        axis.text.x.top = element_text(size = 6, margin = margin(c(0,0,-2,0)))) +
  labs(y = "Rel. abund. (%)")
pdf("InitialFigs/SC_Multipatt_genus.pdf", width = 8, height = 5)
plot_grid(sc.hm.clean, sc.bp.y, sc.l, NULL, nrow = 2, ncol = 2, 
          rel_widths = c(10,2), rel_heights = c(12, 2), align = "hv", axis = "b")
dev.off()

#### __BGC ####
# Add biogeochem. analysis
# Get variables, make corrplot, envfit, ordination with vectors, ordistep, taxa correlations
# Note only BGC for D2 samples in SC

# Get variables
env <- sc$map_loaded %>%
  dplyr::select(CH4_ug_m2_h, CO2_ug_m2_h, CH4_pw_air_ppmv,
                CH4_pot_umol_gdw_h, CO2_pot_umol_gdw_h,
                N2_umol_m2_h,	SOD_umol_m2_h, NO3_umol_m2_h,	NH4_umol_m2_h, SRP_umol_m2_h, DON_umol_m2_h,
                Conductivity_uS_cm, Salinity_ppt_all, DIC_mgL,
                sed_per_C, sed_per_N, sed_CN, sed_per_org, sed_per_inorg, pH)
env_nona <- na.omit(env)

# Corrplot
C <- cor(env_nona)
corrplot(C, method = "number", type = "lower", 
         tl.cex = 0.5, tl.col = "black", number.cex = 0.5)
pdf("InitialFigs/SC_BGC_CH4.pdf", width = 7, height = 5)
meth_corr_by_bgc(env_nona = env_nona)
dev.off()

# Envfit
pcoa <- cmdscale(sc_bc, k = nrow(sc$map_loaded) - 1, eig = T)
set.seed(105)
ef <- envfit(pcoa, env, permutations = 999, na.rm = TRUE)
arrow_factor <- ordiArrowMul(ef)
manual_factor <- 0.2
vec.df <- as.data.frame(ef$vectors$arrows*sqrt(ef$vectors$r)) %>%
mutate(Dim1 = Dim1 * manual_factor,
       Dim2 = Dim2 * manual_factor) %>%
  mutate(variables = rownames(.)) %>%
  filter(ef$vectors$pvals < 0.05) %>%
  filter(variables != "Conductivity_uS_cm") %>%
  mutate(shortnames = c("CO2", "N2", "SOD", "Salinity", "N", "C:N", "pH"))

# Plot with significant vectors
ordiplot(pcoa)
plot(ef, p.max = 0.05, cex = 0.5)
plot(pcoa)
plot(ef, p.max = 0.05, cex = 0.5)
pcoaA1 <- round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[1]*100, digits = 1)
pcoaA2 <- round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[2]*100, digits = 1)
sc$map_loaded$Axis01 <- scores(pcoa)[,1]
sc$map_loaded$Axis02 <- scores(pcoa)[,2]
micro.hulls <- ddply(sc$map_loaded, c("TrtDepth"), find_hull)
g <- ggplot(sc$map_loaded, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = TrtDepth, fill = TrtDepth),
               alpha = 0.1, show.legend = F) +
  geom_point(size = 3, alpha = 1, aes(colour = TrtDepth, shape = TrtDepth),
             show.legend = F) +
  geom_segment(data = vec.df,
               aes(x = 0, xend = Dim1, y = 0, yend = Dim2),
               arrow = arrow(length = unit(0.5, "cm")),
               colour = "gray", alpha = 0.5,
               inherit.aes = FALSE) + 
  geom_text(data = vec.df,
            aes(x = Dim1, y = Dim2, label = shortnames),
            size = 3.5, color = "black") +
  labs(x = paste("PC1: ", pcoaA1, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2, "%", sep = "")) +
  scale_colour_manual(values = c("#440154FF", "#440154FF", "#21908CFF", 
                                 "#21908CFF", "#FDE725FF", "#FDE725FF")) +
  scale_fill_manual(values = c("#440154FF", "#440154FF", "#21908CFF", 
                               "#21908CFF", "#FDE725FF", "#FDE725FF")) +
  scale_shape_manual(values = c(16,17,16,17,16,17)) +
  coord_fixed() +
  theme_bw() +  
  theme(legend.position = "none",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        plot.margin = margin(5, -5, 5, 5, "pt"))
g
leg <- get_legend(ggplot(sc$map_loaded, aes(Axis01, Axis02)) +
                    geom_point(size = 3, alpha = 1, aes(colour = Treatment, shape = Depth)) +
                    theme_bw() +
                    scale_colour_viridis_d(labels = c("Control", "+Freshwater", "+Saltwater"),
                                           guide = guide_legend(reverse = T,
                                                                override.aes = list(shape = 15))) +
                    labs(shape = "Depth (m)"))
pdf("InitialFigs/SC_PCoA_wBGC.pdf", width = 6, height = 4)
plot_grid(g, leg, rel_widths = c(3.5, 1))
dev.off()

# Ordistep
comm_nona <- as.data.frame(t(sc$data_loaded)) %>%
  filter(rownames(.) %in% rownames(env_nona))
mod0 <- rda(comm_nona ~ 1, env_nona)  # Model with intercept only
mod1 <- rda(comm_nona ~ ., env_nona)  # Model with all explanatory variables
mod <- ordistep(mod0, scope = formula(mod1))
mod
mod$anova # SOD and CO2

# Plot correlations for different taxonomic levels at a given % relative abundance threshold
meth_corr_by_taxonomy(input = sc, level = 2, threshold = 0.5)
meth_corr_by_taxonomy(input = sc, level = 3, threshold = 0.5)
meth_corr_by_taxonomy(input = sc, level = 4, threshold = 0.5)
meth_corr_by_taxonomy(input = sc, level = 5, threshold = 0.5)
meth_corr_by_taxonomy(input = sc, level = 6, threshold = 0.5)
meth_corr_by_taxonomy(input = sc, level = 8, threshold = 0.5)
pdf("InitialFigs/SC_Guilds_CH4.pdf", width = 7, height = 5)
meth_corr_by_taxonomy(input = sc, level = 9, threshold = 0)
dev.off()



#### _North Carolina ####
# input_filt <- readRDS("input_filt.rds")
input_filt <- readRDS("input_filt_wBGC.rds")
nc <- filter_data(input_filt,
                  filter_cat = "Estuary",
                  keep_vals = "Alligator")
nc$map_loaded$sampleID <- rownames(nc$map_loaded)
nc <- filter_data(nc,
                  filter_cat = "sampleID",
                  filter_vals = c("TL_nw_d1_DI_ctrl_AF1", "TL_nw_d1_DI_ctrl_AF3", 
                                  "TL_nw_d1_DI_ctrl_AF4", "TL_nw_d1_ASW_noS_BF3",
                                  "TL_nw_d1_ASW_noS_BF4", "TL_nw_d1_ASW_noS_BF5"))
set.seed(530)
nc <- single_rarefy(nc, min(colSums(nc$data_loaded))) # 82350
nc$map_loaded <- nc$map_loaded %>%
  mutate(rich = specnumber(nc$data_loaded, MARGIN = 2),
         shannon = diversity(nc$data_loaded, index = "shannon", MARGIN = 2),
         TrtDepth = paste(nc$map_loaded$Treatment, nc$map_loaded$Depth, sep = ""),
         Depth = as.factor(Depth)) %>%
  mutate_if(is.character, as.factor)

#### __Alpha ####
leveneTest(nc$map_loaded$rich ~ nc$map_loaded$Treatment) # Homogeneous
m <- aov(rich ~ Treatment + Depth, data = nc$map_loaded)
Anova(m, type = "II") # Treatment and Depth
m <- aov(rich ~ Treatment, data = nc$map_loaded)
shapiro.test(m$residuals) # Not normal, but not too bad
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
pdf("InitialFigs/NC_Alpha.pdf", width = 6, height = 3)
ggplot(alpha_long, aes(reorder(Treatment, value, mean), value, 
                       colour = Treatment)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 2, alpha = 0.75, width = 0.2, aes(shape = Depth)) +
  geom_text(data = label_df, aes(Treatment, y, label = str_trim(.group)), 
            size = 4, color = "black") +
  labs(x = "Site", y = "Number of OTUs", shape = "Depth (m)") +
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

#### __Beta ####
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

#### __Taxa ####
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

#### __Simper ####
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


#### __Multipatt ####
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

#### __BGC ####
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



#### _Delaware field ####
# Delaware field (3 sites, Weston et al. 2014)
# Tidal freshwater, oligohaline, mesohaline
# input_filt <- readRDS("input_filt.rds")
input_filt <- readRDS("input_filt_wBGC.rds")
defie <- filter_data(input_filt,
                     filter_cat = "Estuary",
                     keep_vals = "Delaware")
defie <- filter_data(defie,
                     filter_cat = "Experiment",
                     keep_vals = "Soil")
set.seed(530)
defie <- single_rarefy(defie, min(colSums(defie$data_loaded))) # 52468
defie$map_loaded <- defie$map_loaded %>%
  mutate(rich = specnumber(defie$data_loaded, MARGIN = 2),
         shannon = diversity(defie$data_loaded, index = "shannon", MARGIN = 2),
         Treatment = as.factor(Treatment),
         Salt = recode_factor(Treatment, "TFM1_source" = "Freshwater",
                              "TFM2_source" = "Freshwater",
                              "MesoHal_source" = "Mesohaline",
                              "OligoHal_source" = "Oligohaline")) %>%
  mutate(Salt = factor(Salt, levels = c("Freshwater", "Oligohaline", "Mesohaline"))) %>%
  unite("TrtDepth", c("Salt", "Depth"), sep = "", remove = F) %>%
  mutate(Depth = as.factor(Depth)) %>%
  mutate_if(is.character, as.factor)

#### __Alpha ####
leveneTest(defie$map_loaded$rich ~ defie$map_loaded$Salt) # Homogeneous
m <- aov(rich ~ Salt + Depth, data = defie$map_loaded)
Anova(m, type = "II") # Salt, not Depth
m <- aov(rich ~ Salt, data = defie$map_loaded)
shapiro.test(m$residuals) # Normal
summary(m)
t <- emmeans(object = m, specs = "Salt") %>%
  cld(object = ., adjust = "Tukey", Letters = letters, alpha = 0.05) %>%
  mutate(name = "rich",
         y = max(defie$map_loaded$rich)+(max(defie$map_loaded$rich)-min(defie$map_loaded$rich))/20)
leveneTest(defie$map_loaded$shannon ~ defie$map_loaded$Salt) # Homogeneous
m1 <- aov(shannon ~ Salt + Depth, data = defie$map_loaded)
Anova(m1, type = "II") # Salt, not Depth
m1 <- aov(shannon ~ Salt, data = defie$map_loaded)
shapiro.test(m1$residuals) # Normal
summary(m1)
t1 <- emmeans(object = m1, specs = "Salt") %>%
  cld(object = ., adjust = "Tukey", Letters = letters, alpha = 0.05) %>%
  mutate(name = "shannon",
         y = max(defie$map_loaded$shannon)+(max(defie$map_loaded$shannon)-min(defie$map_loaded$shannon))/20)
label_df <- rbind(t, t1)
facet_df <- c("rich" = "(a) Richness",
              "shannon" = "(b) Shannon")
alpha_long <- defie$map_loaded %>%
  pivot_longer(cols = c("rich", "shannon"))
pdf("InitialFigs/DEfie_Alpha.pdf", width = 6, height = 3)
ggplot(alpha_long, aes(reorder(Salt, value, mean), value, 
                       colour = Salt)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 2, alpha = 0.75, width = 0.2, aes(shape = Depth)) +
  geom_text(data = label_df, aes(Salt, y, label = str_trim(.group)), 
            size = 4, color = "black") +
  labs(x = "Site", y = "Number of OTUs", shape = "Depth (m)") +
  scale_colour_viridis_d() +
  guides(colour = "none") +
  facet_wrap(~ name, ncol = 2, scales = "free_y", labeller = as_labeller(facet_df)) +
  theme_bw() +
  theme(legend.position = c(1,0),
        legend.justification = c(1,0),
        legend.background = element_blank(),
        legend.title.align = 0.5,
        axis.title = element_blank(),
        axis.text = element_text(size = 10),
        strip.text = element_text(size = 10))
dev.off()

#### __Beta ####
defie_bc <- calc_dm(defie$data_loaded)
set.seed(1150)
adonis2(defie_bc ~ defie$map_loaded$Salt + defie$map_loaded$Depth) # Salt sig
anova(betadisper(defie_bc, defie$map_loaded$Salt)) # Dispersion homogeneous
anova(betadisper(defie_bc, defie$map_loaded$Depth)) # Dispersion homogeneous
defie_pcoa <- cmdscale(defie_bc, k = nrow(defie$map_loaded) - 1, eig = T)
pcoaA1 <- round((eigenvals(defie_pcoa)/sum(eigenvals(defie_pcoa)))[1]*100, digits = 1)
pcoaA2 <- round((eigenvals(defie_pcoa)/sum(eigenvals(defie_pcoa)))[2]*100, digits = 1)
defie$map_loaded$Axis01 <- scores(defie_pcoa)[,1]
defie$map_loaded$Axis02 <- scores(defie_pcoa)[,2]
micro.hulls <- ddply(defie$map_loaded, c("Salt"), find_hull)
pdf("InitialFigs/DEfie_PCoA.pdf", width = 7, height = 5)
ggplot(defie$map_loaded, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = Salt, fill = Salt),
               alpha = 0.1, show.legend = F) +
  geom_point(size = 3, alpha = 1, aes(colour = Salt, shape = Depth),
             show.legend = T) +
  labs(x = paste("PC1: ", pcoaA1, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2, "%", sep = ""),
       shape = "Depth (m)") +
  scale_fill_viridis_d() +
  scale_colour_viridis_d() +
  theme_bw() +  
  theme(legend.position = "right",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10))
dev.off()

#### __Taxa ####
defie_phyla <- summarize_taxonomy(defie, level = 2, report_higher_tax = F)
plot_ts_heatmap(defie_phyla, defie$map_loaded, 0.01, 'Salt', rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
defie_barsP <- plot_taxa_bars(defie_phyla, defie$map_loaded, "Salt", num_taxa = 12, data_only = T) %>%
  mutate(taxon = fct_rev(taxon))
pdf("InitialFigs/DEfie_Phyla.pdf", width = 7, height = 5)
ggplot(defie_barsP, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Sample", y = "Relative abundance", fill = "Phylum") +
  scale_fill_manual(values = c("grey90", brewer.pal(12, "Paired")[12:1])) +
  scale_y_continuous(expand = c(0.01, 0.01)) + 
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 10))
dev.off()
taxa_summary_by_sample_type(defie_phyla, defie$map_loaded, 'Salt', 0.01, 'KW')

defie_class <- summarize_taxonomy(defie, level = 3, report_higher_tax = F)
plot_ts_heatmap(defie_class, defie$map_loaded, 0.01, 'Salt', rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
defie_barsC <- plot_taxa_bars(defie_class, defie$map_loaded, "Salt", 
                                 num_taxa = 12, data_only = T) %>%
  mutate(taxon = gsub("NA", "Unclassified", taxon)) %>%
  mutate(taxon = fct_relevel(taxon, "Other", after = Inf)) %>%
  mutate(taxon = fct_relevel(taxon, "Unclassified", after = Inf)) %>%
  mutate(taxon = fct_rev(taxon))
ggplot(defie_barsC, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = NULL, y = "Relative abundance", fill = "Class") +
  scale_fill_manual(values = c("grey75", "grey90", brewer.pal(12, "Paired")[11:1])) +
  scale_y_continuous(expand = c(0.01, 0.01)) + 
  theme_classic() +
  theme(axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10))
taxa_summary_by_sample_type(defie_class, defie$map_loaded, 'Salt', 0.01, 'KW')

defie_order <- summarize_taxonomy(defie, level = 4, report_higher_tax = F)
plot_ts_heatmap(defie_order, defie$map_loaded, 0.01, 'Salt', rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
defie_barsO <- plot_taxa_bars(defie_order, defie$map_loaded, "Salt", 
                              num_taxa = 12, data_only = T) %>%
  mutate(taxon = gsub("NA", "Unclassified", taxon)) %>%
  mutate(taxon = fct_relevel(taxon, "Other", after = Inf)) %>%
  mutate(taxon = fct_relevel(taxon, "Unclassified", after = Inf)) %>%
  mutate(taxon = fct_rev(taxon))
ggplot(defie_barsO, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = NULL, y = "Relative abundance", fill = "Order") +
  scale_fill_manual(values = c("grey75", "grey90", brewer.pal(12, "Paired")[11:1])) +
  scale_y_continuous(expand = c(0.01, 0.01)) + 
  theme_classic() +
  theme(axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10))
taxa_summary_by_sample_type(defie_order, defie$map_loaded, 'Salt', 0.01, 'KW')

defie_family <- summarize_taxonomy(defie, level = 5, report_higher_tax = F)
plot_ts_heatmap(defie_family, defie$map_loaded, 0.01, 'Salt', rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
defie_barsF <- plot_taxa_bars(defie_family, defie$map_loaded, "Salt", 
                              num_taxa = 12, data_only = T) %>%
  mutate(taxon = gsub("NA", "Unclassified", taxon)) %>%
  mutate(taxon = fct_relevel(taxon, "Other", after = Inf)) %>%
  mutate(taxon = fct_relevel(taxon, "Unclassified", after = Inf)) %>%
  mutate(taxon = fct_rev(taxon))
ggplot(defie_barsF, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = NULL, y = "Relative abundance", fill = "Family") +
  scale_fill_manual(values = c("grey75", "grey90", brewer.pal(12, "Paired")[11:1])) +
  scale_y_continuous(expand = c(0.01, 0.01)) + 
  theme_classic() +
  theme(axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10))
taxa_summary_by_sample_type(defie_family, defie$map_loaded, 'Salt', 0.01, 'KW')

defie_genus <- summarize_taxonomy(defie, level = 6, report_higher_tax = F)
plot_ts_heatmap(defie_genus, defie$map_loaded, 0.01, 'Salt', rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
defie_barsG <- plot_taxa_bars(defie_genus, defie$map_loaded, "Salt", 
                              num_taxa = 12, data_only = T) %>%
  mutate(taxon = gsub("NA", "Unclassified", taxon)) %>%
  mutate(taxon = fct_relevel(taxon, "Other", after = Inf)) %>%
  mutate(taxon = fct_relevel(taxon, "Unclassified", after = Inf)) %>%
  mutate(taxon = fct_rev(taxon))
ggplot(defie_barsG, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = NULL, y = "Relative abundance", fill = "Genus") +
  scale_fill_manual(values = c("grey75", "grey90", brewer.pal(12, "Paired")[11:1])) +
  scale_y_continuous(expand = c(0.01, 0.01)) + 
  theme_classic() +
  theme(axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10))
taxa_summary_by_sample_type(defie_genus, defie$map_loaded, 'Salt', 0.01, 'KW')

defie_guilds <- summarize_taxonomy(defie, level = 9, report_higher_tax = F)
plot_ts_heatmap(defie_guilds, defie$map_loaded, 0.01, 'Salt', rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
defie_barsGu <- plot_taxa_bars(defie_guilds, defie$map_loaded, "Salt", 
                               num_taxa = 20, data_only = T) %>%
  filter(taxon != "NA") %>%
  droplevels() %>%
  mutate(taxon = factor(taxon,
                        levels = Guild_cols$Guild))
tallest_bar <- defie_barsGu %>%
  group_by(group_by) %>%
  summarise(sum = sum(mean_value))
pdf("InitialFigs/DEfie_Guilds.pdf", width = 7, height = 5)
ggplot(defie_barsGu, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Salt", y = "Relative abundance", fill = "Guild") +
  scale_fill_manual(values = Guild_cols$color) +
  scale_y_continuous(expand = c(max(tallest_bar$sum)/100, max(tallest_bar$sum)/100)) + 
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        axis.title.x = element_blank())
dev.off()
taxa_summary_by_sample_type(defie_guilds, 
                            defie$map_loaded, 
                            type_header = 'Salt', 
                            filter_level = 0.01, 
                            test_type = 'KW')

# Guilds by sample (to show variation)
defie$map_loaded$sampleID <- rownames(defie$map_loaded)
defie_barsGu <- plot_taxa_bars(defie_guilds,
                            defie$map_loaded,
                            "sampleID",
                            num_taxa = 20,
                            data_only = TRUE) %>%
  filter(taxon != "NA") %>%
  droplevels() %>%
  mutate(taxon = factor(taxon,
                        levels = Guild_cols$Guild)) %>%
  left_join(., defie$map_loaded, by = c("group_by" = "sampleID"))
facet_names <- c("0.02" = "Depth 2 cm",
                 "0.12" = "Depth 12 cm",
                 "Freshwater" = "Freshwater",
                 "Oligohaline" = "Oligohaline",
                 "Mesohaline" = "Mesohaline")
tallest_bar <- defie_barsGu %>%
  group_by(group_by) %>%
  summarise(sum = sum(mean_value))
pdf("InitialFigs/DEfie_Guilds_samples.pdf", width = 9, height = 5)
ggplot(defie_barsGu, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Sample", y = "Relative abundance", fill = "Guild") +
  scale_fill_manual(values = Guild_cols$color) +
  scale_y_continuous(expand = c(max(tallest_bar$sum)/100, max(tallest_bar$sum)/100)) + 
  facet_nested(~ Depth + Salt, space = "free", scales = "free_x",
               labeller = as_labeller(facet_names)) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 90, hjust = 1, vjust = 0.5),
        strip.background = element_rect(size = 0.2),
        axis.line.y = element_blank())
dev.off()

#### __Simper ####
defie_sim <- simper(t(defie$data_loaded), 
                 defie$map_loaded$Salt)
defie_s <- summary(defie_sim)
head(defie_s$Freshwater_Mesohaline)
defie_simper_results <- head(defie_s$Freshwater_Mesohaline, n = 20) %>%
  mutate(ASV = rownames(.),
         Comparison = "Freshwater/Mesohaline") %>%
  left_join(., defie$taxonomy_loaded, by = c("ASV" = "taxonomy8")) %>%
  mutate(SaltResponse = ifelse(avb > ava, "Postive", "Negative")) %>%
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
write_xlsx(defie_simper_results, 
           "simper_results_DEfie.xlsx",
           format_headers = F)

#### __Multipatt ####
set.seed(1202)
defie_mp <- multipatt(t(defie$data_loaded), 
                   defie$map_loaded$Salt, 
                   func = "r.g", 
                   control = how(nperm=999))
# None with Q, use P
defie_mp_results <- defie_mp$sign %>%
  mutate(num_sites = rowSums(cbind(s.Freshwater, s.Oligohaline, s.Mesohaline)),
         q.value = qvalue(defie_mp$sign$p.value)$qvalues,
         Group = "NA",
         ASV = rownames(defie_mp$sign)) %>%
  filter(p.value == 0.001) %>%
  filter(stat >= 0.9) %>%
  filter(num_sites <= 2) %>%
  left_join(., defie$taxonomy_loaded, by = c("ASV" = "taxonomy8"))
for (i in 1:nrow(defie_mp_results)) {
  if (defie_mp_results$num_sites[i] == 1 & defie_mp_results$s.Freshwater[i] == 1) {
    defie_mp_results$Group[i] <- "Freshwater"
  }
}
for (i in 1:nrow(defie_mp_results)) {
  if (defie_mp_results$num_sites[i] == 1 & defie_mp_results$s.Oligohaline[i] == 1) {
    defie_mp_results$Group[i] <- "Oligohaline"
  }
}
for (i in 1:nrow(defie_mp_results)) {
  if (defie_mp_results$num_sites[i] == 1 & defie_mp_results$s.Mesohaline[i] == 1) {
    defie_mp_results$Group[i] <- "Mesohaline"
  }
}
table(defie_mp_results$Group)
defie_asv <- summarize_taxonomy(defie, level = 8, report_higher_tax = F)
defie_asv_all <- data.frame("RelAbundance" = round(rowMeans(defie_asv) * 100, digits = 5)) %>%
  mutate(ASV = rownames(.))
defie_mp_corrs <- as.data.frame(defie_mp$str) %>%
  dplyr::select(1:length(levels(defie$map_loaded$Salt))) %>%
  mutate(ASV = rownames(.)) %>%
  filter(ASV %in% defie_mp_results$ASV) %>%
  set_names(c("Freshwater", "Oligohaline", "Mesolhaline", "ASV"))
# Add corrs and taxonomy
defie_mp_results <- defie_mp_results %>%
  filter(Group == "Mesohaline") %>%
  left_join(., defie_asv_all, by = "ASV") %>%
  left_join(., defie_mp_corrs, by = "ASV") %>%
  dplyr::select(taxonomy2, taxonomy6, ASV, Group, names(defie_mp_corrs)[1:3], "RelAbundance") %>%
  mutate(ASV = gsub("ASV", "OTU", ASV)) %>%
  arrange(Group, taxonomy2, taxonomy6) %>%
  mutate(ASV = factor(ASV, levels = ASV)) %>%
  set_names(c("Phylum", "Genus", "ASV", "Group", names(defie_mp_corrs)[1:3], "RelAbundance")) %>%
  mutate(taxon = paste(Phylum, Genus, ASV, sep = "; ")) %>%
  mutate(taxon = factor(taxon, levels = taxon))
# Plot heatmap with abundance barplot
defie.hm.melted <- defie_mp_results %>%
  dplyr::select(taxon, names(defie_mp_corrs)[1:3]) %>%
  melt(., id.vars = c("taxon"))
defie.hm <- ggplot(data = defie.hm.melted, 
                aes(x = factor(taxon), y = variable, fill = value)) + 
  geom_tile() + 
  scale_fill_distiller(name = "Indicator correlation index", palette = "RdBu", 
                       direction = -1, na.value = "transparent", type = "div", limits = c(-1, 1)) +
  scale_x_discrete(breaks = unique(defie.hm.melted$taxon), labels = unique(defie.hm.melted$taxon),
                   limits = rev(levels(defie.hm.melted$taxon))) +
  scale_y_discrete(position = "right") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        legend.title = element_text(size = 10), legend.key.height = unit(0.5,"cm"),
        legend.key.width = unit(1, "cm"), legend.text = element_text(size = 7)) +
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5)) +
  coord_flip()
defie.l <- get_legend(defie.hm)
defie.hm.clean <- defie.hm +
  theme(axis.title.y = element_blank(), axis.text.y = element_text(margin = margin(c(0,-10,0,0)),
                                                                   size = 5),
        axis.ticks.y = element_blank(), axis.ticks.x = element_blank(), 
        panel.grid.major = element_blank(), legend.position="none",
        axis.text.x.top = element_text(size = 6, margin = margin(c(0,0,-2,0))), 
        plot.margin = margin(c(0,-2,0,5))) +
  labs(y = "Group")
defie.bp.y <- ggplot(data = defie_mp_results, aes(x = taxon, y = RelAbundance)) + 
  geom_bar(stat = "identity", fill = "grey") + 
  scale_y_continuous(position = "right") +
  scale_x_discrete(limits = rev(levels(defie_mp_results$taxon))) +
  coord_flip() + 
  theme_minimal() +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position="none", plot.margin = margin(c(0,5,0,-5)),
        axis.text.x.top = element_text(size = 6, margin = margin(c(0,0,-2,0)))) +
  labs(y = "Rel. abund. (%)")
pdf("InitialFigs/DEfie_Multipatt.pdf", width = 8, height = 5)
plot_grid(defie.hm.clean, defie.bp.y, defie.l, NULL, nrow = 2, ncol = 2, 
          rel_widths = c(10,2), rel_heights = c(12, 2), align = "hv", axis = "b")
dev.off()

# Problem is that at OTU level, many OTUs are too rare, not even 1%. 
# So let's see what indicators there are at the genus level.
# First get aggregate taxonomy table
defie_tax <- defie$taxonomy_loaded %>%
  group_by(taxonomy6) %>%
  slice_head(n = 1) %>%
  dplyr::select(-taxonomy7, -taxonomy8)
defie_genus <- summarize_taxonomy(defie, level = 6, relative = F, report_higher_tax = F)
set.seed(1202)
defie_mp <- multipatt(t(defie_genus), 
                   defie$map_loaded$Salt, 
                   func = "r.g", 
                   control = how(nperm=999))
defie_mp_results <- defie_mp$sign %>%
  mutate(num_sites = rowSums(cbind(s.Freshwater, s.Oligohaline, s.Mesohaline)),
         q.value = qvalue(defie_mp$sign$p.value)$qvalues,
         Group = "NA",
         Genus = rownames(defie_mp$sign)) %>%
  filter(q.value < 0.05) %>%
  filter(stat >= 0.7) %>%
  filter(num_sites <= 2) %>%
  left_join(., defie_tax, by = c("Genus" = "taxonomy6"))
for (i in 1:nrow(defie_mp_results)) {
  if (defie_mp_results$num_sites[i] == 1 & defie_mp_results$s.Freshwater[i] == 1) {
    defie_mp_results$Group[i] <- "Freshwater"
  }
}
for (i in 1:nrow(defie_mp_results)) {
  if (defie_mp_results$num_sites[i] == 1 & defie_mp_results$s.Oligohaline[i] == 1) {
    defie_mp_results$Group[i] <- "Oligohaline"
  }
}
for (i in 1:nrow(defie_mp_results)) {
  if (defie_mp_results$num_sites[i] == 1 & defie_mp_results$s.Mesohaline[i] == 1) {
    defie_mp_results$Group[i] <- "Mesohaline"
  }
}
table(defie_mp_results$Group)
defie_genus_all <- data.frame("RelAbundance" = round(rowMeans(defie_genus)/min(colSums(defie$data_loaded)) * 100, 
                                                     digits = 4)) %>%
  mutate(Genus = rownames(.))
defie_mp_corrs <- as.data.frame(defie_mp$str) %>%
  dplyr::select(1:length(levels(defie$map_loaded$Salt))) %>%
  mutate(Genus = rownames(.)) %>%
  filter(Genus %in% defie_mp_results$Genus) %>%
  set_names(c("Freshwater", "Oligohaline", "Mesohaline", "Genus"))
# Add corrs and taxonomy
defie_mp_results <- defie_mp_results %>%
  filter(Group == "Mesohaline") %>%
  left_join(., defie_genus_all, by = "Genus") %>%
  left_join(., defie_mp_corrs, by = "Genus") %>%
  dplyr::select(taxonomy2, Genus, Group, names(defie_mp_corrs)[1:3], "RelAbundance") %>%
  arrange(Group, taxonomy2) %>%
  mutate(Genus = factor(Genus, levels = Genus)) %>%
  set_names(c("Phylum", "Genus", "Group", names(defie_mp_corrs)[1:3], "RelAbundance")) %>%
  mutate(taxon = paste(Phylum, Genus, sep = "; ")) %>%
  mutate(taxon = factor(taxon, levels = taxon))
# Plot heatmap with abundance barplot
defie.hm.melted <- defie_mp_results %>%
  dplyr::select(taxon, names(defie_mp_corrs)[1:3]) %>%
  melt(., id.vars = c("taxon"))
defie.hm <- ggplot(data = defie.hm.melted, 
                aes(x = factor(taxon), y = variable, fill = value)) + 
  geom_tile() + 
  scale_fill_distiller(name = "Indicator correlation index", palette = "RdBu", 
                       direction = -1, na.value = "transparent", type = "div", limits = c(-1, 1)) +
  scale_x_discrete(breaks = unique(defie.hm.melted$taxon), labels = unique(defie.hm.melted$taxon),
                   limits = rev(levels(defie.hm.melted$taxon))) +
  scale_y_discrete(position = "right") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        legend.title = element_text(size = 10), legend.key.height = unit(0.5,"cm"),
        legend.key.width = unit(1, "cm"), legend.text = element_text(size = 7)) +
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5)) +
  coord_flip()
defie.l <- get_legend(defie.hm)
defie.hm.clean <- defie.hm +
  theme(axis.title.y = element_blank(), axis.text.y = element_text(margin = margin(c(0,-5,0,0))),
        axis.ticks.y = element_blank(), axis.ticks.x = element_blank(), 
        panel.grid.major = element_blank(), legend.position="none",
        axis.text.x.top = element_text(size = 6, margin = margin(c(0,0,-2,0))), 
        plot.margin = margin(c(0,-2,0,5))) +
  labs(y = "Group")
defie.bp.y <- ggplot(data = defie_mp_results, aes(x = taxon, y = RelAbundance)) + 
  geom_bar(stat = "identity", fill = "grey") + 
  scale_y_continuous(position = "right") +
  scale_x_discrete(limits = rev(levels(defie_mp_results$taxon))) +
  coord_flip() + 
  theme_minimal() +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position="none", plot.margin = margin(c(0,5,0,-5)),
        axis.text.x.top = element_text(size = 6, margin = margin(c(0,0,-2,0)))) +
  labs(y = "Rel. abund. (%)")
pdf("InitialFigs/DEfie_Multipatt_genus.pdf", width = 8, height = 5)
plot_grid(defie.hm.clean, defie.bp.y, defie.l, NULL, nrow = 2, ncol = 2, 
          rel_widths = c(10,2), rel_heights = c(12, 2), align = "hv", axis = "b")
dev.off()



#### __BGC ####
# Add biogeochem. analysis
# Get variables, make corrplot, envfit, ordination with vectors, ordistep, taxa correlations
# Note, no GHG flux for Delaware field samples

# Get variables
env <- defie$map_loaded %>%
  dplyr::select(Salinity_ppt_all, 
                NH4_mgL, PO4_mgL, Cl_mgL, SO4_mgL,
                Fe_mgL, Porosity, Acetate_mgL, TotalVFA_uM, SR_umol_cm3_d, AMG_umol_cm3_d)
env_nona <- na.omit(env)

# Corrplot
C <- cor(env_nona)
pdf("InitialFigs/DEfie_BGC_corr.pdf", width = 7, height = 5)
corrplot(C, method = "number", type = "lower", 
         tl.cex = 0.5, tl.col = "black", number.cex = 0.5)
dev.off()

# Envfit
pcoa <- cmdscale(defie_bc, k = nrow(defie$map_loaded) - 1, eig = T)
set.seed(100)
ef <- envfit(pcoa, env, permutations = 999, na.rm = TRUE)
ef
arrow_factor <- ordiArrowMul(ef)
manual_factor <- 0.45
vec.df <- as.data.frame(ef$vectors$arrows*sqrt(ef$vectors$r)) %>%
  mutate(Dim1 = Dim1 * manual_factor,
         Dim2 = Dim2 * manual_factor) %>%
  mutate(variables = rownames(.)) %>%
  filter(ef$vectors$pvals < 0.05) %>%
  filter(variables != "Cl_mgL") %>%
  mutate(shortnames = c("Salinity", "NH4", "PO4", "SO4", "Fe", "SR"))

# Plot with significant vectors
ordiplot(pcoa)
plot(ef, p.max = 0.05, cex = 0.5)
pcoaA1 <- round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[1]*100, digits = 1)
pcoaA2 <- round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[2]*100, digits = 1)
defie$map_loaded$Axis01 <- scores(pcoa)[,1]
defie$map_loaded$Axis02 <- scores(pcoa)[,2]
micro.hulls <- ddply(defie$map_loaded, c("Salt"), find_hull)
g <- ggplot(defie$map_loaded, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = Salt, fill = Salt),
               alpha = 0.1, show.legend = F) +
  geom_point(size = 3, alpha = 1, aes(colour = Salt, shape = Depth),
             show.legend = T) +
  labs(x = paste("PC1: ", pcoaA1, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2, "%", sep = ""),
       shape = "Depth (m)",
       colour = "Salinity",
       fill = "Salinity") +
  scale_fill_viridis_d() +
  scale_colour_viridis_d(guide = guide_legend(override.aes = list(shape = 15))) +
  theme_bw() +  
  theme(legend.position = "right",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10)) +
  geom_segment(data = vec.df,
               aes(x = 0, xend = Dim1, y = 0, yend = Dim2),
               arrow = arrow(length = unit(0.5, "cm")),
               colour = "gray", alpha = 0.5,
               inherit.aes = FALSE) + 
  geom_text(data = vec.df,
            aes(x = Dim1, y = Dim2, label = shortnames),
            size = 3.5, color = "black") +
  coord_fixed()
g
pdf("InitialFigs/DEfie_PCoA_wBGC.pdf", width = 6, height = 4)
g
dev.off()

# Ordistep
comm_nona <- as.data.frame(t(defie$data_loaded)) %>%
  filter(rownames(.) %in% rownames(env_nona))
mod0 <- rda(comm_nona ~ 1, env_nona)  # Model with intercept only
mod1 <- rda(comm_nona ~ ., env_nona)  # Model with all explanatory variables
mod <- ordistep(mod0, scope = formula(mod1))
mod
mod$anova # Cl



#### _Delaware inc ####
input_filt_rare <- readRDS("input_filt.rds")
deinc <- filter_data(input_filt_rare,
                    filter_cat = "Estuary",
                    keep_vals = "Delaware")
deinc <- filter_data(deinc,
                    filter_cat = "Experiment",
                    keep_vals = "Soil incubation") # 28 samples

deinc$map_loaded$Depth <- as.factor(deinc$map_loaded$Depth)

# Diagnose outliers and errors
deinc_bc <- calc_dm(deinc$data_loaded)
deinc_pcoa <- cmdscale(deinc_bc, k = nrow(deinc$map_loaded) - 1, eig = T)
eigenvals(deinc_pcoa)/sum(eigenvals(deinc_pcoa)) # 30.3, 17.1 % variation explained
deinc$map_loaded$Axis01 <- scores(deinc_pcoa)[,1]
deinc$map_loaded$Axis02 <- scores(deinc_pcoa)[,2]
micro.hulls <- ddply(deinc$map_loaded, c("Treatment"), find_hull)
ggplotly(ggplot(deinc$map_loaded, aes(Axis01, Axis02)) +
           geom_polygon(data = micro.hulls, 
                        aes(colour = Treatment, fill = Treatment),
                        alpha = 0.1, show.legend = F) +
           geom_point(size = 3, alpha = 1, aes(colour = Treatment, shape = Depth)) +
           labs(x = "PC1: 30.3%", 
                y = "PC2: 17.1%",
                shape = "Depth (m)") +
           scale_colour_viridis_d() +
           scale_fill_viridis_d() +
           theme_bw() +  
           theme(legend.position = c(0,0),
                 legend.justification = c(0,0),
                 legend.background = element_blank(),
                 axis.title = element_text(face = "bold", size = 12), 
                 axis.text = element_text(size = 10),
                 plot.margin = margin(5, 5, 5, 5, "pt")))

# Looks like there are two outliers and potentially 2 with depths mixed up!
# Restart
# input_filt <- readRDS("input_filt.rds")
input_filt <- readRDS("input_filt_wBGC.rds")
deinc <- filter_data(input_filt,
                     filter_cat = "Estuary",
                     keep_vals = "Delaware")
deinc <- filter_data(deinc,
                     filter_cat = "Experiment",
                     keep_vals = "Soil incubation") # 28 samples
deinc$map_loaded$sampleID <- rownames(deinc$map_loaded)
deinc <- filter_data(deinc, 
                     filter_cat = "sampleID", 
                     filter_vals = c("TS_FW_d1_12_2", "TS_FW_d2_12_2"))
for (i in 1:nrow(deinc$map_loaded)) {
  if (deinc$map_loaded$sampleID[i] == "TS_FW_d2_04_1") {
    deinc$map_loaded$Depth[i] <- "0.02"
  }
}
for (i in 1:nrow(deinc$map_loaded)) {
  if (deinc$map_loaded$sampleID[i] == "TS_FW_d1_04_2") {
    deinc$map_loaded$Depth[i] <- "0.12"
  }
}

set.seed(530)
deinc <- single_rarefy(deinc, min(colSums(deinc$data_loaded))) # 31264
deinc$map_loaded <- deinc$map_loaded %>%
  mutate(rich = specnumber(deinc$data_loaded, MARGIN = 2),
         shannon = diversity(deinc$data_loaded, index = "shannon", MARGIN = 2),
         Depth = as.factor(Depth),
         sampleID = rownames(.)) %>%
  mutate(Treatment = gsub("5 ppt ASW", "ASW", Treatment)) %>%
  mutate(Treatment = gsub("wk", "", Treatment)) %>%
  mutate(Treatment = gsub("Freshwater", "Fresh", Treatment)) %>%
  separate(Treatment, into = c("Salt", "Time"), sep = " ", remove = F) %>%
  mutate(Time = replace_na(Time, 0)) %>%
  mutate(Time = as.integer(Time))
for (i in 1:nrow(deinc$map_loaded)) {
  if (deinc$map_loaded$Time[i] == 0) {
    deinc$map_loaded$Treatment[i] <- "Initial"
    deinc$map_loaded$Salt[i] <- "Initial"
  }
}
deinc$map_loaded <- deinc$map_loaded %>%
  mutate_if(is.character, as.factor) %>%
  mutate(Treatment = factor(Treatment,
                          levels = c("Initial", "Fresh 4", "Fresh 7",
                                     "Fresh 12", "ASW 4", "ASW 7", "ASW 12")),
         TrtDepth = paste(deinc$map_loaded$Treatment, deinc$map_loaded$Depth, sep = ""))

#### __Alpha ####
leveneTest(deinc$map_loaded$rich ~ deinc$map_loaded$Treatment) # Homogeneous
m <- aov(rich ~ Treatment + Depth, data = deinc$map_loaded)
Anova(m, type = "II") # Neither
m <- aov(rich ~ Treatment, data = deinc$map_loaded)
shapiro.test(m$residuals) # Normal
summary(m)
t <- emmeans(object = m, specs = "Treatment") %>%
  cld(object = ., adjust = "Tukey", Letters = letters, alpha = 0.05) %>%
  mutate(name = "rich",
         y = max(deinc$map_loaded$rich)+(max(deinc$map_loaded$rich)-min(deinc$map_loaded$rich))/20)
leveneTest(deinc$map_loaded$shannon ~ deinc$map_loaded$Treatment) # Not homogeneous
m1 <- aov(shannon ~ Treatment + Depth, data = deinc$map_loaded)
Anova(m1, type = "II") # Neither
m1 <- aov(shannon ~ Treatment, data = deinc$map_loaded)
shapiro.test(m1$residuals) # Normal
summary(m1)
t1 <- emmeans(object = m1, specs = "Treatment") %>%
  cld(object = ., adjust = "Tukey", Letters = letters, alpha = 0.05) %>%
  mutate(name = "shannon",
         y = max(deinc$map_loaded$shannon)+(max(deinc$map_loaded$shannon)-min(deinc$map_loaded$shannon))/20)
label_df <- rbind(t, t1)
facet_df <- c("rich" = "(a) Richness",
              "shannon" = "(b) Shannon")
alpha_long <- deinc$map_loaded %>%
  pivot_longer(cols = c("rich", "shannon"))
pdf("InitialFigs/DEinc_Alpha.pdf", width = 7, height = 3)
ggplot(alpha_long, aes(Treatment, value, 
                       colour = Salt)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 2, alpha = 0.75, width = 0.2, aes(shape = Depth)) +
  geom_text(data = label_df, aes(Treatment, y, label = str_trim(.group)), 
            size = 4, color = "black") +
  labs(x = "Site", y = "Number of OTUs", shape = "Depth (m)") +
  scale_colour_viridis_d() +
  guides(colour = "none") +
  facet_wrap(~ name, ncol = 2, scales = "free_y", labeller = as_labeller(facet_df)) +
  theme_bw() +
  theme(legend.position = c(0,0),
        legend.justification = c(0,0),
        legend.background = element_blank(),
        legend.title.align = 0.5,
        axis.title = element_blank(),
        axis.text = element_text(size = 8),
        strip.text = element_text(size = 10))
dev.off()

#### __Beta ####
# From here on just look at Initial and week 12
deinc <- filter_data(deinc, 
                     filter_cat = "Time",
                     keep_vals = c("0", "12"))
deinc_bc <- calc_dm(deinc$data_loaded)
deinc_pcoa <- cmdscale(deinc_bc, k = nrow(deinc$map_loaded) - 1, eig = T)
pcoaA1 <- round((eigenvals(deinc_pcoa)/sum(eigenvals(deinc_pcoa)))[1]*100, digits = 1)
pcoaA2 <- round((eigenvals(deinc_pcoa)/sum(eigenvals(deinc_pcoa)))[2]*100, digits = 1)
deinc$map_loaded$Axis01 <- scores(deinc_pcoa)[,1]
deinc$map_loaded$Axis02 <- scores(deinc_pcoa)[,2]
micro.hulls <- ddply(deinc$map_loaded, c("Salt"), find_hull)
pdf("InitialFigs/DEinc_PCoA.pdf", width = 7, height = 5)
ggplot(deinc$map_loaded, aes(Axis01, Axis02)) +
           geom_point(size = 3, alpha = 1, aes(colour = Salt, shape = Depth)) +
           labs(x = paste("PC1: ", pcoaA1, "%", sep = ""), 
                y = paste("PC2: ", pcoaA2, "%", sep = ""),
                shape = "Depth (m)") +
           scale_colour_viridis_d() +
           scale_fill_viridis_d() +
           theme_bw() +  
           theme(legend.position = "right",
                 axis.title = element_text(face = "bold", size = 12), 
                 axis.text = element_text(size = 10),
                 plot.margin = margin(5, 5, 5, 5, "pt"))
dev.off()

set.seed(1150)
adonis2(deinc_bc ~ Salt + Depth, data = deinc$map_loaded) # Depth
anova(betadisper(deinc_bc, deinc$map_loaded$Treatment)) # Dispersion homogeneous
anova(betadisper(deinc_bc, deinc$map_loaded$Depth)) # Dispersion homogeneous
anova(betadisper(deinc_bc, deinc$map_loaded$Time)) # Dispersion homogeneous

#### __Taxa ####
bar_text <- data.frame(group_by = c("Fresh 120.02", "Fresh 120.12"),
                       y = c(1.05, 1.05),
                       label = c("|-----2 cm-----|",
                                 "|-----12 cm-----|"))
deinc_phyla <- summarize_taxonomy(deinc, level = 2, report_higher_tax = F)
plot_ts_heatmap(deinc_phyla, deinc$map_loaded, 0.01, 'TrtDepth', rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
plot_taxa_bars(deinc_phyla, deinc$map_loaded, "TrtDepth", 
                              num_taxa = 12, data_only = F) +
  theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
deinc_barsP <- plot_taxa_bars(deinc_phyla, deinc$map_loaded, "TrtDepth", 
                              num_taxa = 12, data_only = T) %>%
  mutate(group_by = factor(group_by, levels = c("Initial0.02","Fresh 120.02",
                                                "ASW 120.02", "Initial0.12",
                                                "Fresh 120.12", "ASW 120.12"))) %>%
  mutate(taxon = gsub("NA", "Unclassified", taxon)) %>%
  mutate(taxon = fct_relevel(taxon, "Other", after = Inf)) %>%
  mutate(taxon = fct_relevel(taxon, "Unclassified", after = Inf)) %>%
  mutate(taxon = fct_rev(taxon))
pdf("InitialFigs/DEinc_Phyla.pdf", width = 7, height = 5)
ggplot(deinc_barsP, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  geom_text(data = bar_text,
            aes(x = group_by, y = y, label = label),
            inherit.aes = F) +
  labs(x = "Sample", y = "Relative abundance", fill = "Phylum") +
  scale_fill_manual(values = c("grey75", "grey90", brewer.pal(12, "Paired")[1:11])) +
  scale_y_continuous(expand = c(0.01, 0.01)) + 
  scale_x_discrete(labels = c("Initial", "+Fresh", "+ASW", "Initial", "+Fresh", "+ASW")) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10))
dev.off()
taxa_summary_by_sample_type(deinc_phyla, deinc$map_loaded, 'Treatment', 0.01, 'KW')

deinc_class <- summarize_taxonomy(deinc, level = 3, report_higher_tax = F)
plot_ts_heatmap(deinc_class, deinc$map_loaded, 0.01, 'Treatment', rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
deinc_barsC <- plot_taxa_bars(deinc_class, deinc$map_loaded, "TrtDepth", 
                              num_taxa = 12, data_only = T) %>%
  mutate(group_by = factor(group_by, levels = c("Initial0.02","Fresh 120.02",
                                                "ASW 120.02", "Initial0.12",
                                                "Fresh 120.12", "ASW 120.12"))) %>%
  mutate(taxon = gsub("NA", "Unclassified", taxon)) %>%
  mutate(taxon = fct_relevel(taxon, "Other", after = Inf)) %>%
  mutate(taxon = fct_relevel(taxon, "Unclassified", after = Inf)) %>%
  mutate(taxon = fct_rev(taxon))
ggplot(deinc_barsC, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  geom_text(data = bar_text,
            aes(x = group_by, y = y, label = label),
            inherit.aes = F) +
  labs(x = "Sample", y = "Relative abundance", fill = "Class") +
  scale_fill_manual(values = c("grey75", "grey90", brewer.pal(12, "Paired")[1:11])) +
  scale_y_continuous(expand = c(0.01, 0.01)) + 
  scale_x_discrete(labels = c("Initial", "Fresh", "ASW", "Initial", "Fresh", "ASW")) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10))
taxa_summary_by_sample_type(deinc_class, deinc$map_loaded, 'Treatment', 0.01, 'KW')

deinc_order <- summarize_taxonomy(deinc, level = 4, report_higher_tax = F)
plot_ts_heatmap(deinc_order, deinc$map_loaded, 0.01, 'Treatment', rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
deinc_barsO <- plot_taxa_bars(deinc_order, deinc$map_loaded, "TrtDepth", 
                              num_taxa = 12, data_only = T) %>%
  mutate(group_by = factor(group_by, levels = c("Initial0.02","Fresh 120.02",
                                                "ASW 120.02", "Initial0.12",
                                                "Fresh 120.12", "ASW 120.12"))) %>%
  mutate(taxon = gsub("NA", "Unclassified", taxon)) %>%
  mutate(taxon = fct_relevel(taxon, "Other", after = Inf)) %>%
  mutate(taxon = fct_relevel(taxon, "Unclassified", after = Inf)) %>%
  mutate(taxon = fct_rev(taxon))
ggplot(deinc_barsO, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  geom_text(data = bar_text,
            aes(x = group_by, y = y, label = label),
            inherit.aes = F) +
  labs(x = "Sample", y = "Relative abundance", fill = "Order") +
  scale_fill_manual(values = c("grey75", "grey90", brewer.pal(12, "Paired")[1:11])) +
  scale_y_continuous(expand = c(0.01, 0.01)) + 
  scale_x_discrete(labels = c("Initial", "Fresh", "ASW", "Initial", "Fresh", "ASW")) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10))
taxa_summary_by_sample_type(deinc_order, deinc$map_loaded, 'Treatment', 0.01, 'KW')

deinc_family <- summarize_taxonomy(deinc, level = 5, report_higher_tax = F)
plot_ts_heatmap(deinc_family, deinc$map_loaded, 0.01, 'Treatment', rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
deinc_barsF <- plot_taxa_bars(deinc_family, deinc$map_loaded, "TrtDepth", 
                              num_taxa = 12, data_only = T) %>%
  mutate(group_by = factor(group_by, levels = c("Initial0.02","Fresh 120.02",
                                                "ASW 120.02", "Initial0.12",
                                                "Fresh 120.12", "ASW 120.12"))) %>%
  mutate(taxon = gsub("NA", "Unclassified", taxon)) %>%
  mutate(taxon = fct_relevel(taxon, "Other", after = Inf)) %>%
  mutate(taxon = fct_relevel(taxon, "Unclassified", after = Inf)) %>%
  mutate(taxon = fct_rev(taxon))
ggplot(deinc_barsF, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  geom_text(data = bar_text,
            aes(x = group_by, y = y, label = label),
            inherit.aes = F) +
  labs(x = "Sample", y = "Relative abundance", fill = "Family") +
  scale_fill_manual(values = c("grey75", "grey90", brewer.pal(12, "Paired")[1:11])) +
  scale_y_continuous(expand = c(0.01, 0.01)) + 
  scale_x_discrete(labels = c("Initial", "Fresh", "ASW", "Initial", "Fresh", "ASW")) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10))
taxa_summary_by_sample_type(deinc_family, deinc$map_loaded, 'Treatment', 0.01, 'KW')

deinc_genus<- summarize_taxonomy(deinc, level = 6, report_higher_tax = F)
plot_ts_heatmap(deinc_genus, deinc$map_loaded, 0.01, 'Treatment', rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
deinc_barsG <- plot_taxa_bars(deinc_genus, deinc$map_loaded, "TrtDepth", 
                              num_taxa = 12, data_only = T) %>%
  mutate(group_by = factor(group_by, levels = c("Initial0.02","Fresh 120.02",
                                                "ASW 120.02", "Initial0.12",
                                                "Fresh 120.12", "ASW 120.12"))) %>%
  mutate(taxon = gsub("NA", "Unclassified", taxon)) %>%
  mutate(taxon = fct_relevel(taxon, "Other", after = Inf)) %>%
  mutate(taxon = fct_relevel(taxon, "Unclassified", after = Inf)) %>%
  mutate(taxon = fct_rev(taxon))
ggplot(deinc_barsG, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  geom_text(data = bar_text,
            aes(x = group_by, y = y, label = label),
            inherit.aes = F) +
  labs(x = "Sample", y = "Relative abundance", fill = "Genus") +
  scale_fill_manual(values = c("grey75", "grey90", brewer.pal(12, "Paired")[1:11])) +
  scale_y_continuous(expand = c(0.01, 0.01)) + 
  scale_x_discrete(labels = c("Initial", "Fresh", "ASW", "Initial", "Fresh", "ASW")) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10))
taxa_summary_by_sample_type(deinc_genus, deinc$map_loaded, 'Treatment', 0.01, 'KW')

deinc_guilds <- summarize_taxonomy(deinc, level = 9, report_higher_tax = F)
plot_ts_heatmap(deinc_guilds, deinc$map_loaded, 0.01, 'Treatment', rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
deinc_barsGu <- plot_taxa_bars(deinc_guilds, deinc$map_loaded, "TrtDepth", 
                               num_taxa = 20, data_only = T) %>%
  filter(taxon != "NA") %>%
  droplevels() %>%
  mutate(taxon = factor(taxon,
                        levels = Guild_cols$Guild)) %>%
  mutate(group_by = factor(group_by, levels = c("Initial0.02","Fresh 120.02",
                                                "ASW 120.02", "Initial0.12",
                                                "Fresh 120.12", "ASW 120.12")))
tallest_bar <- deinc_barsGu %>%
  group_by(group_by) %>%
  summarise(sum = sum(mean_value))
bar_textGu <- data.frame(group_by = c("Fresh 120.02", "Fresh 120.12"),
                         y = c(max(tallest_bar$sum) + 0.01, max(tallest_bar$sum) + 0.01),
                         label = c("|-----2 cm-----|",
                                   "|-----12 cm-----|"))
pdf("InitialFigs/DEinc_Guilds.pdf", width = 7, height = 5)
ggplot(deinc_barsGu, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  geom_text(data = bar_textGu,
            aes(x = group_by, y = y, label = label),
            inherit.aes = F) +
  labs(x = "TrtDepth", y = "Relative abundance", fill = "Guild") +
  scale_fill_manual(values = Guild_cols$color) +
  scale_y_continuous(expand = c(max(tallest_bar$sum)/100, max(tallest_bar$sum)/100)) + 
  scale_x_discrete(labels = c("Initial", "+Fresh", "+ASW", "Initial", "+Fresh", "+ASW")) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10))
dev.off()
taxa_summary_by_sample_type(deinc_guilds, 
                            deinc$map_loaded, 
                            type_header = 'Treatment', 
                            filter_level = 0.01, 
                            test_type = 'KW')

# Guilds by sample (to show variation)
deinc$map_loaded$sampleID <- rownames(deinc$map_loaded)
deinc_barsGu <- plot_taxa_bars(deinc_guilds,
                               deinc$map_loaded,
                               "sampleID",
                               num_taxa = 20,
                               data_only = TRUE) %>%
  filter(taxon != "NA") %>%
  droplevels() %>%
  mutate(taxon = factor(taxon,
                        levels = Guild_cols$Guild)) %>%
  left_join(., deinc$map_loaded, by = c("group_by" = "sampleID"))
facet_names <- c("0.02" = "Depth 2 cm",
                 "0.12" = "Depth 12 cm",
                 "Initial" = "Initial",
                 "Fresh" = "+Fresh",
                 "ASW" = "+ASW")
tallest_bar <- deinc_barsGu %>%
  group_by(group_by) %>%
  summarise(sum = sum(mean_value))
pdf("InitialFigs/DEinc_Guilds_samples.pdf", width = 9, height = 5)
ggplot(deinc_barsGu, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Sample", y = "Relative abundance", fill = "Guild") +
  scale_fill_manual(values = Guild_cols$color) +
  scale_y_continuous(expand = c(max(tallest_bar$sum)/100, max(tallest_bar$sum)/100)) + 
  facet_nested(~ Depth + Salt, space = "free", scales = "free_x",
               labeller = as_labeller(facet_names)) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 90, hjust = 1, vjust = 0.5),
        strip.background = element_rect(size = 0.2),
        axis.line.y = element_blank())
dev.off()

#### __Simper ####
deinc_sim <- simper(t(deinc$data_loaded), 
                    deinc$map_loaded$TrtDepth)
deinc_s <- summary(deinc_sim)
head(deinc_s$`Fresh 120.02_ASW 120.02`)
head(deinc_s$`Fresh 120.12_ASW 120.12`)
deinc_df1 <- head(deinc_s$`Fresh 120.02_ASW 120.02`, n = 20) %>%
  mutate(Comparison = "ASW/Control 2 cm",
         ASV = rownames(.))
deinc_df2 <- head(deinc_s$`Fresh 120.12_ASW 120.12`, n = 20) %>%
  mutate(Comparison = "ASW/Control 12 cm",
         ASV = rownames(.))
deinc_simper_results <- rbind(deinc_df1, deinc_df2) %>%
  left_join(., deinc$taxonomy_loaded, by = c("ASV" = "taxonomy8")) %>%
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
           "Meadeincontrol" = "avb",
           "CumulativeContribution" = "cumsum")) %>%
  dplyr::select(Comparison, SaltResponse, Domain, Phylum, Class, Order, Family, Genus,
                Species, OTU, MeanSalt, Meadeincontrol, CumulativeContribution)
write_xlsx(deinc_simper_results, 
           "simper_results_DEinc.xlsx",
           format_headers = F)

#### __Multipatt ####
set.seed(1202)
deinc_mp <- multipatt(t(deinc$data_loaded), 
                   deinc$map_loaded$TrtDepth, 
                   func = "r.g", 
                   control = how(nperm=999))
# Note nothing with significant q value, no strong correlations
deinc_mp_results <- deinc_mp$sign %>%
  mutate(num_sites = rowSums(cbind(`s.ASW 120.02`, `s.ASW 120.12`, `s.Fresh 120.02`,
                                   `s.Fresh 120.12`, `s.Initial0.02`, `s.Initial0.12`)),
         q.value = qvalue(deinc_mp$sign$p.value)$qvalues,
         Group = "NA",
         ASV = rownames(deinc_mp$sign)) %>%
  filter(p.value < 0.05) %>%
  filter(stat >= 0.3) %>%
  filter(num_sites < 2) %>%
  left_join(., deinc$taxonomy_loaded, by = c("ASV" = "taxonomy8"))
for (i in 1:nrow(deinc_mp_results)) {
  if (deinc_mp_results$num_sites[i] == 1 & deinc_mp_results$`s.ASW 120.02`[i] == 1) {
    deinc_mp_results$Group[i] <- "ASW 2 cm"
  }
}
for (i in 1:nrow(deinc_mp_results)) {
  if (deinc_mp_results$num_sites[i] == 1 & deinc_mp_results$`s.ASW 120.12`[i] == 1) {
    deinc_mp_results$Group[i] <- "ASW 12 cm"
  }
}
for (i in 1:nrow(deinc_mp_results)) {
  if (deinc_mp_results$num_sites[i] == 1 & deinc_mp_results$`s.Fresh 120.02`[i] == 1) {
    deinc_mp_results$Group[i] <- "Fresh 2 cm"
  }
}
for (i in 1:nrow(deinc_mp_results)) {
  if (deinc_mp_results$num_sites[i] == 1 & deinc_mp_results$`s.Fresh 120.12`[i] == 1) {
    deinc_mp_results$Group[i] <- "Fresh 12 cm"
  }
}
for (i in 1:nrow(deinc_mp_results)) {
  if (deinc_mp_results$num_sites[i] == 1 & deinc_mp_results$s.Initial0.02[i] == 1) {
    deinc_mp_results$Group[i] <- "Initial 2 cm"
  }
}
for (i in 1:nrow(deinc_mp_results)) {
  if (deinc_mp_results$num_sites[i] == 1 & deinc_mp_results$s.Initial0.12[i] == 1) {
    deinc_mp_results$Group[i] <- "Initial 12 cm"
  }
}
table(deinc_mp_results$Group)
deinc_asv <- summarize_taxonomy(deinc, level = 8, report_higher_tax = F)
deinc_asv_all <- data.frame("RelAbundance" = round(rowMeans(deinc_asv) * 100, digits = 4)) %>%
  mutate(ASV = rownames(.))
deinc_mp_corrs <- as.data.frame(deinc_mp$str) %>%
  dplyr::select(1:6) %>%
  mutate(ASV = rownames(.)) %>%
  filter(ASV %in% deinc_mp_results$ASV) %>%
  set_names(c("ASW 2 cm", "ASW 12 cm", "Fresh 2 cm", "Fresh 12 cm", "Initial 2 cm", 
              "Initial 12 cm", "ASV"))
# Add corrs and taxonomy
deinc_mp_results <- deinc_mp_results %>%
  filter(Group == "ASW 2 cm" | Group == "ASW 12 cm") %>%
  left_join(., deinc_asv_all, by = "ASV") %>%
  left_join(., deinc_mp_corrs, by = "ASV") %>%
  dplyr::select(taxonomy2, taxonomy6, ASV, Group, "ASW 2 cm", "ASW 12 cm", 
                "Fresh 2 cm", "Fresh 12 cm", "Initial 2 cm", "Initial 12 cm", "RelAbundance") %>%
  mutate(ASV = gsub("ASV", "OTU", ASV)) %>%
  arrange(Group, taxonomy2, taxonomy6) %>%
  mutate(ASV = factor(ASV, levels = ASV)) %>%
  set_names(c("Phylum", "Genus", "ASV", "Group", "ASW 2 cm", "ASW 12 cm", 
              "Fresh 2 cm", "Fresh 12 cm", "Initial 2 cm", "Initial 12 cm", "RelAbundance")) %>%
  mutate(taxon = paste(Phylum, Genus, ASV, sep = "; ")) %>%
  mutate(taxon = factor(taxon, levels = taxon))
# Plot heatmap with abundance barplot
hm.melted <- deinc_mp_results %>%
  dplyr::select(taxon, "ASW 2 cm", "ASW 12 cm", "Fresh 2 cm", "Fresh 12 cm", "Initial 2 cm", 
                "Initial 12 cm") %>%
  melt(., id.vars = c("taxon"))
hm <- ggplot(data = hm.melted, 
             aes(x = factor(taxon), y = variable, fill = value)) + 
  geom_tile() + 
  scale_fill_distiller(name = "Indicator correlation index", palette = "RdBu", 
                       direction = -1, na.value = "transparent", type = "div", limits = c(-1, 1)) +
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
  theme(axis.title.y = element_blank(), axis.text.y = element_text(margin = margin(c(0,-5,0,0))),
        axis.ticks.y = element_blank(), axis.ticks.x = element_blank(), 
        panel.grid.major = element_blank(), legend.position="none",
        axis.text.x.top = element_text(size = 6, margin = margin(c(0,0,-3,0))), 
        plot.margin = margin(c(0,-2,0,5))) +
  labs(y = "Group")
bp.y <- ggplot(data = deinc_mp_results, aes(x = taxon, y = RelAbundance)) + 
  geom_bar(stat = "identity", fill = "grey") + 
  scale_y_continuous(position = "right") +
  scale_x_discrete(limits = rev(levels(deinc_mp_results$taxon))) +
  coord_flip() + 
  theme_minimal() +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position="none", axis.text.x.top = element_text(size = 6, margin = margin(c(0,0,-3,0))), 
        plot.margin = margin(c(0,-2,0,-5))) +
  labs(y = "Rel. abund. (%)")
pdf("InitialFigs/DEinc_Multipatt.pdf", width = 8, height = 5)
plot_grid(hm.clean, bp.y, l, NULL, nrow = 2, ncol = 2, 
          rel_widths = c(10,2), rel_heights = c(12, 2), align = "hv", axis = "b")
dev.off()



#### __BGC ####
# Add biogeochem. analysis
# Get variables, make corrplot, envfit, ordination with vectors, ordistep, taxa correlations
# Note, no data for Delware inc yet, but get set up.

# Get variables
env <- deinc$map_loaded %>%
  dplyr::select(Salinity_ppt_all, 
                NH4_mgL, PO4_mgL, Cl_mgL, SO4_mgL,
                Fe_mgL, Porosity, Acetate_mgL, TotalVFA_uM, SR_umol_cm3_d, AMG_umol_cm3_d)
env_nona <- na.omit(env)

# Corrplot
C <- cor(env_nona)
pdf("InitialFigs/DEinc_BGC_corr.pdf", width = 7, height = 5)
corrplot(C, method = "number", type = "lower", 
         tl.cex = 0.5, tl.col = "black", number.cex = 0.5)
dev.off()

# Envfit
pcoa <- cmdscale(deinc_bc, k = nrow(deinc$map_loaded) - 1, eig = T)
set.seed(100)
ef <- envfit(pcoa, env, permutations = 999, na.rm = TRUE)
ef
arrow_factor <- ordiArrowMul(ef)
manual_factor <- 0.45
vec.df <- as.data.frame(ef$vectors$arrows*sqrt(ef$vectors$r)) %>%
  mutate(Dim1 = Dim1 * manual_factor,
         Dim2 = Dim2 * manual_factor) %>%
  mutate(variables = rownames(.)) %>%
  filter(ef$vectors$pvals < 0.05) %>%
  filter(variables != "Cl_mgL") %>%
  mutate(shortnames = c("Salinity", "NH4", "PO4", "SO4", "Fe", "SR"))

# Plot with significant vectors
ordiplot(pcoa)
plot(ef, p.max = 0.05, cex = 0.5)
pcoaA1 <- round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[1]*100, digits = 1)
pcoaA2 <- round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[2]*100, digits = 1)
deinc$map_loaded$Axis01 <- scores(pcoa)[,1]
deinc$map_loaded$Axis02 <- scores(pcoa)[,2]
micro.hulls <- ddply(deinc$map_loaded, c("Salt"), find_hull)
g <- ggplot(deinc$map_loaded, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = Salt, fill = Salt),
               alpha = 0.1, show.legend = F) +
  geom_point(size = 3, alpha = 1, aes(colour = Salt, shape = Depth),
             show.legend = T) +
  labs(x = paste("PC1: ", pcoaA1, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2, "%", sep = ""),
       shape = "Depth (m)",
       colour = "Salinity") +
  scale_fill_viridis_d() +
  scale_colour_viridis_d() +
  theme_bw() +  
  theme(legend.position = "right",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10)) +
  geom_segment(data = vec.df,
               aes(x = 0, xend = Dim1, y = 0, yend = Dim2),
               arrow = arrow(length = unit(0.5, "cm")),
               colour = "gray", alpha = 0.5,
               inherit.aes = FALSE) + 
  geom_text(data = vec.df,
            aes(x = Dim1, y = Dim2, label = shortnames),
            size = 3.5, color = "black") +
  coord_fixed()
g
pdf("InitialFigs/DEinc_PCoA_wBGC.pdf", width = 6, height = 4)
g
dev.off()

# Ordistep
comm_nona <- as.data.frame(t(defie$data_loaded)) %>%
  filter(rownames(.) %in% rownames(env_nona))
mod0 <- rda(comm_nona ~ 1, env_nona)  # Model with intercept only
mod1 <- rda(comm_nona ~ ., env_nona)  # Model with all explanatory variables
mod <- ordistep(mod0, scope = formula(mod1))
mod
mod$anova #



#### _Delaware transplant ####
# input_filt <- readRDS("input_filt.rds")
input_filt <- readRDS("input_filt_wBGC.rds")
detra <- filter_data(input_filt,
                     filter_cat = "Estuary",
                     keep_vals = "Delaware")
detra <- filter_data(detra,
                     filter_cat = "Experiment",
                     keep_vals = c("Soil mesocosm"))
set.seed(530)
detra <- single_rarefy(detra, min(colSums(detra$data_loaded))) # 26450
detra$map_loaded <- detra$map_loaded %>%
  mutate(rich = specnumber(detra$data_loaded, MARGIN = 2),
         shannon = diversity(detra$data_loaded, index = "shannon", MARGIN = 2),
         Treatment = factor(Treatment,
                            levels = c("TFM1", "TFM2", "TFM1@TFM2", 
                                       "TFM1@OligoHal", "OligoHal.", "TFM1@MesoHal",
                                       "MesoHal", "TFM1-40cm", "TFM2-40cm", 
                                       "TFM1@TFM2-40cm", "TFM1@OligoHal-40cm", 
                                       "OligoHal-40cm", "TFM1@MesoHal-40cm",
                                       "MesoHal-40cm"))) %>%
  unite("TrtDepth", c("Treatment", "Depth"), sep = "", remove = F) %>%
  separate(Treatment, into = c("Treatment2", "Level"), sep = "-", remove = F) %>%
  mutate(Treatment2 = gsub("\\.", "", Treatment2)) %>%
  replace_na(list(Level = "0cm")) %>%
  mutate(Depth = as.factor(Depth),
         Level = as.factor(Level),
         Treatment2 = factor(Treatment2,
                             levels = c("TFM1", "TFM2", "TFM1@TFM2", "TFM1@OligoHal",
                                        "OligoHal", "TFM1@MesoHal", "MesoHal"))) %>%
  mutate(Treatment3 = paste(Treatment2, Level, sep = "_")) %>%
  mutate_if(is.character, as.factor)

#### __Alpha ####
leveneTest(detra$map_loaded$rich ~ detra$map_loaded$Treatment) # Homogeneous
m <- aov(rich ~ Treatment + Depth, data = detra$map_loaded)
Anova(m, type = "II") # Treatment and Depth
m <- aov(rich ~ Treatment, data = detra$map_loaded)
shapiro.test(m$residuals) # Not normal
summary(m)
TukeyHSD(m)
t <- emmeans(object = m, specs = "Treatment") %>%
  cld(object = ., adjust = "Tukey", Letters = letters, alpha = 0.05) %>%
  mutate(name = "rich",
         y = max(detra$map_loaded$rich)+(max(detra$map_loaded$rich)-min(detra$map_loaded$rich))/20)
leveneTest(detra$map_loaded$shannon ~ detra$map_loaded$Treatment) # Homogeneous
m1 <- aov(shannon ~ Treatment + Depth, data = detra$map_loaded)
Anova(m1, type = "II") # Depth
m1 <- aov(shannon ~ Treatment, data = detra$map_loaded)
shapiro.test(m1$residuals)
summary(m1)
TukeyHSD(m1)
t1 <- emmeans(object = m1, specs = "Treatment") %>%
  cld(object = ., adjust = "Tukey", Letters = letters, alpha = 0.05) %>%
  mutate(name = "shannon",
         y = max(detra$map_loaded$shannon)+(max(detra$map_loaded$shannon)-min(detra$map_loaded$shannon))/20)
label_df <- rbind(t, t1)
facet_df <- c("rich" = "(a) Richness",
              "shannon" = "(b) Shannon")
alpha_long <- detra$map_loaded %>%
  pivot_longer(cols = c("rich", "shannon"))
pdf("InitialFigs/DEtra_Alpha.pdf", width = 7, height = 4)
ggplot(alpha_long, aes(reorder(Treatment, value, mean), value, 
                       colour = Treatment)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 2, alpha = 0.75, width = 0.2, aes(shape = Depth)) +
  geom_text(data = label_df, aes(Treatment, y, label = str_trim(.group)), 
            size = 4, color = "black") +
  labs(x = "Site", y = "Number of OTUs", shape = "Depth (m)") +
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
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        strip.text = element_text(size = 10))
dev.off()

#### __Beta ####
detra_bc <- calc_dm(detra$data_loaded)
set.seed(1150)
adonis2(detra_bc ~ detra$map_loaded$Treatment+detra$map_loaded$Depth+detra$map_loaded$Level) # Both sig
anova(betadisper(detra_bc, detra$map_loaded$Treatment)) # Dispersion homogeneous
anova(betadisper(detra_bc, detra$map_loaded$Depth)) # Dispersion homogeneous
anova(betadisper(detra_bc, detra$map_loaded$Level)) # Dispersion homogeneous
detra_pcoa <- cmdscale(detra_bc, k = nrow(detra$map_loaded) - 1, eig = T)
pcoaA1 <- round((eigenvals(detra_pcoa)/sum(eigenvals(detra_pcoa)))[1]*100, digits = 1)
pcoaA2 <- round((eigenvals(detra_pcoa)/sum(eigenvals(detra_pcoa)))[2]*100, digits = 1)
detra$map_loaded$Axis01 <- scores(detra_pcoa)[,1]
detra$map_loaded$Axis02 <- scores(detra_pcoa)[,2]
micro.hulls <- ddply(detra$map_loaded, c("Treatment"), find_hull)
pdf("InitialFigs/DEtra_PCoA.pdf", width = 7, height = 5)
ggplot(detra$map_loaded, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = Treatment, fill = Treatment),
               alpha = 0.1, show.legend = F) +
  geom_point(size = 3, alpha = 1, aes(colour = Treatment, shape = Depth),
             show.legend = T) +
  labs(x = paste("PC1: ", pcoaA1, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2, "%", sep = "")) +
  scale_fill_viridis_d() +
  scale_colour_viridis_d() +
  theme_bw() +  
  theme(legend.position = "right",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10))
dev.off()

g <- ggplot(detra$map_loaded, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = Treatment, fill = Treatment),
               alpha = 0.1, show.legend = F) +
  geom_point(size = 3, alpha = 1, aes(colour = Treatment, shape = Depth),
             show.legend = T) +
  labs(x = paste("PC1: ", pcoaA1, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2, "%", sep = "")) +
  scale_fill_viridis_d() +
  scale_colour_viridis_d() +
  theme_bw() +  
  theme(legend.position = "right",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10))
ggplotly(g)

#### __Taxa ####
detra_phyla <- summarize_taxonomy(detra, level = 2, report_higher_tax = F)
plot_ts_heatmap(detra_phyla, detra$map_loaded, 0.01, 'Treatment', rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
detra_barsP <- plot_taxa_bars(detra_phyla, 
                              detra$map_loaded, 
                              "Treatment3", 
                              num_taxa = 12, 
                              data_only = T) %>%
  mutate(taxon = fct_rev(taxon)) %>%
  separate(group_by, into = c("group_by", "Level"), sep = "_") %>%
  mutate(group_by = factor(group_by,
                           levels = c("TFM1", "TFM2", "TFM1@TFM2", "TFM1@OligoHal",
                                      "OligoHal", "TFM1@MesoHal", "MesoHal")))
facet_names <- c("0cm" = "Elevation 0 cm", 
                 "40cm" = "Elevation -40 cm") 
pdf("InitialFigs/DEtra_Phyla.pdf", width = 7, height = 5)
ggplot(detra_barsP, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Sample", y = "Relative abundance", fill = "Phylum") +
  scale_fill_manual(values = c("grey90", brewer.pal(12, "Paired")[12:1])) +
  scale_y_continuous(expand = c(0.01, 0.01)) + 
  facet_wrap(~ Level, labeller = as_labeller(facet_names)) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        strip.background = element_rect(size = 0.2),
        axis.line.y = element_blank())
dev.off()
taxa_summary_by_sample_type(detra_phyla, detra$map_loaded, 'Treatment', 0.01, 'KW')

detra_class <- summarize_taxonomy(detra, level = 3, report_higher_tax = F)
plot_ts_heatmap(detra_class, detra$map_loaded, 0.01, 'Treatment', rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
detra_barsC <- plot_taxa_bars(detra_class, detra$map_loaded, "Treatment3", 
                              num_taxa = 12, data_only = T) %>%
  mutate(taxon = gsub("NA", "Unclassified", taxon)) %>%
  mutate(taxon = fct_relevel(taxon, "Other", after = Inf)) %>%
  mutate(taxon = fct_relevel(taxon, "Unclassified", after = Inf)) %>%
  mutate(taxon = fct_rev(taxon)) %>%
  separate(group_by, into = c("group_by", "Level"), sep = "_") %>%
  mutate(group_by = factor(group_by,
                           levels = c("TFM1", "TFM2", "TFM1@TFM2", "TFM1@OligoHal",
                                      "OligoHal", "TFM1@MesoHal", "MesoHal")))
ggplot(detra_barsC, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = NULL, y = "Relative abundance", fill = "Class") +
  scale_fill_manual(values = c("grey75", "grey90", brewer.pal(12, "Paired")[11:1])) +
  scale_y_continuous(expand = c(0.01, 0.01)) + 
  facet_wrap(~ Level, labeller = as_labeller(facet_names)) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        strip.background = element_rect(size = 0.2),
        axis.line.y = element_blank())
taxa_summary_by_sample_type(detra_class, detra$map_loaded, 'Treatment', 0.01, 'KW')

detra_order <- summarize_taxonomy(detra, level = 4, report_higher_tax = F)
plot_ts_heatmap(detra_order, detra$map_loaded, 0.01, 'Treatment', rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
detra_barsO <- plot_taxa_bars(detra_order, detra$map_loaded, "Treatment3", 
                              num_taxa = 12, data_only = T) %>%
  mutate(taxon = gsub("NA", "Unclassified", taxon)) %>%
  mutate(taxon = fct_relevel(taxon, "Other", after = Inf)) %>%
  mutate(taxon = fct_relevel(taxon, "Unclassified", after = Inf)) %>%
  mutate(taxon = fct_rev(taxon)) %>%
  separate(group_by, into = c("group_by", "Level"), sep = "_") %>%
  mutate(group_by = factor(group_by,
                           levels = c("TFM1", "TFM2", "TFM1@TFM2", "TFM1@OligoHal",
                                      "OligoHal", "TFM1@MesoHal", "MesoHal")))
ggplot(detra_barsO, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = NULL, y = "Relative abundance", fill = "Order") +
  scale_fill_manual(values = c("grey75", "grey90", brewer.pal(12, "Paired")[11:1])) +
  scale_y_continuous(expand = c(0.01, 0.01)) + 
  facet_wrap(~ Level, labeller = as_labeller(facet_names)) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        strip.background = element_rect(size = 0.2),
        axis.line.y = element_blank())
taxa_summary_by_sample_type(detra_order, detra$map_loaded, 'Treatment', 0.01, 'KW')

detra_family <- summarize_taxonomy(detra, level = 5, report_higher_tax = F)
plot_ts_heatmap(detra_family, detra$map_loaded, 0.01, 'Treatment', rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
detra_barsF <- plot_taxa_bars(detra_family, detra$map_loaded, "Treatment3", 
                              num_taxa = 12, data_only = T) %>%
  mutate(taxon = gsub("NA", "Unclassified", taxon)) %>%
  mutate(taxon = fct_relevel(taxon, "Other", after = Inf)) %>%
  mutate(taxon = fct_relevel(taxon, "Unclassified", after = Inf)) %>%
  mutate(taxon = fct_rev(taxon)) %>%
  separate(group_by, into = c("group_by", "Level"), sep = "_") %>%
  mutate(group_by = factor(group_by,
                           levels = c("TFM1", "TFM2", "TFM1@TFM2", "TFM1@OligoHal",
                                      "OligoHal", "TFM1@MesoHal", "MesoHal")))
ggplot(detra_barsF, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = NULL, y = "Relative abundance", fill = "Family") +
  scale_fill_manual(values = c("grey75", "grey90", brewer.pal(12, "Paired")[11:1])) +
  scale_y_continuous(expand = c(0.01, 0.01)) + 
  facet_wrap(~ Level, labeller = as_labeller(facet_names)) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        strip.background = element_rect(size = 0.2),
        axis.line.y = element_blank())
taxa_summary_by_sample_type(detra_family, detra$map_loaded, 'Treatment', 0.01, 'KW')

detra_genus <- summarize_taxonomy(detra, level = 6, report_higher_tax = F)
plot_ts_heatmap(detra_genus, detra$map_loaded, 0.01, 'Treatment', rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
detra_barsG <- plot_taxa_bars(detra_genus, detra$map_loaded, "Treatment3", 
                              num_taxa = 12, data_only = T) %>%
  mutate(taxon = gsub("NA", "Unclassified", taxon)) %>%
  mutate(taxon = fct_relevel(taxon, "Other", after = Inf)) %>%
  mutate(taxon = fct_relevel(taxon, "Unclassified", after = Inf)) %>%
  mutate(taxon = fct_rev(taxon)) %>%
  separate(group_by, into = c("group_by", "Level"), sep = "_") %>%
  mutate(group_by = factor(group_by,
                           levels = c("TFM1", "TFM2", "TFM1@TFM2", "TFM1@OligoHal",
                                      "OligoHal", "TFM1@MesoHal", "MesoHal")))
ggplot(detra_barsG, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = NULL, y = "Relative abundance", fill = "Genus") +
  scale_fill_manual(values = c("grey75", "grey90", brewer.pal(12, "Paired")[11:1])) +
  scale_y_continuous(expand = c(0.01, 0.01)) + 
  facet_wrap(~ Level, labeller = as_labeller(facet_names)) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        strip.background = element_rect(size = 0.2),
        axis.line.y = element_blank())
taxa_summary_by_sample_type(detra_genus, detra$map_loaded, 'Treatment', 0.01, 'KW')

detra_guilds <- summarize_taxonomy(detra, level = 9, report_higher_tax = F)
plot_ts_heatmap(detra_guilds, detra$map_loaded, 0.01, 'Treatment', rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
detra_barsGu <- plot_taxa_bars(detra_guilds, detra$map_loaded, "Treatment3", num_taxa = 20, 
                               data_only = T) %>%
  filter(taxon != "NA") %>%
  droplevels() %>%
  mutate(taxon = factor(taxon,
                        levels = Guild_cols$Guild)) %>%
  separate(group_by, into = c("group_by", "Level"), sep = "_") %>%
  mutate(group_by = factor(group_by,
                           levels = c("TFM1", "TFM2", "TFM1@TFM2", "TFM1@OligoHal",
                                      "OligoHal", "TFM1@MesoHal", "MesoHal")))
tallest_bar <- detra_barsGu %>%
  group_by(group_by) %>%
  summarise(sum = sum(mean_value))
pdf("InitialFigs/DEtra_Guilds.pdf", width = 7, height = 6)
ggplot(detra_barsGu, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Treatment", y = "Relative abundance", fill = "Guild") +
  scale_fill_manual(values = Guild_cols$color) +
  scale_y_continuous(expand = c(max(tallest_bar$sum)/100, max(tallest_bar$sum)/100)) + 
  facet_wrap(~ Level, labeller = as_labeller(facet_names)) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        strip.background = element_rect(size = 0.2),
        axis.line.y = element_blank())
dev.off()
taxa_summary_by_sample_type(detra_guilds, 
                            detra$map_loaded, 
                            type_header = 'Treatment', 
                            filter_level = 0.01, 
                            test_type = 'KW')

# Guilds by sample (to show variation)
detra$map_loaded$sampleID <- rownames(detra$map_loaded)
detra_barsGu <- plot_taxa_bars(detra_guilds,
                               detra$map_loaded,
                               "sampleID",
                               num_taxa = 20,
                               data_only = TRUE) %>%
  filter(taxon != "NA") %>%
  droplevels() %>%
  mutate(taxon = factor(taxon,
                        levels = Guild_cols$Guild)) %>%
  left_join(., detra$map_loaded, by = c("group_by" = "sampleID"))
facet_names <- c("0cm" = "Elevation 0 cm", 
                 "40cm" = "Elevation -40 cm", 
                 "0.02" = "Depth 2 cm",
                 "0.12" = "Depth 12 cm",
                 "TFM1" = "TFM1",
                 "TFM2" = "TFM2", 
                 "TFM1@TFM2" = "TFM1@TFM2",
                 "TFM1@OligoHal" = "TFM1@Oligo",
                 "OligoHal" = "Oligo", 
                 "TFM1@MesoHal" = "TFM1@Meso", 
                 "MesoHal" = "Meso")
tallest_bar <- detra_barsGu %>%
  group_by(group_by) %>%
  summarise(sum = sum(mean_value))
pdf("InitialFigs/DEtra_Guilds_samples.pdf", width = 12, height = 6)
ggplot(detra_barsGu, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Sample", y = "Relative abundance", fill = "Guild") +
  scale_fill_manual(values = Guild_cols$color) +
  scale_y_continuous(expand = c(max(tallest_bar$sum)/100, max(tallest_bar$sum)/100)) + 
  facet_nested(~ Level + Depth + Treatment2, space = "free", scales = "free_x",
               labeller = as_labeller(facet_names)) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 8, angle = 90, hjust = 1, vjust = 0.5),
        strip.text = element_text(size = 3),
        strip.background = element_rect(size = 0.2),
        axis.line.y = element_blank())
dev.off()

#### __Simper ####
detra_sim <- simper(t(detra$data_loaded), 
                    detra$map_loaded$Treatment)
detra_s <- summary(detra_sim)
detra_df1 <- head(detra_s$`TFM1@MesoHal_TFM1`, n = 20) %>%
  mutate(ASV = rownames(.),
         Comparison = "Freshwater/Mesohaline")
detra_df2 <- head(detra_s$`TFM1@MesoHal-40cm_TFM1-40cm`, n = 20) %>%
  mutate(ASV = rownames(.),
         Comparison = "Freshwater/Mesohaline 40 cm")
detra_simper_results <- rbind(detra_df1, detra_df2) %>%
  left_join(., detra$taxonomy_loaded, by = c("ASV" = "taxonomy8")) %>%
  mutate(TreatmentResponse = ifelse(avb > ava, "Postive", "Negative")) %>%
  rename("OTU" = "ASV") %>%
  mutate(OTU = gsub("ASV", "OTU", OTU)) %>%
  rename(c("Domain" = "taxonomy1",
           "Phylum" = "taxonomy2",
           "Class" = "taxonomy3",
           "Order" = "taxonomy4",
           "Family" = "taxonomy5",
           "Genus" = "taxonomy6",
           "Species" = "taxonomy7",
           "MeanTreatment" = "ava",
           "MeanControl" = "avb",
           "CumulativeContribution" = "cumsum")) %>%
  dplyr::select(Comparison, TreatmentResponse, Domain, Phylum, Class, Order, Family, Genus,
                Species, OTU, MeanTreatment, MeanControl, CumulativeContribution)
write_xlsx(detra_simper_results, 
           "simper_results_DEtra.xlsx",
           format_headers = F)

#### __Multipatt ####
set.seed(1202)
detra_mp <- multipatt(t(detra$data_loaded), 
                      detra$map_loaded$Treatment, 
                      func = "r.g", 
                      control = how(nperm=999),
                      max.order = 2)
# None with Q, use P
detra_mp_results <- detra_mp$sign %>%
  mutate(num_sites = rowSums(cbind(`s.TFM1@MesoHal`, `s.TFM1@MesoHal-40cm`)),
         q.value = qvalue(detra_mp$sign$p.value)$qvalues,
         Group = "NA",
         ASV = rownames(detra_mp$sign)) %>%
  filter(p.value == 0.001) %>%
  filter(stat >= 0.4) %>%
  # filter(num_sites <= 2) %>%
  left_join(., detra$taxonomy_loaded, by = c("ASV" = "taxonomy8"))
for (i in 1:nrow(detra_mp_results)) {
  if (detra_mp_results$num_sites[i] == 1 & detra_mp_results$`s.TFM1@MesoHal`[i] == 1) {
    detra_mp_results$Group[i] <- "TFM1@MesoHal"
  }
}
for (i in 1:nrow(detra_mp_results)) {
  if (detra_mp_results$num_sites[i] == 1 & detra_mp_results$`s.TFM1@MesoHal-40cm`[i] == 1) {
    detra_mp_results$Group[i] <- "TFM1@MesoHal-40cm"
  }
}
table(detra_mp_results$Group)
detra_asv <- summarize_taxonomy(detra, level = 8, report_higher_tax = F)
detra_asv_all <- data.frame("RelAbundance" = round(rowMeans(detra_asv) * 100, digits = 5)) %>%
  mutate(ASV = rownames(.))
detra_mp_corrs <- as.data.frame(detra_mp$str) %>%
  dplyr::select(`TFM1@MesoHal`, `TFM1@MesoHal-40cm`) %>%
  mutate(ASV = rownames(.)) %>%
  filter(ASV %in% detra_mp_results$ASV) %>%
  set_names(c("TFM1@MesoHal", "TFM1@MesoHal-40cm", "ASV"))
# Add corrs and taxonomy
detra_mp_results <- detra_mp_results %>%
  filter(Group == "TFM1@MesoHal" | Group == "TFM1@MesoHal-40cm") %>%
  left_join(., detra_asv_all, by = "ASV") %>%
  left_join(., detra_mp_corrs, by = "ASV") %>%
  dplyr::select(taxonomy2, taxonomy6, ASV, Group, names(detra_mp_corrs)[1:2], "RelAbundance") %>%
  mutate(ASV = gsub("ASV", "OTU", ASV)) %>%
  arrange(Group, taxonomy2, taxonomy6) %>%
  mutate(ASV = factor(ASV, levels = ASV)) %>%
  set_names(c("Phylum", "Genus", "ASV", "Group", names(detra_mp_corrs)[1:2], "RelAbundance")) %>%
  mutate(taxon = paste(Phylum, Genus, ASV, sep = "; ")) %>%
  mutate(taxon = factor(taxon, levels = taxon))
# Plot heatmap with abundance barplot
detra.hm.melted <- detra_mp_results %>%
  dplyr::select(taxon, names(detra_mp_corrs)[1:2]) %>%
  melt(., id.vars = c("taxon"))
detra.hm <- ggplot(data = detra.hm.melted, 
                   aes(x = factor(taxon), y = variable, fill = value)) + 
  geom_tile() + 
  scale_fill_distiller(name = "Indicator correlation index", palette = "RdBu", 
                       direction = -1, na.value = "transparent", type = "div", limits = c(-1, 1)) +
  scale_x_discrete(breaks = unique(detra.hm.melted$taxon), labels = unique(detra.hm.melted$taxon),
                   limits = rev(levels(detra.hm.melted$taxon))) +
  scale_y_discrete(position = "right") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        legend.title = element_text(size = 10), legend.key.height = unit(0.5,"cm"),
        legend.key.width = unit(1, "cm"), legend.text = element_text(size = 7)) +
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5)) +
  coord_flip()
detra.l <- get_legend(detra.hm)
detra.hm.clean <- detra.hm +
  theme(axis.title.y = element_blank(), axis.text.y = element_text(margin = margin(c(0,-5,0,0))),
        axis.ticks.y = element_blank(), axis.ticks.x = element_blank(), 
        panel.grid.major = element_blank(), legend.position="none",
        axis.text.x.top = element_text(size = 6, margin = margin(c(0,0,-2,0))), 
        plot.margin = margin(c(0,-2,0,5))) +
  labs(y = "Group")
detra.bp.y <- ggplot(data = detra_mp_results, aes(x = taxon, y = RelAbundance)) + 
  geom_bar(stat = "identity", fill = "grey") + 
  scale_y_continuous(position = "right") +
  scale_x_discrete(limits = rev(levels(detra_mp_results$taxon))) +
  coord_flip() + 
  theme_minimal() +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position="none", plot.margin = margin(c(0,5,0,-5)),
        axis.text.x.top = element_text(size = 6, margin = margin(c(0,0,-2,0)))) +
  labs(y = "Rel. abund. (%)")
pdf("InitialFigs/DEtra_Multipatt.pdf", width = 8, height = 5)
plot_grid(detra.hm.clean, detra.bp.y, detra.l, NULL, nrow = 2, ncol = 2, 
          rel_widths = c(10,2), rel_heights = c(12, 2), align = "hv", axis = "b")
dev.off()

# Problem is that at OTU level, OTUs are too rare, not even 1%. 
# So let's see what indicators there are at the genus level.
# First get aggregate taxonomy table
detra_tax <- detra$taxonomy_loaded %>%
  group_by(taxonomy6) %>%
  slice_head(n = 1) %>%
  dplyr::select(-taxonomy7, -taxonomy8)
detra_genus <- summarize_taxonomy(detra, level = 6, relative = F, report_higher_tax = F)
set.seed(1202)
detra_mp <- multipatt(t(detra_genus), 
                      detra$map_loaded$Treatment, 
                      func = "r.g", 
                      control = how(nperm=999),
                      max.order = 2)
detra_mp_results <- detra_mp$sign %>%
  mutate(num_sites = rowSums(cbind(`s.TFM1@MesoHal`, `s.TFM1@MesoHal-40cm`)),
         q.value = qvalue(detra_mp$sign$p.value)$qvalues,
         Group = "NA",
         Genus = rownames(detra_mp$sign)) %>%
  filter(p.value < 0.01) %>%
  filter(stat >= 0.5) %>%
  filter(num_sites <= 2) %>%
  left_join(., detra_tax, by = c("Genus" = "taxonomy6"))
for (i in 1:nrow(detra_mp_results)) {
  if (detra_mp_results$num_sites[i] == 1 & detra_mp_results$`s.TFM1@MesoHal`[i] == 1) {
    detra_mp_results$Group[i] <- "TFM1@MesoHal"
  }
}
for (i in 1:nrow(detra_mp_results)) {
  if (detra_mp_results$num_sites[i] == 1 & detra_mp_results$`s.TFM1@MesoHal-40cm`[i] == 1) {
    detra_mp_results$Group[i] <- "TFM1@MesoHal-40cm"
  }
}
table(detra_mp_results$Group)
detra_genus_all <- data.frame("RelAbundance" = round(rowMeans(detra_genus)/min(colSums(detra$data_loaded)) * 100, 
                                                     digits = 4)) %>%
  mutate(Genus = rownames(.))
detra_mp_corrs <- as.data.frame(detra_mp$str) %>%
  dplyr::select(`TFM1@MesoHal`, `TFM1@MesoHal-40cm`) %>%
  mutate(Genus = rownames(.)) %>%
  filter(Genus %in% detra_mp_results$Genus) %>%
  set_names(c("TFM1@MesoHal", "TFM1@MesoHal-40cm", "Genus"))
# Add corrs and taxonomy
detra_mp_results <- detra_mp_results %>%
  filter(Group == "TFM1@MesoHal" | Group == "TFM1@MesoHal-40cm") %>%
  left_join(., detra_genus_all, by = "Genus") %>%
  left_join(., detra_mp_corrs, by = "Genus") %>%
  dplyr::select(taxonomy2, Genus, Group, names(detra_mp_corrs)[1:2], "RelAbundance") %>%
  arrange(Group, taxonomy2) %>%
  mutate(Genus = factor(Genus, levels = Genus)) %>%
  set_names(c("Phylum", "Genus", "Group", names(detra_mp_corrs)[1:2], "RelAbundance")) %>%
  mutate(taxon = paste(Phylum, Genus, sep = "; ")) %>%
  mutate(taxon = factor(taxon, levels = taxon))
# Plot heatmap with abundance barplot
detra.hm.melted <- detra_mp_results %>%
  dplyr::select(taxon, names(detra_mp_corrs)[1:2]) %>%
  melt(., id.vars = c("taxon"))
detra.hm <- ggplot(data = detra.hm.melted, 
                   aes(x = factor(taxon), y = variable, fill = value)) + 
  geom_tile() + 
  scale_fill_distiller(name = "Indicator correlation index", palette = "RdBu", 
                       direction = -1, na.value = "transparent", type = "div", limits = c(-1, 1)) +
  scale_x_discrete(breaks = unique(detra.hm.melted$taxon), labels = unique(detra.hm.melted$taxon),
                   limits = rev(levels(detra.hm.melted$taxon))) +
  scale_y_discrete(position = "right") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        legend.title = element_text(size = 10), legend.key.height = unit(0.5,"cm"),
        legend.key.width = unit(1, "cm"), legend.text = element_text(size = 7)) +
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5)) +
  coord_flip()
detra.l <- get_legend(detra.hm)
detra.hm.clean <- detra.hm +
  theme(axis.title.y = element_blank(), axis.text.y = element_text(margin = margin(c(0,-5,0,0))),
        axis.ticks.y = element_blank(), axis.ticks.x = element_blank(), 
        panel.grid.major = element_blank(), legend.position="none",
        axis.text.x.top = element_text(size = 6, margin = margin(c(0,0,-2,0))), 
        plot.margin = margin(c(0,-2,0,5))) +
  labs(y = "Group")
detra.bp.y <- ggplot(data = detra_mp_results, aes(x = taxon, y = RelAbundance)) + 
  geom_bar(stat = "identity", fill = "grey") + 
  scale_y_continuous(position = "right") +
  scale_x_discrete(limits = rev(levels(detra_mp_results$taxon))) +
  coord_flip() + 
  theme_minimal() +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position="none", plot.margin = margin(c(0,5,0,-5)),
        axis.text.x.top = element_text(size = 6, margin = margin(c(0,0,-2,0)))) +
  labs(y = "Rel. abund. (%)")
pdf("InitialFigs/DEtra_Multipatt_genus.pdf", width = 8, height = 5)
plot_grid(detra.hm.clean, detra.bp.y, detra.l, NULL, nrow = 2, ncol = 2, 
          rel_widths = c(10,2), rel_heights = c(12, 2), align = "hv", axis = "b")
dev.off()



#### __BGC ####
# Add biogeochem. analysis
# Get variables, make corrplot, envfit, ordination with vectors, ordistep, taxa correlations

# Get variables
env <- detra$map_loaded %>%
  dplyr::select(CH4_ug_m2_h, N2O_ug_m2_h, 
                Salinity_ppt_all, 
                NH4_mgL, PO4_mgL, Cl_mgL, SO4_mgL,
                Fe_mgL, Porosity, Acetate_mgL, TotalVFA_uM, SR_umol_cm3_d, AMG_umol_cm3_d)
env_nona <- na.omit(env)

# Corrplot
C <- cor(env_nona)
corrplot(C, method = "number", type = "lower", 
         tl.cex = 0.5, tl.col = "black", number.cex = 0.5)
pdf("InitialFigs/DEtra_BGC_CH4.pdf", width = 7, height = 5)
meth_corr_by_bgc(env_nona = env_nona)
dev.off()

# Envfit
pcoa <- cmdscale(detra_bc, k = nrow(detra$map_loaded) - 1, eig = T)
set.seed(100)
ef <- envfit(pcoa, env, permutations = 999, na.rm = TRUE)
ef
arrow_factor <- ordiArrowMul(ef)
manual_factor <- 0.5
vec.df <- as.data.frame(ef$vectors$arrows*sqrt(ef$vectors$r)) %>%
  mutate(Dim1 = Dim1 * manual_factor,
         Dim2 = Dim2 * manual_factor) %>%
  mutate(variables = rownames(.)) %>%
  filter(ef$vectors$pvals < 0.05) %>%
  filter(variables != "Cl_mgL") %>%
  mutate(shortnames = c("Salinity", "SO4", "Porosity", "SR"))

# Plot with significant vectors
ordiplot(pcoa)
plot(ef, p.max = 0.05, cex = 0.5)
pcoaA1 <- round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[1]*100, digits = 1)
pcoaA2 <- round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[2]*100, digits = 1)
detra$map_loaded$Axis01 <- scores(pcoa)[,1]
detra$map_loaded$Axis02 <- scores(pcoa)[,2]
micro.hulls <- ddply(detra$map_loaded, c("Treatment"), find_hull)
g <- ggplot(detra$map_loaded, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = Treatment, fill = Treatment),
               alpha = 0.1, show.legend = F) +
  geom_point(size = 3, alpha = 1, aes(colour = Treatment, shape = Depth)) +
  geom_segment(data = vec.df,
               aes(x = 0, xend = Dim1, y = 0, yend = Dim2),
               arrow = arrow(length = unit(0.5, "cm")),
               colour = "gray", alpha = 0.5,
               inherit.aes = FALSE) + 
  geom_text(data = vec.df,
            aes(x = Dim1, y = Dim2, label = shortnames),
            size = 3.5, color = "black") +
  labs(x = paste("PC1: ", pcoaA1, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2, "%", sep = ""),
       shape = "Depth (m)") +
  scale_colour_viridis_d() +
  scale_fill_viridis_d() +
  guides(colour = guide_legend(override.aes = list(shape = 15),
                               order = 1)) +
  coord_fixed() +
  theme_bw() +  
  theme(legend.position = "right",
        legend.background = element_blank(),
        legend.key.size = unit(0.4, "cm"),
        legend.margin = margin(2,2,2,2, "pt"),
        legend.spacing.y = unit(0, "cm"),
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        plot.margin = margin(5, 5, 5, 5, "pt"))
g
pdf("InitialFigs/DEtra_PCoA_wBGC.pdf", width = 6, height = 4)
g
dev.off()
ggplotly(g)

# Ordistep
comm_nona <- as.data.frame(t(detra$data_loaded)) %>%
  filter(rownames(.) %in% rownames(env_nona))
mod0 <- rda(comm_nona ~ 1, env_nona)  # Model with intercept only
mod1 <- rda(comm_nona ~ ., env_nona)  # Model with all explanatory variables
mod <- ordistep(mod0, scope = formula(mod1))
mod
mod$anova # Salinity

# Plot correlations for different taxonomic levels at a given % relative abundance threshold
pdf("InitialFigs/DEtra_Phyla_CH4.pdf", width = 7, height = 5)
meth_corr_by_taxonomy(input = detra, level = 2, threshold = 0.5)
dev.off()
meth_corr_by_taxonomy(input = detra, level = 3, threshold = 0.5)
meth_corr_by_taxonomy(input = detra, level = 4, threshold = 0.5)
meth_corr_by_taxonomy(input = detra, level = 5, threshold = 0.5)
meth_corr_by_taxonomy(input = detra, level = 6, threshold = 0.5)
meth_corr_by_taxonomy(input = detra, level = 8, threshold = 0.5)
pdf("InitialFigs/DEtra_Guilds_CH4.pdf", width = 7, height = 5)
meth_corr_by_taxonomy(input = detra, level = 9, threshold = 0)
dev.off()



#### ...................................... ####
#### 6. SF ####
# input_filt <- readRDS("input_filt.rds")
input_filt <- readRDS("input_filt_comb_wBGC.rds")
sf <- filter_data(input_filt,
                  filter_cat = "Estuary",
                  keep_vals = "SF")
sort(colSums(sf$data_loaded))
mean(colSums(sf$data_loaded))
se(colSums(sf$data_loaded))
# Rarefy at 65222
# Drop Sandmound_TuleB_D1 (1296 reads) and Muzzi_PWB_D2 (5287 reads)
set.seed(530)
sf <- single_rarefy(sf, 65222)
sf$map_loaded <- sf$map_loaded %>%
  mutate(rich = specnumber(sf$data_loaded, MARGIN = 2),
         shannon = diversity(sf$data_loaded, index = "shannon", MARGIN = 2),
         Salt = NA,
         Depth = as.factor(Depth)) %>%
  mutate_if(is.character, as.factor)

# Add column for Salt
for (i in 1:nrow(sf$map_loaded)) {
  if (sf$map_loaded$Site[i] == "Sandmound" | sf$map_loaded$Site[i] == "West Pond") {
    sf$map_loaded$Salt[i] <- "Freshwater"
  }
}
for (i in 1:nrow(sf$map_loaded)) {
  if (sf$map_loaded$Site[i] == "Mayberry" | sf$map_loaded$Site[i] == "Browns") {
    sf$map_loaded$Salt[i] <- "Oligohaline"
  }
}
for (i in 1:nrow(sf$map_loaded)) {
  if (sf$map_loaded$Site[i] == "Joice" | sf$map_loaded$Site[i] == "Rush Ranch") {
    sf$map_loaded$Salt[i] <- "Mesohaline"
  }
}
for (i in 1:nrow(sf$map_loaded)) {
  if (sf$map_loaded$Site[i] == "Goodyear" | sf$map_loaded$Site[i] == "White Slough" |
      sf$map_loaded$Site[i] == "Tolay" | sf$map_loaded$Site[i] == "China Camp" |
      sf$map_loaded$Site[i] == "Muzzi") {
    sf$map_loaded$Salt[i] <- "Polyhaline"
  }
}
sf$map_loaded <- sf$map_loaded %>%
  mutate(Salt = factor(Salt,
                       levels = c("Freshwater", "Oligohaline", 
                                  "Mesohaline", "Polyhaline")))

#### __Alpha ####
leveneTest(sf$map_loaded$rich ~ sf$map_loaded$Salt) # Homogeneous
m <- aov(rich ~ Salt + Depth, data = sf$map_loaded)
Anova(m, type = "II") # Salt, not Depth
m <- aov(rich ~ Salt, data = sf$map_loaded)
shapiro.test(m$residuals) # Normal
summary(m)
t <- emmeans(object = m, specs = "Salt") %>%
  cld(object = ., adjust = "Tukey", Letters = letters, alpha = 0.05) %>%
  mutate(name = "rich",
         y = max(sf$map_loaded$rich)+(max(sf$map_loaded$rich)-min(sf$map_loaded$rich))/20)
leveneTest(sf$map_loaded$shannon ~ sf$map_loaded$Salt) # Homogeneous
m1 <- aov(shannon ~ Salt + Depth, data = sf$map_loaded)
Anova(m1, type = "II") # Salt, not Depth
m1 <- aov(shannon ~ Salt, data = sf$map_loaded)
shapiro.test(m1$residuals) # Normal
summary(m1)
t1 <- emmeans(object = m1, specs = "Salt") %>%
  cld(object = ., adjust = "Tukey", Letters = letters, alpha = 0.05) %>%
  mutate(name = "shannon",
         y = max(sf$map_loaded$shannon)+(max(sf$map_loaded$shannon)-min(sf$map_loaded$shannon))/20)
label_df <- rbind(t, t1)
facet_df <- c("rich" = "(a) Richness",
              "shannon" = "(b) Shannon")
alpha_long <- sf$map_loaded %>%
  pivot_longer(cols = c("rich", "shannon"))
pdf("InitialFigs/SF_Alpha.pdf", width = 6, height = 3)
ggplot(alpha_long, aes(reorder(Salt, value, mean), value, 
                       colour = Salt)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 2, alpha = 0.75, width = 0.2, aes(shape = Depth)) +
  geom_text(data = label_df, aes(Salt, y, label = str_trim(.group)), 
            size = 4, color = "black") +
  labs(x = "Site", y = "Number of OTUs", shape = "Depth (m)") +
  scale_colour_viridis_d() +
  guides(colour = "none") +
  facet_wrap(~ name, ncol = 2, scales = "free_y", labeller = as_labeller(facet_df)) +
  theme_bw() +
  theme(legend.position = c(1,0),
        legend.justification = c(1,0),
        legend.background = element_blank(),
        legend.title.align = 0.5,
        axis.title = element_blank(),
        axis.text = element_text(size = 10),
        strip.text = element_text(size = 10))
dev.off()

#### __Beta ####
sf_bc <- calc_dm(sf$data_loaded)
set.seed(1150)
adonis2(sf_bc ~ sf$map_loaded$Salt + sf$map_loaded$Depth) # Salt and Depth sig
anova(betadisper(sf_bc, sf$map_loaded$Salt)) # Dispersion not homogeneous
anova(betadisper(sf_bc, sf$map_loaded$Depth)) # Dispersion homogeneous
sf_pcoa <- cmdscale(sf_bc, k = nrow(sf$map_loaded) - 1, eig = T)
pcoaA1 <- round((eigenvals(sf_pcoa)/sum(eigenvals(sf_pcoa)))[1]*100, digits = 1)
pcoaA2 <- round((eigenvals(sf_pcoa)/sum(eigenvals(sf_pcoa)))[2]*100, digits = 1)
sf$map_loaded$Axis01 <- scores(sf_pcoa)[,1]
sf$map_loaded$Axis02 <- scores(sf_pcoa)[,2]
micro.hulls <- ddply(sf$map_loaded, c("Salt"), find_hull)
pdf("InitialFigs/SF_PCoA.pdf", width = 7, height = 5)
ggplot(sf$map_loaded, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = Salt, fill = Salt),
               alpha = 0.1, show.legend = F) +
  geom_point(size = 3, alpha = 1, aes(colour = Salt, shape = Depth),
             show.legend = T) +
  labs(x = paste("PC1: ", pcoaA1, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2, "%", sep = ""),
       shape = "Depth (cm)") +
  scale_fill_viridis_d() +
  scale_colour_viridis_d() +
  theme_bw() +  
  theme(legend.position = "right",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10))
dev.off()

#### __Taxa ####
sf_phyla <- summarize_taxonomy(sf, level = 2, report_higher_tax = F)
plot_ts_heatmap(sf_phyla, sf$map_loaded, 0.01, 'Salt', rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
sf_barsP <- plot_taxa_bars(sf_phyla, sf$map_loaded, "Salt", num_taxa = 12, data_only = T) %>%
  mutate(taxon = fct_rev(taxon))
pdf("InitialFigs/SF_Phyla.pdf", width = 7, height = 5)
ggplot(sf_barsP, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Sample", y = "Relative abundance", fill = "Phylum") +
  scale_fill_manual(values = c("grey90", brewer.pal(12, "Paired")[12:1])) +
  scale_y_continuous(expand = c(0.01, 0.01)) + 
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 10))
dev.off()
taxa_summary_by_sample_type(sf_phyla, sf$map_loaded, 'Salt', 0.01, 'KW')

sf_class <- summarize_taxonomy(sf, level = 3, report_higher_tax = F)
plot_ts_heatmap(sf_class, sf$map_loaded, 0.01, 'Salt', rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
sf_barsC <- plot_taxa_bars(sf_class, sf$map_loaded, "Salt", 
                              num_taxa = 12, data_only = T) %>%
  mutate(taxon = gsub("NA", "Unclassified", taxon)) %>%
  mutate(taxon = fct_relevel(taxon, "Other", after = Inf)) %>%
  mutate(taxon = fct_relevel(taxon, "Unclassified", after = Inf)) %>%
  mutate(taxon = fct_rev(taxon))
ggplot(sf_barsC, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = NULL, y = "Relative abundance", fill = "Class") +
  scale_fill_manual(values = c("grey75", "grey90", brewer.pal(12, "Paired")[11:1])) +
  scale_y_continuous(expand = c(0.01, 0.01)) + 
  theme_classic() +
  theme(axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10))
taxa_summary_by_sample_type(sf_class, sf$map_loaded, 'Salt', 0.01, 'KW')

sf_order <- summarize_taxonomy(sf, level = 4, report_higher_tax = F)
plot_ts_heatmap(sf_order, sf$map_loaded, 0.01, 'Salt', rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
sf_barsO <- plot_taxa_bars(sf_order, sf$map_loaded, "Salt", 
                              num_taxa = 12, data_only = T) %>%
  mutate(taxon = gsub("NA", "Unclassified", taxon)) %>%
  mutate(taxon = fct_relevel(taxon, "Other", after = Inf)) %>%
  mutate(taxon = fct_relevel(taxon, "Unclassified", after = Inf)) %>%
  mutate(taxon = fct_rev(taxon))
ggplot(sf_barsO, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = NULL, y = "Relative abundance", fill = "Order") +
  scale_fill_manual(values = c("grey75", "grey90", brewer.pal(12, "Paired")[11:1])) +
  scale_y_continuous(expand = c(0.01, 0.01)) + 
  theme_classic() +
  theme(axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10))
taxa_summary_by_sample_type(sf_order, sf$map_loaded, 'Salt', 0.01, 'KW')

sf_family <- summarize_taxonomy(sf, level = 5, report_higher_tax = F)
plot_ts_heatmap(sf_family, sf$map_loaded, 0.01, 'Salt', rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
sf_barsF <- plot_taxa_bars(sf_family, sf$map_loaded, "Salt", 
                              num_taxa = 12, data_only = T) %>%
  mutate(taxon = gsub("NA", "Unclassified", taxon)) %>%
  mutate(taxon = fct_relevel(taxon, "Other", after = Inf)) %>%
  mutate(taxon = fct_relevel(taxon, "Unclassified", after = Inf)) %>%
  mutate(taxon = fct_rev(taxon))
ggplot(sf_barsF, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = NULL, y = "Relative abundance", fill = "Family") +
  scale_fill_manual(values = c("grey75", "grey90", brewer.pal(12, "Paired")[11:1])) +
  scale_y_continuous(expand = c(0.01, 0.01)) + 
  theme_classic() +
  theme(axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10))
taxa_summary_by_sample_type(sf_family, sf$map_loaded, 'Salt', 0.01, 'KW')

sf_genus <- summarize_taxonomy(sf, level = 6, report_higher_tax = F)
plot_ts_heatmap(sf_genus, sf$map_loaded, 0.01, 'Salt', rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
sf_barsG <- plot_taxa_bars(sf_genus, sf$map_loaded, "Salt", 
                              num_taxa = 12, data_only = T) %>%
  mutate(taxon = gsub("NA", "Unclassified", taxon)) %>%
  mutate(taxon = fct_relevel(taxon, "Other", after = Inf)) %>%
  mutate(taxon = fct_relevel(taxon, "Unclassified", after = Inf)) %>%
  mutate(taxon = fct_rev(taxon))
ggplot(sf_barsG, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = NULL, y = "Relative abundance", fill = "Genus") +
  scale_fill_manual(values = c("grey75", "grey90", brewer.pal(12, "Paired")[11:1])) +
  scale_y_continuous(expand = c(0.01, 0.01)) + 
  theme_classic() +
  theme(axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10))
taxa_summary_by_sample_type(sf_genus, sf$map_loaded, 'Salt', 0.01, 'KW')

sf_guilds <- summarize_taxonomy(sf, level = 9, report_higher_tax = F)
plot_ts_heatmap(sf_guilds, sf$map_loaded, 0.01, 'Salt', rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
sf_barsGu <- plot_taxa_bars(sf_guilds, sf$map_loaded, "Salt", 
                               num_taxa = 20, data_only = T) %>%
  filter(taxon != "NA") %>%
  droplevels() %>%
  mutate(taxon = factor(taxon,
                        levels = Guild_cols$Guild))
tallest_bar <- sf_barsGu %>%
  group_by(group_by) %>%
  summarise(sum = sum(mean_value))
pdf("InitialFigs/SF_Guilds.pdf", width = 7, height = 5)
ggplot(sf_barsGu, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Salt", y = "Relative abundance", fill = "Guild") +
  scale_fill_manual(values = Guild_cols$color) +
  scale_y_continuous(expand = c(max(tallest_bar$sum)/100, max(tallest_bar$sum)/100)) + 
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        axis.title.x = element_blank())
dev.off()
taxa_summary_by_sample_type(sf_guilds, 
                            sf$map_loaded, 
                            type_header = 'Salt', 
                            filter_level = 0.01, 
                            test_type = 'KW')

# Guilds by sample (to show variation)
sf$map_loaded$sampleID <- rownames(sf$map_loaded)
sf_barsGu <- plot_taxa_bars(sf_guilds,
                               sf$map_loaded,
                               "sampleID",
                               num_taxa = 20,
                               data_only = TRUE) %>%
  filter(taxon != "NA") %>%
  droplevels() %>%
  mutate(taxon = factor(taxon,
                        levels = Guild_cols$Guild)) %>%
  left_join(., sf$map_loaded, by = c("group_by" = "sampleID"))
facet_names <- c(" 0-5" = "Depth 5 cm",
                 " 5-15" = "Depth 15 cm",
                 "Freshwater" = "Freshwater",
                 "Oligohaline" = "Oligohaline",
                 "Mesohaline" = "Mesohaline",
                 "Polyhaline" = "Polyhaline")
tallest_bar <- sf_barsGu %>%
  group_by(group_by) %>%
  summarise(sum = sum(mean_value))
pdf("InitialFigs/SF_Guilds_samples.pdf", width = 9, height = 5)
ggplot(sf_barsGu, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Sample", y = "Relative abundance", fill = "Guild") +
  scale_fill_manual(values = Guild_cols$color) +
  scale_y_continuous(expand = c(max(tallest_bar$sum)/100, max(tallest_bar$sum)/100)) + 
  facet_nested(~ Depth + Salt, space = "free", scales = "free_x",
               labeller = as_labeller(facet_names)) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 90, hjust = 1, vjust = 0.5),
        strip.background = element_rect(size = 0.2),
        axis.line.y = element_blank())
dev.off()

#### __Simper ####
sf_sim <- simper(t(sf$data_loaded), 
                    sf$map_loaded$Salt)
sf_s <- summary(sf_sim)
head(sf_s$Freshwater_Mesohaline)
sf_simper_results <- head(sf_s$Freshwater_Mesohaline, n = 20) %>%
  mutate(ASV = rownames(.),
         Comparison = "Freshwater/Mesohaline") %>%
  left_join(., sf$taxonomy_loaded, by = c("ASV" = "taxonomy8")) %>%
  mutate(SaltResponse = ifelse(avb > ava, "Postive", "Negative")) %>%
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
write_xlsx(sf_simper_results, 
           "simper_results_sf.xlsx",
           format_headers = F)

#### __Multipatt ####
set.seed(1202)
sf_mp <- multipatt(t(sf$data_loaded), 
                      sf$map_loaded$Salt, 
                      func = "r.g", 
                      control = how(nperm=999))
# None with Q, use P
sf_mp_results <- sf_mp$sign %>%
  mutate(num_sites = rowSums(cbind(s.Freshwater, s.Oligohaline, s.Mesohaline)),
         q.value = qvalue(sf_mp$sign$p.value)$qvalues,
         Group = "NA",
         ASV = rownames(sf_mp$sign)) %>%
  filter(p.value == 0.001) %>%
  filter(stat >= 0.9) %>%
  filter(num_sites <= 2) %>%
  left_join(., sf$taxonomy_loaded, by = c("ASV" = "taxonomy8"))
for (i in 1:nrow(sf_mp_results)) {
  if (sf_mp_results$num_sites[i] == 1 & sf_mp_results$s.Freshwater[i] == 1) {
    sf_mp_results$Group[i] <- "Freshwater"
  }
}
for (i in 1:nrow(sf_mp_results)) {
  if (sf_mp_results$num_sites[i] == 1 & sf_mp_results$s.Oligohaline[i] == 1) {
    sf_mp_results$Group[i] <- "Oligohaline"
  }
}
for (i in 1:nrow(sf_mp_results)) {
  if (sf_mp_results$num_sites[i] == 1 & sf_mp_results$s.Mesohaline[i] == 1) {
    sf_mp_results$Group[i] <- "Mesohaline"
  }
}
table(sf_mp_results$Group)
sf_asv <- summarize_taxonomy(sf, level = 8, report_higher_tax = F)
sf_asv_all <- data.frame("RelAbundance" = round(rowMeans(sf_asv) * 100, digits = 5)) %>%
  mutate(ASV = rownames(.))
sf_mp_corrs <- as.data.frame(sf_mp$str) %>%
  dplyr::select(1:length(levels(sf$map_loaded$Salt))) %>%
  mutate(ASV = rownames(.)) %>%
  filter(ASV %in% sf_mp_results$ASV) %>%
  set_names(c("Freshwater", "Oligohaline", "Mesolhaline", "ASV"))
# Add corrs and taxonomy
sf_mp_results <- sf_mp_results %>%
  filter(Group == "Mesohaline") %>%
  left_join(., sf_asv_all, by = "ASV") %>%
  left_join(., sf_mp_corrs, by = "ASV") %>%
  dplyr::select(taxonomy2, taxonomy6, ASV, Group, names(sf_mp_corrs)[1:3], "RelAbundance") %>%
  mutate(ASV = gsub("ASV", "OTU", ASV)) %>%
  arrange(Group, taxonomy2, taxonomy6) %>%
  mutate(ASV = factor(ASV, levels = ASV)) %>%
  set_names(c("Phylum", "Genus", "ASV", "Group", names(sf_mp_corrs)[1:3], "RelAbundance")) %>%
  mutate(taxon = paste(Phylum, Genus, ASV, sep = "; ")) %>%
  mutate(taxon = factor(taxon, levels = taxon))
# Plot heatmap with abundance barplot
sf.hm.melted <- sf_mp_results %>%
  dplyr::select(taxon, names(sf_mp_corrs)[1:3]) %>%
  melt(., id.vars = c("taxon"))
sf.hm <- ggplot(data = sf.hm.melted, 
                   aes(x = factor(taxon), y = variable, fill = value)) + 
  geom_tile() + 
  scale_fill_distiller(name = "Indicator correlation index", palette = "RdBu", 
                       direction = -1, na.value = "transparent", type = "div", limits = c(-1, 1)) +
  scale_x_discrete(breaks = unique(sf.hm.melted$taxon), labels = unique(sf.hm.melted$taxon),
                   limits = rev(levels(sf.hm.melted$taxon))) +
  scale_y_discrete(position = "right") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        legend.title = element_text(size = 10), legend.key.height = unit(0.5,"cm"),
        legend.key.width = unit(1, "cm"), legend.text = element_text(size = 7)) +
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5)) +
  coord_flip()
sf.l <- get_legend(sf.hm)
sf.hm.clean <- sf.hm +
  theme(axis.title.y = element_blank(), axis.text.y = element_text(margin = margin(c(0,-10,0,0)),
                                                                   size = 5),
        axis.ticks.y = element_blank(), axis.ticks.x = element_blank(), 
        panel.grid.major = element_blank(), legend.position="none",
        axis.text.x.top = element_text(size = 6, margin = margin(c(0,0,-2,0))), 
        plot.margin = margin(c(0,-2,0,5))) +
  labs(y = "Group")
sf.bp.y <- ggplot(data = sf_mp_results, aes(x = taxon, y = RelAbundance)) + 
  geom_bar(stat = "identity", fill = "grey") + 
  scale_y_continuous(position = "right") +
  scale_x_discrete(limits = rev(levels(sf_mp_results$taxon))) +
  coord_flip() + 
  theme_minimal() +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position="none", plot.margin = margin(c(0,5,0,-5)),
        axis.text.x.top = element_text(size = 6, margin = margin(c(0,0,-2,0)))) +
  labs(y = "Rel. abund. (%)")
pdf("InitialFigs/sf_Multipatt.pdf", width = 8, height = 5)
plot_grid(sf.hm.clean, sf.bp.y, sf.l, NULL, nrow = 2, ncol = 2, 
          rel_widths = c(10,2), rel_heights = c(12, 2), align = "hv", axis = "b")
dev.off()

# Problem is that at OTU level, many OTUs are too rare, not even 1%. 
# So let's see what indicators there are at the genus level.
# First get aggregate taxonomy table
sf_tax <- sf$taxonomy_loaded %>%
  group_by(taxonomy6) %>%
  slice_head(n = 1) %>%
  dplyr::select(-taxonomy7, -taxonomy8)
sf_genus <- summarize_taxonomy(sf, level = 6, relative = F, report_higher_tax = F)
set.seed(1202)
sf_mp <- multipatt(t(sf_genus), 
                      sf$map_loaded$Salt, 
                      func = "r.g", 
                      control = how(nperm=999))
sf_mp_results <- sf_mp$sign %>%
  mutate(num_sites = rowSums(cbind(s.Freshwater, s.Oligohaline, s.Mesohaline)),
         q.value = qvalue(sf_mp$sign$p.value)$qvalues,
         Group = "NA",
         Genus = rownames(sf_mp$sign)) %>%
  filter(q.value < 0.05) %>%
  filter(stat >= 0.7) %>%
  filter(num_sites <= 2) %>%
  left_join(., sf_tax, by = c("Genus" = "taxonomy6"))
for (i in 1:nrow(sf_mp_results)) {
  if (sf_mp_results$num_sites[i] == 1 & sf_mp_results$s.Freshwater[i] == 1) {
    sf_mp_results$Group[i] <- "Freshwater"
  }
}
for (i in 1:nrow(sf_mp_results)) {
  if (sf_mp_results$num_sites[i] == 1 & sf_mp_results$s.Oligohaline[i] == 1) {
    sf_mp_results$Group[i] <- "Oligohaline"
  }
}
for (i in 1:nrow(sf_mp_results)) {
  if (sf_mp_results$num_sites[i] == 1 & sf_mp_results$s.Mesohaline[i] == 1) {
    sf_mp_results$Group[i] <- "Mesohaline"
  }
}
table(sf_mp_results$Group)
sf_genus_all <- data.frame("RelAbundance" = round(rowMeans(sf_genus)/min(colSums(sf$data_loaded)) * 100, 
                                                     digits = 4)) %>%
  mutate(Genus = rownames(.))
sf_mp_corrs <- as.data.frame(sf_mp$str) %>%
  dplyr::select(1:length(levels(sf$map_loaded$Salt))) %>%
  mutate(Genus = rownames(.)) %>%
  filter(Genus %in% sf_mp_results$Genus) %>%
  set_names(c("Freshwater", "Oligohaline", "Mesohaline", "Genus"))
# Add corrs and taxonomy
sf_mp_results <- sf_mp_results %>%
  filter(Group == "Mesohaline") %>%
  left_join(., sf_genus_all, by = "Genus") %>%
  left_join(., sf_mp_corrs, by = "Genus") %>%
  dplyr::select(taxonomy2, Genus, Group, names(sf_mp_corrs)[1:3], "RelAbundance") %>%
  arrange(Group, taxonomy2) %>%
  mutate(Genus = factor(Genus, levels = Genus)) %>%
  set_names(c("Phylum", "Genus", "Group", names(sf_mp_corrs)[1:3], "RelAbundance")) %>%
  mutate(taxon = paste(Phylum, Genus, sep = "; ")) %>%
  mutate(taxon = factor(taxon, levels = taxon))
# Plot heatmap with abundance barplot
sf.hm.melted <- sf_mp_results %>%
  dplyr::select(taxon, names(sf_mp_corrs)[1:3]) %>%
  melt(., id.vars = c("taxon"))
sf.hm <- ggplot(data = sf.hm.melted, 
                   aes(x = factor(taxon), y = variable, fill = value)) + 
  geom_tile() + 
  scale_fill_distiller(name = "Indicator correlation index", palette = "RdBu", 
                       direction = -1, na.value = "transparent", type = "div", limits = c(-1, 1)) +
  scale_x_discrete(breaks = unique(sf.hm.melted$taxon), labels = unique(sf.hm.melted$taxon),
                   limits = rev(levels(sf.hm.melted$taxon))) +
  scale_y_discrete(position = "right") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        legend.title = element_text(size = 10), legend.key.height = unit(0.5,"cm"),
        legend.key.width = unit(1, "cm"), legend.text = element_text(size = 7)) +
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5)) +
  coord_flip()
sf.l <- get_legend(sf.hm)
sf.hm.clean <- sf.hm +
  theme(axis.title.y = element_blank(), axis.text.y = element_text(margin = margin(c(0,-5,0,0))),
        axis.ticks.y = element_blank(), axis.ticks.x = element_blank(), 
        panel.grid.major = element_blank(), legend.position="none",
        axis.text.x.top = element_text(size = 6, margin = margin(c(0,0,-2,0))), 
        plot.margin = margin(c(0,-2,0,5))) +
  labs(y = "Group")
sf.bp.y <- ggplot(data = sf_mp_results, aes(x = taxon, y = RelAbundance)) + 
  geom_bar(stat = "identity", fill = "grey") + 
  scale_y_continuous(position = "right") +
  scale_x_discrete(limits = rev(levels(sf_mp_results$taxon))) +
  coord_flip() + 
  theme_minimal() +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position="none", plot.margin = margin(c(0,5,0,-5)),
        axis.text.x.top = element_text(size = 6, margin = margin(c(0,0,-2,0)))) +
  labs(y = "Rel. abund. (%)")
pdf("InitialFigs/sf_Multipatt_genus.pdf", width = 8, height = 5)
plot_grid(sf.hm.clean, sf.bp.y, sf.l, NULL, nrow = 2, ncol = 2, 
          rel_widths = c(10,2), rel_heights = c(12, 2), align = "hv", axis = "b")
dev.off()



#### __BGC ####
# Add biogeochem. analysis
# Get variables, make corrplot, envfit, ordination with vectors, ordistep, taxa correlations
# Note, porewater chemistry is only on D2 samples (5-15 cm depth)

# Get variables
env <- sf$map_loaded %>%
  dplyr::select(CH4_ug_m2_h, CO2_ug_m2_h, 
                Salinity_ppt_all, SO4_mgL, DOC_mgL, Fe_mgL, Mn_mgL, Cu_mgL, Zn_mgL,
                sed_pH, sed_NH4_mgL, sed_NO3_mgL, sed_PO4_mgL, sed_Cl_mgL, sed_SO4_mgL,
                sed_per_C, sed_per_N, sed_CN, sed_Bulk_dens,
                sed_Fe_mgL, sed_Mn_mgL, sed_Cu_mgL, sed_Zn_mgL)
env_nona <- na.omit(env)

# Corrplot
C <- cor(env_nona)
corrplot(C, method = "number", type = "lower", 
         tl.cex = 0.5, tl.col = "black", number.cex = 0.5)
pdf("InitialFigs/SF_BGC_corr.pdf", width = 7, height = 5)
meth_corr_by_bgc(env_nona = env_nona)
dev.off()

# Envfit
pcoa <- cmdscale(sf_bc, k = nrow(sf$map_loaded) - 1, eig = T)
set.seed(100)
ef <- envfit(pcoa, env, permutations = 999, na.rm = TRUE)
ef
arrow_factor <- ordiArrowMul(ef)
manual_factor <- 0.5
vec.df <- as.data.frame(ef$vectors$arrows*sqrt(ef$vectors$r)) %>%
  mutate(Dim1 = Dim1 * manual_factor,
         Dim2 = Dim2 * manual_factor) %>%
  mutate(variables = rownames(.)) %>%
  filter(ef$vectors$pvals < 0.05) %>%
  filter(variables != "sed_Cl_mgL") %>%
  filter(variables != "sed_SO4_mgL") %>%
  mutate(shortnames = c("CH4", "CO2", "Salinity", "SO4", "DOC", "Cu",
                        "Zn", "sed_pH", "sed_NH4", "C", "N", "C:N",
                        "BD", "sed_Fe", "sed_Mn", "sed_Cu", "sed_Zn"))
# Note, deleted sediment CL and SO4 as they were very similar to porewater

# Plot with significant vectors
ordiplot(pcoa)
plot(ef, p.max = 0.05, cex = 0.5)
pcoaA1 <- round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[1]*100, digits = 1)
pcoaA2 <- round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[2]*100, digits = 1)
sf$map_loaded$Axis01 <- scores(pcoa)[,1]
sf$map_loaded$Axis02 <- scores(pcoa)[,2]
micro.hulls <- ddply(sf$map_loaded, c("Salt"), find_hull)
g <- ggplot(sf$map_loaded, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = Salt, fill = Salt),
               alpha = 0.1, show.legend = F) +
  geom_point(size = 3, alpha = 0.8, aes(colour = Salt, shape = Depth),
             show.legend = T) +
  labs(x = paste("PC1: ", pcoaA1, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2, "%", sep = ""),
       shape = "Depth (cm)",
       colour = "Salinity",
       fill = "Salinity") +
  scale_fill_viridis_d() +
  scale_colour_viridis_d(guide = guide_legend(override.aes = list(shape = 15))) +
  theme_bw() +  
  theme(legend.position = "right",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        legend.key.size = unit(0.4, "cm"),
        legend.margin = margin(2,2,2,-2, "pt"),
        legend.spacing = unit(0, "cm")) +
  geom_segment(data = vec.df,
               aes(x = 0, xend = Dim1, y = 0, yend = Dim2),
               arrow = arrow(length = unit(0.5, "cm")),
               colour = "gray", alpha = 0.5,
               inherit.aes = FALSE) + 
  geom_text(data = subset(vec.df, shortnames != "SO4"),
            aes(x = Dim1, y = Dim2, label = shortnames),
            size = 3, color = "black") +
  geom_text(data = subset(vec.df, shortnames == "SO4"),
            aes(x = Dim1, y = Dim2 - 0.015, label = shortnames),
            size = 3, color = "black") +
  coord_fixed()
g
pdf("InitialFigs/SF_PCoA_wBGC.pdf", width = 6, height = 4)
g
dev.off()
ggplotly(g)

# Ordistep
comm_nona <- as.data.frame(t(sf$data_loaded)) %>%
  filter(rownames(.) %in% rownames(env_nona))
mod0 <- rda(comm_nona ~ 1, env_nona)  # Model with intercept only
mod1 <- rda(comm_nona ~ ., env_nona)  # Model with all explanatory variables
mod <- ordistep(mod0, scope = formula(mod1))
mod
mod$anova #

# Plot correlations for different taxonomic levels at a given % relative abundance threshold
pdf("InitialFigs/SF_Phyla_CH4.pdf", width = 7, height = 5)
meth_corr_by_taxonomy(input = detra, level = 2, threshold = 0.5)
dev.off()
meth_corr_by_taxonomy(input = detra, level = 3, threshold = 0.5)
meth_corr_by_taxonomy(input = detra, level = 4, threshold = 0.5)
meth_corr_by_taxonomy(input = detra, level = 5, threshold = 0.5)
meth_corr_by_taxonomy(input = detra, level = 6, threshold = 0.5)
meth_corr_by_taxonomy(input = detra, level = 8, threshold = 0.5)
pdf("InitialFigs/SF_Guilds_CH4.pdf", width = 7, height = 5)
meth_corr_by_taxonomy(input = detra, level = 9, threshold = 0)
dev.off()



#### ...................................... ####
#### 7. Comparison Overview ####
# input_filt_rare <- readRDS("input_filt_rare_comb.rds")
input_filt_rare <- readRDS("input_filt_rare_comb_wBGC.rds")
input_filt_rare$map_loaded$Estuary <- factor(input_filt_rare$map_loaded$Estuary,
                                             levels = c("Waccamaw", "Alligator",
                                                        "Delaware", "SF"))
input_filt_rare$map_loaded$Depth2 <- recode_factor(input_filt_rare$map_loaded$Depth,
                                                   " 0-5" = "0-5",
                                                   " 5-15" = "5-15",
                                                   "0.02" = "0-5",
                                                   "0.025" = "0-5",
                                                   "0.1" = "5-15",
                                                   "0.12" = "5-15",
                                                   "0.125" = "5-15")
input_filt_rare_abund <- filter_taxa_from_input(input_filt_rare,
                                                filter_thresh = 0.05) # 94247 taxa removed



#### _Alpha ####
# OTU Richness 
leveneTest(input_filt_rare$map_loaded$rich ~ input_filt_rare$map_loaded$Estuary)
# Variance not homogeneous (p < 0.05)
m <- aov(input_filt_rare$map_loaded$rich ~ input_filt_rare$map_loaded$Estuary)
shapiro.test(m$residuals)
# Residuals not normally distributed (p < 0.05)
summary(m)

pdf("InitialFigs/Comb_All_Richness.pdf", width = 8, height = 4)
ggplot(input_filt_rare$map_loaded, aes(reorder(Info, rich, mean), rich, colour = Estuary)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 1, alpha = 0.75, width = 0.2) +
  labs(x = "Site", y = "Number of OTUs", colour = "Estuary") +
  scale_colour_viridis_d() +
  facet_grid(~ Estuary, scales = "free_x", space = "free") +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
        strip.text = element_text(size = 6.5),
        plot.margin = margin(0.1, 0.1, 0.1, 0.25, "cm"),
        panel.spacing.x = unit(0.05, "cm"))
dev.off()



# Shannon
leveneTest(input_filt_rare$map_loaded$shannon ~ input_filt_rare$map_loaded$Estuary)
# Variance homogeneous (p > 0.05)
m1 <- aov(input_filt_rare$map_loaded$shannon ~ input_filt_rare$map_loaded$Estuary)
shapiro.test(m1$residuals)
# Residuals not normally distributed (p < 0.05)
summary(m1)

pdf("InitialFigs/Comb_All_Shannon.pdf", width = 8, height = 4)
ggplot(input_filt_rare$map_loaded, aes(reorder(Info, shannon, mean), shannon, colour = Estuary)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 1, alpha = 0.75, width = 0.2) +
  labs(x = "Site", y = "Shannon diversity", colour = "Estuary") +
  scale_colour_viridis_d() +
  facet_grid(~ Estuary, scales = "free_x", space = "free") +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
        strip.text = element_text(size = 6.5),
        plot.margin = margin(0.1, 0.1, 0.1, 0.6, "cm"),
        panel.spacing.x = unit(0.05, "cm"))
dev.off()



#### _Beta  ####
bc <- calc_dm(input_filt_rare$data_loaded)
pcoa <- cmdscale(bc, k = nrow(input_filt_rare$map_loaded) - 1, eig = T)
pcoaA1 <- round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[1]*100, digits = 1)
pcoaA2 <- round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[2]*100, digits = 1)
input_filt_rare$map_loaded$Axis01 <- scores(pcoa)[,1]
input_filt_rare$map_loaded$Axis02 <- scores(pcoa)[,2]

# Explore the broad estuary comparison
micro.hulls <- ddply(input_filt_rare$map_loaded, c("Estuary"), find_hull)
pdf("InitialFigs/Comb_All_PCoA.pdf", width = 7, height = 5)
ggplot(input_filt_rare$map_loaded, aes(Axis01, Axis02, colour = Estuary)) +
  geom_polygon(data = micro.hulls, aes(colour = Estuary, fill = Estuary),
               alpha = 0.1, show.legend = F) +
  geom_point(size = 2, alpha = 0.5) +
  labs(x = paste("PC1: ", pcoaA1, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2, "%", sep = ""),
       colour = "Estuary") +
  scale_colour_viridis_d() +
  scale_fill_viridis_d() +
  theme_bw() +  
  theme(legend.position = c(1,1),
        legend.justification = c(1,1),
        legend.background = element_blank(),
        axis.title = element_text(face = "bold", size = 16), 
        axis.text = element_text(size = 14))
dev.off()

# Interactive, with sample or site info
g1 <- ggplot(input_filt_rare$map_loaded, aes(Axis01, Axis02, shape = Estuary, colour = Info)) +
  geom_point(size = 2, alpha = 0.5) +
  labs(x = paste("PC1: ", pcoaA1, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2, "%", sep = ""),
       shape = "Estuary",
       colour = "Info") +
  scale_colour_viridis_d() +
  scale_fill_viridis_d() +
  theme_bw() +  
  theme(legend.position = "right",
        axis.title = element_text(face = "bold", size = 16), 
        axis.text = element_text(size = 14))
g1
ggplotly(g1)

# Try with filtering out rare taxa. Use the 0.05% cutoff of Wyatt
# Overall pretty similar
bc_abund <- calc_dm(input_filt_rare_abund$data_loaded)
pcoa_abund <- cmdscale(bc_abund, k = nrow(input_filt_rare_abund$map_loaded) - 1, eig = T)
pcoaA1 <- round((eigenvals(pcoa_abund)/sum(eigenvals(pcoa_abund)))[1]*100, digits = 1)
pcoaA2 <- round((eigenvals(pcoa_abund)/sum(eigenvals(pcoa_abund)))[2]*100, digits = 1)
input_filt_rare_abund$map_loaded$Axis01 <- scores(pcoa_abund)[,1]
input_filt_rare_abund$map_loaded$Axis02 <- scores(pcoa_abund)[,2]
micro.hulls <- ddply(input_filt_rare_abund$map_loaded, c("Estuary"), find_hull)
pdf("InitialFigs/Comb_All_PCoA_0.05.pdf", width = 7, height = 5)
ggplot(input_filt_rare_abund$map_loaded, aes(Axis01, Axis02, colour = Estuary)) +
  geom_polygon(data = micro.hulls, aes(colour = Estuary, fill = Estuary),
               alpha = 0.1, show.legend = F) +
  geom_point(size = 2, alpha = 0.5) +
  labs(x = paste("PC1: ", pcoaA1, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2, "%", sep = ""),
       colour = "Estuary") +
  scale_colour_viridis_d() +
  scale_fill_viridis_d() +
  theme_bw() +  
  theme(legend.position = c(1,1),
        legend.justification = c(1,1),
        legend.background = element_blank(),
        axis.title = element_text(face = "bold", size = 16), 
        axis.text = element_text(size = 14))
dev.off()

# Stats
set.seed(1150)
adonis2(bc ~ input_filt_rare$map_loaded$Estuary)
set.seed(1150)
adonis2(bc_abund ~ input_filt_rare_abund$map_loaded$Estuary)

anova(betadisper(bc, input_filt_rare$map_loaded$Estuary)) # Dispersion not homogeneous
anova(betadisper(bc_abund, input_filt_rare_abund$map_loaded$Estuary)) # Dispersion not homogeneous

# Hierarchical clustering - color by estuary
# Use dendextend https://cran.r-project.org/web/packages/dendextend/vignettes/dendextend.html
clust <- hclust(bc, method = "ward.D2")
den <- as.dendrogram(clust)
colors_to_use <- as.numeric(as.factor(input_filt_rare$map_loaded$Estuary))
colors_to_use <- colors_to_use[order.dendrogram(den)]
colors_to_use
clusters <- colors_to_use
colors_to_use <- as.character(colors_to_use) %>%
  str_replace_all(., "4", "colorfour") %>%
  str_replace_all(., "3", "colorthree") %>%
  str_replace_all(., "2", "colortwo") %>%
  str_replace_all(., "1", "colorone") %>%
  str_replace_all(., "colorfour", viridis_pal()(4)[4]) %>%
  str_replace_all(., "colorthree", viridis_pal()(4)[3]) %>%
  str_replace_all(., "colortwo", viridis_pal()(4)[2]) %>%
  str_replace_all(., "colorone", viridis_pal()(4)[1])
labels_colors(den) <- colors_to_use
labels_colors(den)
colors_to_use
clusters
clusters <- as.character(clusters) %>%
  str_replace_all(., "4", "colorfour") %>%
  str_replace_all(., "3", "colorthree") %>%
  str_replace_all(., "2", "colortwo") %>%
  str_replace_all(., "1", "colorone") %>%
  str_replace_all(., "colorfour", "1") %>%
  str_replace_all(., "colorthree", "2") %>%
  str_replace_all(., "colortwo", "4") %>%
  str_replace_all(., "colorone", "3") %>%
  as.numeric()
den <- branches_attr_by_clusters(dend = den, 
                                 clusters = clusters, 
                                 values = colors_to_use,
                                 attr = "col")
labels_cex(den) <- 0.1
set(den, "branches_lwd" = c(0.1,0.1,0.1,0.1))
den <- hang.dendrogram(den, hang_height = 0.1)
pdf("InitialFigs/Comb_All_Clust.pdf", height = 8, width = 8)
par(oma = c(0, 0, 0, 0),
    mar = c(3,3,3,7))
plot(den, 
     horiz =  TRUE,  nodePar = list(cex = 0.007))
title("Hierarchical Clustering Dendrogram\n(Ward.D2)", adj = 0.5, line = -0.5)
legend(x = 6.4, y = 355,
       legend = levels(input_filt_rare$map_loaded$Estuary), 
       fill = viridis_pal()(4))
dev.off()



#### _Taxa ####
#### _Indicators ####
sim <- simper(t(input_filt_rare_abund$data_loaded), 
              input_filt_rare_abund$map_loaded$Estuary)
s <- summary(sim)
s1 <- as.data.frame(head(s$Delaware_Waccamaw, n = 10)) %>%
  mutate(Comparison = "Delaware_Waccamaw",
         ASV = rownames(.))
s2 <- as.data.frame(head(s$Delaware_Alligator, n = 10)) %>%
  mutate(Comparison = "Delaware_Alligator",
         ASV = rownames(.))
s3 <- as.data.frame(head(s$Delaware_SF, n = 10)) %>%
  mutate(Comparison = "Delaware_SF",
         ASV = rownames(.))
s4 <- as.data.frame(head(s$Waccamaw_Alligator, n = 10)) %>%
  mutate(Comparison = "Waccamaw_Alligator",
         ASV = rownames(.))
s5 <- as.data.frame(head(s$Waccamaw_SF, n = 10)) %>%
  mutate(Comparison = "Waccamaw_SF",
         ASV = rownames(.))
s6 <- as.data.frame(head(s$Alligator_SF, n = 10)) %>%
  mutate(Comparison = "Alligator_SF",
         ASV = rownames(.))
simper_results <- rbind(s1, s2, s3, s4, s5, s6) %>%
  left_join(., input_filt_rare_abund$taxonomy_loaded, by = c("ASV" = "taxonomy8"))
saveRDS(simper_results, "simper_results.RDS")
simper_results <- readRDS("simper_results.RDS")
length(unique(simper_results$ASV)) # 27 unique ASVs
write_xlsx(simper_results, 
           "simper_results_Comb.xlsx",
           format_headers = F)

# MULTIPATT (list ASVs associated with each group)
# Don't plot here because later do at Field Control/Lab/Field Exp. level
set.seed(1202)
mp <- multipatt(t(input_filt_rare_abund$data_loaded), 
                input_filt_rare_abund$map_loaded$Estuary, 
                func = "IndVal.g", 
                control = how(nperm=999))
summary(mp)
View(mp$sign)
multipatt_results <- mp$sign %>%
  filter(p.value == 0.001) %>%
  mutate(num_sites = rowSums(cbind(s.Waccamaw, s.Alligator, s.Delaware, s.SF))) %>%
  filter(num_sites == 1) %>%
  mutate(Group = "NA")
for (i in 1:nrow(multipatt_results)) {
  if (multipatt_results$s.Waccamaw[i] == 1) {
    multipatt_results$Group[i] <- "Waccamaw"
  }
}
for (i in 1:nrow(multipatt_results)) {
  if (multipatt_results$s.Alligator[i] == 1) {
    multipatt_results$Group[i] <- "Alligator"
  }
}
for (i in 1:nrow(multipatt_results)) {
  if (multipatt_results$s.Delaware[i] == 1) {
    multipatt_results$Group[i] <- "Delaware"
  }
}
for (i in 1:nrow(multipatt_results)) {
  if (multipatt_results$s.SF[i] == 1) {
    multipatt_results$Group[i] <- "SF"
  }
}
table(multipatt_results$Group)

#### _Domain ####
tax_sum_domain <- summarize_taxonomy(input_filt_rare, level = 1, report_higher_tax = F)
plot_ts_heatmap(tax_sum_domain, 
                input_filt_rare$map_loaded, 
                0.01, 
                'Estuary',
                rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
bars <- plot_taxa_bars(tax_sum_domain,
                       input_filt_rare$map_loaded,
                       "Estuary",
                       num_taxa = 10,
                       data_only = TRUE)
ggplot(bars, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Estuary", y = "Relative abundance", fill = "Domain") +
  scale_fill_manual(values = c("red", "blue")) +
  scale_y_continuous(expand = c(0.01, 0.01)) + 
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10))
taxa_summary_by_sample_type(tax_sum_domain, 
                            input_filt_rare$map_loaded, 
                            type_header = 'Estuary', 
                            filter_level = 0.01, 
                            test_type = 'KW')

#### _Phylum ####
tax_sum_phyla <- summarize_taxonomy(input_filt_rare, level = 2, report_higher_tax = F)
plot_ts_heatmap(tax_sum_phyla, 
                input_filt_rare$map_loaded, 
                0.001, 
                'Estuary',
                rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
bars <- plot_taxa_bars(tax_sum_phyla,
                       input_filt_rare$map_loaded,
                       "Estuary",
                       num_taxa = 12,
                       data_only = TRUE) %>%
  mutate(taxon = fct_rev(taxon))
pdf("InitialFigs/Comb_All_Phyla.pdf", width = 7, height = 5)
ggplot(bars, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Estuary", y = "Relative abundance", fill = "Phylum") +
  scale_fill_manual(values = c("grey90", brewer.pal(12, "Paired")[12:1])) +
  scale_y_continuous(expand = c(0.01, 0.01)) + 
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 10))
dev.off()
taxa_summary_by_sample_type(tax_sum_phyla, 
                            input_filt_rare$map_loaded, 
                            type_header = 'Estuary', 
                            filter_level = 0.01, 
                            test_type = 'KW')

# Show all phyla
bars <- plot_taxa_bars(tax_sum_phyla,
                       input_filt_rare$map_loaded,
                       "Estuary",
                       num_taxa = nrow(tax_sum_phyla),
                       data_only = TRUE)
nb.cols <- nrow(tax_sum_phyla)
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)
pdf("InitialFigs/Comb_All_Phyla_All.pdf", width = 8, height = 6)
ggplot(bars, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Estuary", y = "Relative abundance", fill = "Phylum") +
  scale_fill_manual(values = mycolors) +
  scale_y_continuous(expand = c(0.01, 0.01)) + 
  guides(fill = guide_legend(ncol = 2)) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 8),
        legend.key.size = unit(0.1, "cm"),
        legend.text = element_text(size = 8))
dev.off()

# Show all samples
input_filt_rare$map_loaded$sampleID <- rownames(input_filt_rare$map_loaded)
bars <- plot_taxa_bars(tax_sum_phyla,
                       input_filt_rare$map_loaded,
                       "sampleID",
                       num_taxa = 12,
                       data_only = TRUE) %>%
  mutate(taxon = fct_rev(taxon)) %>%
  left_join(., input_filt_rare$map_loaded, by = c("group_by" = "sampleID"))
pdf("InitialFigs/Comb_All_Phyla_samples.pdf", width = 12, height = 6)
ggplot(bars, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Sample", y = "Relative abundance", fill = "Phylum") +
  scale_fill_manual(values = c("grey90", brewer.pal(12, "Paired")[12:1])) +
  scale_y_continuous(expand = c(max(tallest_bar$sum)/100, max(tallest_bar$sum)/100)) + 
  facet_nested(~ Estuary, space = "free", scales = "free_x") +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 2, angle = 90, hjust = 1, vjust = 0.5),
        strip.text = element_text(size = 6),
        strip.background = element_rect(size = 0.2),
        axis.line.y = element_blank())
dev.off()

# Look at archaeal phyla
tax_sum_phyla_ar <- summarize_taxonomy(input_filt_rare, level = 2, report_higher_tax = T)
tax_sum_phyla_ar <- tax_sum_phyla_ar[grep("Archaea", rownames(tax_sum_phyla_ar)),]
bars_ar <- plot_taxa_bars(tax_sum_phyla_ar,
                          input_filt_rare$map_loaded,
                          "Estuary",
                          num_taxa = nrow(tax_sum_phyla_ar),
                          data_only = TRUE)
nb.cols <- nrow(tax_sum_phyla_ar)
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)
ggplot(bars_ar, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Estuary", y = "Relative abundance", fill = "Phylum") +
  scale_fill_manual(values = mycolors) +
  scale_y_continuous(expand = c(0.01, 0.01)) + 
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 10))

#### _Class ####
tax_sum_class <- summarize_taxonomy(input_filt_rare, level = 3, report_higher_tax = F)
plot_ts_heatmap(tax_sum_class, 
                input_filt_rare$map_loaded, 
                0.01, 
                'Estuary',
                rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
bars <- plot_taxa_bars(tax_sum_class,
                       input_filt_rare$map_loaded,
                       "Estuary",
                       num_taxa = 12, data_only = T) %>%
  mutate(taxon = gsub("NA", "Unclassified", taxon)) %>%
  mutate(taxon = fct_relevel(taxon, "Other", after = Inf)) %>%
  mutate(taxon = fct_relevel(taxon, "Unclassified", after = Inf)) %>%
  mutate(taxon = fct_rev(taxon))
ggplot(bars, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Estuary", y = "Relative abundance", fill = "Class") +
  scale_fill_manual(values = c("grey75", "grey90", brewer.pal(12, "Paired")[11:1])) +
  scale_y_continuous(expand = c(0.01, 0.01)) + 
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 10))
taxa_summary_by_sample_type(tax_sum_class, 
                            input_filt_rare$map_loaded, 
                            type_header = 'Estuary', 
                            filter_level = 0.01, 
                            test_type = 'KW')

# Check sulfur reducers (not all, just quick Desulo check)
tax_sum_class_su <- summarize_taxonomy(input_filt_rare, level = 3, report_higher_tax = T)
tax_sum_class_su <- tax_sum_class_su[grep("Desulfo", rownames(tax_sum_class_su)),]
bars_su <- plot_taxa_bars(tax_sum_class_su,
                          input_filt_rare$map_loaded,
                          "Estuary",
                          num_taxa = nrow(tax_sum_class_su),
                          data_only = TRUE) %>%
  mutate(taxon = substring(taxon, 11))
nb.cols <- nrow(tax_sum_class_su)
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)
tallest_bar <- bars_su %>%
  group_by(group_by) %>%
  summarise(sum = sum(mean_value))
ggplot(bars_su, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Estuary", y = "Relative abundance", fill = "Class") +
  scale_fill_manual(values = mycolors) +
  scale_y_continuous(expand = c(max(tallest_bar$sum)/100, max(tallest_bar$sum)/100)) + 
  theme_classic() +
  theme(axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10),
        axis.title.x = element_blank())

#### _Order ####
tax_sum_order <- summarize_taxonomy(input_filt_rare, level = 4, report_higher_tax = F)
plot_ts_heatmap(tax_sum_order, 
                input_filt_rare$map_loaded, 
                0.01, 
                'Estuary',
                rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
bars <- plot_taxa_bars(tax_sum_order,
                       input_filt_rare$map_loaded,
                       "Estuary",
                       num_taxa = 12, data_only = T) %>%
  mutate(taxon = gsub("NA", "Unclassified", taxon)) %>%
  mutate(taxon = fct_relevel(taxon, "Other", after = Inf)) %>%
  mutate(taxon = fct_relevel(taxon, "Unclassified", after = Inf)) %>%
  mutate(taxon = fct_rev(taxon))
ggplot(bars, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Estuary", y = "Relative abundance", fill = "Order") +
  scale_fill_manual(values = c("grey75", "grey90", brewer.pal(12, "Paired")[11:1])) +
  scale_y_continuous(expand = c(0.01, 0.01)) + 
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 10))
taxa_summary_by_sample_type(tax_sum_order, 
                            input_filt_rare$map_loaded, 
                            type_header = 'Estuary', 
                            filter_level = 0.01, 
                            test_type = 'KW')

#### _Family ####
tax_sum_families <- summarize_taxonomy(input_filt_rare, level = 5, report_higher_tax = FALSE)
plot_ts_heatmap(tax_sum_families, 
                input_filt_rare$map_loaded, 
                0.01, 
                'Estuary',
                rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
bars <- plot_taxa_bars(tax_sum_families,
                       input_filt_rare$map_loaded,
                       "Estuary",
                       num_taxa = 12, data_only = T) %>%
  mutate(taxon = gsub("NA", "Unclassified", taxon)) %>%
  mutate(taxon = fct_relevel(taxon, "Other", after = Inf)) %>%
  mutate(taxon = fct_relevel(taxon, "Unclassified", after = Inf)) %>%
  mutate(taxon = fct_rev(taxon))
ggplot(bars, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Estuary", y = "Relative abundance", fill = "Family") +
  scale_fill_manual(values = c("grey75", "grey90", brewer.pal(12, "Paired")[11:1])) +
  scale_y_continuous(expand = c(0.01, 0.01)) + 
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 10))
taxa_summary_by_sample_type(tax_sum_families, 
                            input_filt_rare$map_loaded, 
                            type_header = 'Estuary', 
                            filter_level = 0.01, 
                            test_type = 'KW')

# Look at methanogens (careful to do this accurately!)
tax_sum_families_meth <- summarize_taxonomy(input_filt_rare, level = 5, 
                                            report_higher_tax = F)
tax_sum_families_meth <- tax_sum_families_meth[grep("(Methano|Methermicoccaceae|
                                                    Syntrophoarchaeaceae)", 
                                                    rownames(tax_sum_families_meth)),]
tax_sum_families_meth <- tax_sum_families_meth[!grepl("Methanoperedenaceae", 
                                                      rownames(tax_sum_families_meth)),]
bars_meth <- plot_taxa_bars(tax_sum_families_meth,
                            input_filt_rare$map_loaded,
                            "Estuary",
                            num_taxa = nrow(tax_sum_families_meth),
                            data_only = TRUE)
nb.cols <- nrow(tax_sum_families_meth)
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)
tallest_bar <- bars_meth %>%
  group_by(group_by) %>%
  summarise(sum = sum(mean_value))
ggplot(bars_meth, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Estuary", y = "Relative abundance", fill = "Family") +
  scale_fill_manual(values = mycolors) +
  scale_y_continuous(expand = c(max(tallest_bar$sum)/100, max(tallest_bar$sum)/100)) + 
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 10))

#### _Genus ####
tax_sum_genera <- summarize_taxonomy(input_filt_rare, level = 6, report_higher_tax = TRUE)
plot_ts_heatmap(tax_sum_genera, 
                input_filt_rare$map_loaded, 
                0.01, 
                'Estuary',
                rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
bars <- plot_taxa_bars(tax_sum_genera,
                       input_filt_rare$map_loaded,
                       "Estuary",
                       num_taxa = 10,
                       data_only = TRUE)
ggplot(bars, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Estuary", y = "Relative abundance", fill = "Genus") +
  scale_fill_manual(values = c("grey75", "grey90", brewer.pal(12, "Paired")[11:1])) +
  scale_y_continuous(expand = c(0.01, 0.01)) + 
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 10))
taxa_summary_by_sample_type(tax_sum_genera, 
                            input_filt_rare$map_loaded, 
                            type_header = 'Estuary', 
                            filter_level = 0.01, 
                            test_type = 'KW')

#### _Guilds ####
tax_sum_guilds <- summarize_taxonomy(input_filt_rare, level = 9, report_higher_tax = F)
plot_ts_heatmap(tax_sum_guilds, 
                input_filt_rare$map_loaded, 
                0, 
                'Estuary',
                rev_taxa = T,
                remove_other = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
bars <- plot_taxa_bars(tax_sum_guilds,
                       input_filt_rare$map_loaded,
                       "Estuary",
                       num_taxa = 20,
                       data_only = TRUE) %>%
  filter(taxon != "NA") %>%
  droplevels() %>%
  mutate(taxon = factor(taxon,
                        levels = Guild_cols$Guild))
tallest_bar <- bars %>%
  group_by(group_by) %>%
  summarise(sum = sum(mean_value))
pdf("InitialFigs/Comb_All_Guilds.pdf", width = 7, height = 5)
ggplot(bars, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Estuary", y = "Relative abundance", fill = "Guild") +
  scale_fill_manual(values = Guild_cols$color) +
  scale_y_continuous(expand = c(max(tallest_bar$sum)/100, max(tallest_bar$sum)/100)) + 
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10))
dev.off()
taxa_summary_by_sample_type(tax_sum_guilds, 
                            input_filt_rare$map_loaded, 
                            type_header = 'Estuary', 
                            filter_level = 0.01, 
                            test_type = 'KW')

# Guilds by sample (to show variation)
input_filt_rare$map_loaded$sampleID <- rownames(input_filt_rare$map_loaded)
bars <- plot_taxa_bars(tax_sum_guilds,
                       input_filt_rare$map_loaded,
                       "sampleID",
                       num_taxa = 20,
                       data_only = TRUE) %>%
  filter(taxon != "NA") %>%
  droplevels() %>%
  mutate(taxon = factor(taxon,
                        levels = Guild_cols$Guild)) %>%
  left_join(., input_filt_rare$map_loaded, by = c("group_by" = "sampleID"))
tallest_bar <- bars %>%
  group_by(group_by) %>%
  summarise(sum = sum(mean_value))
pdf("InitialFigs/Comb_All_Guilds_samples.pdf", width = 12, height = 6)
ggplot(bars, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Sample", y = "Relative abundance", fill = "Guild") +
  scale_fill_manual(values = Guild_cols$color) +
  scale_y_continuous(expand = c(max(tallest_bar$sum)/100, max(tallest_bar$sum)/100)) + 
  facet_nested(~ Estuary, space = "free", scales = "free_x") +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 2, angle = 90, hjust = 1, vjust = 0.5),
        strip.text = element_text(size = 6),
        strip.background = element_rect(size = 0.2),
        axis.line.y = element_blank())
dev.off()

#### _Venn ####
phy <- summarize_taxonomy(input_filt_rare, level = 2, report_higher_tax = F)
cla <- summarize_taxonomy(input_filt_rare, level = 3, report_higher_tax = F)
ord <- summarize_taxonomy(input_filt_rare, level = 4, report_higher_tax = F)
fam <- summarize_taxonomy(input_filt_rare, level = 5, report_higher_tax = F)
gen <- summarize_taxonomy(input_filt_rare, level = 6, report_higher_tax = F)

input_phylum <- input_filt_rare
input_phylum$data_loaded <- phy
input_class <- input_filt_rare
input_class$data_loaded <- cla
input_order <- input_filt_rare
input_order$data_loaded <- ord
input_family <- input_filt_rare
input_family$data_loaded <- fam
input_genus <- input_filt_rare
input_genus$data_loaded <- gen

pdf("InitialFigs/Comb_All_VennOTU.pdf", width = 7, height = 5)
plot_venn_diagram(input_filt_rare,
                  "Estuary",
                  0.00000000000000001)
dev.off()

pdf("InitialFigs/Comb_All_VennAll.pdf", width = 15, height = 6.5)
plot_grid(plot_venn_diagram(input_phylum, "Estuary", 0.00000000000000001),
          plot_venn_diagram(input_class, "Estuary", 0.00000000000000001),
          plot_venn_diagram(input_order, "Estuary", 0.00000000000000001),
          plot_venn_diagram(input_family, "Estuary", 0.00000000000000001),
          plot_venn_diagram(input_genus, "Estuary", 0.00000000000000001),
          plot_venn_diagram(input_filt_rare, "Estuary", 0.00000000000000001),
          labels = c("(a) Phylum", "(b) Class", "(c) Order", 
                     "(d) Family", "(e) Genus", "(f) OTU"))
dev.off()

#### _BGC ####
# Add biogeochem. analysis
# Get variables, make corrplot, envfit, ordination with vectors, ordistep, taxa correlations
# Note, porewater chemistry is only on D2 samples (5-15 cm depth)

# Get variables
env <- input_filt_rare$map_loaded %>%
  dplyr::select(CH4_ug_m2_h, CO2_ug_m2_h, Salinity_ppt_all, pH)
env_nona <- na.omit(env) # n = 208

# Corrplot
C <- cor(env_nona)
corrplot(C, method = "number", type = "lower", 
         tl.cex = 0.5, tl.col = "black", number.cex = 0.5)
pdf("InitialFigs/Comb_All_BGC_corr.pdf", width = 7, height = 5)
meth_corr_by_bgc(env_nona = env_nona)
dev.off()

# Envfit
pcoa <- cmdscale(bc, k = nrow(input_filt_rare$map_loaded) - 1, eig = T)
set.seed(100)
ef <- envfit(pcoa, env, permutations = 999, na.rm = TRUE)
ef
arrow_factor <- ordiArrowMul(ef)
manual_factor <- 0.36
vec.df <- as.data.frame(ef$vectors$arrows*sqrt(ef$vectors$r)) %>%
  mutate(Dim1 = Dim1 * manual_factor,
         Dim2 = Dim2 * manual_factor) %>%
  mutate(variables = rownames(.)) %>%
  filter(ef$vectors$pvals < 0.05) %>%
  mutate(shortnames = c("CH4", "CO2", "Salinity", "pH"))

# Plot with significant vectors
ordiplot(pcoa)
plot(ef, p.max = 0.05, cex = 0.5)
pcoaA1 <- round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[1]*100, digits = 1)
pcoaA2 <- round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[2]*100, digits = 1)
input_filt_rare$map_loaded$Axis01 <- scores(pcoa)[,1]
input_filt_rare$map_loaded$Axis02 <- scores(pcoa)[,2]
micro.hulls <- ddply(input_filt_rare$map_loaded, c("Estuary"), find_hull)
g <- ggplot(input_filt_rare$map_loaded, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = Estuary, fill = Estuary),
               alpha = 0.1, show.legend = F) +
  geom_point(size = 3, alpha = 0.5, aes(colour = Estuary, shape = Depth2),
             show.legend = T) +
  labs(x = paste("PC1: ", pcoaA1, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2, "%", sep = ""),
       shape = "Depth (cm)",
       colour = "Estuary",
       fill = "Estuary") +
  scale_fill_viridis_d(guide = guide_legend(override.aes = list(alpha = 1))) +
  scale_colour_viridis_d(guide = guide_legend(override.aes = list(alpha = 1, shape = 15))) +
  guides(shape = guide_legend(override.aes = list(alpha = 1))) +
  theme_bw() +  
  theme(legend.position = "right",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10)) +
  geom_segment(data = vec.df,
               aes(x = 0, xend = Dim1, y = 0, yend = Dim2),
               arrow = arrow(length = unit(0.5, "cm")),
               colour = "gray", alpha = 0.5,
               inherit.aes = FALSE) + 
  geom_text(data = vec.df,
            aes(x = Dim1, y = Dim2, label = shortnames),
            size = 3, color = "black") +
  coord_fixed()
g
pdf("InitialFigs/Comb_All_PCoA_wBGC.pdf", width = 6, height = 4)
g
dev.off()

# Ordistep
comm_nona <- as.data.frame(t(input_filt_rare$data_loaded)) %>%
  filter(rownames(.) %in% rownames(env_nona))
mod0 <- rda(comm_nona ~ 1, env_nona)  # Model with intercept only
mod1 <- rda(comm_nona ~ ., env_nona)  # Model with all explanatory variables
mod <- ordistep(mod0, scope = formula(mod1))
mod
mod$anova #

# Plot correlations for different taxonomic levels at a given % relative abundance threshold
pdf("InitialFigs/Comb_All_Phyla_CH4.pdf", width = 7, height = 5)
meth_corr_by_taxonomy(input = input_filt_rare, level = 2, threshold = 0.5)
dev.off()
meth_corr_by_taxonomy(input = input_filt_rare, level = 3, threshold = 0.5)
meth_corr_by_taxonomy(input = input_filt_rare, level = 4, threshold = 0.5)
meth_corr_by_taxonomy(input = input_filt_rare, level = 5, threshold = 0.5)
meth_corr_by_taxonomy(input = input_filt_rare, level = 6, threshold = 0.5)
meth_corr_by_taxonomy(input = input_filt_rare, level = 8, threshold = 0.5) 
# ASV4, ASV 3 pos., ASV2 neg.
pdf("InitialFigs/Comb_All_Guilds_CH4.pdf", width = 7, height = 5)
meth_corr_by_taxonomy(input = input_filt_rare, level = 9, threshold = 0)
dev.off()



#### ...................................... ####
#### 8. Comparison Field Control ####
# Only look at field soils (no manipulations or incubations)
# Classify as fresh/oligo/meso/poly for broad salinity grouping

#### _Setup ####
# input_filt_rare <- readRDS("input_filt_rare_comb.rds")
input_filt_rare <- readRDS("input_filt_rare_comb_wBGC.rds")
input_filt_rare$map_loaded <- input_filt_rare$map_loaded %>%
  mutate(Estuary = factor(Estuary,
                          levels = c("Waccamaw", "Alligator", "Delaware", "SF")),
         sampleID = rownames(.),
         Field = "NA",
         Salt = "NA",
         Depth = gsub("0-5", 0.05, Depth)) %>%
  mutate(Depth = gsub("5-15", 0.15, Depth)) %>%
  mutate(Depth = as.numeric(Depth)) %>%
  mutate(Depth = ifelse(Depth <= 0.05, "< 5 cm", "5 - 15 cm"))

# Add column for if unmanipulated field sample
for (i in 1:nrow(input_filt_rare$map_loaded)) {
  if (input_filt_rare$map_loaded$Estuary[i] == "SF") {
    input_filt_rare$map_loaded$Field[i] <- "Field"
  }
}
for (i in 1:nrow(input_filt_rare$map_loaded)) {
  if (input_filt_rare$map_loaded$Estuary[i] == "Delaware" & 
      input_filt_rare$map_loaded$Site[i] == "Soil") {
    input_filt_rare$map_loaded$Field[i] <- "Field"
  }
}
for (i in 1:nrow(input_filt_rare$map_loaded)) {
  if (input_filt_rare$map_loaded$Estuary[i] == "Alligator" & 
      input_filt_rare$map_loaded$Site[i] == "Soil") {
    input_filt_rare$map_loaded$Field[i] <- "Field"
  }
}
for (i in 1:nrow(input_filt_rare$map_loaded)) {
  if (input_filt_rare$map_loaded$Estuary[i] == "Waccamaw" & 
      input_filt_rare$map_loaded$Detail[i] == "Control") {
    input_filt_rare$map_loaded$Field[i] <- "Field"
  }
}
field <- filter_data(input_filt_rare,
                     filter_cat = "Field",
                     keep_vals = "Field")

# Add column for Salt
for (i in 1:nrow(field$map_loaded)) {
  if (field$map_loaded$Site[i] == "Sandmound" | field$map_loaded$Site[i] == "West Pond") {
    field$map_loaded$Salt[i] <- "Freshwater"
  }
}
for (i in 1:nrow(field$map_loaded)) {
  if (field$map_loaded$Site[i] == "Mayberry" | field$map_loaded$Site[i] == "Browns") {
    field$map_loaded$Salt[i] <- "Oligohaline"
  }
}
for (i in 1:nrow(field$map_loaded)) {
  if (field$map_loaded$Site[i] == "Joice" | field$map_loaded$Site[i] == "Rush Ranch") {
    field$map_loaded$Salt[i] <- "Mesohaline"
  }
}
for (i in 1:nrow(field$map_loaded)) {
  if (field$map_loaded$Site[i] == "Goodyear" | field$map_loaded$Site[i] == "White Slough" |
      field$map_loaded$Site[i] == "Tolay" | field$map_loaded$Site[i] == "China Camp" |
      field$map_loaded$Site[i] == "Muzzi") {
    field$map_loaded$Salt[i] <- "Polyhaline"
  }
}
for (i in 1:nrow(field$map_loaded)) {
  if (field$map_loaded$Estuary[i] == "Waccamaw") {
    field$map_loaded$Salt[i] <- "Freshwater"
  }
}
for (i in 1:nrow(field$map_loaded)) {
  if (field$map_loaded$Estuary[i] == "Alligator") {
    field$map_loaded$Salt[i] <- "Freshwater"
  }
}
for (i in 1:nrow(field$map_loaded)) {
  if (field$map_loaded$Info[i] == "OligoHal_source") {
    field$map_loaded$Salt[i] <- "Oligohaline"
  }
}
for (i in 1:nrow(field$map_loaded)) {
  if (field$map_loaded$Info[i] == "MesoHal_source") {
    field$map_loaded$Salt[i] <- "Mesohaline"
  }
}
for (i in 1:nrow(field$map_loaded)) {
  if (field$map_loaded$Info[i] == "TFM1_source" | field$map_loaded$Info[i] == "TFM2_source") {
    field$map_loaded$Salt[i] <- "Freshwater"
  }
}
field$map_loaded <- field$map_loaded %>%
  unite("EstSalt", c(Estuary, Salt), sep = "_", remove = F) %>%
  mutate(EstSalt = factor(EstSalt,
                          levels = c("SF_Freshwater", "Alligator_Freshwater",
                                     "Delaware_Freshwater", "Waccamaw_Freshwater",
                                     "SF_Oligohaline", "Delaware_Oligohaline",
                                     "SF_Mesohaline", "Delaware_Mesohaline",
                                     "SF_Polyhaline")),
         Salt = factor(Salt,
                       levels = c("Freshwater", "Oligohaline", 
                                  "Mesohaline", "Polyhaline")))

#### _Alpha ####
leveneTest(field$map_loaded$rich ~ field$map_loaded$Salt) # Almost homogeneous
m <- aov(rich ~ Estuary + Salt + Depth, data = field$map_loaded)
Anova(m, type = "II") # Estuary, Salt sig. Depth marginal.
m <- aov(rich ~ Salt, data = field$map_loaded)
shapiro.test(m$residuals) # Normal
summary(m)
t <- emmeans(object = m, specs = "Salt") %>%
  cld(object = ., adjust = "Tukey", Letters = letters, alpha = 0.05) %>%
  mutate(name = "rich",
         y = max(field$map_loaded$rich)+(max(field$map_loaded$rich)-min(field$map_loaded$rich))/20)
leveneTest(field$map_loaded$shannon ~ field$map_loaded$Salt) # Homogeneous
m1 <- aov(shannon ~ Estuary + Salt + Depth, data = field$map_loaded)
Anova(m1, type = "II") # All sig
m1 <- aov(shannon ~ Salt, data = field$map_loaded)
shapiro.test(m1$residuals) # Not normal
summary(m1)
t1 <- emmeans(object = m1, specs = "Salt") %>%
  cld(object = ., adjust = "Tukey", Letters = letters, alpha = 0.05) %>%
  mutate(name = "shannon",
         y = max(field$map_loaded$shannon)+(max(field$map_loaded$shannon)-min(field$map_loaded$shannon))/20)
label_df <- rbind(t, t1)
facet_df <- c("rich" = "(a) Richness",
              "shannon" = "(b) Shannon")
alpha_long <- field$map_loaded %>%
  pivot_longer(cols = c("rich", "shannon"))
pdf("InitialFigs/Comb_Control_Alpha.pdf", width = 8, height = 3)
ggplot(alpha_long, aes(reorder(Salt, value, mean), value)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 2, alpha = 0.75, width = 0.2, aes(shape = Depth, colour = Estuary)) +
  geom_text(data = label_df, aes(Salt, y, label = str_trim(.group)), 
            size = 4, color = "black") +
  labs(x = "Site", y = "Number of OTUs") +
  scale_colour_viridis_d() +
  facet_wrap(~ name, ncol = 2, scales = "free_y", labeller = as_labeller(facet_df)) +
  theme_bw() +
  theme(legend.position = "right",
        axis.title = element_blank(),
        axis.text = element_text(size = 10),
        strip.text = element_text(size = 10))
dev.off()

#### _Beta ####
field_bc <- calc_dm(field$data_loaded)
set.seed(1150)
adonis2(field_bc ~ field$map_loaded$Estuary+field$map_loaded$Salt+field$map_loaded$Depth) # All
anova(betadisper(field_bc, field$map_loaded$Estuary)) # Dispersion not homogeneous
anova(betadisper(field_bc, field$map_loaded$Salt)) # Dispersion not homogeneous
anova(betadisper(field_bc, field$map_loaded$Depth)) # Dispersion homogeneous
field_pcoa <- cmdscale(field_bc, k = nrow(field$map_loaded) - 1, eig = T)
pcoaA1 <- round((eigenvals(field_pcoa)/sum(eigenvals(field_pcoa)))[1]*100, digits = 1)
pcoaA2 <- round((eigenvals(field_pcoa)/sum(eigenvals(field_pcoa)))[2]*100, digits = 1)
field$map_loaded$Axis01 <- scores(field_pcoa)[,1]
field$map_loaded$Axis02 <- scores(field_pcoa)[,2]
micro.hulls <- ddply(field$map_loaded, c("Salt"), find_hull)
pdf("InitialFigs/Comb_Control_PCoA.pdf", width = 7, height = 5)
ggplot(field$map_loaded, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = Salt, fill = Salt),
               alpha = 0.1, show.legend = F) +
  geom_point(size = 3, alpha = 1, aes(colour = Salt, shape = Estuary)) +
  labs(x = paste("PC1: ", pcoaA1, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2, "%", sep = ""),
       colour = "Salinity") +
  scale_fill_viridis_d() +
  scale_colour_viridis_d(guide = guide_legend(override.aes = list(shape = 15))) +
  scale_shape_manual(values = c(16, 17, 18, 3)) +
  theme_bw() +  
  theme(legend.position = "right",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10))
dev.off()

#### _BC Among Within ####
# Plot BC dissimilarity across site and salt combos
# Is Fresh vs. Poly at same site more different than Fresh vs Fresh at different sites?
# Compare SF_Fresh to the other three fresh and to SF_Oligo, SF_Meso, SF_Poly

# Convert from dist object to matrix object
bac_bray_mat <- as.matrix(field_bc)
# Remove duplicates and diagonal by setting upper triangle and diagonal to NA
bac_bray_mat[upper.tri(bac_bray_mat, diag = TRUE)] <- NA
# Convert from matrix object to dataframe
bac_bray_df <- as.data.frame(bac_bray_mat)
# Make a sample ID column from the row names.The other sample in each pairwise comparison will be the variable column after melting
bac_bray_df$sampleID <- rownames(bac_bray_df)
# Melt with reshape2 package
bac_bray_df_long <- melt(bac_bray_df, id.vars = "sampleID")
# Get rid of the NA's (duplicates and diagonal)
bac_bray_df_long <- na.omit(bac_bray_df_long)
# Note, the length of this dataframe should now equal (n*(n-1))/2
nrow(bac_bray_df_long) == (nrow(field$map_loaded)*(nrow(field$map_loaded)-1))/2 # Good!
# Make sampleID a factor
bac_bray_df_long$sampleID <- as.factor(bac_bray_df_long$sampleID)
# Now add EstSalt, matching to col1 and col2
# This will give the estuary and salinity for each of the two samples in each pairwise comparison
EstSalt <- dplyr::select(field$map_loaded, sampleID, EstSalt)
bac_bray_df_long <- inner_join(bac_bray_df_long, EstSalt, 
                               by = c("sampleID" = "sampleID"))
bac_bray_df_long <- inner_join(bac_bray_df_long, EstSalt, 
                               by = c("variable" = "sampleID"))
# Make new column for comparison and filter to comparisons of interest
bac_bray_df_long$comparison <- "NA"
for (i in 1:nrow(bac_bray_df_long)) {
  if (bac_bray_df_long$EstSalt.x[i] == "SF_Freshwater" &
      bac_bray_df_long$EstSalt.y[i] == "Delaware_Freshwater") {
    bac_bray_df_long$comparison[i] <- "Freshwater, Between Sites"
  }
  if (bac_bray_df_long$EstSalt.x[i] == "SF_Freshwater" &
      bac_bray_df_long$EstSalt.y[i] == "Alligator_Freshwater") {
    bac_bray_df_long$comparison[i] <- "Freshwater, Between Sites"
  }
  if (bac_bray_df_long$EstSalt.x[i] == "SF_Freshwater" &
      bac_bray_df_long$EstSalt.y[i] == "Waccamaw_Freshwater") {
    bac_bray_df_long$comparison[i] <- "Freshwater, Between Sites"
  }
  if (bac_bray_df_long$EstSalt.x[i] == "Delaware_Freshwater" &
      bac_bray_df_long$EstSalt.y[i] == "Alligator_Freshwater") {
    bac_bray_df_long$comparison[i] <- "Freshwater, Between Sites"
  }
  if (bac_bray_df_long$EstSalt.x[i] == "Delaware_Freshwater" &
      bac_bray_df_long$EstSalt.y[i] == "Waccamaw_Freshwater") {
    bac_bray_df_long$comparison[i] <- "Freshwater, Between Sites"
  }
  if (bac_bray_df_long$EstSalt.x[i] == "Alligator_Freshwater" &
      bac_bray_df_long$EstSalt.y[i] == "Waccamaw_Freshwater") {
    bac_bray_df_long$comparison[i] <- "Freshwater, Between Sites"
  }
  if (bac_bray_df_long$EstSalt.x[i] == "Delaware_Freshwater" &
      bac_bray_df_long$EstSalt.y[i] == "Waccamaw_Freshwater") {
    bac_bray_df_long$comparison[i] <- "Freshwater, Between Sites"
  }
  if (bac_bray_df_long$EstSalt.x[i] == "SF_Freshwater" &
      bac_bray_df_long$EstSalt.y[i] == "SF_Oligohaline") {
    bac_bray_df_long$comparison[i] <- "SF, Freshwater/Oligohaline"
  }
  if (bac_bray_df_long$EstSalt.x[i] == "SF_Freshwater" &
      bac_bray_df_long$EstSalt.y[i] == "SF_Mesohaline") {
    bac_bray_df_long$comparison[i] <- "SF, Freshwater/Mesohaline"
  }
  if (bac_bray_df_long$EstSalt.x[i] == "SF_Freshwater" &
      bac_bray_df_long$EstSalt.y[i] == "SF_Polyhaline") {
    bac_bray_df_long$comparison[i] <- "SF, Freshwater/Polyhaline"
  }
}
# Make the new column a factor
bac_bray_df_long$comparison <- as.factor(bac_bray_df_long$comparison)
# Filter unwanted comparisons
bac_bray_df_long <- subset(bac_bray_df_long, comparison != "NA")
# Check the sample sizes
table(bac_bray_df_long$comparison)
# ANOVA
m <- aov(value ~ comparison, data = bac_bray_df_long)
summary(m)
t <- emmeans(object = m, specs = "comparison") %>%
  cld(object = ., adjust = "Tukey", Letters = letters, alpha = 0.05) %>%
  mutate(name = "value",
         y = max(bac_bray_df_long$value)+(max(bac_bray_df_long$value)-min(bac_bray_df_long$value))/20)
bc_plot <- ggplot(data = bac_bray_df_long, aes(comparison, value, colour = comparison)) +
  geom_jitter(size = 0.75, alpha = 0.2, width = 0.3) +
  geom_boxplot(outlier.shape = NA, colour = "black", fill = NA) +
  geom_text(data = t, aes(comparison, y, label = str_trim(.group)), 
            size = 4, color = "black") +
  labs(x = "Site Comparison",
       y = "Bray-Curtis Dissimilarity") +
  scale_color_brewer(palette = "Set2") +
  ylim(0.55, 1.05) +
  theme_bw() +
  theme(legend.position = "none",
        axis.title = element_text(face="bold", size = 12),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
bc_plot

# Interactive
ggplotly(bc_plot)

# With density plot
detach(package:dendextend, unload = TRUE)
yplot <- ggdensity(bac_bray_df_long, "value", fill = "comparison", 
                   palette = "Set2", size = 0.25) +
  rotate() + 
  clean_theme() + 
  rremove("legend") + 
  xlim(0.55, 1.05)
yplot
plot_final <- insert_yaxis_grob(bc_plot, yplot, position = "right")
pdf("InitialFigs/Comb_Control_BC.pdf", width = 6, height = 5)
ggdraw(plot_final)
dev.off()

#### _Taxa ####
field_phyla <- summarize_taxonomy(field, level = 2, report_higher_tax = F)
plot_ts_heatmap(field_phyla, field$map_loaded, 0.01, 'EstSalt', rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
field_barsP <- plot_taxa_bars(field_phyla, field$map_loaded, "EstSalt", 
                              num_taxa = 12, data_only = T) %>%
  mutate(taxon = fct_rev(taxon))
pdf("InitialFigs/Comb_Control_Phyla.pdf", width = 7, height = 5)
ggplot(field_barsP, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Sample", y = "Relative abundance", fill = "Phylum") +
  scale_fill_manual(values = c("grey90", brewer.pal(12, "Paired")[12:1])) +
  scale_y_continuous(expand = c(0.01, 0.01)) + 
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
dev.off()
taxa_summary_by_sample_type(field_phyla, field$map_loaded, 'EstSalt', 0.01, 'KW')

# Show all samples
field$map_loaded$sampleID <- rownames(field$map_loaded)
bars <- plot_taxa_bars(field_phyla,
                       field$map_loaded,
                       "sampleID",
                       num_taxa = 12,
                       data_only = TRUE) %>%
  mutate(taxon = fct_rev(taxon)) %>%
  left_join(., field$map_loaded, by = c("group_by" = "sampleID"))
facet_names <- c("Freshwater" = "Freshwater", 
                 "Oligohaline" = "Oligohaline", 
                 "Mesohaline" = "Mesohaline",
                 "Polyhaline" = "Polyhaline",
                 "Waccamaw" = "SC",
                 "Alligator" = "NC", 
                 "Delaware" = "DE",
                 "SF" = "SF")
pdf("InitialFigs/Comb_Control_Phyla_samples.pdf", width = 12, height = 6)
ggplot(bars, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Sample", y = "Relative abundance", fill = "Phylum") +
  scale_fill_manual(values = c("grey90", brewer.pal(12, "Paired")[12:1])) +
  scale_y_continuous(expand = c(0.01, 0.01)) + 
  facet_nested(~ Salt + Estuary, space = "free", scales = "free_x",
               labeller = as_labeller(facet_names)) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 2, angle = 90, hjust = 1, vjust = 0.5),
        strip.text = element_text(size = 6),
        strip.background = element_rect(size = 0.2),
        axis.line.y = element_blank())
dev.off()

field_class <- summarize_taxonomy(field, level = 3, report_higher_tax = F)
plot_ts_heatmap(field_class, field$map_loaded, 0.01, 'EstSalt', rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
field_barsC <- plot_taxa_bars(field_class, field$map_loaded, "EstSalt", 
                              num_taxa = 12, data_only = T) %>%
  mutate(taxon = gsub("NA", "Unclassified", taxon)) %>%
  mutate(taxon = fct_relevel(taxon, "Other", after = Inf)) %>%
  mutate(taxon = fct_relevel(taxon, "Unclassified", after = Inf)) %>%
  mutate(taxon = fct_rev(taxon))
pdf("InitialFigs/Comb_Control_Class.pdf", width = 7, height = 5)
ggplot(field_barsC, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = NULL, y = "Relative abundance", fill = "Class") +
  scale_fill_manual(values = c("grey75", "grey90", brewer.pal(12, "Paired")[11:1])) +
  scale_y_continuous(expand = c(0.01, 0.01)) + 
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
dev.off()
taxa_summary_by_sample_type(field_class, field$map_loaded, 'EstSalt', 0.01, 'KW')

field_order <- summarize_taxonomy(field, level = 4, report_higher_tax = F)
plot_ts_heatmap(field_order, field$map_loaded, 0.01, 'EstSalt', rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
field_barsO <- plot_taxa_bars(field_order, field$map_loaded, "EstSalt", 
                              num_taxa = 12, data_only = T) %>%
  mutate(taxon = gsub("NA", "Unclassified", taxon)) %>%
  mutate(taxon = fct_relevel(taxon, "Other", after = Inf)) %>%
  mutate(taxon = fct_relevel(taxon, "Unclassified", after = Inf)) %>%
  mutate(taxon = fct_rev(taxon))
pdf("InitialFigs/Comb_Control_Order.pdf", width = 7, height = 5)
ggplot(field_barsO, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = NULL, y = "Relative abundance", fill = "Order") +
  scale_fill_manual(values = c("grey75", "grey90", brewer.pal(12, "Paired")[11:1])) +
  scale_y_continuous(expand = c(0.01, 0.01)) + 
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
dev.off()
taxa_summary_by_sample_type(field_order, field$map_loaded, 'EstSalt', 0.01, 'KW')

field_family <- summarize_taxonomy(field, level = 5, report_higher_tax = F)
plot_ts_heatmap(field_family, field$map_loaded, 0.01, 'EstSalt', rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
field_barsF <- plot_taxa_bars(field_family, field$map_loaded, "EstSalt", 
                              num_taxa = 12, data_only = T) %>%
  mutate(taxon = gsub("NA", "Unclassified", taxon)) %>%
  mutate(taxon = fct_relevel(taxon, "Other", after = Inf)) %>%
  mutate(taxon = fct_relevel(taxon, "Unclassified", after = Inf)) %>%
  mutate(taxon = fct_rev(taxon))
pdf("InitialFigs/Comb_Control_Family.pdf", width = 7, height = 5)
ggplot(field_barsF, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = NULL, y = "Relative abundance", fill = "Family") +
  scale_fill_manual(values = c("grey75", "grey90", brewer.pal(12, "Paired")[11:1])) +
  scale_y_continuous(expand = c(0.01, 0.01)) + 
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
dev.off()
taxa_summary_by_sample_type(field_family, field$map_loaded, 'EstSalt', 0.01, 'KW')

field_genus <- summarize_taxonomy(field, level = 6, report_higher_tax = F)
plot_ts_heatmap(field_genus, field$map_loaded, 0.005, 'EstSalt', rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
field_barsG <- plot_taxa_bars(field_genus, field$map_loaded, "EstSalt", 
                              num_taxa = 12, data_only = T) %>%
  mutate(taxon = gsub("NA", "Unclassified", taxon)) %>%
  mutate(taxon = fct_relevel(taxon, "Other", after = Inf)) %>%
  mutate(taxon = fct_relevel(taxon, "Unclassified", after = Inf)) %>%
  mutate(taxon = fct_rev(taxon))
pdf("InitialFigs/Comb_Control_Genus.pdf", width = 7, height = 5)
ggplot(field_barsG, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = NULL, y = "Relative abundance", fill = "Genus") +
  scale_fill_manual(values = c("grey75", "grey90", brewer.pal(12, "Paired")[11:1])) +
  scale_y_continuous(expand = c(0.01, 0.01)) + 
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
dev.off()
taxa_summary_by_sample_type(field_genus, field$map_loaded, 'EstSalt', 0.01, 'KW')

field_guilds <- summarize_taxonomy(field, level = 9, report_higher_tax = F)
plot_ts_heatmap(field_guilds, field$map_loaded, 0, 'EstSalt', rev_taxa = T, remove_other = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
field_barsGu <- plot_taxa_bars(field_guilds, field$map_loaded, "EstSalt",
                               num_taxa = 20, data_only = T) %>%
  filter(taxon != "NA") %>%
  droplevels() %>%
  mutate(taxon = factor(taxon,
                        levels = Guild_cols$Guild))
tallest_bar <- field_barsGu %>%
  group_by(group_by) %>%
  summarise(sum = sum(mean_value))
pdf("InitialFigs/Comb_Control_Guilds.pdf", width = 7, height = 6)
ggplot(field_barsGu, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Estuary", y = "Relative abundance", fill = "Guild") +
  scale_fill_manual(values = Guild_cols$color) +
  scale_y_continuous(expand = c(max(tallest_bar$sum)/100, max(tallest_bar$sum)/100)) + 
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
dev.off()
taxa_summary_by_sample_type(field_guilds, field$map_loaded, 'EstSalt', 0.01, 'KW')

# Guilds by sample (to show variation)
field$map_loaded$sampleID <- rownames(field$map_loaded)
field_barsGu <- plot_taxa_bars(field_guilds,
                               field$map_loaded,
                               "sampleID",
                               num_taxa = 20,
                               data_only = TRUE) %>%
  filter(taxon != "NA") %>%
  droplevels() %>%
  mutate(taxon = factor(taxon,
                        levels = Guild_cols$Guild)) %>%
  left_join(., field$map_loaded, by = c("group_by" = "sampleID"))
facet_names <- c("Freshwater" = "Freshwater", 
                 "Oligohaline" = "Oligohaline", 
                 "Mesohaline" = "Mesohaline",
                 "Polyhaline" = "Polyhaline",
                 "Waccamaw" = "SC",
                 "Alligator" = "NC", 
                 "Delaware" = "DE",
                 "SF" = "SF")
tallest_bar <- field_barsGu %>%
  group_by(group_by) %>%
  summarise(sum = sum(mean_value))
pdf("InitialFigs/Comb_Control_Guilds_samples.pdf", width = 12, height = 6)
ggplot(field_barsGu, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Sample", y = "Relative abundance", fill = "Guild") +
  scale_fill_manual(values = Guild_cols$color) +
  scale_y_continuous(expand = c(max(tallest_bar$sum)/100, max(tallest_bar$sum)/100)) + 
  facet_nested(~ Salt + Estuary, space = "free", scales = "free_x",
               labeller = as_labeller(facet_names)) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 2, angle = 90, hjust = 1, vjust = 0.5),
        strip.text = element_text(size = 6),
        strip.background = element_rect(size = 0.2),
        axis.line.y = element_blank())
dev.off()

#### _Methano ####
# Look at methanogens (careful to do this accurately!)
tax_sum_families_meth <- summarize_taxonomy(field, level = 5, 
                                            report_higher_tax = F)
tax_sum_families_meth <- tax_sum_families_meth[grep("(Methano|Methermicoccaceae|
                                                    Syntrophoarchaeaceae)", 
                                                    rownames(tax_sum_families_meth)),]
tax_sum_families_meth <- tax_sum_families_meth[!grepl("Methanoperedenaceae", 
                                                      rownames(tax_sum_families_meth)),]
field_barsMethano <- plot_taxa_bars(tax_sum_families_meth, 
                                    field$map_loaded, 
                                    "EstSalt",
                                    num_taxa = nrow(tax_sum_families_meth), 
                                    data_only = T)
nb.cols <- nrow(tax_sum_families_meth)
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)
tallest_bar <- field_barsMethano %>%
  group_by(group_by) %>%
  summarise(sum = sum(mean_value))
pdf("InitialFigs/Comb_Control_Methano.pdf", width = 7, height = 5)
ggplot(field_barsMethano, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = NULL, y = "Relative abundance", fill = "Family") +
  scale_fill_manual(values = mycolors) +
  scale_y_continuous(expand = c(max(tallest_bar$sum)/100, max(tallest_bar$sum)/100)) + 
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
dev.off()

#### _Desulfo ####
# Subset taxa to SRB and SRB_syn guilds
field_srb <- filter_taxa_from_input(field,
                                    taxa_to_keep = c("SRB", "SRB_syn"),
                                    at_spec_level = 9)
desulfo_wTax <- summarize_taxonomy(field_srb, level = 3, 
                                   report_higher_tax = T, relative = FALSE)
rownames(desulfo_wTax) <- substring(rownames(desulfo_wTax), 11)
field_barsDesulfo <- plot_taxa_bars(desulfo_wTax, field$map_loaded, "EstSalt", 
                                    num_taxa = nrow(desulfo_wTax), data_only = T) %>%
  mutate(mean_value = mean_value/26429)
nb.cols <- nrow(desulfo_wTax)
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)
tallest_bar <- field_barsDesulfo %>%
  group_by(group_by) %>%
  summarise(sum = sum(mean_value))
pdf("InitialFigs/Comb_Control_Desulfo.pdf", width = 7, height = 7)
ggplot(field_barsDesulfo, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = NULL, y = "Relative abundance", fill = "Class") +
  scale_fill_manual(values = mycolors) +
  scale_y_continuous(expand = c(max(tallest_bar$sum)/100, max(tallest_bar$sum)/100)) + 
  guides(fill = guide_legend(ncol = 1)) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        legend.key.size = unit(0.25, "cm"),
        legend.text = element_text(size = 10))
dev.off()

#### _Simper ####
field_sim <- simper(t(field$data_loaded), 
                    field$map_loaded$Salt)
field_s <- summary(field_sim)
head(field_s$Freshwater_Oligohaline)
head(field_s$Freshwater_Mesohaline)
head(field_s$Freshwater_Polyhaline)
head(field_s$Oligohaline_Mesohaline)
head(field_s$Oligohaline_Polyhaline)
head(field_s$Mesohaline_Polyhaline)
field_df1 <- head(field_s$Freshwater_Oligohaline, n = 20) %>%
  mutate(ASV = rownames(.),
         Comparison = "Freshwater/Oligohaline")
field_df2 <- head(field_s$Freshwater_Mesohaline, n = 20) %>%
  mutate(ASV = rownames(.),
         Comparison = "Freshwater/Mesohaline")
field_df3 <- head(field_s$Freshwater_Polyhaline, n = 20) %>%
  mutate(ASV = rownames(.),
         Comparison = "Freshwater/Polyhaline")
field_simper_results <- rbind(field_df1, field_df2, field_df3) %>%
  mutate(ASV = rownames(.),
         Comparison = "Freshwater/Mesohaline") %>%
  left_join(., field$taxonomy_loaded, by = c("ASV" = "taxonomy8")) %>%
  mutate(SaltResponse = ifelse(avb > ava, "Increase", "Decrease")) %>%
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
write_xlsx(field_simper_results, 
           "simper_results_comb_control.xlsx",
           format_headers = F)

#### _Multipatt ####
set.seed(1202)
field_mp <- multipatt(t(field$data_loaded), 
                      field$map_loaded$Salt, 
                      func = "r.g", 
                      control = how(nperm=999),
                      max.order = 1)
field_mp_results <- field_mp$sign %>%
  mutate(q.value = qvalue(field_mp$sign$p.value)$qvalues,
         Group = "NA",
         ASV = rownames(field_mp$sign)) %>%
  filter(q.value < 0.05) %>%
  filter(stat >= 0.65) %>%
  left_join(., field$taxonomy_loaded, by = c("ASV" = "taxonomy8"))
for (i in 1:nrow(field_mp_results)) {
  if (field_mp_results$s.Freshwater[i] == 1) {
    field_mp_results$Group[i] <- "Freshwater"
  }
}
for (i in 1:nrow(field_mp_results)) {
  if (field_mp_results$s.Oligohaline[i] == 1) {
    field_mp_results$Group[i] <- "Oligohaline"
  }
}
for (i in 1:nrow(field_mp_results)) {
  if (field_mp_results$s.Mesohaline[i] == 1) {
    field_mp_results$Group[i] <- "Mesohaline"
  }
}
for (i in 1:nrow(field_mp_results)) {
  if (field_mp_results$s.Polyhaline[i] == 1) {
    field_mp_results$Group[i] <- "Polyhaline"
  }
}
table(field_mp_results$Group)
field_asv <- summarize_taxonomy(field, level = 8, report_higher_tax = F)
field_asv_all <- data.frame("RelAbundance" = round(rowMeans(field_asv) * 100, digits = 5)) %>%
  mutate(ASV = rownames(.))
field_mp_corrs <- as.data.frame(field_mp$str) %>%
  dplyr::select(1:length(levels(field$map_loaded$Salt))) %>%
  mutate(ASV = rownames(.)) %>%
  filter(ASV %in% field_mp_results$ASV) %>%
  set_names(c("Freshwater", "Oligohaline", "Mesolhaline", "Polyhaline", "ASV"))
# Add corrs and taxonomy
field_mp_results <- field_mp_results %>%
  left_join(., field_asv_all, by = "ASV") %>%
  left_join(., field_mp_corrs, by = "ASV") %>%
  dplyr::select(taxonomy2, taxonomy6, ASV, Group, names(field_mp_corrs)[1:4], "RelAbundance") %>%
  mutate(ASV = gsub("ASV", "OTU", ASV)) %>%
  arrange(Group, taxonomy2, taxonomy6) %>%
  mutate(ASV = factor(ASV, levels = ASV)) %>%
  set_names(c("Phylum", "Genus", "ASV", "Group", names(field_mp_corrs)[1:4], "RelAbundance")) %>%
  mutate(taxon = paste(Phylum, Genus, ASV, sep = "; ")) %>%
  mutate(taxon = factor(taxon, levels = taxon))
# Plot heatmap with abundance barplot
field.hm.melted <- field_mp_results %>%
  dplyr::select(taxon, names(field_mp_corrs)[1:4]) %>%
  melt(., id.vars = c("taxon"))
field.hm <- ggplot(data = field.hm.melted, 
                   aes(x = factor(taxon), y = variable, fill = value)) + 
  geom_tile() + 
  scale_fill_distiller(name = "Indicator correlation index", palette = "RdBu", direction = -1, na.value = "transparent", type = "div", limits = c(-1, 1)) +
  scale_x_discrete(breaks = unique(field.hm.melted$taxon), labels = unique(field.hm.melted$taxon),
                   limits = rev(levels(field.hm.melted$taxon))) +
  scale_y_discrete(position = "right") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        legend.title = element_text(size = 10), legend.key.height = unit(0.5,"cm"),
        legend.key.width = unit(1, "cm"), legend.text = element_text(size = 7)) +
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5)) +
  coord_flip()
field.l <- get_legend(field.hm)
field.hm.clean <- field.hm +
  theme(axis.title.y = element_blank(), axis.text.y = element_text(margin = margin(c(0,-5,0,0)),
                                                                   size = 6),
        axis.ticks.y = element_blank(), axis.ticks.x = element_blank(), 
        panel.grid.major = element_blank(), legend.position="none",
        axis.text.x.top = element_text(size = 6, margin = margin(c(0,0,-2,0))), 
        plot.margin = margin(c(0,-2,0,5))) +
  labs(y = "Group")
field.bp.y <- ggplot(data = field_mp_results, aes(x = taxon, y = RelAbundance)) + 
  geom_bar(stat = "identity", fill = "grey") + 
  scale_y_continuous(position = "right") +
  scale_x_discrete(limits = rev(levels(field_mp_results$taxon))) +
  coord_flip() + 
  theme_minimal() +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position="none", plot.margin = margin(c(0,5,0,-5)),
        axis.text.x.top = element_text(size = 6, margin = margin(c(0,0,-2,0)))) +
  labs(y = "Rel. abund. (%)")
pdf("InitialFigs/Comb_Control_Multipatt.pdf", width = 8, height = 6)
plot_grid(field.hm.clean, field.bp.y, field.l, NULL, nrow = 2, ncol = 2, 
          rel_widths = c(10,2), rel_heights = c(12, 2), align = "hv", axis = "b")
dev.off()

# Problem is that at OTU level, many OTUs are too rare, not even 1%. So let's see what indicators there are at the genus level.
# First get aggregate taxonomy table
field_tax <- field$taxonomy_loaded %>%
  group_by(taxonomy6) %>%
  slice_head(n = 1) %>%
  dplyr::select(-taxonomy7, -taxonomy8)
field_genus <- summarize_taxonomy(field, level = 6, relative = F, report_higher_tax = F)
set.seed(1202)
field_mp <- multipatt(t(field_genus), 
                      field$map_loaded$Salt, 
                      func = "r.g", 
                      control = how(nperm=999),
                      max.order = 1)
field_mp_results <- field_mp$sign %>%
  mutate(q.value = qvalue(field_mp$sign$p.value)$qvalues,
         Group = "NA",
         Genus = rownames(field_mp$sign)) %>%
  filter(q.value < 0.05) %>%
  filter(stat >= 0.57) %>%
  left_join(., field_tax, by = c("Genus" = "taxonomy6"))
for (i in 1:nrow(field_mp_results)) {
  if (field_mp_results$s.Freshwater[i] == 1) {
    field_mp_results$Group[i] <- "Freshwater"
  }
}
for (i in 1:nrow(field_mp_results)) {
  if (field_mp_results$s.Oligohaline[i] == 1) {
    field_mp_results$Group[i] <- "Oligohaline"
  }
}
for (i in 1:nrow(field_mp_results)) {
  if (field_mp_results$s.Mesohaline[i] == 1) {
    field_mp_results$Group[i] <- "Mesohaline"
  }
}
for (i in 1:nrow(field_mp_results)) {
  if (field_mp_results$s.Polyhaline[i] == 1) {
    field_mp_results$Group[i] <- "Polyhaline"
  }
}
table(field_mp_results$Group)
field_genus_all <- data.frame("RelAbundance" = round(rowMeans(field_genus)/min(colSums(field$data_loaded)) * 100, digits = 4)) %>%
  mutate(Genus = rownames(.))
field_mp_corrs <- as.data.frame(field_mp$str) %>%
  dplyr::select(1:length(levels(field$map_loaded$Salt))) %>%
  mutate(Genus = rownames(.)) %>%
  filter(Genus %in% field_mp_results$Genus) %>%
  set_names(c("Freshwater", "Oligohaline", "Mesohaline", "Polyhaline", "Genus"))
# Add corrs and taxonomy
field_mp_results <- field_mp_results %>%
  left_join(., field_genus_all, by = "Genus") %>%
  left_join(., field_mp_corrs, by = "Genus") %>%
  dplyr::select(taxonomy2, Genus, Group, names(field_mp_corrs)[1:4], "RelAbundance") %>%
  arrange(Group, taxonomy2) %>%
  mutate(Genus = factor(Genus, levels = Genus)) %>%
  set_names(c("Phylum", "Genus", "Group", names(field_mp_corrs)[1:4], "RelAbundance")) %>%
  mutate(taxon = paste(Phylum, Genus, sep = "; ")) %>%
  mutate(taxon = factor(taxon, levels = taxon))
# Plot heatmap with abundance barplot
field.hm.melted <- field_mp_results %>%
  dplyr::select(taxon, names(field_mp_corrs)[1:4]) %>%
  melt(., id.vars = c("taxon"))
field.hm <- ggplot(data = field.hm.melted, 
                   aes(x = factor(taxon), y = variable, fill = value)) + 
  geom_tile() + 
  scale_fill_distiller(name = "Indicator correlation index", palette = "RdBu", direction = -1, na.value = "transparent", type = "div", limits = c(-1, 1)) +
  scale_x_discrete(breaks = unique(field.hm.melted$taxon), labels = unique(field.hm.melted$taxon),
                   limits = rev(levels(field.hm.melted$taxon))) +
  scale_y_discrete(position = "right") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        legend.title = element_text(size = 10), legend.key.height = unit(0.5,"cm"),
        legend.key.width = unit(1, "cm"), legend.text = element_text(size = 7)) +
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5)) +
  coord_flip()
field.l <- get_legend(field.hm)
field.hm.clean <- field.hm +
  theme(axis.title.y = element_blank(), axis.text.y = element_text(margin = margin(c(0,-5,0,0)),
                                                                   size = 6),
        axis.ticks.y = element_blank(), axis.ticks.x = element_blank(), 
        panel.grid.major = element_blank(), legend.position="none",
        axis.text.x.top = element_text(size = 6, margin = margin(c(0,0,-2,0))), 
        plot.margin = margin(c(0,-2,0,5))) +
  labs(y = "Group")
field.bp.y <- ggplot(data = field_mp_results, aes(x = taxon, y = RelAbundance)) + 
  geom_bar(stat = "identity", fill = "grey") + 
  scale_y_continuous(position = "right") +
  scale_x_discrete(limits = rev(levels(field_mp_results$taxon))) +
  coord_flip() + 
  theme_minimal() +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position="none", plot.margin = margin(c(0,5,0,-5)),
        axis.text.x.top = element_text(size = 6, margin = margin(c(0,0,-2,0)))) +
  labs(y = "Rel. abund. (%)")
pdf("InitialFigs/Comb_Control_Multipatt_genus.pdf", width = 8, height = 6)
plot_grid(field.hm.clean, field.bp.y, field.l, NULL, nrow = 2, ncol = 2, 
          rel_widths = c(10,2), rel_heights = c(12, 2), align = "hv", axis = "b")
dev.off()

#### _BGC ####
# Add biogeochem. analysis
# Get variables, make corrplot, envfit, ordination with vectors, ordistep, taxa correlations
# Note, porewater chemistry is only on D2 samples (5-15 cm depth)

# Get variables
env <- field$map_loaded %>%
  dplyr::select(CH4_ug_m2_h, CO2_ug_m2_h, Salinity_ppt_all, pH)
env_nona <- na.omit(env) # n = 169

# Corrplot
C <- cor(env_nona)
corrplot(C, method = "number", type = "lower", 
         tl.cex = 0.5, tl.col = "black", number.cex = 0.5)
pdf("InitialFigs/Comb_Control_BGC_corr.pdf", width = 7, height = 5)
meth_corr_by_bgc(env_nona = env_nona)
dev.off()

# Envfit
pcoa <- cmdscale(field_bc, k = nrow(field$map_loaded) - 1, eig = T)
set.seed(100)
ef <- envfit(pcoa, env, permutations = 999, na.rm = TRUE)
ef
arrow_factor <- ordiArrowMul(ef)
manual_factor <- 0.6
vec.df <- as.data.frame(ef$vectors$arrows*sqrt(ef$vectors$r)) %>%
  mutate(Dim1 = Dim1 * manual_factor,
         Dim2 = Dim2 * manual_factor) %>%
  mutate(variables = rownames(.)) %>%
  filter(ef$vectors$pvals < 0.05) %>%
  mutate(shortnames = c("CH4", "CO2", "Salinity", "pH"))

# Plot with significant vectors
ordiplot(pcoa)
plot(ef, p.max = 0.05, cex = 0.5)
pcoaA1 <- round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[1]*100, digits = 1)
pcoaA2 <- round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[2]*100, digits = 1)
field$map_loaded$Axis01 <- scores(pcoa)[,1]
field$map_loaded$Axis02 <- scores(pcoa)[,2]
micro.hulls <- ddply(field$map_loaded, c("Salt"), find_hull)
g <- ggplot(field$map_loaded, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = Salt, fill = Salt),
               alpha = 0.1, show.legend = F) +
  geom_point(size = 3, alpha = 0.75, aes(colour = Salt, shape = Estuary),
             show.legend = T) +
  labs(x = paste("PC1: ", pcoaA1, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2, "%", sep = ""),
       shape = "Depth (cm)",
       colour = "Salinity",
       fill = "Salinity") +
  scale_fill_viridis_d(guide = guide_legend(override.aes = list(alpha = 1))) +
  scale_colour_viridis_d(guide = guide_legend(override.aes = list(alpha = 1, shape = 15))) +
  guides(shape = guide_legend(override.aes = list(alpha = 1))) +
  theme_bw() +  
  theme(legend.position = "right",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10)) +
  geom_segment(data = vec.df,
               aes(x = 0, xend = Dim1, y = 0, yend = Dim2),
               arrow = arrow(length = unit(0.5, "cm")),
               colour = "gray", alpha = 0.5,
               inherit.aes = FALSE) + 
  geom_text(data = vec.df,
            aes(x = Dim1, y = Dim2, label = shortnames),
            size = 3, color = "black") +
  coord_fixed()
g
pdf("InitialFigs/Comb_Control_PCoA_wBGC.pdf", width = 7, height = 4)
g
dev.off()

# Ordistep
comm_nona <- as.data.frame(t(field$data_loaded)) %>%
  filter(rownames(.) %in% rownames(env_nona))
mod0 <- rda(comm_nona ~ 1, env_nona)  # Model with intercept only
mod1 <- rda(comm_nona ~ ., env_nona)  # Model with all explanatory variables
mod <- ordistep(mod0, scope = formula(mod1))
mod
mod$anova #

# Plot correlations for different taxonomic levels at a given % relative abundance threshold
pdf("InitialFigs/Comb_Control_Phyla_CH4.pdf", width = 7, height = 5)
meth_corr_by_taxonomy(input = field, level = 2, threshold = 0.5)
dev.off()
meth_corr_by_taxonomy(input = field, level = 3, threshold = 0.5)
meth_corr_by_taxonomy(input = field, level = 4, threshold = 0.5)
meth_corr_by_taxonomy(input = field, level = 5, threshold = 0.5)
meth_corr_by_taxonomy(input = field, level = 6, threshold = 0.5)
meth_corr_by_taxonomy(input = field, level = 8, threshold = 0.5) 
# ASV 5, ASV 6 neg
pdf("InitialFigs/Comb_Control_Guilds_CH4.pdf", width = 7, height = 5)
meth_corr_by_taxonomy(input = field, level = 9, threshold = 0)
dev.off()



#### 9. Comparison Lab Inc ####
# Only look at incubations
# Only look at control and + ASW
# Only look at final time point
# Delaware and North Carolina
# Pull out methanogens and sulfate reducers
#### _Setup ####
# input_filt_rare <- readRDS("input_filt_rare_comb.rds")
input_filt_rare <- readRDS("input_filt_rare_comb_wBGC.rds")
input_filt_rare$map_loaded <- input_filt_rare$map_loaded %>%
  mutate(Estuary = factor(Estuary,
                          levels = c("Waccamaw", "Alligator", "Delaware", "SF")),
         sampleID = rownames(.),
         Field = "NA",
         Salt = "NA",
         Depth = gsub("0-5", 0.05, Depth)) %>%
  mutate(Depth = gsub("5-15", 0.15, Depth)) %>%
  mutate(Depth = as.numeric(Depth)) %>%
  mutate(Depth = ifelse(Depth <= 0.05, "< 5 cm", "5 - 15 cm"))

# Add column for if unmanipulated field sample
for (i in 1:nrow(input_filt_rare$map_loaded)) {
  if (input_filt_rare$map_loaded$Estuary[i] == "SF") {
    input_filt_rare$map_loaded$Field[i] <- "Field"
  }
}
for (i in 1:nrow(input_filt_rare$map_loaded)) {
  if (input_filt_rare$map_loaded$Estuary[i] == "Delaware" & 
      input_filt_rare$map_loaded$Site[i] == "Soil") {
    input_filt_rare$map_loaded$Field[i] <- "Field"
  }
}
for (i in 1:nrow(input_filt_rare$map_loaded)) {
  if (input_filt_rare$map_loaded$Estuary[i] == "Alligator" & 
      input_filt_rare$map_loaded$Site[i] == "Soil") {
    input_filt_rare$map_loaded$Field[i] <- "Field"
  }
}
for (i in 1:nrow(input_filt_rare$map_loaded)) {
  if (input_filt_rare$map_loaded$Estuary[i] == "Waccamaw" & 
      input_filt_rare$map_loaded$Detail[i] == "Control") {
    input_filt_rare$map_loaded$Field[i] <- "Field"
  }
}
lab <- filter_data(input_filt_rare,
                   filter_cat = "Field",
                   keep_vals = "NA") # 141 samples
table(lab$map_loaded$Estuary)

# Now need to get only final time point and only incubations
# Note that there are 2 incubation experiments
# Both used 5 ppt ASW and went 12 weeks
lab <- filter_data(lab,
                   filter_cat = "Site",
                   keep_vals = c("Soil incubation", "Soil Incubation")) # 67 samples
table(lab$map_loaded$Estuary)

# Note all "Alligator" are 3 months, so get "12 wk" for Delaware
table(lab$map_loaded$Detail)

lab <- filter_data(lab,
                   filter_cat = "Detail",
                   keep_vals = c("5 ppt ASW 12wk",
                                 "Freshwater 12wk",
                                 "5ppt ASW",
                                 "DI_ctrl")) # 26 samples
table(lab$map_loaded$Estuary)

# Make Treatment Column
lab$map_loaded <- lab$map_loaded %>%
  mutate_if(is.character, as.factor) %>%
  mutate(Salt = recode_factor(Detail,
                              "Freshwater 12wk" = "Control",
                              "DI_ctrl" = "Control",
                              "5 ppt ASW 12wk" = "+ASW",
                              "5ppt ASW" = "+ASW"))

# Need combined estuary/salt factor
lab$map_loaded$EstSalt <- paste(lab$map_loaded$Estuary,
                                lab$map_loaded$Salt,
                                sep = "_")

#### _Alpha ####
leveneTest(lab$map_loaded$rich ~ lab$map_loaded$Salt) # Homogeneous
m <- aov(rich ~ Estuary + Salt + Depth, data = lab$map_loaded)
Anova(m, type = "II") # Estuary sig, salt and depth not sig.
shapiro.test(m$residuals) # Almost normal

leveneTest(lab$map_loaded$shannon ~ lab$map_loaded$Salt) # Homogenoeus
m1 <- aov(shannon ~ Estuary + Salt + Depth, data = lab$map_loaded)
Anova(m1, type = "II") # Estuary sig, salt and depth not sig.
shapiro.test(m1$residuals) # Not normal
facet_df <- c("rich" = "(a) Richness",
              "shannon" = "(b) Shannon")
alpha_long <- lab$map_loaded %>%
  pivot_longer(cols = c("rich", "shannon"))
pdf("InitialFigs/Comb_Lab_Alpha.pdf", width = 6, height = 3)
ggplot(alpha_long, aes(reorder(Salt, value, mean), value)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 3, alpha = 1, width = 0.2, aes(shape = Depth, colour = Estuary)) +
  labs(x = "Site", y = "Number of OTUs") +
  scale_colour_viridis_d() +
  facet_wrap(~ name, ncol = 2, scales = "free_y", labeller = as_labeller(facet_df)) +
  theme_bw() +
  theme(legend.position = "right",
        axis.title = element_blank(),
        axis.text = element_text(size = 10),
        strip.text = element_text(size = 10))
dev.off()



#### _Beta ####
lab_bc <- calc_dm(lab$data_loaded)
set.seed(1150)
adonis2(lab_bc ~ lab$map_loaded$Estuary+lab$map_loaded$Salt+lab$map_loaded$Depth) # All
adonis2(lab_bc ~ lab$map_loaded$Depth+lab$map_loaded$Salt+lab$map_loaded$Estuary)
anova(betadisper(lab_bc, lab$map_loaded$Estuary)) # Dispersion homogeneous
anova(betadisper(lab_bc, lab$map_loaded$Salt)) # Dispersion homogeneous
anova(betadisper(lab_bc, lab$map_loaded$Depth)) # Dispersion homogeneous
lab_pcoa <- cmdscale(lab_bc, k = nrow(lab$map_loaded) - 1, eig = T)
pcoaA1 <- round((eigenvals(lab_pcoa)/sum(eigenvals(lab_pcoa)))[1]*100, digits = 1)
pcoaA2 <- round((eigenvals(lab_pcoa)/sum(eigenvals(lab_pcoa)))[2]*100, digits = 1)
lab$map_loaded$Axis01 <- scores(lab_pcoa)[,1]
lab$map_loaded$Axis02 <- scores(lab_pcoa)[,2]
micro.hulls <- ddply(lab$map_loaded, c("Salt"), find_hull)
p_leg <- ggplot(lab$map_loaded, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = Salt, fill = Salt),
               alpha = 0.1, show.legend = F) +
  geom_point(size = 3, alpha = 1, aes(colour = Salt, shape = Estuary),
             show.legend = T) +
  labs(x = paste("PC1: ", pcoaA1, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2, "%", sep = ""),
       colour = "Salt") +
  scale_fill_manual(values = c("blue", "red")) +
  scale_colour_manual(values = c("blue", "red")) +
  guides(colour = guide_legend(override.aes = list(alpha = 1,
                                                   shape = 15)),
         fill = "none") +
  theme_bw() +  
  theme(legend.position = "right",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10))
p_leg
leg <- get_legend(p_leg)
micro.hulls <- ddply(lab$map_loaded, c("EstSalt"), find_hull)
p_noleg <- ggplot(lab$map_loaded, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = EstSalt, fill = EstSalt),
               alpha = 0.1, show.legend = F) +
  geom_point(size = 3, alpha = 1, aes(colour = EstSalt, shape = Estuary),
             show.legend = T) +
  labs(x = paste("PC1: ", pcoaA1, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2, "%", sep = ""),
       colour = "Salt") +
  scale_fill_manual(values = c("red", "blue", "red", "blue")) +
  scale_colour_manual(values = c("red", "blue", "red", "blue")) +
  theme_bw() +  
  theme(legend.position = "none",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10))
p_noleg

pdf("InitialFigs/Comb_Lab_PCoA.pdf", width = 7, height = 5)
plot_grid(p_noleg, leg, rel_widths = c(5,1))
dev.off()


#### _Taxa ####
lab_phyla <- summarize_taxonomy(lab, level = 2, report_higher_tax = F)
plot_ts_heatmap(lab_phyla, lab$map_loaded, 0.01, 'EstSalt', rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
lab_barsP <- plot_taxa_bars(lab_phyla, lab$map_loaded, "EstSalt", 
                              num_taxa = 12, data_only = T) %>%
  mutate(taxon = fct_rev(taxon)) %>%
  separate(group_by, into = c("Estuary", "Salt"), sep = "_") %>%
  mutate(Salt = fct_rev(Salt))
pdf("InitialFigs/Comb_Lab_Phyla.pdf", width = 7, height = 5)
ggplot(lab_barsP, aes(Salt, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Sample", y = "Relative abundance", fill = "Phylum") +
  scale_fill_manual(values = c("grey90", brewer.pal(12, "Paired")[12:1])) +
  scale_y_continuous(expand = c(0.01, 0.01)) + 
  facet_wrap(~ Estuary) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        strip.background = element_rect(size = 0.2),
        axis.line.y = element_blank())
dev.off()
taxa_summary_by_sample_type(lab_phyla, lab$map_loaded, 'EstSalt', 0.01, 'KW')

# Show all samples
lab$map_loaded$sampleID <- rownames(lab$map_loaded)
lab_barsP <- plot_taxa_bars(lab_phyla,
                            lab$map_loaded,
                            "sampleID",
                            num_taxa = 12,
                            data_only = TRUE) %>%
  mutate(taxon = fct_rev(taxon)) %>%
  left_join(., lab$map_loaded, by = c("group_by" = "sampleID"))
pdf("InitialFigs/Comb_Lab_Phyla_samples.pdf", width = 8, height = 6)
ggplot(lab_barsP, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Sample", y = "Relative abundance", fill = "Phylum") +
  scale_fill_manual(values = c("grey90", brewer.pal(12, "Paired")[12:1])) +
  scale_y_continuous(expand = c(0.01, 0.01)) + 
  facet_nested(~ Estuary + Salt, space = "free", scales = "free_x") +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 4, angle = 90, hjust = 1, vjust = 0.5),
        strip.text = element_text(size = 6),
        strip.background = element_rect(size = 0.2),
        axis.line.y = element_blank())
dev.off()

lab_class <- summarize_taxonomy(lab, level = 3, report_higher_tax = F)
plot_ts_heatmap(lab_class, lab$map_loaded, 0.01, 'EstSalt', rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
lab_barsC <- plot_taxa_bars(lab_class, lab$map_loaded, "EstSalt", 
                              num_taxa = 12, data_only = T) %>%
  mutate(taxon = gsub("NA", "Unclassified", taxon)) %>%
  mutate(taxon = fct_relevel(taxon, "Other", after = Inf)) %>%
  mutate(taxon = fct_relevel(taxon, "Unclassified", after = Inf)) %>%
  mutate(taxon = fct_rev(taxon)) %>%
  separate(group_by, into = c("Estuary", "Salt"), sep = "_") %>%
  mutate(Salt = fct_rev(Salt))
pdf("InitialFigs/Comb_Lab_Class.pdf", width = 7, height = 5)
ggplot(lab_barsC, aes(Salt, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = NULL, y = "Relative abundance", fill = "Class") +
  scale_fill_manual(values = c("grey75", "grey90", brewer.pal(12, "Paired")[11:1])) +
  scale_y_continuous(expand = c(0.01, 0.01)) + 
  facet_wrap(~ Estuary) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        strip.background = element_rect(size = 0.2),
        axis.line.y = element_blank())
dev.off()
taxa_summary_by_sample_type(lab_class, lab$map_loaded, 'EstSalt', 0.01, 'KW')

lab_order <- summarize_taxonomy(lab, level = 4, report_higher_tax = F)
plot_ts_heatmap(lab_order, lab$map_loaded, 0.01, 'EstSalt', rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
lab_barsO <- plot_taxa_bars(lab_order, lab$map_loaded, "EstSalt", 
                              num_taxa = 12, data_only = T) %>%
  mutate(taxon = gsub("NA", "Unclassified", taxon)) %>%
  mutate(taxon = fct_relevel(taxon, "Other", after = Inf)) %>%
  mutate(taxon = fct_relevel(taxon, "Unclassified", after = Inf)) %>%
  mutate(taxon = fct_rev(taxon)) %>%
  separate(group_by, into = c("Estuary", "Salt"), sep = "_") %>%
  mutate(Salt = fct_rev(Salt))
pdf("InitialFigs/Comb_Lab_Order.pdf", width = 7, height = 5)
ggplot(lab_barsO, aes(Salt, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = NULL, y = "Relative abundance", fill = "Order") +
  scale_fill_manual(values = c("grey75", "grey90", brewer.pal(12, "Paired")[11:1])) +
  scale_y_continuous(expand = c(0.01, 0.01)) + 
  facet_wrap(~ Estuary) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        strip.background = element_rect(size = 0.2),
        axis.line.y = element_blank())
dev.off()
taxa_summary_by_sample_type(lab_order, lab$map_loaded, 'EstSalt', 0.01, 'KW')

lab_family <- summarize_taxonomy(lab, level = 5, report_higher_tax = F)
plot_ts_heatmap(lab_family, lab$map_loaded, 0.01, 'EstSalt', rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
lab_barsF <- plot_taxa_bars(lab_family, lab$map_loaded, "EstSalt", 
                              num_taxa = 12, data_only = T) %>%
  mutate(taxon = gsub("NA", "Unclassified", taxon)) %>%
  mutate(taxon = fct_relevel(taxon, "Other", after = Inf)) %>%
  mutate(taxon = fct_relevel(taxon, "Unclassified", after = Inf)) %>%
  mutate(taxon = fct_rev(taxon)) %>%
  separate(group_by, into = c("Estuary", "Salt"), sep = "_") %>%
  mutate(Salt = fct_rev(Salt))
pdf("InitialFigs/Comb_Lab_Family.pdf", width = 7, height = 5)
ggplot(lab_barsF, aes(Salt, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = NULL, y = "Relative abundance", fill = "Family") +
  scale_fill_manual(values = c("grey75", "grey90", brewer.pal(12, "Paired")[11:1])) +
  scale_y_continuous(expand = c(0.01, 0.01)) + 
  facet_wrap(~ Estuary) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        strip.background = element_rect(size = 0.2),
        axis.line.y = element_blank())
dev.off()
taxa_summary_by_sample_type(lab_family, lab$map_loaded, 'EstSalt', 0.01, 'KW')

lab_genus <- summarize_taxonomy(lab, level = 6, report_higher_tax = F)
plot_ts_heatmap(lab_genus, lab$map_loaded, 0.005, 'EstSalt', rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
lab_barsG <- plot_taxa_bars(lab_genus, lab$map_loaded, "EstSalt", 
                              num_taxa = 12, data_only = T) %>%
  mutate(taxon = gsub("NA", "Unclassified", taxon)) %>%
  mutate(taxon = fct_relevel(taxon, "Other", after = Inf)) %>%
  mutate(taxon = fct_relevel(taxon, "Unclassified", after = Inf)) %>%
  mutate(taxon = fct_rev(taxon)) %>%
  separate(group_by, into = c("Estuary", "Salt"), sep = "_") %>%
  mutate(Salt = fct_rev(Salt))
pdf("InitialFigs/Comb_Lab_Genus.pdf", width = 7, height = 5)
ggplot(lab_barsG, aes(Salt, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = NULL, y = "Relative abundance", fill = "Genus") +
  scale_fill_manual(values = c("grey75", "grey90", brewer.pal(12, "Paired")[11:1])) +
  scale_y_continuous(expand = c(0.01, 0.01)) + 
  facet_wrap(~ Estuary) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        strip.background = element_rect(size = 0.2),
        axis.line.y = element_blank())
dev.off()
taxa_summary_by_sample_type(lab_genus, lab$map_loaded, 'EstSalt', 0.01, 'KW')

lab_guilds <- summarize_taxonomy(lab, level = 9, report_higher_tax = F)
plot_ts_heatmap(lab_guilds, lab$map_loaded, 0, 'EstSalt', rev_taxa = T, remove_other = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
lab_barsGu <- plot_taxa_bars(lab_guilds, lab$map_loaded, "EstSalt",
                               num_taxa = 20, data_only = T) %>%
  filter(taxon != "NA") %>%
  droplevels() %>%
  mutate(taxon = factor(taxon,
                        levels = Guild_cols$Guild)) %>%
  separate(group_by, into = c("Estuary", "Salt"), sep = "_") %>%
  mutate(Salt = fct_rev(Salt))
pdf("InitialFigs/Comb_Lab_Guilds.pdf", width = 7, height = 5)
ggplot(lab_barsGu, aes(Salt, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Estuary", y = "Relative abundance", fill = "Guild") +
  scale_fill_manual(values = Guild_cols$color) +
  scale_y_continuous(expand = c(max(tallest_bar$sum)/100, max(tallest_bar$sum)/100)) + 
  facet_wrap(~ Estuary) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        strip.background = element_rect(size = 0.2),
        axis.line.y = element_blank())
dev.off()
taxa_summary_by_sample_type(lab_guilds, lab$map_loaded, 'EstSalt', 0.01, 'KW')

# Show all samples
lab_barsGu <- plot_taxa_bars(lab_guilds,
                             lab$map_loaded,
                             "sampleID",
                             num_taxa = 20,
                             data_only = TRUE) %>%
  filter(taxon != "NA") %>%
  droplevels() %>%
  mutate(taxon = factor(taxon,
                        levels = Guild_cols$Guild)) %>%
  left_join(., lab$map_loaded, by = c("group_by" = "sampleID"))
tallest_bar <- lab_barsGu %>%
  group_by(group_by) %>%
  summarise(sum = sum(mean_value))
pdf("InitialFigs/Comb_Lab_Guilds_samples.pdf", width = 8, height = 6)
ggplot(lab_barsGu, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Sample", y = "Relative abundance", fill = "Guild") +
  scale_fill_manual(values = Guild_cols$color) +
  scale_y_continuous(expand = c(max(tallest_bar$sum)/100, max(tallest_bar$sum)/100)) + 
  facet_nested(~ Estuary + Salt, space = "free", scales = "free_x") +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 4, angle = 90, hjust = 1, vjust = 0.5),
        strip.text = element_text(size = 6),
        strip.background = element_rect(size = 0.2),
        axis.line.y = element_blank())
dev.off()

#### _Methano ####
# Look at methanogens (careful to do this accurately!)
tax_sum_families_meth <- summarize_taxonomy(lab, level = 5, 
                                            report_higher_tax = F)
tax_sum_families_meth <- tax_sum_families_meth[grep("(Methano|Methermicoccaceae|
                                                    Syntrophoarchaeaceae)", 
                                                    rownames(tax_sum_families_meth)),]
tax_sum_families_meth <- tax_sum_families_meth[!grepl("Methanoperedenaceae", 
                                                      rownames(tax_sum_families_meth)),]
lab_barsMethano <- plot_taxa_bars(tax_sum_families_meth, 
                                  lab$map_loaded, 
                                  "EstSalt",
                                  num_taxa = nrow(tax_sum_families_meth), 
                                  data_only = T) %>%
  separate(group_by, into = c("Estuary", "Salt"), sep = "_") %>%
  mutate(Salt = fct_rev(Salt))
nb.cols <- nrow(tax_sum_families_meth)
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)
tallest_bar <- lab_barsMethano %>%
  group_by(Estuary, Salt) %>%
  summarise(sum = sum(mean_value))
pdf("InitialFigs/Comb_Lab_Methano.pdf", width = 7, height = 5)
ggplot(lab_barsMethano, aes(Salt, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = NULL, y = "Relative abundance", fill = "Family") +
  scale_fill_manual(values = mycolors) +
  scale_y_continuous(expand = c(max(tallest_bar$sum)/100, max(tallest_bar$sum)/100)) + 
  facet_wrap(~ Estuary) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        strip.background = element_rect(size = 0.2),
        axis.line.y = element_blank())
dev.off()

# Show all samples
lab_barsMethano <- plot_taxa_bars(tax_sum_families_meth,
                                  lab$map_loaded,
                                  "sampleID",
                                  num_taxa = nrow(tax_sum_families_meth),
                                  data_only = TRUE) %>%
  left_join(., lab$map_loaded, by = c("group_by" = "sampleID"))
tallest_bar <- lab_barsMethano %>%
  group_by(group_by) %>%
  summarise(sum = sum(mean_value))
pdf("InitialFigs/Comb_Lab_Methano_samples.pdf", width = 8, height = 6)
ggplot(lab_barsMethano, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Sample", y = "Relative abundance", fill = "Family") +
  scale_fill_manual(values = mycolors) +
  scale_y_continuous(expand = c(max(tallest_bar$sum)/100, max(tallest_bar$sum)/100)) +  
  facet_nested(~ Estuary + Salt, space = "free", scales = "free_x") +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 4, angle = 90, hjust = 1, vjust = 0.5),
        strip.text = element_text(size = 6),
        strip.background = element_rect(size = 0.2),
        axis.line.y = element_blank())
dev.off()

# Delve into Methanobacteriaceae result
# Check methanogen family raw counts
tax_sum_families_meth <- summarize_taxonomy(lab, level = 5, 
                                            report_higher_tax = F,
                                            relative = F)
tax_sum_families_meth <- tax_sum_families_meth[grep("(Methano|Methermicoccaceae|
                                                    Syntrophoarchaeaceae)", 
                                                    rownames(tax_sum_families_meth)),]
tax_sum_families_meth <- tax_sum_families_meth[!grepl("Methanoperedenaceae", 
                                                      rownames(tax_sum_families_meth)),]

# Check Methanobacteriaceae genera
lab_Methanobacteriaceae <- filter_taxa_from_input(lab,
                                                  taxa_to_keep = c("Methanobacteriaceae"),
                                                  at_spec_level = 5)
Methanobacteriaceae_wTax <- summarize_taxonomy(lab_Methanobacteriaceae, 
                                               level = 6, 
                                               report_higher_tax = T, 
                                               relative = FALSE)
# There are two genera (Methanobacterium and Methanosphaera) and NA
rownames(Methanobacteriaceae_wTax) <- substring(rownames(Methanobacteriaceae_wTax), 62)
nb.cols <- nrow(Methanobacteriaceae_wTax)
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)
lab_barsMethanobacteriaceae <- plot_taxa_bars(Methanobacteriaceae_wTax, 
                                              lab$map_loaded, 
                                              "EstSalt", 
                                              num_taxa = nrow(Methanobacteriaceae_wTax), 
                                              data_only = T) %>%
  mutate(mean_value = mean_value/26429) %>%
  separate(group_by, into = c("Estuary", "Salt"), sep = "_") %>%
  mutate(Salt = fct_rev(Salt))
tallest_bar <- lab_barsMethanobacteriaceae %>%
  group_by(Estuary, Salt) %>%
  summarise(sum = sum(mean_value))
pdf("InitialFigs/Comb_Lab_Methanobacteriaceae_Genera.pdf", width = 7, height = 5)
ggplot(lab_barsMethanobacteriaceae, aes(Salt, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = NULL, y = "Relative abundance", fill = "Genus") +
  scale_fill_manual(values = mycolors) +
  scale_y_continuous(expand = c(max(tallest_bar$sum)/100, max(tallest_bar$sum)/100)) + 
  facet_wrap(~ Estuary) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        strip.background = element_rect(size = 0.2),
        axis.line.y = element_blank())
dev.off()

# Now also check OTU level
Methanobacteriaceae_wTax <- summarize_taxonomy(lab_Methanobacteriaceae, 
                                               level = 8, 
                                               report_higher_tax = T, 
                                               relative = FALSE)
rownames(Methanobacteriaceae_wTax) <- substring(rownames(Methanobacteriaceae_wTax), 62)
# Replace ASV with OTU in figure to avoid confusion
rownames(Methanobacteriaceae_wTax) <- str_replace(rownames(Methanobacteriaceae_wTax), "ASV", "OTU")
nb.cols <- nrow(Methanobacteriaceae_wTax)
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)
lab_barsMethanobacteriaceae <- plot_taxa_bars(Methanobacteriaceae_wTax, 
                                              lab$map_loaded, 
                                              "EstSalt", 
                                              num_taxa = nrow(Methanobacteriaceae_wTax), 
                                              data_only = T) %>%
  mutate(mean_value = mean_value/26429) %>%
  separate(group_by, into = c("Estuary", "Salt"), sep = "_") %>%
  mutate(Salt = fct_rev(Salt))
tallest_bar <- lab_barsMethanobacteriaceae %>%
  group_by(Estuary, Salt) %>%
  summarise(sum = sum(mean_value))
pdf("InitialFigs/Comb_Lab_Methanobacteriaceae_OTUs.pdf", width = 7, height = 5)
ggplot(lab_barsMethanobacteriaceae, aes(Salt, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = NULL, y = "Relative abundance", fill = "Family; Genus; Species; OTU ID") +
  scale_fill_manual(values = mycolors) +
  scale_y_continuous(expand = c(max(tallest_bar$sum)/100, max(tallest_bar$sum)/100)) + 
  facet_wrap(~ Estuary) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        strip.background = element_rect(size = 0.2),
        axis.line.y = element_blank(),
        legend.key.size = unit(0.5, "cm"),
        legend.text = element_text(size = 10))
dev.off()


#### _Desulfo ####
# Subset taxa to SRB and SRB_syn guilds
lab_srb <- filter_taxa_from_input(lab,
                                  taxa_to_keep = c("SRB", "SRB_syn"),
                                  at_spec_level = 9)
desulfo_wTax <- summarize_taxonomy(lab_srb, level = 3, 
                                   report_higher_tax = T, relative = FALSE)
rownames(desulfo_wTax) <- substring(rownames(desulfo_wTax), 11)
nb.cols <- nrow(desulfo_wTax)
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)
lab_barsDesulfo <- plot_taxa_bars(desulfo_wTax, lab$map_loaded, "EstSalt", 
                                  num_taxa = nrow(desulfo_wTax), data_only = T) %>%
  mutate(mean_value = mean_value/26429) %>%
  separate(group_by, into = c("Estuary", "Salt"), sep = "_") %>%
  mutate(Salt = fct_rev(Salt))
tallest_bar <- lab_barsDesulfo %>%
  group_by(Estuary, Salt) %>%
  summarise(sum = sum(mean_value))
pdf("InitialFigs/Comb_Lab_Desulfo.pdf", width = 7, height = 5)
ggplot(lab_barsDesulfo, aes(Salt, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = NULL, y = "Relative abundance", fill = "Class") +
  scale_fill_manual(values = mycolors) +
  scale_y_continuous(expand = c(max(tallest_bar$sum)/100, max(tallest_bar$sum)/100)) + 
  facet_wrap(~ Estuary) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        strip.background = element_rect(size = 0.2),
        axis.line.y = element_blank())
dev.off()

#### _Simper ####
lab_sim <- simper(t(lab$data_loaded), 
                    lab$map_loaded$Salt)
lab_s <- summary(lab_sim)
head(lab_s$`Control_+ASW`)
lab_df1 <- head(lab_s$`Control_+ASW`, n = 20) %>%
  mutate(ASV = rownames(.),
         Comparison = "Control/+ASW") %>%
  left_join(., lab$taxonomy_loaded, by = c("ASV" = "taxonomy8")) %>%
  mutate(SaltResponse = ifelse(avb > ava, "Increase", "Decrease")) %>%
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
write_xlsx(lab_df1, 
           "simper_results_comb_lab.xlsx",
           format_headers = F)

#### _Multipatt ####
set.seed(1202)
lab_mp <- multipatt(t(lab$data_loaded), 
                      lab$map_loaded$Salt, 
                      func = "r.g", 
                      control = how(nperm=999),
                      max.order = 1)
lab_mp_results <- lab_mp$sign %>%
  mutate(q.value = qvalue(lab_mp$sign$p.value)$qvalues,
         Group = "NA",
         ASV = rownames(lab_mp$sign)) %>%
  filter(p.value < 0.01) %>%
  filter(stat >= 0.5) %>%
  left_join(., lab$taxonomy_loaded, by = c("ASV" = "taxonomy8"))
for (i in 1:nrow(lab_mp_results)) {
  if (lab_mp_results$s.Control[i] == 1) {
    lab_mp_results$Group[i] <- "Control"
  }
}
for (i in 1:nrow(lab_mp_results)) {
  if (lab_mp_results$`s.+ASW`[i] == 1) {
    lab_mp_results$Group[i] <- "+ASW"
  }
}
table(lab_mp_results$Group)
lab_asv <- summarize_taxonomy(lab, level = 8, report_higher_tax = F)
lab_asv_all <- data.frame("RelAbundance" = round(rowMeans(lab_asv) * 100, digits = 5)) %>%
  mutate(ASV = rownames(.))
lab_mp_corrs <- as.data.frame(lab_mp$str) %>%
  dplyr::select(1:length(levels(lab$map_loaded$Salt))) %>%
  mutate(ASV = rownames(.)) %>%
  filter(ASV %in% lab_mp_results$ASV) %>%
  set_names(c("Control", "+ASW", "ASV"))
# Add corrs and taxonomy
lab_mp_results <- lab_mp_results %>%
  left_join(., lab_asv_all, by = "ASV") %>%
  left_join(., lab_mp_corrs, by = "ASV") %>%
  dplyr::select(taxonomy2, taxonomy6, ASV, Group, names(lab_mp_corrs)[1:2], "RelAbundance") %>%
  mutate(ASV = gsub("ASV", "OTU", ASV)) %>%
  arrange(Group, taxonomy2, taxonomy6) %>%
  mutate(ASV = factor(ASV, levels = ASV)) %>%
  set_names(c("Phylum", "Genus", "ASV", "Group", names(lab_mp_corrs)[1:2], "RelAbundance")) %>%
  mutate(taxon = paste(Phylum, Genus, ASV, sep = "; ")) %>%
  mutate(taxon = factor(taxon, levels = taxon))
# Plot heatmap with abundance barplot
lab.hm.melted <- lab_mp_results %>%
  dplyr::select(taxon, names(lab_mp_corrs)[1:2]) %>%
  melt(., id.vars = c("taxon"))
lab.hm <- ggplot(data = lab.hm.melted, 
                   aes(x = factor(taxon), y = variable, fill = value)) + 
  geom_tile() + 
  scale_fill_distiller(name = "Indicator correlation index", palette = "RdBu", direction = -1, na.value = "transparent", type = "div", limits = c(-1, 1)) +
  scale_x_discrete(breaks = unique(lab.hm.melted$taxon), labels = unique(lab.hm.melted$taxon),
                   limits = rev(levels(lab.hm.melted$taxon))) +
  scale_y_discrete(position = "right") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        legend.title = element_text(size = 10), legend.key.height = unit(0.5,"cm"),
        legend.key.width = unit(1, "cm"), legend.text = element_text(size = 7)) +
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5)) +
  coord_flip()
lab.l <- get_legend(lab.hm)
lab.hm.clean <- lab.hm +
  theme(axis.title.y = element_blank(), axis.text.y = element_text(margin = margin(c(0,-20,0,0)), 
                                                                   size = 2),
        axis.ticks.y = element_blank(), axis.ticks.x = element_blank(), 
        panel.grid.major = element_blank(), legend.position="none",
        axis.text.x.top = element_text(size = 6, margin = margin(c(0,0,-2,0))), 
        plot.margin = margin(c(0,-2,0,5))) +
  labs(y = "Group")
lab.bp.y <- ggplot(data = lab_mp_results, aes(x = taxon, y = RelAbundance)) + 
  geom_bar(stat = "identity", fill = "grey") + 
  scale_y_continuous(position = "right") +
  scale_x_discrete(limits = rev(levels(lab_mp_results$taxon))) +
  coord_flip() + 
  theme_minimal() +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position="none", plot.margin = margin(c(0,5,0,-5)),
        axis.text.x.top = element_text(size = 6, margin = margin(c(0,0,-2,0)))) +
  labs(y = "Rel. abund. (%)")
pdf("InitialFigs/Comb_Lab_Multipatt.pdf", width = 8, height = 10)
plot_grid(lab.hm.clean, lab.bp.y, lab.l, NULL, nrow = 2, ncol = 2, 
          rel_widths = c(10,2), rel_heights = c(12, 2), align = "hv", axis = "b")
dev.off()

#### _BGC ####
# Add biogeochem. analysis
# Get variables, make corrplot, envfit, ordination with vectors, ordistep, taxa correlations
# Note, porewater chemistry is only on D2 samples (5-15 cm depth)

# Get variables
env <- input_filt_rare$map_loaded %>%
  dplyr::select(CH4_ug_m2_h, CO2_ug_m2_h, Salinity_ppt_all, pH)
env_nona <- na.omit(env) # n = 208

# Corrplot
C <- cor(env_nona)
corrplot(C, method = "number", type = "lower", 
         tl.cex = 0.5, tl.col = "black", number.cex = 0.5)
pdf("InitialFigs/Comb_All_BGC_corr.pdf", width = 7, height = 5)
meth_corr_by_bgc(env_nona = env_nona)
dev.off()

# Envfit
pcoa <- cmdscale(bc, k = nrow(input_filt_rare$map_loaded) - 1, eig = T)
set.seed(100)
ef <- envfit(pcoa, env, permutations = 999, na.rm = TRUE)
ef
arrow_factor <- ordiArrowMul(ef)
manual_factor <- 0.36
vec.df <- as.data.frame(ef$vectors$arrows*sqrt(ef$vectors$r)) %>%
  mutate(Dim1 = Dim1 * manual_factor,
         Dim2 = Dim2 * manual_factor) %>%
  mutate(variables = rownames(.)) %>%
  filter(ef$vectors$pvals < 0.05) %>%
  mutate(shortnames = c("CH4", "CO2", "Salinity", "pH"))

# Plot with significant vectors
ordiplot(pcoa)
plot(ef, p.max = 0.05, cex = 0.5)
pcoaA1 <- round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[1]*100, digits = 1)
pcoaA2 <- round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[2]*100, digits = 1)
input_filt_rare$map_loaded$Axis01 <- scores(pcoa)[,1]
input_filt_rare$map_loaded$Axis02 <- scores(pcoa)[,2]
micro.hulls <- ddply(input_filt_rare$map_loaded, c("Estuary"), find_hull)
g <- ggplot(input_filt_rare$map_loaded, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = Estuary, fill = Estuary),
               alpha = 0.1, show.legend = F) +
  geom_point(size = 3, alpha = 0.5, aes(colour = Estuary, shape = Depth2),
             show.legend = T) +
  labs(x = paste("PC1: ", pcoaA1, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2, "%", sep = ""),
       shape = "Depth (cm)",
       colour = "Estuary",
       fill = "Estuary") +
  scale_fill_viridis_d(guide = guide_legend(override.aes = list(alpha = 1))) +
  scale_colour_viridis_d(guide = guide_legend(override.aes = list(alpha = 1, shape = 15))) +
  guides(shape = guide_legend(override.aes = list(alpha = 1))) +
  theme_bw() +  
  theme(legend.position = "right",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10)) +
  geom_segment(data = vec.df,
               aes(x = 0, xend = Dim1, y = 0, yend = Dim2),
               arrow = arrow(length = unit(0.5, "cm")),
               colour = "gray", alpha = 0.5,
               inherit.aes = FALSE) + 
  geom_text(data = vec.df,
            aes(x = Dim1, y = Dim2, label = shortnames),
            size = 3, color = "black") +
  coord_fixed()
g
pdf("InitialFigs/Comb_All_PCoA_wBGC.pdf", width = 6, height = 4)
g
dev.off()

# Ordistep
comm_nona <- as.data.frame(t(input_filt_rare$data_loaded)) %>%
  filter(rownames(.) %in% rownames(env_nona))
mod0 <- rda(comm_nona ~ 1, env_nona)  # Model with intercept only
mod1 <- rda(comm_nona ~ ., env_nona)  # Model with all explanatory variables
mod <- ordistep(mod0, scope = formula(mod1))
mod
mod$anova #

# Plot correlations for different taxonomic levels at a given % relative abundance threshold
pdf("InitialFigs/Comb_All_Phyla_CH4.pdf", width = 7, height = 5)
meth_corr_by_taxonomy(input = input_filt_rare, level = 2, threshold = 0.5)
dev.off()
meth_corr_by_taxonomy(input = input_filt_rare, level = 3, threshold = 0.5)
meth_corr_by_taxonomy(input = input_filt_rare, level = 4, threshold = 0.5)
meth_corr_by_taxonomy(input = input_filt_rare, level = 5, threshold = 0.5)
meth_corr_by_taxonomy(input = input_filt_rare, level = 6, threshold = 0.5)
meth_corr_by_taxonomy(input = input_filt_rare, level = 8, threshold = 0.5) 
# ASV4, ASV 3 pos., ASV2 neg.
pdf("InitialFigs/Comb_All_Guilds_CH4.pdf", width = 7, height = 5)
meth_corr_by_taxonomy(input = input_filt_rare, level = 9, threshold = 0)
dev.off()



#### 10. Comparison Field Exp ####
# Only look at field manipulations
# South Carolina and Delaware
# Delaware just look at transplant to oligo which is similar to seawater addition
# Delaware just look at surface depth, similar to SC depth
# "Soil Field plots", and "Soil mesocosm"
#### _Setup ####
# input_filt_rare <- readRDS("input_filt_rare_comb.rds")
input_filt_rare <- readRDS("input_filt_rare_comb_wBGC.rds")
input_filt_rare$map_loaded <- input_filt_rare$map_loaded %>%
  mutate(Estuary = factor(Estuary,
                          levels = c("Waccamaw", "Alligator", "Delaware", "SF")),
         sampleID = rownames(.),
         Field = "NA",
         Salt = "NA",
         Depth = gsub("0-5", 0.05, Depth)) %>%
  mutate(Depth = gsub("5-15", 0.15, Depth)) %>%
  mutate(Depth = as.numeric(Depth)) %>%
  mutate(Depth = ifelse(Depth <= 0.05, "< 5 cm", "5 - 15 cm"))

exp <- filter_data(input_filt_rare,
                   filter_cat = "Site",
                   keep_vals = c("Soil Field plots", "Soil mesocosm")) # 84 samples
table(exp$map_loaded$Estuary)

# Get surface transplants to TFM and OligoHal from Delaware
# Get freshwater amended and saltwater amended from SC
exp <- filter_data(exp,
                   filter_cat = "Detail",
                   keep_vals = c("TFM1@TFM2",
                                 "TFM1@OligoHal",
                                 "Freshwater amended",
                                 "Saltwater amended")) # 28 samples
table(exp$map_loaded$Estuary) # 20 SC, 8 Del

# Make Treatment Column
exp$map_loaded <- exp$map_loaded %>%
  mutate_if(is.character, as.factor) %>%
  mutate(Salt = recode_factor(Detail,
                              "TFM1@TFM2" = "Control",
                              "Freshwater amended" = "Control",
                              "TFM1@OligoHal" = "+ASW",
                              "Saltwater amended" = "+ASW"))

# Need combined estuary/salt factor
exp$map_loaded$EstSalt <- paste(exp$map_loaded$Estuary,
                                exp$map_loaded$Salt,
                                sep = "_")

#### _Alpha ####
leveneTest(exp$map_loaded$rich ~ exp$map_loaded$Salt) # Homogeneous
m <- aov(rich ~ Estuary + Salt + Depth, data = exp$map_loaded)
Anova(m, type = "II") # Salt sig, estuary and depth not sig
shapiro.test(m$residuals) # Normal

leveneTest(exp$map_loaded$shannon ~ exp$map_loaded$Salt) # Homogenoeus
m1 <- aov(shannon ~ Estuary + Salt + Depth, data = exp$map_loaded)
Anova(m1, type = "II") # Salt sig, estuary and depth not sig
shapiro.test(m1$residuals) # Not normal
facet_df <- c("rich" = "(a) Richness",
              "shannon" = "(b) Shannon")
alpha_long <- exp$map_loaded %>%
  pivot_longer(cols = c("rich", "shannon"))
pdf("InitialFigs/Comb_Exp_Alpha.pdf", width = 6, height = 3)
ggplot(alpha_long, aes(reorder(Salt, value, mean), value)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 3, alpha = 1, width = 0.2, aes(shape = Depth, colour = Estuary)) +
  labs(x = "Site", y = "Number of OTUs") +
  scale_colour_viridis_d() +
  facet_wrap(~ name, ncol = 2, scales = "free_y", labeller = as_labeller(facet_df)) +
  theme_bw() +
  theme(legend.position = "right",
        axis.title = element_blank(),
        axis.text = element_text(size = 10),
        strip.text = element_text(size = 10))
dev.off()

#### _Beta ####
exp_bc <- calc_dm(exp$data_loaded)
set.seed(1150)
adonis2(exp_bc ~ Estuary + Salt + Depth, data = exp$map_loaded) # All
adonis2(exp_bc ~ Depth + Salt + Estuary, data = exp$map_loaded) # All
anova(betadisper(exp_bc, exp$map_loaded$Estuary)) # Dispersion homogeneous
anova(betadisper(exp_bc, exp$map_loaded$Salt)) # Dispersion homogeneous
anova(betadisper(exp_bc, exp$map_loaded$Depth)) # Dispersion homogeneous
exp_pcoa <- cmdscale(exp_bc, k = nrow(exp$map_loaded) - 1, eig = T)
pcoaA1 <- round((eigenvals(exp_pcoa)/sum(eigenvals(exp_pcoa)))[1]*100, digits = 1)
pcoaA2 <- round((eigenvals(exp_pcoa)/sum(eigenvals(exp_pcoa)))[2]*100, digits = 1)
exp$map_loaded$Axis01 <- scores(exp_pcoa)[,1]
exp$map_loaded$Axis02 <- scores(exp_pcoa)[,2]
micro.hulls <- ddply(exp$map_loaded, c("Salt"), find_hull)
p_leg2 <- ggplot(exp$map_loaded, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = Salt, fill = Salt),
               alpha = 0.1, show.legend = F) +
  geom_point(size = 3, alpha = 1, aes(colour = Salt, shape = Estuary),
             show.legend = T) +
  labs(x = paste("PC1: ", pcoaA1, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2, "%", sep = ""),
       colour = "Salt") +
  scale_fill_manual(values = c("blue", "red")) +
  scale_colour_manual(values = c("blue", "red")) +
  guides(colour = guide_legend(override.aes = list(alpha = 1,
                                                   shape = 15)),
         fill = "none") +
  theme_bw() +  
  theme(legend.position = "right",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10))
leg2 <- get_legend(p_leg2)
micro.hulls <- ddply(exp$map_loaded, c("EstSalt"), find_hull)
p_noleg2 <- ggplot(exp$map_loaded, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = EstSalt, fill = EstSalt),
               alpha = 0.1, show.legend = F) +
  geom_point(size = 3, alpha = 1, aes(colour = EstSalt, shape = Estuary),
             show.legend = T) +
  labs(x = paste("PC1: ", pcoaA1, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2, "%", sep = ""),
       colour = "Salt") +
  scale_fill_manual(values = c("red", "blue", "red", "blue")) +
  scale_colour_manual(values = c("red", "blue", "red", "blue")) +
  guides(colour = guide_legend(override.aes = list(alpha = 1,
                                                   shape = 15)),
         fill = "none") +
  theme_bw() +  
  theme(legend.position = "none",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10))

pdf("InitialFigs/Comb_Exp_PCoA.pdf", width = 7, height = 5)
plot_grid(p_noleg2, leg2, rel_widths = c(5,1))
dev.off()

#### _Taxa ####
exp_phyla <- summarize_taxonomy(exp, level = 2, report_higher_tax = F)
plot_ts_heatmap(exp_phyla, exp$map_loaded, 0.01, 'EstSalt', rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
exp_barsP <- plot_taxa_bars(exp_phyla, exp$map_loaded, "EstSalt", 
                            num_taxa = 12, data_only = T) %>%
  mutate(taxon = fct_rev(taxon)) %>%
  separate(group_by, into = c("Estuary", "Salt"), sep = "_") %>%
  mutate(Salt = fct_rev(Salt))
pdf("InitialFigs/Comb_Exp_Phyla.pdf", width = 7, height = 5)
ggplot(exp_barsP, aes(Salt, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Sample", y = "Relative abundance", fill = "Phylum") +
  scale_fill_manual(values = c("grey90", brewer.pal(12, "Paired")[12:1])) +
  scale_y_continuous(expand = c(0.01, 0.01)) + 
  facet_wrap(~ Estuary) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        strip.background = element_rect(size = 0.2),
        axis.line.y = element_blank())
dev.off()
taxa_summary_by_sample_type(exp_phyla, exp$map_loaded, 'EstSalt', 0.01, 'KW')

exp_class <- summarize_taxonomy(exp, level = 3, report_higher_tax = F)
plot_ts_heatmap(exp_class, exp$map_loaded, 0.01, 'EstSalt', rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
exp_barsC <- plot_taxa_bars(exp_class, exp$map_loaded, "EstSalt", 
                            num_taxa = 12, data_only = T) %>%
  mutate(taxon = gsub("NA", "Unclassified", taxon)) %>%
  mutate(taxon = fct_relevel(taxon, "Other", after = Inf)) %>%
  mutate(taxon = fct_relevel(taxon, "Unclassified", after = Inf)) %>%
  mutate(taxon = fct_rev(taxon)) %>%
  separate(group_by, into = c("Estuary", "Salt"), sep = "_") %>%
  mutate(Salt = fct_rev(Salt))
pdf("InitialFigs/Comb_Exp_Class.pdf", width = 7, height = 5)
ggplot(exp_barsC, aes(Salt, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = NULL, y = "Relative abundance", fill = "Class") +
  scale_fill_manual(values = c("grey75", "grey90", brewer.pal(12, "Paired")[11:1])) +
  scale_y_continuous(expand = c(0.01, 0.01)) + 
  facet_wrap(~ Estuary) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        strip.background = element_rect(size = 0.2),
        axis.line.y = element_blank())
dev.off()
taxa_summary_by_sample_type(exp_class, exp$map_loaded, 'EstSalt', 0.01, 'KW')

exp_order <- summarize_taxonomy(exp, level = 4, report_higher_tax = F)
plot_ts_heatmap(exp_order, exp$map_loaded, 0.01, 'EstSalt', rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
exp_barsO <- plot_taxa_bars(exp_order, exp$map_loaded, "EstSalt", 
                            num_taxa = 12, data_only = T) %>%
  mutate(taxon = gsub("NA", "Unclassified", taxon)) %>%
  mutate(taxon = fct_relevel(taxon, "Other", after = Inf)) %>%
  mutate(taxon = fct_relevel(taxon, "Unclassified", after = Inf)) %>%
  mutate(taxon = fct_rev(taxon)) %>%
  separate(group_by, into = c("Estuary", "Salt"), sep = "_") %>%
  mutate(Salt = fct_rev(Salt))
pdf("InitialFigs/Comb_Exp_Order.pdf", width = 7, height = 5)
ggplot(exp_barsO, aes(Salt, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = NULL, y = "Relative abundance", fill = "Order") +
  scale_fill_manual(values = c("grey75", "grey90", brewer.pal(12, "Paired")[11:1])) +
  scale_y_continuous(expand = c(0.01, 0.01)) + 
  facet_wrap(~ Estuary) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        strip.background = element_rect(size = 0.2),
        axis.line.y = element_blank())
dev.off()
taxa_summary_by_sample_type(exp_order, exp$map_loaded, 'EstSalt', 0.01, 'KW')

exp_family <- summarize_taxonomy(exp, level = 5, report_higher_tax = F)
plot_ts_heatmap(exp_family, exp$map_loaded, 0.01, 'EstSalt', rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
exp_barsF <- plot_taxa_bars(exp_family, exp$map_loaded, "EstSalt", 
                            num_taxa = 12, data_only = T) %>%
  mutate(taxon = gsub("NA", "Unclassified", taxon)) %>%
  mutate(taxon = fct_relevel(taxon, "Other", after = Inf)) %>%
  mutate(taxon = fct_relevel(taxon, "Unclassified", after = Inf)) %>%
  mutate(taxon = fct_rev(taxon)) %>%
  separate(group_by, into = c("Estuary", "Salt"), sep = "_") %>%
  mutate(Salt = fct_rev(Salt))
pdf("InitialFigs/Comb_Exp_Family.pdf", width = 7, height = 5)
ggplot(exp_barsF, aes(Salt, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = NULL, y = "Relative abundance", fill = "Family") +
  scale_fill_manual(values = c("grey75", "grey90", brewer.pal(12, "Paired")[11:1])) +
  scale_y_continuous(expand = c(0.01, 0.01)) + 
  facet_wrap(~ Estuary) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        strip.background = element_rect(size = 0.2),
        axis.line.y = element_blank())
dev.off()
taxa_summary_by_sample_type(exp_family, exp$map_loaded, 'EstSalt', 0.01, 'KW')

exp_genus <- summarize_taxonomy(exp, level = 6, report_higher_tax = F)
plot_ts_heatmap(exp_genus, exp$map_loaded, 0.005, 'EstSalt', rev_taxa = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
exp_barsG <- plot_taxa_bars(exp_genus, exp$map_loaded, "EstSalt", 
                            num_taxa = 12, data_only = T) %>%
  mutate(taxon = gsub("NA", "Unclassified", taxon)) %>%
  mutate(taxon = fct_relevel(taxon, "Other", after = Inf)) %>%
  mutate(taxon = fct_relevel(taxon, "Unclassified", after = Inf)) %>%
  mutate(taxon = fct_rev(taxon)) %>%
  separate(group_by, into = c("Estuary", "Salt"), sep = "_") %>%
  mutate(Salt = fct_rev(Salt))
pdf("InitialFigs/Comb_Exp_Genus.pdf", width = 7, height = 5)
ggplot(exp_barsG, aes(Salt, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = NULL, y = "Relative abundance", fill = "Genus") +
  scale_fill_manual(values = c("grey75", "grey90", brewer.pal(12, "Paired")[11:1])) +
  scale_y_continuous(expand = c(0.01, 0.01)) + 
  facet_wrap(~ Estuary) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        strip.background = element_rect(size = 0.2),
        axis.line.y = element_blank())
dev.off()
taxa_summary_by_sample_type(exp_genus, exp$map_loaded, 'EstSalt', 0.01, 'KW')

exp_guilds <- summarize_taxonomy(exp, level = 9, report_higher_tax = F)
plot_ts_heatmap(exp_guilds, exp$map_loaded, 0, 'EstSalt', rev_taxa = T, remove_other = T) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5))
exp_barsGu <- plot_taxa_bars(exp_guilds, exp$map_loaded, "EstSalt",
                             num_taxa = 20, data_only = T) %>%
  filter(taxon != "NA") %>%
  droplevels() %>%
  mutate(taxon = factor(taxon,
                        levels = Guild_cols$Guild)) %>%
  separate(group_by, into = c("Estuary", "Salt"), sep = "_") %>%
  mutate(Salt = fct_rev(Salt))
tallest_bar <- exp_barsGu %>%
  group_by(Estuary, Salt) %>%
  summarise(sum = sum(mean_value))
pdf("InitialFigs/Comb_Exp_Guilds.pdf", width = 7, height = 5)
ggplot(exp_barsGu, aes(Salt, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = "Estuary", y = "Relative abundance", fill = "Guild") +
  scale_fill_manual(values = Guild_cols$color) +
  scale_y_continuous(expand = c(max(tallest_bar$sum)/100, max(tallest_bar$sum)/100)) + 
  facet_wrap(~ Estuary) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        strip.background = element_rect(size = 0.2),
        axis.line.y = element_blank())
dev.off()
taxa_summary_by_sample_type(exp_guilds, exp$map_loaded, 'EstSalt', 0.01, 'KW')

#### _Methano ####
# Look at methanogens (careful to do this accurately!)
tax_sum_families_meth <- summarize_taxonomy(exp, level = 5, 
                                            report_higher_tax = F)
tax_sum_families_meth <- tax_sum_families_meth[grep("(Methano|Methermicoccaceae|
                                                    Syntrophoarchaeaceae)", 
                                                    rownames(tax_sum_families_meth)),]
tax_sum_families_meth <- tax_sum_families_meth[!grepl("Methanoperedenaceae", 
                                                      rownames(tax_sum_families_meth)),]
exp_barsMethano <- plot_taxa_bars(tax_sum_families_meth, 
                                  exp$map_loaded, 
                                  "EstSalt",
                                  num_taxa = nrow(tax_sum_families_meth), 
                                  data_only = T) %>%
  separate(group_by, into = c("Estuary", "Salt"), sep = "_") %>%
  mutate(Salt = fct_rev(Salt))
nb.cols <- nrow(tax_sum_families_meth)
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)
tallest_bar <- exp_barsMethano %>%
  group_by(Estuary, Salt) %>%
  summarise(sum = sum(mean_value))
pdf("InitialFigs/Comb_Exp_Methano.pdf", width = 7, height = 5)
ggplot(exp_barsMethano, aes(Salt, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = NULL, y = "Relative abundance", fill = "Family") +
  scale_fill_manual(values = mycolors) +
  scale_y_continuous(expand = c(max(tallest_bar$sum)/100, max(tallest_bar$sum)/100)) + 
  facet_wrap(~ Estuary) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        strip.background = element_rect(size = 0.2),
        axis.line.y = element_blank())
dev.off()

#### _Desulfo ####
# Subset taxa to SRB and SRB_syn guilds
exp_srb <- filter_taxa_from_input(exp,
                                  taxa_to_keep = c("SRB", "SRB_syn"),
                                  at_spec_level = 9)
desulfo_wTax <- summarize_taxonomy(exp_srb, level = 3, 
                                   report_higher_tax = T, relative = FALSE)
rownames(desulfo_wTax) <- substring(rownames(desulfo_wTax), 11)
nb.cols <- nrow(desulfo_wTax)
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)
exp_barsDesulfo <- plot_taxa_bars(desulfo_wTax, exp$map_loaded, "EstSalt", 
                                  num_taxa = nrow(desulfo_wTax), data_only = T) %>%
  mutate(mean_value = mean_value/26429) %>%
  separate(group_by, into = c("Estuary", "Salt"), sep = "_") %>%
  mutate(Salt = fct_rev(Salt))
tallest_bar <- exp_barsDesulfo %>%
  group_by(Estuary, Salt) %>%
  summarise(sum = sum(mean_value))
pdf("InitialFigs/Comb_Exp_Desulfo.pdf", width = 7, height = 5)
ggplot(exp_barsDesulfo, aes(Salt, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, size = 0.25) +
  labs(x = NULL, y = "Relative abundance", fill = "Class") +
  scale_fill_manual(values = mycolors) +
  scale_y_continuous(expand = c(max(tallest_bar$sum)/100, max(tallest_bar$sum)/100)) + 
  facet_wrap(~ Estuary) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        strip.background = element_rect(size = 0.2),
        axis.line.y = element_blank(),
        legend.key.size = unit(0.25, "cm"),
        legend.text = element_text(size = 10))
dev.off()

#### _Simper ####
exp_sim <- simper(t(exp$data_loaded), 
                  exp$map_loaded$Salt)
exp_s <- summary(exp_sim)
head(exp_s$`Control_+ASW`)
exp_df1 <- head(exp_s$`Control_+ASW`, n = 20) %>%
  mutate(ASV = rownames(.),
         Comparison = "Control/+ASW") %>%
  left_join(., exp$taxonomy_loaded, by = c("ASV" = "taxonomy8")) %>%
  mutate(SaltResponse = ifelse(avb > ava, "Increase", "Decrease")) %>%
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
write_xlsx(exp_df1, 
           "simper_results_comb_exp.xlsx",
           format_headers = F)

#### _Multipatt ####
set.seed(1202)
exp_mp <- multipatt(t(exp$data_loaded), 
                    exp$map_loaded$Salt, 
                    func = "r.g", 
                    control = how(nperm=999),
                    max.order = 1)
exp_mp_results <- exp_mp$sign %>%
  mutate(q.value = qvalue(exp_mp$sign$p.value)$qvalues,
         Group = "NA",
         ASV = rownames(exp_mp$sign)) %>%
  filter(p.value < 0.01) %>%
  filter(stat >= 0.55) %>%
  left_join(., exp$taxonomy_loaded, by = c("ASV" = "taxonomy8"))
for (i in 1:nrow(exp_mp_results)) {
  if (exp_mp_results$s.Control[i] == 1) {
    exp_mp_results$Group[i] <- "Control"
  }
}
for (i in 1:nrow(exp_mp_results)) {
  if (exp_mp_results$`s.+ASW`[i] == 1) {
    exp_mp_results$Group[i] <- "+ASW"
  }
}
table(exp_mp_results$Group)
exp_asv <- summarize_taxonomy(exp, level = 8, report_higher_tax = F)
exp_asv_all <- data.frame("RelAbundance" = round(rowMeans(exp_asv) * 100, digits = 5)) %>%
  mutate(ASV = rownames(.))
exp_mp_corrs <- as.data.frame(exp_mp$str) %>%
  dplyr::select(1:length(levels(exp$map_loaded$Salt))) %>%
  mutate(ASV = rownames(.)) %>%
  filter(ASV %in% exp_mp_results$ASV) %>%
  set_names(c("Control", "+ASW", "ASV"))
# Add corrs and taxonomy
exp_mp_results <- exp_mp_results %>%
  left_join(., exp_asv_all, by = "ASV") %>%
  left_join(., exp_mp_corrs, by = "ASV") %>%
  dplyr::select(taxonomy2, taxonomy6, ASV, Group, names(exp_mp_corrs)[1:2], "RelAbundance") %>%
  mutate(ASV = gsub("ASV", "OTU", ASV)) %>%
  arrange(Group, taxonomy2, taxonomy6) %>%
  mutate(ASV = factor(ASV, levels = ASV)) %>%
  set_names(c("Phylum", "Genus", "ASV", "Group", names(exp_mp_corrs)[1:2], "RelAbundance")) %>%
  mutate(taxon = paste(Phylum, Genus, ASV, sep = "; ")) %>%
  mutate(taxon = factor(taxon, levels = taxon))
# Plot heatmap with abundance barplot
exp.hm.melted <- exp_mp_results %>%
  dplyr::select(taxon, names(exp_mp_corrs)[1:2]) %>%
  melt(., id.vars = c("taxon"))
exp.hm <- ggplot(data = exp.hm.melted, 
                 aes(x = factor(taxon), y = variable, fill = value)) + 
  geom_tile() + 
  scale_fill_distiller(name = "Indicator correlation index", palette = "RdBu", direction = -1, na.value = "transparent", type = "div", limits = c(-1, 1)) +
  scale_x_discrete(breaks = unique(exp.hm.melted$taxon), labels = unique(exp.hm.melted$taxon),
                   limits = rev(levels(exp.hm.melted$taxon))) +
  scale_y_discrete(position = "right") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        legend.title = element_text(size = 10), legend.key.height = unit(0.5,"cm"),
        legend.key.width = unit(1, "cm"), legend.text = element_text(size = 7)) +
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5)) +
  coord_flip()
exp.l <- get_legend(exp.hm)
exp.hm.clean <- exp.hm +
  theme(axis.title.y = element_blank(), axis.text.y = element_text(margin = margin(c(0,-5,0,0)), size = 5),
        axis.ticks.y = element_blank(), axis.ticks.x = element_blank(), 
        panel.grid.major = element_blank(), legend.position="none",
        axis.text.x.top = element_text(size = 6, margin = margin(c(0,0,-2,0))), 
        plot.margin = margin(c(0,-2,0,5))) +
  labs(y = "Group")
exp.bp.y <- ggplot(data = exp_mp_results, aes(x = taxon, y = RelAbundance)) + 
  geom_bar(stat = "identity", fill = "grey") + 
  scale_y_continuous(position = "right") +
  scale_x_discrete(limits = rev(levels(exp_mp_results$taxon))) +
  coord_flip() + 
  theme_minimal() +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position="none", plot.margin = margin(c(0,5,0,-5)),
        axis.text.x.top = element_text(size = 6, margin = margin(c(0,0,-2,0)))) +
  labs(y = "Rel. abund. (%)")
pdf("InitialFigs/Comb_Exp_Multipatt.pdf", width = 8, height = 10)
plot_grid(exp.hm.clean, exp.bp.y, exp.l, NULL, nrow = 2, ncol = 2, 
          rel_widths = c(10,2), rel_heights = c(12, 2), align = "hv", axis = "b")
dev.off()

#### _BGC ####
# Add biogeochem. analysis
# Get variables, make corrplot, envfit, ordination with vectors, ordistep, taxa correlations
# Note, porewater chemistry is only on D2 samples (5-15 cm depth)

# Get variables
env <- exp$map_loaded %>%
  dplyr::select(CH4_ug_m2_h, Salinity_ppt_all)
env_nona <- na.omit(env) # n = 

# Corrplot
C <- cor(env_nona)
corrplot(C, method = "number", type = "lower", 
         tl.cex = 0.5, tl.col = "black", number.cex = 0.5)
pdf("InitialFigs/Comb_Exp_BGC_corr.pdf", width = 7, height = 5)
meth_corr_by_bgc(env_nona = env_nona)
dev.off()

# Envfit
pcoa <- cmdscale(exp_bc, k = nrow(exp$map_loaded) - 1, eig = T)
set.seed(100)
ef <- envfit(pcoa, env, permutations = 999, na.rm = TRUE)
ef
arrow_factor <- ordiArrowMul(ef)
manual_factor <- 0.5
vec.df <- as.data.frame(ef$vectors$arrows*sqrt(ef$vectors$r)) %>%
  mutate(Dim1 = Dim1 * manual_factor,
         Dim2 = Dim2 * manual_factor) %>%
  mutate(variables = rownames(.)) %>%
  filter(ef$vectors$pvals < 0.05) %>%
  mutate(shortnames = c("CH4"))

# Plot with significant vectors
ordiplot(pcoa)
plot(ef, p.max = 0.05, cex = 0.5)
pcoaA1 <- round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[1]*100, digits = 1)
pcoaA2 <- round((eigenvals(pcoa)/sum(eigenvals(pcoa)))[2]*100, digits = 1)
exp$map_loaded$Axis01 <- scores(pcoa)[,1]
exp$map_loaded$Axis02 <- scores(pcoa)[,2]
micro.hulls <- ddply(exp$map_loaded, c("EstSalt"), find_hull)
g <- ggplot(exp$map_loaded, aes(Axis01, Axis02)) +
  geom_polygon(data = micro.hulls, 
               aes(colour = EstSalt, fill = EstSalt),
               alpha = 0.1, show.legend = F) +
  geom_point(size = 3, alpha = 0.5, aes(colour = EstSalt, shape = Estuary),
             show.legend = T) +
  labs(x = paste("PC1: ", pcoaA1, "%", sep = ""), 
       y = paste("PC2: ", pcoaA2, "%", sep = ""),
       shape = "Depth (cm)",
       colour = "Salt") +
  scale_fill_manual(values = c("red", "blue", "red", "blue")) +
  scale_colour_manual(values = c("red", "blue", "red", "blue")) +
  guides(colour = guide_legend(override.aes = list(alpha = 1,
                                                   shape = 15)),
         fill = "none") +
  theme_bw() +  
  theme(legend.position = "none",
        axis.title = element_text(face = "bold", size = 12), 
        axis.text = element_text(size = 10)) +
  geom_segment(data = vec.df,
               aes(x = 0, xend = Dim1, y = 0, yend = Dim2),
               arrow = arrow(length = unit(0.5, "cm")),
               colour = "gray", alpha = 0.5,
               inherit.aes = FALSE) + 
  geom_text(data = vec.df,
            aes(x = Dim1, y = Dim2, label = shortnames),
            size = 4, color = "black") +
  coord_fixed()
plot_grid(g, leg2, rel_widths = c(5,1))
pdf("InitialFigs/Comb_Exp_PCoA_wBGC.pdf", width = 6, height = 4)
plot_grid(g, leg2, rel_widths = c(5,1))
dev.off()

# Ordistep
comm_nona <- as.data.frame(t(exp$data_loaded)) %>%
  filter(rownames(.) %in% rownames(env_nona))
mod0 <- rda(comm_nona ~ 1, env_nona)  # Model with intercept only
mod1 <- rda(comm_nona ~ ., env_nona)  # Model with all explanatory variables
mod <- ordistep(mod0, scope = formula(mod1))
mod
mod$anova # Variables aren't better than null model...

# Plot correlations for different taxonomic levels at a given % relative abundance threshold
pdf("InitialFigs/Comb_All_Phyla_CH4.pdf", width = 7, height = 5)
meth_corr_by_taxonomy(input = exp, level = 2, threshold = 0.5)
dev.off()
meth_corr_by_taxonomy(input = exp, level = 3, threshold = 0.5)
meth_corr_by_taxonomy(input = exp, level = 4, threshold = 0.5)
meth_corr_by_taxonomy(input = exp, level = 5, threshold = 0.5)
meth_corr_by_taxonomy(input = exp, level = 6, threshold = 0.5)
meth_corr_by_taxonomy(input = exp, level = 8, threshold = 0.5) 
# ASV4, ASV 3 pos., ASV2 neg.
pdf("InitialFigs/Comb_All_Guilds_CH4.pdf", width = 7, height = 5)
meth_corr_by_taxonomy(input = exp, level = 9, threshold = 0)
dev.off()



#### End Script ####