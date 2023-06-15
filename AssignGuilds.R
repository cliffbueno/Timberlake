#################################################################################################
#################################################################################################
# OTU subsetting modules v.0.1
# Split from plotting heatmaps 3.4

#################################################################################################
#################################################################################################
#### OUTLINE ####

### 1) Filter OTU table by sample subset
##   a) OTU sample subsetting Function
#        Delta_OTU <- Samp_subset_OTU(otu_V, Delta_sites)

### 2) Filter corrRank OTUs by correlation cutoff
##   a) Corr Ranks filtered --  For heatmaps
#        corrRank_gU <- unique(corrRank_OTU[1:(ncol(corrRank_OTU)-3)]) # Unique genus data
                                                                       # Drops OTU id, "Taxon" corr[last 3 cols]
#        Corr_filt_r <- unique(corrRanks_filt(otu_V, corrRank_gU, "Genus", 0.5))

##   b) Corr filtered OTUs  -- for barplots
#        CorrFilt_OTU <- corrFilt_OTU(otu_V, CH4corrRanksOTU, "Genus", 0.5)

### 3) Get OTU functional "GUILDS" for BGC cycling
##   a) Function to make grep output DF with OTU #'s
#        methano <- Subs_to_DF(methano)

##   b) Function to retrive 16S Guilds using grep
#        Guild_OTUs <- Get_16S_Guilds(otu_V)

##   c) Generate color pallete for Guilds
#        Made manually, saved as fixed

#################################################################################################
# Import packagaes
suppressMessages(library(Hmisc))
library(ggplot2)
library(RColorBrewer)
library(reshape2)
#suppressMessages(library(cowplot))
#theme_set(theme_grey())

#################################################################################################
## Import testing Data

# Read import module for testing convenience here, still a bit of mucking around in pre-processing
#source("Import_SalOTU_dat_Plot_test_v0.1.R")

# Reimport color legend from barplots
#Taxons3d <-read.table("Test_taxonomy_color_pallete.txt", sep="\t", header = T)
# head(Taxons3d)

#################################################################################################
#################################################################################################

#################################################################################################
#### 1) Sample Subsetting

### Generate separate Delta sites metadata
# Separate Delta Sites (oligo and FW)
#levels(Meta_iTag$SALTgroup)
#Meta_iTag_FW <- Meta_iTag[Meta_iTag$SALTgroup =="FW",]
#Meta_iTag_Oligo <- Meta_iTag[Meta_iTag$SALTgroup =="Oligo",]

#Meta_iTag_Delta <-rbind(Meta_iTag_FW, Meta_iTag_Oligo)
# Meta_iTag_Delta

#Delta_sites <-data.frame(Meta_iTag_Delta[,"Sample"])                 # DF of Delta sites
#colnames(Delta_sites) [1] <- "Sample"                                # Rename column "Sample"
#Delta_sites["Samp_Index"] <- seq(1:nrow(Delta_sites))                # Make sample index for reordering post-merges
# Delta_sites

#####################################
### a) OTU sample subsetting function

# Takes as input an OTU Table and Sample subset, which contains a sample column and index column
# TODO: ABSTRACT VARS IN FUNCTION, ADDING DEFAULTS

# get OTU sample subset function
Samp_subset_OTU = function(otu_V, Samp_subset) { 
# Separate data and taxonomy columns
    Count_cols <- sapply(otu_V, is.numeric)                               # are columns numeric?
    Tax_cols <- otu_V[!Count_cols]                                        # Get Taxonomy columns  dim(Tax_cols)  # head(Tax_cols)
    otu_d <- data.frame(otu_V[Count_cols])                                # Get count data (Sample) columns

    # Subset data columns by Sample subset, using merge on transp.        # Passing column list didn't work, why?
    otu_dT <-data.frame(t(otu_V[Count_cols]))                             # transp. counts and DF 
    otu_dT["Sample"] <- row.names(otu_dT)                                 # Make sample column 
    otu_SubsetT <- merge(Samp_subset, otu_dT)                             # Merge with subset of samples 
    otu_SubsetT <- otu_SubsetT[order(otu_SubsetT["Samp_Index"]),]         # Sort by Sample Index

    # Retranspose and cbind Taxonomy data
    row.names(otu_SubsetT) <- otu_SubsetT[,"Sample"]                      # replace rownames w. Sample names (pre-transp.)
    otu_SubsetT <- otu_SubsetT[3:ncol(otu_SubsetT)]                       # Drop first two cols (Sample, Samp_Index from merge)
    otu_Subset_d <-data.frame(t(otu_SubsetT))                             # Transpose to OTU table               #dim(otu_Delta_d) #head(otu_Delta_d)
    otu_Subset <- cbind(otu_Subset_d, Tax_cols)                           # add back Taxonomy data
    return(otu_Subset)
}

# Subset Delta samples from OTU table
#Delta_OTU <- Samp_subset_OTU(otu_V, Delta_sites)
#head(Delta_OTU)



#################################################################################################
#### 2) Filter corrRank object by correlation cutoff ####

#########################################
### a) Corr Ranks filtered --for heatmaps

# Takes a corrRanks table of correlations, cuts rows with r < r_cut in rank_var
# returns only Taxonomy data and Corr across ranks

corrRanks_filt = function(otu_V, corrRanks, rank_var, r_cut){          # Simplified from v.3.1

    #corrRanks <- CH4corrRanksOTU_Delta 
    #chunk_var <- "Genus"
    rank_r <- paste0(rank_var,"_r")
    abs_rank_r <-paste0("abs_", rank_r)
    #r_cut <- 0.5

    # Get abs val of ranks
    corrRanks[abs_rank_r] <-abs(corrRanks[,rank_r])

    # Get corr filtered ranks 
    corrRanks_F <- corrRanks[corrRanks[,abs_rank_r] > r_cut,]         
    corrRanks_F <- corrRanks_F[,1:(ncol(corrRanks_F)-1)]
    row.names(corrRanks_F) <-  make.unique(as.character(corrRanks_F[,rank_var]))  # Rename rows
    return(corrRanks_F)
    
    
    # Select column of interest, make DF
    #r_Filt_rank <-data.frame(corrRanks_F[, rank_var])
    #colnames(r_Filt_rank)[1] <- rank_var
    #r_Filt_rank[rank_r]<-corrRanks_F[,rank_r]
    #return(r_Filt_rank)
    return(corrRanks_F)
}

# prepare data 
#corrRank_OTU <- CH4corrRanksOTU_Delta
#corrRank_gU <-unique(corrRank_OTU[1:(ncol(corrRank_OTU)-3)])                     # Drop OTU data, Taxonomy corr...   # remove last 3 cols: OTU and Tax corr

# Use corrRanks_Filt function
#Corr_filt_r <- unique(corrRanks_filt(otu_V, corrRank_gU, "Genus", 0.5))
#dim(Corr_filt_r); head(Corr_filt_r)


#########################################
### b) Corr filtered OTUs --for barplots

# Takes OTU table and corrRanks, returns OTUs > corr cutoff r value

# Function to get corr. filtered OTU table from corrRanks object
corrFilt_OTU = function(otu_V, corrRanks, chunk_var, r_cut){

    # Filter corrRanks by abs(chunk_var) > r_cut 
    chunk_r <- paste0(chunk_var,"_r")                               # Get colname tax rank corr. var., here [chunk] 
    abs_chunk_r <-paste0("abs_", chunk_r)                           # Get colname for abs of corr(chunk) 
    corrRanks[abs_chunk_r] <-abs(corrRanks[,chunk_r])               # Get abs val of chunk
    corrRanks_F <- corrRanks[corrRanks[,abs_chunk_r] > r_cut,]      # Get corr filtered corrRanks data 

    # Make DF of corr. filtered Taxa, corr. r values, OTU_ID
    r_Filt_rank <-data.frame(corrRanks_F[,chunk_var])               # Taxa names 
    colnames(r_Filt_rank)[1] <- chunk_var                           # Label Taxa col
    r_Filt_rank[chunk_r]<-corrRanks_F[,chunk_r]                     # Get correlation r values
    r_Filt_rank["OTU"] <- corrRanks_F[,"OTU"]                       # Add OTU column for merge below...      # dim(r_Filt_rank); head(r_Filt_rank)

    # Merge with abundance table 
    otu_Vo <- otu_V                                                 # copy otu table
    otu_Vo["OTU"]<-row.names(otu_Vo)                                # add otu_IDs as col.
    r_Filt_OTU <- merge(otu_Vo, r_Filt_rank)                        # Merge with abundance 
    return(r_Filt_OTU)                                                                                  # dim(r_Filt_OTU); head(r_Filt_OTU)
}


# Test OTU_corr_filt function
#CorrFilt_OTU <- corrFilt_OTU(otu_V, CH4corrRanksOTU, "Genus", 0.5)
#dim(CorrFilt_OTU); head(CorrFilt_OTU)

#CorrFilt_OTU_Delta <- corrFilt_OTU(Delta_OTU, CH4corrRanksOTU_Delta, "Genus", 0.5)
# dim(CorrFilt_OTU_Delta); head(CorrFilt_OTU_Delta)


#################################################################################################
#### 3) Get OTU functional "Guilds" for BGC cycling ####

###################################################
### a) Function to make grep output DF with OTU #'s
# Needed to parse many greps below

## function for grep results to strip Taxa and make factors
Subs_to_DF = function(otu_Sub, Tax_sort="Consensus.lineage"){
    Count_cols <- sapply(otu_Sub, is.numeric)                               # are columns numeric?
    Tax_cols <- otu_Sub[!Count_cols]                                        # Get Taxonomy columns  dim(Tax_cols)  # head(Tax_cols)
    Tax_cols[] <- lapply(Tax_cols, factor)                                   # Convert each col to factor (needed after grep)                      

    #Tax_Sub[] <- lapply(Tax_cols, factor)                                   # Convert each col to factor (needed after grep)                      
    Tax_Sub <-data.frame(Tax_cols)                                           # Make data frame, otu rownames dropped

    #Tax_Sub <-data.frame(Tax_Sub[])                                           # Make data frame, otu rownames dropped
    #Tax_Sub["OTU"]<-row.names(Tax_cols)                                     # Replace OTU names as column
    Tax_Sub <-Tax_Sub[order(Tax_Sub[,Tax_sort]),]                           # sort by Tax_sort
    return(Tax_Sub)
}

##################################################
### b) Function to retrieve 16S Guilds using grep
# approximation of BGC function using organism names x Review Papers on function

Get_16S_Guilds = function(otu_V){

    #############################################
    # Methanogens

    methano <- otu_V[grepl("Methano", otu_V$Consensus.lineage),]                        # Get methanogens 
    methano <- Subs_to_DF(methano)                                                      # Fxn to get Taxon. DF with factors   
    # 2022 CHANGE HERE
    methano["Guild"] <- ifelse(methano$Order == 'Methanosarciniales', "CH4_ac", "CH4_H2")
    # Note here did not separate the strictly acetoclastic families Methanosaetaceae and Methanotrichaceae from other mixotrophic Methanosarcina(s)
    # dim(methano); levels(methano$Genus); # unique(methano)  # methano

    #############################################
    # Methylotrophs 

    methylo <- otu_V[grepl("Methylo", otu_V$Consensus.lineage),]                                 # Get methylo from OTU table    
    methylo <- Subs_to_DF(methylo)                                                               # use function to get Taxon. DF with factors   

    # Use ifelse to make new column of methylo GUILDS
    # MOB IIa based on Kneif et al. 2015 Front. Microb, other fine differences ignored (e.g. IIb, Ia-c)
    methylo["Guild"] <- ifelse(methylo$Class=='Gammaproteobacteria', "MOB_I", 
                            #  ifelse(methylo$Family=='Methylocystaceae', "MOB_IIa",            # MOB_IIa called before...
                            ifelse(methylo$Family=='Beijerinckiaceae', "MOB_IIa",            # MOB_IIa called before...
                            
                            ### f__Beijerinckiaceae? or
                            # o__Rhizobiales;
                            
                            ifelse(methylo$Class=='Alphaproteobacteria', "MOB_II",           # all MOB II
                            ifelse(methylo$Class=='Betaproteobacteria', "MeOB", NA ))))      # "Methyoltrophs"

    methylo["Guild"]<-factor(methylo$Guild, levels = c("MOB_I","MOB_II", "MOB_IIa", "MeOB"))     # Reorder Guild factor
    methylo <- methylo[order(methylo$Guild),]                                                    # Sort by Guild
    #dim(methylo); levels(methylo$Genus); levels(methylo$Guild); unique(methylo); #methylo  #levels(MOB_II$Class);

    #############################################
    # ANME  -- weren't included in total, no counts
    ANME <- otu_V[grepl("ANME", otu_V$Consensus.lineage),]                        # Get methylo from OTU table    
    ANME <- Subs_to_DF(ANME)                                                      # use function to get Taxon. DF with factors 
    #ANME

    #############################################
    # Nitrifiers

    # Get ammonia oxidizers and nitrite oxidizers                                 # Slightly trickier than other groups
    nitros <- otu_V[grepl("Nitros", otu_V$Consensus.lineage),]   

    # Get only AOA, AOB (nitroso)
    nitroso <- nitros[grepl("Nitroso", nitros$Consensus.lineage),]       
    nitroso <- Subs_to_DF(nitroso)                                                # use function to get Taxon. DF with factors      
    nitroso["Guild"] <- ifelse(nitroso$Kingdom=='Archaea', "AOA", "AOB")          # Label AOA vs. AOB by Kingdom 
    # dim(nitroso); levels(nitroso$Genus);  #unique(nitroso)

    # Get NOB as opposite set from AOA/AOB, drop SRB
    nitros2 <- nitros[!grepl("Nitroso", nitros$Consensus.lineage),]       
    nitros2 <- nitros2[!grepl("desulfo", nitros2$Consensus.lineage),]             # Drop sulfate reducing Nitrospira (Not MOB), ambigs
    # nitros2 <- nitros2[!grepl("NitrospiralesOR", nitros2$Consensus.lineage),]   # Could drop ambig (Order, could be desulf)
                                                                              # But appears at least one is (+)~CH4
    nitros2 <- Subs_to_DF(nitros2)                                                # use function to get Taxon. DF with factors      
    nitros2["Guild"] <- "NOB"                                                     # dim(nitros2); levels(nitros2$Family);   #unique(nitros2) #nitros2# [1:40,]
    
    #############################################
    # ANAMMOX
    anamox_list <-"(Kuenenia|Anammoxoglobus|Scalindua|Brocadia|Jettenia)"
    anamox <- otu_V[grepl(anamox_list, otu_V$Consensus.lineage),] 
    anamox <- Subs_to_DF(anamox)
    anamox["Guild"] <-"Anamx"                                                     # levels(anamox$Genus); unique(anamox)# anamox; 


    #############################################
    # Sulfate Reducers
    
    ##### Split out syntrophs SEPARATELY
    # several non SRB, non-"Desulfo" syntrophs; eg. several ;o__Clostridiales;f__Syntrophomonadaceae; Delta D;f__Syntrophobacteraceae
    # SIVLA updated list
    Syntroph_list<-"(Desulfobacterales|Desulfatiglandales|Desulfobaccales|Desulfomonilales)"
    Syntroph <- otu_V[grepl(Syntroph_list, otu_V$Consensus.lineage),] 
    Syntroph <- Subs_to_DF(Syntroph)
    Syntroph["Guild"] <-"SRB_syn"                                             # dim(Syntroph); levels(Syntroph$Genus); #unique(Syntroph)
    
    # Sourced from SRB reviews: Plugge et al. 2011 Front Microb; Muller et al. 2015 ISMEJ
    Desulf <- otu_V[grepl("Desulf", otu_V$Consensus.lineage),]               # Get Desulfo --ALL

    # Drop non SRB "Desulf"
    Desulf <- Desulf[!grepl("Nitros", Desulf$Consensus.lineage),]            # Drops Nitrospinaceae, in Order: Desulfobacterales 
    Desulf <-Desulf[!grepl(Syntroph_list, Desulf$Consensus.lineage),]        # Drops the above Syntrophs 
    Desulf <-Desulf[!grepl("Geobacter",Desulf$Consensus.lineage),]           # Drops Geobacter (own group), according to Plugge et al. 2011 Front. Microb.
    Desulf <- Subs_to_DF(Desulf)                                             # dim(Desulf); levels(Desulf$Genus); # unique(Desulf)

    ##### Get misc. SRB (not Desulf), based on Trees in Muller et al. 2015 ISME
    SRB_xtra<-"(Caldiserica|Thermanaeromonas|Sporomusaceae|Carboxydothermus|Pelotomaculum|
                Moorella |Ammonifex|Acetonema|Thermosinus|Thermanaeromonas|Carboxydothermus|Caldiserica|
                Gordonibacter|Thermodesulfobium|Thermodesulfovibrio|Magnetobacterium)" 
    Desulf2 <- otu_V[grepl(SRB_xtra, otu_V$Consensus.lineage),]               # grep SRB_xtra 
    Desulf2 <- Subs_to_DF(Desulf2)                                            # dim(Desulf2); levels(Desulf2$Genus); #unique(Desulf2)# Desulf2; 

    # Sulfate reducing ARCHAEA, from Muller et al. 2015 ISME
    SRA_list<-"(Archeoglobus|Pyrobaculum|Vulcanisaeta|Caldirvirga)"           # Sulfate red. archaea
    SRA <- otu_V[grepl(SRA_list, otu_V$Consensus.lineage),]                   # Get SRA
    SRA <- Subs_to_DF(SRA)                                                    # levels(SRA$Genus); unique(SRA)# SRA

    ## Combine SRB lists and label (note no SRA here so passing over)
    SRB <- rbind(Desulf, Desulf2)
    SRB["Guild"]<-"SRB"
    #dim(SRB); levels(SRB$Genus); #SRB



    #############################################
    # Sulfur Oxidizers
    # While many are Thio, add rest from Muller et al. 2015 ISMEJ to same list
    SOxB_list<-"(Thio|Allochromatium|Marichromatium|Halochromatium|Alkalilimnicola|Halorhodospira|Ruthia|Vesicomyosocius|Sedimenticola
            Sideroxydans|Sulfuricella|Riegeria|Azospirillum|Rhodomicrobium|Magnetospirillum|Magnetococcus|
            Chlorobium|Chlorobaculum|Prosthecochloris)"
    # 1st line are Gamma prot, Fig. 3;     # 2nd line are B and A prot.;      #  3rd line Chlorobi

    SOxB <- otu_V[grepl(SOxB_list, otu_V$Consensus.lineage),] 
    SOxB <- Subs_to_DF(SOxB)
    SOxB["Guild"] <- "SOxB"                                                   # dim(SOxB); levels(SOxB$Genus); #unique(SOxB)# Desulf2; 
                                                                          # unsure if should add all Chromatiales...substantially increases OTU number
    #############################################
    # Iron Reducers
    #  Get FeRB, probably too limited of a list but taxonomy poorly known.  # Below list from Wikipedia...
    FeRB_list<-"(Geobacter|Shewanella|Thermoanaerobacter|Deferribacter|Geothrix|Albidiferax|Deferrisomatota)"     
    FeRB <- otu_V[grepl(FeRB_list, otu_V$Consensus.lineage),] 
    FeRB <- Subs_to_DF(FeRB)
    FeRB["Guild"] <- "FeRB"
    #dim(FeRB); levels(FeRB$Genus); unique(FeRB)# FeRB; 

    #############################################
    # Iron Oxidizers
    # From Kato et al. 2015 Front Microb, Illbert et al. 2013 Biochemica et Biophysica Acta, recheck against Hedrich 2011 Microbiol.
    FeOB_list<-"(Ferro|Ferro|Lepto|Metallo|Mariprofundus|Gallionella|Sideroxydans|Acidithiobacillus|Rhodopseudomonas|Sulfolobus)"                                 
    FeOB <- otu_V[grepl(FeOB_list, otu_V$Consensus.lineage),] 
    FeOB <- Subs_to_DF(FeOB)
    FeOB["Guild"] <- "FeOB"
    #dim(FeOB); levels(FeOB$Genus); unique(FeOB) # FeOB


    #############################################
    #############################################
    # Combine lists from each guild category
    Guild_fxn_taxa <- rbind(methano, methylo, nitroso, nitros2, anamox, SOxB, Syntroph, SRB, FeOB, FeRB)

    # Make Guilds factor 
    Guilds <- as.character(unique(Guild_fxn_taxa$"Guild"))                              # List Guilds in order of appearance
    Guild_fxn_taxa["Guild"]<-factor(Guild_fxn_taxa$"Guild", levels = Guilds)            # Guilds as factor, in order of appearance

    # Sort by Guild, Consensus.lineage
    Guild_fxn_taxa <-Guild_fxn_taxa[order(Guild_fxn_taxa$Guild, Guild_fxn_taxa['Consensus.lineage']),]
    #dim(Guild_fxn_taxa); levels(Guild_fxn_taxa$Guild); #unique(Guild_fxn_taxa); #head(Guild_fxn_taxa);

    # Add OTUs
    Guild_fxn_taxa["OTU"] <-row.names(Guild_fxn_taxa)

    # Get OTUs as column and extract OTU - Guild mapping for module export
    Guild_otu_cols <-c("Guild", "OTU")
    Guild_OTUs <- Guild_fxn_taxa[,Guild_otu_cols]
    # dim(Guild_OTUs); levels(Guild_OTUs$Guild); head(Guild_OTUs)
return(Guild_OTUs)
}


#### Test get guilds function
#Guild_OTUs <- Get_16S_Guilds(otu_V)
#head(Guild_OTUs)

# note that this only returns a DF of OTUs and Guilds, for OTUs that have guild assignments



##################################################
### b) Alternative version: Function to retrieve 16S Guilds using grep, with 4 different methanogen guilds
# Making alternate function since if this is used, much downstream will need to be updated like the colors (because more guilds)
# approximation of BGC function using organism names x Review Papers on function
# Update to separate some methanogens and ANME, see Chadwick et al. 2022
# List of taxa in Chadwick:
# Methanovorans 
# Methanogasteraceae/Methanogaster
# Methanocomedenaceae/Methanocomedens:Methanomarinus
# Methanophagales/Methanophagaceae/Methanophaga:Methanoalium,
# List of Chadwick taxa in SILVA 138.1:
# Methanoperedenaceae/Methanoperedens, Syntrophoarchaeaceae
# Fixed typo Archaeoglobus instead of Archeoglobus

Get_16S_Guilds_alt = function(otu_V){

    #############################################
    # Methanogens (remove ANME "methano"'s in subsequent section)

    methano <- otu_V[grepl("(Methano|Methermicoccaceae)", otu_V$Consensus.lineage),]                    # Get methanogens 
    methano <- Subs_to_DF(methano)                                                      # Fxn to get Taxon. DF with factors   
    # 2022 CHANGE HERE
    methano["Guild"] <- ifelse(methano$Family == 'Methanosarcinaceae' | methano$Family == 'Methanobacteriaceae', "CH4_mix",
                               ifelse(methano$Family == 'Methanosaetaceae', 'CH4_ac',
                                     ifelse(methano$Family == 'Methanotrichaceae', 'CH4_ac', 
                                           ifelse(methano$Order == 'Methanomassiliicoccales' | methano$Order == 'Methanonatronarchaeales', 'CH4_me', 
                                                  ifelse(methano$Family == 'Methermicoccaceae', 'CH4_me','CH4_H2')))))
    methano <- methano[!grepl("Syntrophoarchaeaceae", methano$Consensus.lineage),] # drop
    methano <- methano[!grepl("Methanoperedenaceae", methano$Consensus.lineage),] # drop
    methano <- methano[!grepl("ANME", methano$Consensus.lineage),] # drop
    
    # Separating the strictly acetoclastic families Methanosaetaceae and Methanotrichaceae from other mixotrophic Methanosarcinaceae
    # dim(methano); levels(methano$Genus); # unique(methano)  # methano

    #############################################
    # Methylotrophs 

    methylo <- otu_V[grepl("Methylo", otu_V$Consensus.lineage),]                                 # Get methylo from OTU table    
    methylo <- Subs_to_DF(methylo)                                                               # use function to get Taxon. DF with factors   

    # Use ifelse to make new column of methylo GUILDS
    # MOB IIa based on Kneif et al. 2015 Front. Microb, other fine differences ignored (e.g. IIb, Ia-c)
    methylo["Guild"] <- ifelse(methylo$Class=='Gammaproteobacteria', "MOB_I", 
                            #  ifelse(methylo$Family=='Methylocystaceae', "MOB_IIa",            # MOB_IIa called before...
                            ifelse(methylo$Family=='Beijerinckiaceae', "MOB_IIa",            # MOB_IIa called before...
                            
                            ### f__Beijerinckiaceae? or
                            # o__Rhizobiales;
                            
                            ifelse(methylo$Class=='Alphaproteobacteria', "MOB_II",           # all MOB II
                            ifelse(methylo$Class=='Betaproteobacteria', "MeOB", NA ))))      # "Methyoltrophs"

    methylo["Guild"]<-factor(methylo$Guild, levels = c("MOB_I","MOB_II", "MOB_IIa", "MeOB"))     # Reorder Guild factor
    methylo <- methylo[order(methylo$Guild),]                                                    # Sort by Guild
    #dim(methylo); levels(methylo$Genus); levels(methylo$Guild); unique(methylo); #methylo  #levels(MOB_II$Class);

    #############################################
    # ANME
    if (nrow(otu_V[grepl("(ANME|Methanoperedenaceae|Syntrophoarchaeaceae)", otu_V$Consensus.lineage),]) > 0) {
    ANME <- otu_V[grepl("(ANME|Methanoperedenaceae|Syntrophoarchaeaceae)", otu_V$Consensus.lineage),] # Get ANME from OTU table    
    ANME <- Subs_to_DF(ANME)                                                      # use function to get Taxon. DF with factors 
    ANME["Guild"] <- "ANME"
    #ANME
    }

    #############################################
    # Nitrifiers

    # Get ammonia oxidizers and nitrite oxidizers                                 # Slightly trickier than other groups
    nitros <- otu_V[grepl("Nitros", otu_V$Consensus.lineage),]   

    # Get only AOA, AOB (nitroso)
    nitroso <- nitros[grepl("Nitroso", nitros$Consensus.lineage),]       
    nitroso <- Subs_to_DF(nitroso)                                                # use function to get Taxon. DF with factors      
    nitroso["Guild"] <- ifelse(nitroso$Kingdom=='Archaea', "AOA", "AOB")          # Label AOA vs. AOB by Kingdom 
    # dim(nitroso); levels(nitroso$Genus);  #unique(nitroso)

    # Get NOB as opposite set from AOA/AOB, drop SRB
    nitros2 <- nitros[!grepl("Nitroso", nitros$Consensus.lineage),]       
    nitros2 <- nitros2[!grepl("desulfo", nitros2$Consensus.lineage),]             # Drop sulfate reducing Nitrospira (Not MOB), ambigs
    # nitros2 <- nitros2[!grepl("NitrospiralesOR", nitros2$Consensus.lineage),]   # Could drop ambig (Order, could be desulf)
                                                                              # But appears at least one is (+)~CH4
    nitros2 <- Subs_to_DF(nitros2)                                                # use function to get Taxon. DF with factors      
    nitros2["Guild"] <- "NOB"                                                     # dim(nitros2); levels(nitros2$Family);   #unique(nitros2) #nitros2# [1:40,]
    
    #############################################
    # ANAMMOX
    anamox_list <-"(Kuenenia|Anammoxoglobus|Scalindua|Brocadia|Jettenia)"
    if (nrow(otu_V[grepl(anamox_list, otu_V$Consensus.lineage),]) > 0) {
    anamox <- otu_V[grepl(anamox_list, otu_V$Consensus.lineage),] 
    anamox <- Subs_to_DF(anamox)
    anamox["Guild"] <-"Anamx"                                                     # levels(anamox$Genus); unique(anamox)# anamox; 
    }

    #############################################
    # Sulfate Reducers
    
    # Sourced from SRB reviews: Plugge et al. 2011 Front Microb; Muller et al. 2015 ISMEJ
    Desulf <- otu_V[grepl("Desulf", otu_V$Consensus.lineage),]               # Get Desulfo --ALL

    # Drop non SRB "Desulf"
    Desulf <- Desulf[!grepl("Nitros", Desulf$Consensus.lineage),]            # Drops Nitrospinaceae, in Order: Desulfobacterales 
    Desulf <-Desulf[!grepl("Syntroph",Desulf$Consensus.lineage),]            # Drops Syntrophs (own group), according to Plugge et al. 2011 Front. Microb.
    Desulf <-Desulf[!grepl("Geobacter",Desulf$Consensus.lineage),]           # Drops Geobacter (own group), according to Plugge et al. 2011 Front. Microb.
    Desulf <- Subs_to_DF(Desulf)                                             # dim(Desulf); levels(Desulf$Genus); # unique(Desulf)

    ##### Get misc. SRB (not Desulf), based on Trees in Muller et al. 2015 ISME
    SRB_xtra<-"(Caldiserica|Thermanaeromonas|Sporomusaceae|Carboxydothermus|Pelotomaculum|
                Moorella |Ammonifex|Acetonema|Thermosinus|Thermanaeromonas|Carboxydothermus|Caldiserica|
                Gordonibacter|Thermodesulfobium|Thermodesulfovibrio|Magnetobacterium)" 
    Desulf2 <- otu_V[grepl(SRB_xtra, otu_V$Consensus.lineage),]               # grep SRB_xtra 
    Desulf2 <- Subs_to_DF(Desulf2)                                            # dim(Desulf2); levels(Desulf2$Genus); #unique(Desulf2)# Desulf2; 

    # Sulfate reducing ARCHAEA, from Muller et al. 2015 ISME
    SRA_list<-"(Archaeoglobus|Pyrobaculum|Vulcanisaeta|Caldirvirga)"           # Sulfate red. archaea
    SRA <- otu_V[grepl(SRA_list, otu_V$Consensus.lineage),]                   # Get SRA
    SRA <- Subs_to_DF(SRA)                                                    # levels(SRA$Genus); unique(SRA)# SRA

    ## Combine SRB lists and label (note no SRA here so passing over)
    SRB <- rbind(Desulf, Desulf2)
    SRB["Guild"]<-"SRB"
    #dim(SRB); levels(SRB$Genus); #SRB

    ##### Split out syntrophs SEPARATELY
    # several non SRB, non-"Desulfo" syntrophs; eg. several ;o__Clostridiales;f__Syntrophomonadaceae; Delta D;f__Syntrophobacteraceae
    # remove Syntrophoarchaeaceae here so it retains ANME assignment
    Syntroph <- otu_V[grepl("Syntroph", otu_V$Consensus.lineage),]
    Syntroph <- Syntroph[!grepl("Syntrophoarchaeaceae", Syntroph$Consensus.lineage),] # drop
    Syntroph <- Subs_to_DF(Syntroph)
    Syntroph["Guild"] <-"SRB_syn"                                             # dim(Syntroph); levels(Syntroph$Genus); #unique(Syntroph)

    #############################################
    # Sulfur Oxidizers
    # While many are Thio, add rest from Muller et al. 2015 ISMEJ to same list
    SOxB_list<-"(Thio|Allochromatium|Marichromatium|Halochromatium|Alkalilimnicola|Halorhodospira|Ruthia|Vesicomyosocius|Sedimenticola
            Sideroxydans|Sulfuricella|Riegeria|Azospirillum|Rhodomicrobium|Magnetospirillum|Magnetococcus|
            Chlorobium|Chlorobaculum|Prosthecochloris)"
    # 1st line are Gamma prot, Fig. 3;     # 2nd line are B and A prot.;      #  3rd line Chlorobi

    SOxB <- otu_V[grepl(SOxB_list, otu_V$Consensus.lineage),] 
    SOxB <- Subs_to_DF(SOxB)
    SOxB["Guild"] <- "SOxB"                                                   # dim(SOxB); levels(SOxB$Genus); #unique(SOxB)# Desulf2; 
                                                                          # unsure if should add all Chromatiales...substantially increases OTU number
    #############################################
    # Iron Reducers
    #  Get FeRB, probably too limited of a list but taxonomy poorly known.  # Below list from Wikipedia...
    FeRB_list<-"(Geobacter|Shewanella|Thermoanaerobacter|Deferribacter|Geothrix|Albidiferax|Deferrisomatota)"     
    FeRB <- otu_V[grepl(FeRB_list, otu_V$Consensus.lineage),] 
    FeRB <- Subs_to_DF(FeRB)
    FeRB["Guild"] <- "FeRB"
    #dim(FeRB); levels(FeRB$Genus); unique(FeRB)# FeRB; 

    #############################################
    # Iron Oxidizers
    # From Kato et al. 2015 Front Microb, Illbert et al. 2013 Biochemica et Biophysica Acta, recheck against Hedrich 2011 Microbiol.
    FeOB_list<-"(Ferro|Ferro|Lepto|Metallo|Mariprofundus|Gallionella|Sideroxydans|Acidithiobacillus|Rhodopseudomonas|Sulfolobus)"                                 
    FeOB <- otu_V[grepl(FeOB_list, otu_V$Consensus.lineage),] 
    FeOB <- Subs_to_DF(FeOB)
    FeOB["Guild"] <- "FeOB"
    #dim(FeOB); levels(FeOB$Genus); unique(FeOB) # FeOB


    #############################################
    #############################################
    # Combine lists from each guild category
    Guild_fxn_taxa <- rbind(methano, methylo, ANME, nitroso, nitros2, anamox, SOxB, Syntroph, SRB, FeOB, FeRB)

    # Make Guilds factor 
    Guilds <- as.character(unique(Guild_fxn_taxa$"Guild"))                              # List Guilds in order of appearance
    Guild_fxn_taxa["Guild"]<-factor(Guild_fxn_taxa$"Guild", levels = Guilds)            # Guilds as factor, in order of appearance

    # Sort by Guild, Consensus.lineage
    Guild_fxn_taxa <-Guild_fxn_taxa[order(Guild_fxn_taxa$Guild, Guild_fxn_taxa['Consensus.lineage']),]
    #dim(Guild_fxn_taxa); levels(Guild_fxn_taxa$Guild); #unique(Guild_fxn_taxa); #head(Guild_fxn_taxa);

    # Add OTUs
    Guild_fxn_taxa["OTU"] <-row.names(Guild_fxn_taxa)

    # Get OTUs as column and extract OTU - Guild mapping for module export
    Guild_otu_cols <-c("Guild", "OTU")
    Guild_OTUs <- Guild_fxn_taxa[,Guild_otu_cols]
    # dim(Guild_OTUs); levels(Guild_OTUs$Guild); head(Guild_OTUs)
return(Guild_OTUs)
}


#### Test get guilds function
#Guild_OTUs <- Get_16S_Guilds(otu_V)
#head(Guild_OTUs)

# note that this only returns a DF of OTUs and Guilds, for OTUs that have guild assignments



#############################################
### c) Generate color pallete for Guilds  -- Deac
# Here manually


# REDO using HEX, following -- note colors listed, then their hex
# http://research.stowers.org/mcm/efg/R/Color/Chart/ColorChart.pdf

#reds <- c("tomato3","tomato4")                                      
reds <- c("#CD4F39","#8B3626")                                             # For Methanos 
#blues <-c("steelblue1", "steelblue3", "steelblue4", "slateblue1")          
blues <-c("#63B8FF", "#4F94CD", "#36648B", "#836FFF")                      # For MOB   
#greens <-c("lawngreen", "green3", "green4", "springgreen4")               
greens <-c("#7CFC00", "#00CD00", "#008B00", "#00FF7F")                     # For AOA,B etc.
#purps <- c("orchid3", "maroon3", "magenta4")                          
purps <- c("#CD69C9", "#CD2990", "#8B008B")                                # For SRB...
#orbrn <- c("darkorgange3", "saddlebrown")                          
orbrn <- c("#CD6600", "#8B4513")                                           # For FeOB/RB

Guild_pal0 <- c(reds, blues, greens, purps, orbrn)                         # Combine guild palletes 
# Guild_pal0

# Get list of Guilds
# Guilds <- levels(Guild_OTUs$Guild)

#### Fixed to keep color pallete running after test data commented out above ####
#dim(Guild_OTUs); levels(Guild_OTUs$Guild); head(Guild_OTUs)
Guilds <- c("CH4_H2", "CH4_ac", "MOB_I", "MOB_II", "MOB_IIa", "MeOB", "AOA", "AOB", "NOB", "Anamx", "SOxB", "SRB_syn", "SRB", "FeOB", "FeRB")

# Make palette with index
color <-as.character(Guild_pal0)                                                       # pal as char
Guild_pal <-cbind(Guilds, color)                                                       # Join Guilds, colors
                         
Guild_colors <-data.frame(Guild_pal)                                                   # Make DF
Guild_colors["G_index"]<-seq(1:nrow(Guild_pal))                                        # Make sorting index
Guild_colors["Guilds"] <- reorder(Guild_colors[,"Guilds"], Guild_colors[,"G_index"])   # Reorder Guild levels
Guild_colors["color"] <- reorder(Guild_colors[,"color"], Guild_colors[,"G_index"])     # Reorder color levels
colnames(Guild_colors)[1]<-"Guild"                                                     # Rename to Guild for merges          
# levels(Guild_colors$Guild); #levels(Guild_colors$color); Guild_colors

# Guild_colors
# WRITE TABLE NEEDED...
# write.table(Guild_colors, "Guild_color_palette.txt", sep="\t")
# OR wrap in function

## Show color key
Guild <- levels(Guild_colors$Guild)
Guild_counts <- c(rep(1, length(Guild)))
colors <- as.character(Guild_colors$color)

#options(repr.plot.width=1.5, repr.plot.height=6) 
#barplot(Guilds_counts, col=rev(Guild_pal), names.arg=rev(Guilds), horiz=TRUE, las=2)

#options(repr.plot.width=6, repr.plot.height=2)
#barplot(Guild_counts, col=colors, names.arg=Guild, las=2)                               # las = 2 is for rotating text labels 90 degrees



