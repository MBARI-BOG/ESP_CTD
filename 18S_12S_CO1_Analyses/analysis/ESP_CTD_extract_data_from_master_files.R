# Subset relevant samples from master files for ESP/CTD analysis

library(tidyverse)
library(phyloseq)
library(here)
library(readxl)

# Load data
##### COI #####

### Prep metadata
# Load master metadata
metadata_COI <- read.csv("/Users/markusmin/Documents/MBARI-2167/local/eDNA_meta_analysis_local_data/banzai_dada2/COI/for_MBON/raw_data/072720/COI_master_metadata.csv")

# Import 16S metadata (it doesn't matter which metadata, they're all the same samples) that contains information on matching samples
ESP_CTD_16S_metadata <- read.csv("/Users/markusmin/Documents/MBARI-2167/local/ESP_to_CTD_local_data/ESP_CTD_Data/banzai_dada2/16S/ESP_CTD_16S_metadata.csv",row.names = 1)

# Take out only the environmental (matched) samples
ESP_CTD_16S_envt_metadata <- subset(ESP_CTD_16S_metadata, !(ESP_CTD_16S_metadata$ESP_CTD_ID %in% c("")))

ESP_CTD_16S_envt_metadata$ESP_CTD_ID <- as.character(ESP_CTD_16S_envt_metadata$ESP_CTD_ID)
ESP_CTD_16S_envt_metadata %>% mutate(.,matching_ID = ifelse(SAMPLING_cruise == "CN18S",paste0("CN18S",as.numeric(sub(".{3}","",ESP_CTD_ID))),
                                                            ifelse(SAMPLING_cruise == "CN18F",paste0("CN18F",as.numeric(sub(".{3}","",ESP_CTD_ID))),"unknown"))) -> ESP_CTD_16S_envt_metadata
ESP_CTD_16S_envt_metadata$matching_ID

# Change sample name so it matches
ESP_CTD_16S_envt_metadata$sample_name <- as.character(ESP_CTD_16S_envt_metadata$sample_name)
# change "-" to "_"
ESP_CTD_16S_envt_metadata$sample_name <- str_replace_all(ESP_CTD_16S_envt_metadata$sample_name,"-","_")

# Join the two together
metadata_COI_ESP_CTD <- subset(metadata_COI, libraryID %in% c("FF","GG","Y"))
# Drop the plate suffix
metadata_COI_ESP_CTD$sample_name <- str_replace_all(metadata_COI_ESP_CTD$sample_name,"_FF","")
metadata_COI_ESP_CTD$sample_name <- str_replace_all(metadata_COI_ESP_CTD$sample_name,"_GG","")
metadata_COI_ESP_CTD$sample_name <- str_replace_all(metadata_COI_ESP_CTD$sample_name,"_Y","")

setdiff(ESP_CTD_16S_envt_metadata$sample_name,metadata_COI_ESP_CTD$sample_name)

# Now, let's extract the samples we want, based on the matching ESP/CTD sample, plus the blanks
metadata_COI_ESP_CTD_subset <- subset(metadata_COI_ESP_CTD, sample_name %in% ESP_CTD_16S_envt_metadata$sample_name | !(sample_type %in% c("positive","environmental") & !(str_detect(sample_name,"BP"))))
# Drop those weird ESP blanks - DONT! These are important
# metadata_COI_ESP_CTD_subset <- subset(metadata_COI_ESP_CTD_subset, !(sample_name %in% c("CN18FESPkoa_SC59","CN18FESPkoa_SC1","CN18S_koa3G_SC59","CN18S_koa3G_SC58","CN18S_koa3G_SC57","CN18S_koa3G_SC31","CN18S_koa3G_SC30","CN18S_koa3G_SC29")))
# # metadata_COI_ESP_CTD_subset$sample_name


## add in the metadata that we want about matching IDs
# Change sample name so it matches
ESP_CTD_16S_metadata$sample_name <- as.character(ESP_CTD_16S_metadata$sample_name)
# change "-" to "_"
ESP_CTD_16S_metadata$sample_name <- str_replace_all(ESP_CTD_16S_metadata$sample_name,"-","_")
# For pre/post blanks, change to include an underscore and "koa"
ESP_CTD_16S_metadata$sample_name <- str_replace_all(ESP_CTD_16S_metadata$sample_name,"CN18FB","CN18FESPkoa_B")
# change "pst" to "post
ESP_CTD_16S_metadata$sample_name <- str_replace_all(ESP_CTD_16S_metadata$sample_name,"pst","post")

# Add the matchingID column
ESP_CTD_16S_metadata$ESP_CTD_ID <- as.character(ESP_CTD_16S_metadata$ESP_CTD_ID)
ESP_CTD_16S_metadata %>% mutate(.,matching_ID = ifelse(!(ESP_CTD_ID %in% c("")), paste0(SAMPLING_cruise, as.numeric(sub(".{3}","",ESP_CTD_ID))),"")) -> ESP_CTD_16S_metadata

# Take out only the columns of interest from the 16S metadata
ESP_CTD_16S_metadata <- ESP_CTD_16S_metadata[,c("sample_name","sample_type","CTD_or_ESP","ESP_CTD","ESP_CTD_ID","matching_ID")]

# setdiff(metadata_COI_ESP_CTD$sample_name,ESP_CTD_16S_metadata$sample_name)
metadata_COI_ESP_CTD_subset <- left_join(metadata_COI_ESP_CTD_subset,ESP_CTD_16S_metadata, by = "sample_name")

# Add the plate suffix back in so the phyloseq object will merge
metadata_COI_ESP_CTD_subset %>% mutate(.,sample_name = paste0(sample_name,"_",libraryID)) -> metadata_COI_ESP_CTD_subset

# Drop all of the metadata columns that we don't want
metadata_COI_ESP_CTD_subset <- metadata_COI_ESP_CTD_subset[c("sample_name","SAMPLING_cruise","depth","library","libraryID","locus","original_name","samp_collection_device","sample_type.x","sample_type.y","CTD_or_ESP","ESP_CTD","ESP_CTD_ID","matching_ID")]

# Update column to identify ESP and CTD samples
metadata_COI_ESP_CTD_subset$CTD_or_ESP <- as.character(metadata_COI_ESP_CTD_subset$CTD_or_ESP)
metadata_COI_ESP_CTD_subset %>% mutate(.,CTD_or_ESP = ifelse(CTD_or_ESP %in% c(""),gsub('[0-9]+', '', ESP_CTD_ID),CTD_or_ESP)) -> metadata_COI_ESP_CTD_subset
metadata_COI_ESP_CTD_subset %>% mutate(.,CTD_or_ESP = ifelse(CTD_or_ESP == "ESP ","ESP",CTD_or_ESP)) -> metadata_COI_ESP_CTD_subset


### Make some tweaks to tax table
tax_table_COI <- read.csv("/Users/markusmin/Documents/MBARI-2167/local/eDNA_meta_analysis_local_data/banzai_dada2/COI/for_MBON/raw_data/072720/COI_master_tax_table.csv")

tax_table_COI %>% mutate_all(as.character) -> tax_table_COI

subset(tax_table_COI,Phylum == "unknown")

# Fill in the Phylum level for some common taxa
tax_table_COI %>% mutate(.,Phylum = ifelse(Class == "Haptophyta","Haptophyta",
                                           ifelse(Class == "Oomycetes","Oomycota",
                                                  ifelse(Class %in% c("Phaeophyceae","Chrysophyceae","Dictyochophyceae"),"Ochrophyta",Phylum)))) -> tax_table_COI

tax_table_COI <- tibble::column_to_rownames(tax_table_COI,"ASV")


### Create phyloseq object

# Load in master files
asv_table_COI <- read.csv("/Users/markusmin/Documents/MBARI-2167/local/eDNA_meta_analysis_local_data/banzai_dada2/COI/for_MBON/raw_data/072720/COI_master_ASV_table.csv",row.names = 1)
# DESTROY X!!!
destroyX = function(es) {
  f = es
  for (col in c(1:ncol(f))){ #for each column in dataframe
    if (startsWith(colnames(f)[col], "X") == TRUE)  { #if starts with 'X' ..
      colnames(f)[col] <- substr(colnames(f)[col], 2, 100) #get rid of it
    }
  }
  assign(deparse(substitute(es)), f, inherits = TRUE) #assign corrected data to original name
}
asv_table_COI <- destroyX(asv_table_COI)


# Change metadata column to rownames
metadata_COI_ESP_CTD_subset %>% tibble::column_to_rownames(.,"sample_name") -> metadata_COI_ESP_CTD_subset_for_phyloseq

setdiff(rownames(metadata_COI_ESP_CTD_subset_for_phyloseq),colnames(asv_table_COI))

ESP_CTD_COI_phyloseq <- merge_phyloseq(otu_table(asv_table_COI,taxa_are_rows = TRUE),tax_table(as.matrix(tax_table_COI)),sample_data(metadata_COI_ESP_CTD_subset_for_phyloseq))
ESP_CTD_COI_phyloseq <- prune_taxa(taxa_sums(ESP_CTD_COI_phyloseq) > 0, ESP_CTD_COI_phyloseq)
# ESP_CTD_COI_phyloseq <- prune_samples(sample_data(ESP_CTD_COI_phyloseq)$CTD_or_ESP %in% c("ESP","CTD"),ESP_CTD_COI_phyloseq)

sample_names(ESP_CTD_COI_phyloseq)

# Export three datafiles (asv_table, tax_table, and metadata) for other scripts.
ESP_CTD_COI_tax_table <- as.data.frame(tax_table(ESP_CTD_COI_phyloseq))
write.csv(ESP_CTD_COI_tax_table, here("data", "COI", "ESP_CTD_COI_tax_table.csv"))

ESP_CTD_COI_asv_table <- as.data.frame(otu_table(ESP_CTD_COI_phyloseq))
write.csv(ESP_CTD_COI_asv_table, here("data", "COI", "ESP_CTD_COI_asv_table.csv"))

ESP_CTD_COI_metadata <- as.data.frame(sample_data(ESP_CTD_COI_phyloseq))
write.csv(ESP_CTD_COI_metadata, here("data", "COI", "ESP_CTD_COI_metadata.csv"))


##### COI: Subset and export direct comparison samples #####

# Load the sample names from Katie's file, edit to match those in master files
read_excel(here("data", "Katie_ESP_CTD_matched_samples_metadata.xlsx"), sheet = "CN17S") %>% 
  subset(., Type %in% c("combination_depth", "bench_ESP")) %>% 
  mutate(COI_name = paste0(sample_name, "_X")) -> match_metadata

# Subset the samples from the master metadata, move sample names to row names
subset(metadata_COI, sample_name %in% match_metadata$COI_name) -> matched_metadata_COI
rownames(matched_metadata_COI) <- NULL
matched_metadata_COI %>% 
  column_to_rownames("sample_name") -> matched_metadata_COI

# Join as phyloseq
matched_COI_phyloseq <- merge_phyloseq(otu_table(asv_table_COI,taxa_are_rows = TRUE),tax_table(as.matrix(tax_table_COI)),sample_data(matched_metadata_COI))
matched_COI_phyloseq <- prune_taxa(taxa_sums(matched_COI_phyloseq) > 0, matched_COI_phyloseq)

# Export three datafiles (asv_table, tax_table, and metadata) for other scripts.
matched_ESP_CTD_COI_tax_table <- as.data.frame(tax_table(matched_COI_phyloseq))
write.csv(matched_ESP_CTD_COI_tax_table, here("data", "COI", "direct_comparisons", "directcomp_ESP_CTD_COI_tax_table.csv"))

matched_ESP_CTD_COI_asv_table <- as.data.frame(otu_table(matched_COI_phyloseq))
write.csv(matched_ESP_CTD_COI_asv_table, here("data", "COI", "direct_comparisons", "directcomp_ESP_CTD_COI_asv_table.csv"))

matched_ESP_CTD_COI_metadata <- as.data.frame(sample_data(matched_COI_phyloseq))
write.csv(matched_ESP_CTD_COI_metadata, here("data", "COI", "direct_comparisons", "directcomp_ESP_CTD_COI_metadata.csv"))

##### 18S #####

### Prep metadata
# Load master metadata
metadata_18S <- read.csv("/Users/markusmin/Documents/MBARI-2167/local/eDNA_meta_analysis_local_data/banzai_dada2/18S/for_MBON/raw_data/072720/18S_master_metadata.csv")

# Import 16S metadata (it doesn't matter which metadata, they're all the same samples) that contains information on matching samples
ESP_CTD_16S_metadata <- read.csv("/Users/markusmin/Documents/MBARI-2167/local/ESP_to_CTD_local_data/ESP_CTD_Data/banzai_dada2/16S/ESP_CTD_16S_metadata.csv",row.names = 1)

# Take out only the environmental (matched) samples
ESP_CTD_16S_envt_metadata <- subset(ESP_CTD_16S_metadata, !(ESP_CTD_16S_metadata$ESP_CTD_ID %in% c("")))

ESP_CTD_16S_envt_metadata$ESP_CTD_ID <- as.character(ESP_CTD_16S_envt_metadata$ESP_CTD_ID)
ESP_CTD_16S_envt_metadata %>% mutate(.,matching_ID = ifelse(SAMPLING_cruise == "CN18S",paste0("CN18S",as.numeric(sub(".{3}","",ESP_CTD_ID))),
                                                            ifelse(SAMPLING_cruise == "CN18F",paste0("CN18F",as.numeric(sub(".{3}","",ESP_CTD_ID))),"unknown"))) -> ESP_CTD_16S_envt_metadata
ESP_CTD_16S_envt_metadata$matching_ID

# Change sample name so it matches
ESP_CTD_16S_envt_metadata$sample_name <- as.character(ESP_CTD_16S_envt_metadata$sample_name)
# change "-" to "_"
ESP_CTD_16S_envt_metadata$sample_name <- str_replace_all(ESP_CTD_16S_envt_metadata$sample_name,"-","_")

# Join the two together
metadata_18S_ESP_CTD <- subset(metadata_18S, libraryID %in% c("AA","CC","HH"))
# Drop the plate suffix
metadata_18S_ESP_CTD$sample_name <- str_replace_all(metadata_18S_ESP_CTD$sample_name,"_AA","")
metadata_18S_ESP_CTD$sample_name <- str_replace_all(metadata_18S_ESP_CTD$sample_name,"_CC","")
metadata_18S_ESP_CTD$sample_name <- str_replace_all(metadata_18S_ESP_CTD$sample_name,"_HH","")

setdiff(ESP_CTD_16S_envt_metadata$sample_name,metadata_18S_ESP_CTD$sample_name)

# Now, let's extract the samples we want, based on the matching ESP/CTD sample, plus the blanks
metadata_18S_ESP_CTD_subset <- subset(metadata_18S_ESP_CTD, sample_name %in% ESP_CTD_16S_envt_metadata$sample_name | !(sample_type %in% c("positive","environmental") & !(str_detect(sample_name,"BP"))))
# Drop those weird ESP blanks - DON'T - these are important!
# metadata_18S_ESP_CTD_subset <- subset(metadata_18S_ESP_CTD_subset, !(sample_name %in% c("CN18FESPkoa_SC59","CN18FESPkoa_SC1","CN18S_koa3G_SC59","CN18S_koa3G_SC58","CN18S_koa3G_SC57","CN18S_koa3G_SC31","CN18S_koa3G_SC30","CN18S_koa3G_SC29")))
# metadata_18S_ESP_CTD_subset$sample_name


## add in the metadata that we want about matching IDs
# Change sample name so it matches
ESP_CTD_16S_metadata$sample_name <- as.character(ESP_CTD_16S_metadata$sample_name)
# change "-" to "_"
ESP_CTD_16S_metadata$sample_name <- str_replace_all(ESP_CTD_16S_metadata$sample_name,"-","_")
# For pre/post blanks, change to include an underscore and "koa"
ESP_CTD_16S_metadata$sample_name <- str_replace_all(ESP_CTD_16S_metadata$sample_name,"CN18FB","CN18FESPkoa_B")
# change "pst" to "post
ESP_CTD_16S_metadata$sample_name <- str_replace_all(ESP_CTD_16S_metadata$sample_name,"pst","post")

# Add the matchingID column
ESP_CTD_16S_metadata$ESP_CTD_ID <- as.character(ESP_CTD_16S_metadata$ESP_CTD_ID)
ESP_CTD_16S_metadata %>% mutate(.,matching_ID = ifelse(!(ESP_CTD_ID %in% c("")), paste0(SAMPLING_cruise, as.numeric(sub(".{3}","",ESP_CTD_ID))),"")) -> ESP_CTD_16S_metadata

# Take out only the columns of interest from the 16S metadata
ESP_CTD_16S_metadata <- ESP_CTD_16S_metadata[,c("sample_name","sample_type","CTD_or_ESP","ESP_CTD","ESP_CTD_ID","matching_ID")]

setdiff(metadata_18S_ESP_CTD_subset$sample_name,ESP_CTD_16S_metadata$sample_name)
metadata_18S_ESP_CTD_subset <- left_join(metadata_18S_ESP_CTD_subset,ESP_CTD_16S_metadata, by = "sample_name")

# Add the plate suffix back in so the phyloseq object will merge
metadata_18S_ESP_CTD_subset %>% mutate(.,sample_name = paste0(sample_name,"_",libraryID)) -> metadata_18S_ESP_CTD_subset

# Drop all of the metadata columns that we don't want
metadata_18S_ESP_CTD_subset <- metadata_18S_ESP_CTD_subset[c("sample_name","SAMPLING_cruise","depth","library","libraryID","locus","original_name","samp_collection_device","sample_type.x","sample_type.y","CTD_or_ESP","ESP_CTD","ESP_CTD_ID","matching_ID")]

# Update column to identify ESP and CTD samples
metadata_18S_ESP_CTD_subset$CTD_or_ESP <- as.character(metadata_18S_ESP_CTD_subset$CTD_or_ESP)
metadata_18S_ESP_CTD_subset %>% mutate(.,CTD_or_ESP = ifelse(CTD_or_ESP %in% c(""),gsub('[0-9]+', '', ESP_CTD_ID),CTD_or_ESP)) -> metadata_18S_ESP_CTD_subset
metadata_18S_ESP_CTD_subset %>% mutate(.,CTD_or_ESP = ifelse(CTD_or_ESP == "ESP ","ESP",CTD_or_ESP)) -> metadata_18S_ESP_CTD_subset

### Make some tweaks to tax table

tax_table_18S <- read.csv("/Users/markusmin/Documents/MBARI-2167/local/eDNA_meta_analysis_local_data/banzai_dada2/18S/for_MBON/raw_data/072720/18S_master_tax_table.csv")

tax_table_18S %>% mutate_all(as.character) -> tax_table_18S

subset(tax_table_18S,Phylum == "unknown")

# Fill in the Phylum level for some common taxa
tax_table_18S %>% mutate(.,Phylum = ifelse(Class == "Haptophyta","Haptophyta",
                                           ifelse(Class == "Dinophyceae","Dinoflagellata",
                                                  ifelse(Class %in% c("Phaeophyceae","Chrysophyceae","Dictyochophyceae"),"Ochrophyta",
                                                         ifelse(Order == "Diplonemea","Euglenozoa",
                                                                Phylum))))) -> tax_table_18S

tax_table_18S <- tibble::column_to_rownames(tax_table_18S,"ASV")


### Create phyloseq object

# Load in master files
asv_table_18S <- read.csv("/Users/markusmin/Documents/MBARI-2167/local/eDNA_meta_analysis_local_data/banzai_dada2/18S/for_MBON/raw_data/072720/18S_master_ASV_table.csv",row.names = 1)
# DESTROY X!!!
destroyX = function(es) {
  f = es
  for (col in c(1:ncol(f))){ #for each column in dataframe
    if (startsWith(colnames(f)[col], "X") == TRUE)  { #if starts with 'X' ..
      colnames(f)[col] <- substr(colnames(f)[col], 2, 100) #get rid of it
    }
  }
  assign(deparse(substitute(es)), f, inherits = TRUE) #assign corrected data to original name
}
asv_table_18S <- destroyX(asv_table_18S)

# Change metadata column to rownames
metadata_18S_ESP_CTD_subset %>% tibble::column_to_rownames(.,"sample_name") -> metadata_18S_ESP_CTD_subset_for_phyloseq

setdiff(rownames(metadata_18S_ESP_CTD_subset_for_phyloseq),colnames(asv_table_18S))

ESP_CTD_18S_phyloseq <- merge_phyloseq(otu_table(asv_table_18S,taxa_are_rows = TRUE),tax_table(as.matrix(tax_table_18S)),sample_data(metadata_18S_ESP_CTD_subset_for_phyloseq))
ESP_CTD_18S_phyloseq <- prune_taxa(taxa_sums(ESP_CTD_18S_phyloseq) > 0, ESP_CTD_18S_phyloseq)
# ESP_CTD_18S_phyloseq <- prune_samples(sample_data(ESP_CTD_18S_phyloseq)$CTD_or_ESP %in% c("ESP","CTD"),ESP_CTD_18S_phyloseq)

# sample_names(ESP_CTD_18S_phyloseq)


# Export three datafiles (asv_table, tax_table, and metadata) for other scripts.
ESP_CTD_18S_tax_table <- as.data.frame(tax_table(ESP_CTD_18S_phyloseq))
write.csv(ESP_CTD_18S_tax_table, here("data", "18S", "ESP_CTD_18S_tax_table.csv"))

ESP_CTD_18S_asv_table <- as.data.frame(otu_table(ESP_CTD_18S_phyloseq))
write.csv(ESP_CTD_18S_asv_table, here("data", "18S", "ESP_CTD_18S_asv_table.csv"))

ESP_CTD_18S_metadata <- as.data.frame(sample_data(ESP_CTD_18S_phyloseq))
write.csv(ESP_CTD_18S_metadata, here("data", "18S", "ESP_CTD_18S_metadata.csv"))




##### 12S #####

### Prep metadata

# Load master metadata
metadata_12S <- read.csv("/Users/markusmin/Documents/MBARI-2167/local/eDNA_meta_analysis_local_data/banzai_dada2/12S/for_MBON/raw_data/072720/12S_master_metadata.csv")

# Import 16S metadata (it doesn't matter which metadata, they're all the same samples) that contains information on matching samples
ESP_CTD_16S_metadata <- read.csv("/Users/markusmin/Documents/MBARI-2167/local/ESP_to_CTD_local_data/ESP_CTD_Data/banzai_dada2/16S/ESP_CTD_16S_metadata.csv",row.names = 1)

# Take out only the environmental (matched) samples
ESP_CTD_16S_envt_metadata <- subset(ESP_CTD_16S_metadata, !(ESP_CTD_16S_metadata$ESP_CTD_ID %in% c("")))

ESP_CTD_16S_envt_metadata$ESP_CTD_ID <- as.character(ESP_CTD_16S_envt_metadata$ESP_CTD_ID)
ESP_CTD_16S_envt_metadata %>% mutate(.,matching_ID = ifelse(SAMPLING_cruise == "CN18S",paste0("CN18S",as.numeric(sub(".{3}","",ESP_CTD_ID))),
                                                            ifelse(SAMPLING_cruise == "CN18F",paste0("CN18F",as.numeric(sub(".{3}","",ESP_CTD_ID))),"unknown"))) -> ESP_CTD_16S_envt_metadata
ESP_CTD_16S_envt_metadata$matching_ID

# Change sample name so it matches
ESP_CTD_16S_envt_metadata$sample_name <- as.character(ESP_CTD_16S_envt_metadata$sample_name)
# change "-" to "_"
ESP_CTD_16S_envt_metadata$sample_name <- str_replace_all(ESP_CTD_16S_envt_metadata$sample_name,"-","_")
# For the KOA samples on plate NN, change name
ESP_CTD_16S_envt_metadata$sample_name <- str_replace_all(ESP_CTD_16S_envt_metadata$sample_name,"CN18SESPkoa_SC","CN18SKOA_ESP")
# CN18SESPkoa_SC49 to CN18SKOA_ESP49

# Take out our three plates of interest
metadata_12S_ESP_CTD <- subset(metadata_12S, libraryID %in% c("QQ","OO","NN"))

# Remove the four duplicated samples
metadata_12S_ESP_CTD <- subset(metadata_12S_ESP_CTD, !(sample_name %in% c("CN18Sc14_8_eDNA_td_QQ","CN18Sc15_8_eDNA_td_QQ","CN18Sc18_2_eDNA_td_QQ","CN18Sc27_2_eDNA_td_QQ")))
# Drop the classic replicates as well
metadata_12S_ESP_CTD <- subset(metadata_12S_ESP_CTD, !(sample_name %in% c("CN18Sc14_8_eDNA_QQ","CN18Sc15_8_eDNA_QQ","CN18Sc18_2_eDNA_QQ","CN18Sc27_2_eDNA_QQ")))

# metadata_12S_ESP_CTD_subset$sample_name

# Drop the plate suffix
metadata_12S_ESP_CTD$sample_name <- str_replace_all(metadata_12S_ESP_CTD$sample_name,"_QQ","")
metadata_12S_ESP_CTD$sample_name <- str_replace_all(metadata_12S_ESP_CTD$sample_name,"_OO","")
metadata_12S_ESP_CTD$sample_name <- str_replace_all(metadata_12S_ESP_CTD$sample_name,"_NN","")

# Now, let's extract the samples we want, based on the matching ESP/CTD sample, plus the blanks
metadata_12S_ESP_CTD_subset <- subset(metadata_12S_ESP_CTD, sample_name %in% ESP_CTD_16S_envt_metadata$sample_name | sample_name == "CN18FESPkoa_SC1" | !(sample_type %in% c("positive","environmental") & !(str_detect(sample_name,"BP"))))
# Drop those weird ESP blanks - DON'T - these are important!
# metadata_12S_ESP_CTD_subset <- subset(metadata_12S_ESP_CTD_subset, !(sample_name %in% c("CN18FESPkoa_SC59","CN18FESPkoa_SC1","CN18S_koa3G_SC59","CN18S_koa3G_SC58","CN18S_koa3G_SC57","CN18S_koa3G_SC31","CN18S_koa3G_SC30","CN18S_koa3G_SC29")))


## add in the metadata that we want about matching IDs
# Change sample name so it matches
ESP_CTD_16S_metadata$sample_name <- as.character(ESP_CTD_16S_metadata$sample_name)
# change "-" to "_"
ESP_CTD_16S_metadata$sample_name <- str_replace_all(ESP_CTD_16S_metadata$sample_name,"-","_")
# For pre/post blanks, change to include an underscore and "koa"
ESP_CTD_16S_metadata$sample_name <- str_replace_all(ESP_CTD_16S_metadata$sample_name,"CN18FB","CN18FESPkoa_B")
# change "pst" to "post
ESP_CTD_16S_metadata$sample_name <- str_replace_all(ESP_CTD_16S_metadata$sample_name,"pst","post")
# For the KOA samples on plate NN, change name
ESP_CTD_16S_metadata$sample_name <- str_replace_all(ESP_CTD_16S_metadata$sample_name,"CN18SESPkoa_SC","CN18SKOA_ESP")

# Add the matchingID column
ESP_CTD_16S_metadata$ESP_CTD_ID <- as.character(ESP_CTD_16S_metadata$ESP_CTD_ID)
ESP_CTD_16S_metadata %>% mutate(.,matching_ID = ifelse(!(ESP_CTD_ID %in% c("")), paste0(SAMPLING_cruise, as.numeric(sub(".{3}","",ESP_CTD_ID))),"")) -> ESP_CTD_16S_metadata

# Take out only the columns of interest from the 16S metadata
ESP_CTD_16S_metadata <- ESP_CTD_16S_metadata[,c("sample_name","sample_type","CTD_or_ESP","ESP_CTD","ESP_CTD_ID","matching_ID")]

setdiff(metadata_12S_ESP_CTD_subset$sample_name,ESP_CTD_16S_metadata$sample_name)
metadata_12S_ESP_CTD_subset <- left_join(metadata_12S_ESP_CTD_subset,ESP_CTD_16S_metadata, by = "sample_name")

# Add the plate suffix back in so the phyloseq object will merge
metadata_12S_ESP_CTD_subset %>% mutate(.,sample_name = paste0(sample_name,"_",libraryID)) -> metadata_12S_ESP_CTD_subset

# Drop all of the metadata columns that we don't want
metadata_12S_ESP_CTD_subset <- metadata_12S_ESP_CTD_subset[c("sample_name","SAMPLING_cruise","depth","library","libraryID","locus","original_name","samp_collection_device","sample_type.x","sample_type.y","CTD_or_ESP","ESP_CTD","ESP_CTD_ID","matching_ID")]

# Update column to identify ESP and CTD samples
metadata_12S_ESP_CTD_subset$CTD_or_ESP <- as.character(metadata_12S_ESP_CTD_subset$CTD_or_ESP)
metadata_12S_ESP_CTD_subset %>% mutate(.,CTD_or_ESP = ifelse(CTD_or_ESP %in% c(""),gsub('[0-9]+', '', ESP_CTD_ID),CTD_or_ESP)) -> metadata_12S_ESP_CTD_subset
metadata_12S_ESP_CTD_subset %>% mutate(.,CTD_or_ESP = ifelse(CTD_or_ESP == "ESP ","ESP",CTD_or_ESP)) -> metadata_12S_ESP_CTD_subset

sort(subset(metadata_12S_ESP_CTD_subset, !(is.na(matching_ID) | matching_ID %in% c("")))$matching_ID)
subset(metadata_12S_ESP_CTD_subset, is.na(matching_ID) | matching_ID %in% c(""))

subset(metadata_12S_ESP_CTD_subset,CTD_or_ESP == "CTD")
subset(metadata_12S_ESP_CTD_subset,CTD_or_ESP == "ESP")

metadata_12S_ESP_CTD_subset

# ESP: CN18F 4, 5, 21, 22


### Create phyloseq object

# Load in master files
asv_table_12S <- read.csv("/Users/markusmin/Documents/MBARI-2167/local/eDNA_meta_analysis_local_data/banzai_dada2/12S/for_MBON/raw_data/072720/12S_master_ASV_table.csv",row.names = 1)
# DESTROY X!!!
destroyX = function(es) {
  f = es
  for (col in c(1:ncol(f))){ #for each column in dataframe
    if (startsWith(colnames(f)[col], "X") == TRUE)  { #if starts with 'X' ..
      colnames(f)[col] <- substr(colnames(f)[col], 2, 100) #get rid of it
    }
  }
  assign(deparse(substitute(es)), f, inherits = TRUE) #assign corrected data to original name
}
asv_table_12S <- destroyX(asv_table_12S)
tax_table_12S <- read.csv("/Users/markusmin/Documents/MBARI-2167/local/eDNA_meta_analysis_local_data/banzai_dada2/12S/for_MBON/raw_data/072720/12S_master_tax_table.csv",row.names = 1)


# Change metadata column to rownames
metadata_12S_ESP_CTD_subset %>% tibble::column_to_rownames(.,"sample_name") -> metadata_12S_ESP_CTD_subset_for_phyloseq

setdiff(rownames(metadata_12S_ESP_CTD_subset_for_phyloseq),colnames(asv_table_12S))

ESP_CTD_12S_phyloseq <- merge_phyloseq(otu_table(asv_table_12S,taxa_are_rows = TRUE),tax_table(as.matrix(tax_table_12S)),sample_data(metadata_12S_ESP_CTD_subset_for_phyloseq))
ESP_CTD_12S_phyloseq <- prune_taxa(taxa_sums(ESP_CTD_12S_phyloseq) > 0, ESP_CTD_12S_phyloseq)
# ESP_CTD_12S_phyloseq <- prune_samples(sample_data(ESP_CTD_12S_phyloseq)$CTD_or_ESP %in% c("ESP","CTD"),ESP_CTD_12S_phyloseq)

# sample_names(ESP_CTD_12S_phyloseq)
# sum(sample_sums(ESP_CTD_12S_phyloseq))
# sort(sample_names(ESP_CTD_12S_phyloseq))

# Take out only chordates
ESP_CTD_12S_phyloseq_chordates <- subset_taxa(ESP_CTD_12S_phyloseq,Phylum == "Chordata")

# Export three datafiles (asv_table, tax_table, and metadata) for other scripts.
# Export the non-subset datafiles
ESP_CTD_12S_tax_table <- as.data.frame(tax_table(ESP_CTD_12S_phyloseq))
write.csv(ESP_CTD_12S_tax_table, here("data", "12S", "ESP_CTD_12S_tax_table.csv"))

ESP_CTD_12S_asv_table <- as.data.frame(otu_table(ESP_CTD_12S_phyloseq))
write.csv(ESP_CTD_12S_asv_table, here("data", "12S", "ESP_CTD_12S_asv_table.csv"))

ESP_CTD_12S_metadata <- as.data.frame(sample_data(ESP_CTD_12S_phyloseq))
write.csv(ESP_CTD_12S_metadata, here("data", "12S", "ESP_CTD_12S_metadata.csv"))


##### 16S #####

### Prep metadata

# Load master metadata
metadata_16S <- read.csv("/Users/markusmin/Documents/MBARI-2167/local/ESP_to_CTD_local_data/ESP_CTD_Data/banzai_dada2/ESP_CTD_16S/MB_20200625_1614_16S_analysis_metadata.csv")

metadata_16S_ESP_CTD_subset <- metadata_16S

# Add the matchingID column
metadata_16S_ESP_CTD_subset$ESP_CTD_ID <- as.character(metadata_16S_ESP_CTD_subset$ESP_CTD_ID)
metadata_16S_ESP_CTD_subset %>% mutate(.,matching_ID = ifelse(!(ESP_CTD_ID %in% c("")), paste0(SAMPLING_cruise, as.numeric(sub(".{3}","",ESP_CTD_ID))),"")) -> metadata_16S_ESP_CTD_subset

# The one sample that's missing is the one where the fasta files were empty

# Update column to identify ESP and CTD samples
metadata_16S_ESP_CTD_subset$CTD_or_ESP <- as.character(metadata_16S_ESP_CTD_subset$CTD_or_ESP)
metadata_16S_ESP_CTD_subset %>% mutate(.,CTD_or_ESP = ifelse(CTD_or_ESP %in% c(""),gsub('[0-9]+', '', ESP_CTD_ID),CTD_or_ESP)) -> metadata_16S_ESP_CTD_subset
metadata_16S_ESP_CTD_subset %>% mutate(.,CTD_or_ESP = ifelse(CTD_or_ESP == "ESP ","ESP",CTD_or_ESP)) -> metadata_16S_ESP_CTD_subset

sort(subset(metadata_16S_ESP_CTD_subset, !(is.na(matching_ID) | matching_ID %in% c("")))$matching_ID)
subset(metadata_16S_ESP_CTD_subset, is.na(matching_ID) | matching_ID %in% c(""))

subset(metadata_16S_ESP_CTD_subset,CTD_or_ESP == "CTD")
subset(metadata_16S_ESP_CTD_subset,CTD_or_ESP == "ESP")

# Change metadata names so that they match
metadata_16S_ESP_CTD_subset$sample_name <- str_replace_all(metadata_16S_ESP_CTD_subset$sample_name,"-","_")

# ESP: CN18F 4, 5, 21, 22

write.csv(metadata_16S_ESP_CTD_subset, here("data", "16S", "ESP_CTD_16S_metadata.csv"))




