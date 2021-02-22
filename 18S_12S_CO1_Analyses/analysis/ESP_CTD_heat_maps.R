## Heat map figures

library(dplyr)
library(ggplot2)
library(knitr)
library(shadowtext)
library(compositions)
library(zCompositions)
library(phyloseq)
library(reshape2)
library(ggplot2)
library(stringr)
library(forcats)
library(dichromat)
library(Rmisc)
library(superheat)
library(viridis)
library(cowplot)
library(tidyverse)
library(here)
library(ggpubr)
library(gridExtra)
library(egg)

## ----Load data----------------------------------------------------------------

#### COI
asv_table_path_COI <- here::here("data", "COI", "ESP_CTD_COI_asv_table.csv")
tax_table_path_COI <- here::here("data", "COI", "ESP_CTD_COI_tax_table.csv")
metadata_path_COI <- here::here("data", "COI", "ESP_CTD_COI_metadata.csv")

ESP_CTD_COI_phyloseq <- merge_phyloseq(otu_table(read.csv(asv_table_path_COI, row.names = 1),taxa_are_rows = TRUE),
                                       tax_table(as.matrix(read.csv(tax_table_path_COI, row.names = 1))),
                                       sample_data(read.csv(metadata_path_COI, row.names = 1)))

#### 18S
asv_table_path_18S <- here::here("data", "18S", "ESP_CTD_18S_asv_table.csv")
tax_table_path_18S <- here::here("data", "18S", "ESP_CTD_18S_tax_table.csv")
metadata_path_18S <- here::here("data", "18S", "ESP_CTD_18S_metadata.csv")

ESP_CTD_18S_phyloseq <- merge_phyloseq(otu_table(read.csv(asv_table_path_18S, row.names = 1),taxa_are_rows = TRUE),
                                       tax_table(as.matrix(read.csv(tax_table_path_18S, row.names = 1))),
                                       sample_data(read.csv(metadata_path_18S, row.names = 1)))

#### 12S
asv_table_path_12S <- here::here("data", "12S", "ESP_CTD_12S_asv_table.csv")
tax_table_path_12S <- here::here("data", "12S", "ESP_CTD_12S_tax_table.csv")
metadata_path_12S <- here::here("data", "12S", "ESP_CTD_12S_metadata.csv")

ESP_CTD_12S_phyloseq <- merge_phyloseq(otu_table(read.csv(asv_table_path_12S, row.names = 1),taxa_are_rows = TRUE),
                                       tax_table(as.matrix(read.csv(tax_table_path_12S, row.names = 1))),
                                       sample_data(read.csv(metadata_path_12S, row.names = 1)))

### 16S
# Read in combined taxonomy + ASV table
asv_tax_combined_path_16S <- here::here("data", "16S", "ESP_CTD_16S_dada2_table_with_tax.csv")
asv_tax_combined_16S <- read.csv(asv_tax_combined_path_16S)
# Split off asv table
dplyr::select(asv_tax_combined_16S, -Taxon) %>% 
  column_to_rownames("FeatureID") -> asv_table_16S
# Split off taxonomy, separate into fields per rank
tax_table_16S <- dplyr::select(asv_tax_combined_16S, FeatureID, Taxon)
tax_table_16S$Taxon <- gsub("[a-zA-Z]__","", tax_table_16S$Taxon)
tax_table_16S %>% separate(., Taxon, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";") %>% 
  column_to_rownames("FeatureID") -> tax_table_16S

metadata_path_16S <- here::here("data", "16S", "ESP_CTD_16S_metadata.csv")

ESP_CTD_16S_phyloseq <- merge_phyloseq(otu_table(asv_table_16S,taxa_are_rows = TRUE),
                                       tax_table(as.matrix(tax_table_16S)),
                                       sample_data(read.csv(metadata_path_16S, row.names = 1)))

# Set figure directory
fig_dir <- here::here("figures", "heatmaps")


## ----Plot COI ----------------------------------------------------------------

### Prep data

# Subset environmental samples for comparison
ESP_CTD_COI_envt_phyloseq <- prune_samples(sample_data(ESP_CTD_COI_phyloseq)$CTD_or_ESP %in% c("CTD","ESP"),ESP_CTD_COI_phyloseq)

# Merge at the phylum level
ESP_CTD_COI_envt_phyloseq_phylum <- tax_glom(ESP_CTD_COI_envt_phyloseq,taxrank=rank_names(ESP_CTD_COI_envt_phyloseq)[2], NArm = FALSE)

# Extract tax table
ESP_CTD_COI_phylum_tax_table <- as.data.frame(as(tax_table(ESP_CTD_COI_envt_phyloseq_phylum),"matrix"))
ESP_CTD_COI_phylum_tax_table %>% tibble::rownames_to_column("ASV") -> ESP_CTD_COI_phylum_tax_table
ESP_CTD_COI_phylum_tax_table <- ESP_CTD_COI_phylum_tax_table[,c("ASV","Kingdom","Phylum")]

# Extract asv table
ESP_CTD_COI_phylum_asv_table <- as.data.frame(as(otu_table(ESP_CTD_COI_envt_phyloseq_phylum),"matrix"))
ESP_CTD_COI_phylum_asv_table %>% tibble::rownames_to_column("ASV") -> ESP_CTD_COI_phylum_asv_table

# Convert all asv counts to relative abundance
ESP_CTD_COI_phylum_asv_table <- tibble::column_to_rownames(ESP_CTD_COI_phylum_asv_table, "ASV")
ESP_CTD_COI_phylum_asv_table_matrix <- as.matrix(ESP_CTD_COI_phylum_asv_table)
ESP_CTD_COI_phylum_sample_sums <- as.vector(colSums(ESP_CTD_COI_phylum_asv_table))
t(t(ESP_CTD_COI_phylum_asv_table_matrix)/ESP_CTD_COI_phylum_sample_sums) -> ESP_CTD_COI_phylum_normalized_asv_table
ESP_CTD_COI_phylum_normalized_asv_table <- as.data.frame(ESP_CTD_COI_phylum_normalized_asv_table)
ESP_CTD_COI_phylum_normalized_asv_table <- rownames_to_column(ESP_CTD_COI_phylum_normalized_asv_table,"ASV")

# Replace ASVs with the phylum it belongs to in asv table

ESP_CTD_COI_phylum_asv_table_combined <- left_join(ESP_CTD_COI_phylum_normalized_asv_table,ESP_CTD_COI_phylum_tax_table,by = "ASV")
ESP_CTD_COI_phylum_asv_table_combined <- dplyr::select(ESP_CTD_COI_phylum_asv_table_combined,-c(Kingdom,ASV))
# combine no_hit, unknown, and unassigned together
ESP_CTD_COI_phylum_asv_table_combined$Phylum <- as.character(ESP_CTD_COI_phylum_asv_table_combined$Phylum)
ESP_CTD_COI_phylum_asv_table_combined %>% mutate(., Phylum = ifelse(Phylum %in% c("no_hit","unassigned","unknown"),"unassigned",Phylum)) -> ESP_CTD_COI_phylum_asv_table_combined
# add the values for all "unassigned" of the same phyla together
ESP_CTD_COI_phylum_asv_table_combined %>% group_by(Phylum) %>% summarise_all(funs(sum)) -> ESP_CTD_COI_phylum_asv_table_combined
# remove index row names
rownames(ESP_CTD_COI_phylum_asv_table_combined) <- NULL
ESP_CTD_COI_phylum_asv_table_combined <- tibble::column_to_rownames(ESP_CTD_COI_phylum_asv_table_combined,"Phylum")

# Transpose
ESP_CTD_COI_phylum_asv_table_combined <- as.data.frame(t(ESP_CTD_COI_phylum_asv_table_combined))
ESP_CTD_COI_phylum_asv_table_combined <- tibble::rownames_to_column(ESP_CTD_COI_phylum_asv_table_combined,"sample_name")

# Determine which phyla are the most abundant
COI_phyla_sums <- as.data.frame(colSums(dplyr::select(ESP_CTD_COI_phylum_asv_table_combined,-sample_name)))
COI_phyla_sums <- tibble::rownames_to_column(COI_phyla_sums,"Phylum")
colnames(COI_phyla_sums)[2] <- "total_reads"
COI_phyla_sums <- COI_phyla_sums[order(-COI_phyla_sums$total_reads),]
top_10_COI_phyla <- subset(COI_phyla_sums, Phylum != "unassigned")$Phylum[1:10] 

# Extract metadata
ESP_CTD_COI_envt_metadata <- as.data.frame(as(sample_data(ESP_CTD_COI_envt_phyloseq_phylum),"matrix"))
ESP_CTD_COI_envt_metadata <- tibble::rownames_to_column(ESP_CTD_COI_envt_metadata,"sample_name")
ESP_CTD_COI_envt_metadata %>% mutate(.,name = paste0(matching_ID,"_",CTD_or_ESP)) -> ESP_CTD_COI_envt_metadata # add a new name that combines the matching ID and sample type
ESP_CTD_COI_envt_metadata_subset <- dplyr::select(ESP_CTD_COI_envt_metadata,c(sample_name,CTD_or_ESP,matching_ID,name)) # subset only the data that we want

# Join together
ESP_CTD_COI_phylum_superheat_df <- left_join(ESP_CTD_COI_phylum_asv_table_combined,ESP_CTD_COI_envt_metadata_subset,by = "sample_name")

ESP_CTD_COI_phylum_superheat_df <- tibble::column_to_rownames(ESP_CTD_COI_phylum_superheat_df,"name")
# Subset only the top 10 phyla

ESP_CTD_COI_phylum_superheat_df_top10phyla <- ESP_CTD_COI_phylum_superheat_df[,top_10_COI_phyla]
ESP_CTD_COI_phylum_superheat_df_top10phyla <- as.data.frame(t(ESP_CTD_COI_phylum_superheat_df_top10phyla))


# Determine the sample names of the ESP and CTD samples
COI_CTD_sample_names <- subset(ESP_CTD_COI_envt_metadata,CTD_or_ESP == "CTD")$name
COI_ESP_sample_names <- subset(ESP_CTD_COI_envt_metadata,CTD_or_ESP == "ESP")$name


### Plot in ggplot

#### Prep data

# Convert to long form
gathercols <- colnames(ESP_CTD_COI_phylum_superheat_df_top10phyla)
ESP_CTD_COI_phylum_superheat_df_top10phyla_v2 <- tibble::rownames_to_column(ESP_CTD_COI_phylum_superheat_df_top10phyla,"Phylum")

ESP_CTD_COI_phylum_superheat_df_top10phyla_gather <- gather(ESP_CTD_COI_phylum_superheat_df_top10phyla_v2, key = "sample", value = "reads", gathercols, factor_key=TRUE)

# Make phyla into a factor (for plotting)
ESP_CTD_COI_phylum_superheat_df_top10phyla_gather$Phylum <- as.factor(ESP_CTD_COI_phylum_superheat_df_top10phyla_gather$Phylum)

# Split into ESP and CTD DFs
CTD_COI_phylum_superheat_df_top10phyla_gather <- subset(ESP_CTD_COI_phylum_superheat_df_top10phyla_gather, sample %in% COI_CTD_sample_names)
ESP_COI_phylum_superheat_df_top10phyla_gather <- subset(ESP_CTD_COI_phylum_superheat_df_top10phyla_gather, sample %in% COI_ESP_sample_names)

# Reorder the CTD samples so that they're in the same order as the ESP samples
# rev(as.vector(unique(CTD_COI_phylum_superheat_df_top10phyla_gather$sample)))
CTD_COI_phylum_superheat_df_top10phyla_gather$sample <- factor(CTD_COI_phylum_superheat_df_top10phyla_gather$sample, levels = c("CN18F1_CTD", "CN18F2_CTD", "CN18F3_CTD", "CN18F4_CTD", "CN18F5_CTD", "CN18F6_CTD", "CN18F7_CTD", "CN18F8_CTD", "CN18F9_CTD", "CN18F10_CTD", "CN18F11_CTD", "CN18F12_CTD", "CN18F13_CTD", "CN18F14_CTD", "CN18F15_CTD", "CN18F16_CTD", "CN18F17_CTD", "CN18F18_CTD", "CN18F19_CTD", "CN18F20_CTD", "CN18F21_CTD", "CN18F22_CTD", "CN18S1_CTD", "CN18S2_CTD", "CN18S3_CTD", "CN18S4_CTD", "CN18S5_CTD", "CN18S6_CTD", "CN18S7_CTD", "CN18S8_CTD", "CN18S9_CTD", "CN18S10_CTD"))


# Change zeros to NAs for plotting purposes
CTD_COI_phylum_superheat_df_top10phyla_gather[CTD_COI_phylum_superheat_df_top10phyla_gather$reads == 0,]$reads <- NA
ESP_COI_phylum_superheat_df_top10phyla_gather[ESP_COI_phylum_superheat_df_top10phyla_gather$reads == 0,]$reads <- NA

#### Make the heatmap - two panel

CTD_COI_phylum_heatmap_gg <- ggplot(CTD_COI_phylum_superheat_df_top10phyla_gather, aes(x = sample, y = Phylum, fill= reads)) + 
  geom_tile(color = "white")+
  scale_fill_viridis(discrete = FALSE, breaks = c(0.001,0.25, 0.5, 0.75, 1), limits = c(0,1), na.value = "gray50")+
  scale_y_discrete(limits = rev(levels(CTD_COI_phylum_superheat_df_top10phyla_gather$Phylum)))+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 15),
        axis.title = element_blank(),
        axis.ticks.length = unit(0,"cm"),
        legend.position = "none",
        plot.title = element_text(size = 15, face = "bold",hjust = 0.5),
        panel.background = element_rect(fill = "white", color = "white"))

ESP_COI_phylum_heatmap_gg <- ggplot(ESP_COI_phylum_superheat_df_top10phyla_gather, aes(x = sample, y = Phylum, fill= reads)) + 
  geom_tile(color = "white")+
  scale_fill_viridis(discrete = FALSE, breaks = c(0.001,0.25, 0.5, 0.75, 1), limits = c(0,1), na.value = "gray50")+
  labs(fill = "Relative\nabundance")+
  scale_y_discrete(limits = rev(levels(CTD_COI_phylum_superheat_df_top10phyla_gather$Phylum)))+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_blank(),
        axis.ticks.length = unit(0,"cm"),
        legend.position = "none",
        plot.title = element_text(size = 15, face = "bold",hjust = 0.5),
        panel.background = element_rect(fill = "white", color = "white"),
        # plot.margin = margin(0, 1.5, 0, 0, "cm"))
        plot.margin = margin(1, 0, 0, 0, "cm"))

gA <- ggplotGrob(CTD_COI_phylum_heatmap_gg)
gB <- ggplotGrob(ESP_COI_phylum_heatmap_gg)
ESP_CTD_COI_gg_heatmap <- arrangeGrob(cbind(gA, gB))


ggsave(paste0(fig_dir, "/ESP_CTD_COI_gg_heatmap.png"), ESP_CTD_COI_gg_heatmap, width = 8.5, height = 2.5)






## ----Plot 18S ----------------------------------------------------------------

### Prep data

# Subset environmental samples for comparison
ESP_CTD_18S_envt_phyloseq <- prune_samples(sample_data(ESP_CTD_18S_phyloseq)$CTD_or_ESP %in% c("CTD","ESP"),ESP_CTD_18S_phyloseq)

# Merge at the phylum level
ESP_CTD_18S_envt_phyloseq_phylum <- tax_glom(ESP_CTD_18S_envt_phyloseq,taxrank=rank_names(ESP_CTD_18S_envt_phyloseq)[2], NArm = FALSE)

# Extract tax table
ESP_CTD_18S_phylum_tax_table <- as.data.frame(as(tax_table(ESP_CTD_18S_envt_phyloseq_phylum),"matrix"))
ESP_CTD_18S_phylum_tax_table %>% tibble::rownames_to_column("ASV") -> ESP_CTD_18S_phylum_tax_table
ESP_CTD_18S_phylum_tax_table <- ESP_CTD_18S_phylum_tax_table[,c("ASV","Kingdom","Phylum")]

# Extract asv table
ESP_CTD_18S_phylum_asv_table <- as.data.frame(as(otu_table(ESP_CTD_18S_envt_phyloseq_phylum),"matrix"))
ESP_CTD_18S_phylum_asv_table %>% tibble::rownames_to_column("ASV") -> ESP_CTD_18S_phylum_asv_table

# Convert all asv counts to relative abundance
ESP_CTD_18S_phylum_asv_table <- tibble::column_to_rownames(ESP_CTD_18S_phylum_asv_table, "ASV")
ESP_CTD_18S_phylum_asv_table_matrix <- as.matrix(ESP_CTD_18S_phylum_asv_table)
ESP_CTD_18S_phylum_sample_sums <- as.vector(colSums(ESP_CTD_18S_phylum_asv_table))
t(t(ESP_CTD_18S_phylum_asv_table_matrix)/ESP_CTD_18S_phylum_sample_sums) -> ESP_CTD_18S_phylum_normalized_asv_table
ESP_CTD_18S_phylum_normalized_asv_table <- as.data.frame(ESP_CTD_18S_phylum_normalized_asv_table)
ESP_CTD_18S_phylum_normalized_asv_table <- rownames_to_column(ESP_CTD_18S_phylum_normalized_asv_table,"ASV")


# Replace ASVs with the phylum it belongs to in asv table
ESP_CTD_18S_phylum_asv_table_combined <- left_join(ESP_CTD_18S_phylum_normalized_asv_table,ESP_CTD_18S_phylum_tax_table,by = "ASV")
ESP_CTD_18S_phylum_asv_table_combined <- dplyr::select(ESP_CTD_18S_phylum_asv_table_combined,-c(Kingdom,ASV))
# combine no_hit, unknown, and unassigned together
ESP_CTD_18S_phylum_asv_table_combined$Phylum <- as.character(ESP_CTD_18S_phylum_asv_table_combined$Phylum)
ESP_CTD_18S_phylum_asv_table_combined %>% mutate(., Phylum = ifelse(Phylum %in% c("no_hit","unassigned","unknown"),"unassigned",Phylum)) -> ESP_CTD_18S_phylum_asv_table_combined
# add the values for all "unassigned" of the same phyla together
ESP_CTD_18S_phylum_asv_table_combined %>% group_by(Phylum) %>% summarise_all(funs(sum)) -> ESP_CTD_18S_phylum_asv_table_combined
# remove index row names
rownames(ESP_CTD_18S_phylum_asv_table_combined) <- NULL
ESP_CTD_18S_phylum_asv_table_combined <- tibble::column_to_rownames(ESP_CTD_18S_phylum_asv_table_combined,"Phylum")

# Transpose
ESP_CTD_18S_phylum_asv_table_combined <- as.data.frame(t(ESP_CTD_18S_phylum_asv_table_combined))
ESP_CTD_18S_phylum_asv_table_combined <- tibble::rownames_to_column(ESP_CTD_18S_phylum_asv_table_combined,"sample_name")

# Determine which phyla are the most abundant
phyla_sums_18S <- as.data.frame(colSums(dplyr::select(ESP_CTD_18S_phylum_asv_table_combined,-sample_name)))
phyla_sums_18S <- tibble::rownames_to_column(phyla_sums_18S,"Phylum")
colnames(phyla_sums_18S)[2] <- "total_reads"
phyla_sums_18S <- phyla_sums_18S[order(-phyla_sums_18S$total_reads),]
top_10_18S_phyla <- subset(phyla_sums_18S, Phylum != "unassigned")$Phylum[1:10] 

# Extract metadata
ESP_CTD_18S_envt_metadata <- as.data.frame(as(sample_data(ESP_CTD_18S_envt_phyloseq_phylum),"matrix"))
ESP_CTD_18S_envt_metadata <- tibble::rownames_to_column(ESP_CTD_18S_envt_metadata,"sample_name")
ESP_CTD_18S_envt_metadata %>% mutate(.,name = paste0(matching_ID,"_",CTD_or_ESP)) -> ESP_CTD_18S_envt_metadata # add a new name that combines the matching ID and sample type
ESP_CTD_18S_envt_metadata_subset <- dplyr::select(ESP_CTD_18S_envt_metadata,c(sample_name,CTD_or_ESP,matching_ID,name)) # subset only the data that we want

# Join together
ESP_CTD_18S_phylum_superheat_df <- left_join(ESP_CTD_18S_phylum_asv_table_combined,ESP_CTD_18S_envt_metadata_subset,by = "sample_name")

ESP_CTD_18S_phylum_superheat_df <- tibble::column_to_rownames(ESP_CTD_18S_phylum_superheat_df,"name")
# Subset only the top 10 phyla

ESP_CTD_18S_phylum_superheat_df_top10phyla <- ESP_CTD_18S_phylum_superheat_df[,top_10_18S_phyla]
ESP_CTD_18S_phylum_superheat_df_top10phyla <- as.data.frame(t(ESP_CTD_18S_phylum_superheat_df_top10phyla))


# Determine the sample names of the ESP and CTD samples
CTD_sample_names_18S <- subset(ESP_CTD_18S_envt_metadata,CTD_or_ESP == "CTD")$name
ESP_sample_names_18S <- subset(ESP_CTD_18S_envt_metadata,CTD_or_ESP == "ESP")$name


### Plot in ggplot

#### Prep data

# Convert to long form
gathercols <- colnames(ESP_CTD_18S_phylum_superheat_df_top10phyla)
ESP_CTD_18S_phylum_superheat_df_top10phyla_v2 <- tibble::rownames_to_column(ESP_CTD_18S_phylum_superheat_df_top10phyla,"Phylum")

ESP_CTD_18S_phylum_superheat_df_top10phyla_gather <- gather(ESP_CTD_18S_phylum_superheat_df_top10phyla_v2, key = "sample", value = "reads", gathercols, factor_key=TRUE)

# Make phyla into a factor (for plotting)
ESP_CTD_18S_phylum_superheat_df_top10phyla_gather$Phylum <- as.factor(ESP_CTD_18S_phylum_superheat_df_top10phyla_gather$Phylum)

# Split into ESP and CTD DFs
CTD_18S_phylum_superheat_df_top10phyla_gather <- subset(ESP_CTD_18S_phylum_superheat_df_top10phyla_gather, sample %in% CTD_sample_names_18S)
ESP_18S_phylum_superheat_df_top10phyla_gather <- subset(ESP_CTD_18S_phylum_superheat_df_top10phyla_gather, sample %in% ESP_sample_names_18S)

# Reorder the CTD samples so that they're in the same order as the ESP samples
# rev(as.vector(unique(CTD_18S_phylum_superheat_df_top10phyla_gather$sample)))
CTD_18S_phylum_superheat_df_top10phyla_gather$sample <- factor(CTD_18S_phylum_superheat_df_top10phyla_gather$sample, levels = c("CN18F1_CTD", "CN18F2_CTD", "CN18F3_CTD", "CN18F4_CTD", "CN18F5_CTD", "CN18F6_CTD", "CN18F7_CTD", "CN18F8_CTD", "CN18F9_CTD", "CN18F10_CTD", "CN18F11_CTD", "CN18F12_CTD", "CN18F13_CTD", "CN18F14_CTD", "CN18F15_CTD", "CN18F16_CTD", "CN18F17_CTD", "CN18F18_CTD", "CN18F19_CTD", "CN18F20_CTD", "CN18F21_CTD", "CN18F22_CTD", "CN18S1_CTD", "CN18S2_CTD", "CN18S3_CTD", "CN18S4_CTD", "CN18S5_CTD", "CN18S6_CTD", "CN18S7_CTD", "CN18S8_CTD", "CN18S9_CTD", "CN18S10_CTD"))

# Change zeros to NAs for plotting purposes
CTD_18S_phylum_superheat_df_top10phyla_gather[CTD_18S_phylum_superheat_df_top10phyla_gather$reads == 0,]$reads <- NA
ESP_18S_phylum_superheat_df_top10phyla_gather[ESP_18S_phylum_superheat_df_top10phyla_gather$reads == 0,]$reads <- NA

#### Make the heatmap - two panel

CTD_18S_phylum_heatmap_gg <- ggplot(CTD_18S_phylum_superheat_df_top10phyla_gather, aes(x = sample, y = Phylum, fill= reads)) + 
  geom_tile(color = "white")+
  scale_fill_viridis(discrete = FALSE, breaks = c(0.001,0.25, 0.5, 0.75, 1), limits = c(0,1), na.value = "gray50")+
  scale_y_discrete(limits = rev(levels(CTD_18S_phylum_superheat_df_top10phyla_gather$Phylum)))+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 15),
        axis.title = element_blank(),
        axis.ticks.length = unit(0,"cm"),
        legend.position = "none",
        plot.title = element_text(size = 15, face = "bold",hjust = 0.5),
        panel.background = element_rect(fill = "white", color = "white"))

ESP_18S_phylum_heatmap_gg <- ggplot(ESP_18S_phylum_superheat_df_top10phyla_gather, aes(x = sample, y = Phylum, fill= reads)) + 
  geom_tile(color = "white")+
  scale_fill_viridis(discrete = FALSE, breaks = c(0.001,0.25, 0.5, 0.75, 1), limits = c(0,1), na.value = "gray50")+
  labs(fill = "Relative\nabundance")+
  scale_y_discrete(limits = rev(levels(CTD_18S_phylum_superheat_df_top10phyla_gather$Phylum)))+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_blank(),
        axis.ticks.length = unit(0,"cm"),
        legend.position = "none",
        plot.title = element_text(size = 15, face = "bold",hjust = 0.5),
        panel.background = element_rect(fill = "white", color = "white"),
        # plot.margin = margin(0, 1.5, 0, 0, "cm"))
        plot.margin = margin(1, 0, 0, 0, "cm"))

gA <- ggplotGrob(CTD_18S_phylum_heatmap_gg)
gB <- ggplotGrob(ESP_18S_phylum_heatmap_gg)
ESP_CTD_18S_gg_heatmap <- arrangeGrob(cbind(gA, gB))


ggsave(paste0(fig_dir, "/ESP_CTD_18S_gg_heatmap.png"), ESP_CTD_18S_gg_heatmap, width = 8.5, height = 2.5)



## ----Plot 12S ----------------------------------------------------------------

### Prep data

# Subset environmental samples for comparison
ESP_CTD_12S_envt_phyloseq <- prune_samples(sample_data(ESP_CTD_12S_phyloseq)$CTD_or_ESP %in% c("CTD","ESP"),ESP_CTD_12S_phyloseq)

# Merge at the family level
ESP_CTD_12S_envt_phyloseq_family <- tax_glom(ESP_CTD_12S_envt_phyloseq,taxrank=rank_names(ESP_CTD_12S_envt_phyloseq)[5], NArm = FALSE)

# Extract tax table
ESP_CTD_12S_family_tax_table <- as.data.frame(as(tax_table(ESP_CTD_12S_envt_phyloseq_family),"matrix"))
ESP_CTD_12S_family_tax_table %>% tibble::rownames_to_column("ASV") -> ESP_CTD_12S_family_tax_table
ESP_CTD_12S_family_tax_table <- ESP_CTD_12S_family_tax_table[,c("ASV","Kingdom","Family")]

# Extract asv table
ESP_CTD_12S_family_asv_table <- as.data.frame(as(otu_table(ESP_CTD_12S_envt_phyloseq_family),"matrix"))
ESP_CTD_12S_family_asv_table %>% tibble::rownames_to_column("ASV") -> ESP_CTD_12S_family_asv_table

# Convert all asv counts to relative abundance
ESP_CTD_12S_family_asv_table <- tibble::column_to_rownames(ESP_CTD_12S_family_asv_table, "ASV")
ESP_CTD_12S_family_asv_table_matrix <- as.matrix(ESP_CTD_12S_family_asv_table)
ESP_CTD_12S_family_sample_sums <- as.vector(colSums(ESP_CTD_12S_family_asv_table))
t(t(ESP_CTD_12S_family_asv_table_matrix)/ESP_CTD_12S_family_sample_sums) -> ESP_CTD_12S_family_normalized_asv_table
ESP_CTD_12S_family_normalized_asv_table <- as.data.frame(ESP_CTD_12S_family_normalized_asv_table)
ESP_CTD_12S_family_normalized_asv_table <- rownames_to_column(ESP_CTD_12S_family_normalized_asv_table,"ASV")


# Replace ASVs with the family it belongs to in asv table
ESP_CTD_12S_family_asv_table_combined <- left_join(ESP_CTD_12S_family_normalized_asv_table,ESP_CTD_12S_family_tax_table,by = "ASV")
ESP_CTD_12S_family_asv_table_combined <- dplyr::select(ESP_CTD_12S_family_asv_table_combined,-c(Kingdom,ASV))
# combine no_hit, unknown, and unassigned together
ESP_CTD_12S_family_asv_table_combined$Family <- as.character(ESP_CTD_12S_family_asv_table_combined$Family)
ESP_CTD_12S_family_asv_table_combined %>% mutate(., Family = ifelse(Family %in% c("no_hit","unassigned","unknown"),"unassigned",Family)) -> ESP_CTD_12S_family_asv_table_combined
# add the values for all "unassigned" of the same family together
ESP_CTD_12S_family_asv_table_combined %>% group_by(Family) %>% summarise_all(funs(sum)) -> ESP_CTD_12S_family_asv_table_combined
# remove index row names
rownames(ESP_CTD_12S_family_asv_table_combined) <- NULL
ESP_CTD_12S_family_asv_table_combined <- tibble::column_to_rownames(ESP_CTD_12S_family_asv_table_combined,"Family")

# Transpose
ESP_CTD_12S_family_asv_table_combined <- as.data.frame(t(ESP_CTD_12S_family_asv_table_combined))
ESP_CTD_12S_family_asv_table_combined <- tibble::rownames_to_column(ESP_CTD_12S_family_asv_table_combined,"sample_name")

# Determine which family are the most abundant
family_sums_12S <- as.data.frame(colSums(dplyr::select(ESP_CTD_12S_family_asv_table_combined,-sample_name)))
family_sums_12S <- tibble::rownames_to_column(family_sums_12S,"Family")
colnames(family_sums_12S)[2] <- "total_reads"
family_sums_12S <- family_sums_12S[order(-family_sums_12S$total_reads),]
top_10_12S_family <- subset(family_sums_12S, Family != "unassigned")$Family[1:10] 

# Extract metadata
ESP_CTD_12S_envt_metadata <- as.data.frame(as(sample_data(ESP_CTD_12S_envt_phyloseq_family),"matrix"))
ESP_CTD_12S_envt_metadata <- tibble::rownames_to_column(ESP_CTD_12S_envt_metadata,"sample_name")
ESP_CTD_12S_envt_metadata %>% mutate(.,name = paste0(matching_ID,"_",CTD_or_ESP)) -> ESP_CTD_12S_envt_metadata # add a new name that combines the matching ID and sample type
ESP_CTD_12S_envt_metadata_subset <- dplyr::select(ESP_CTD_12S_envt_metadata,c(sample_name,CTD_or_ESP,matching_ID,name)) # subset only the data that we want

# Join together
ESP_CTD_12S_family_superheat_df <- left_join(ESP_CTD_12S_family_asv_table_combined,ESP_CTD_12S_envt_metadata_subset,by = "sample_name")

ESP_CTD_12S_family_superheat_df <- tibble::column_to_rownames(ESP_CTD_12S_family_superheat_df,"name")
# Subset only the top 10 family

ESP_CTD_12S_family_superheat_df_top10family <- ESP_CTD_12S_family_superheat_df[,top_10_12S_family]
ESP_CTD_12S_family_superheat_df_top10family <- as.data.frame(t(ESP_CTD_12S_family_superheat_df_top10family))


# Determine the sample names of the ESP and CTD samples
CTD_sample_names_12S <- subset(ESP_CTD_12S_envt_metadata,CTD_or_ESP == "CTD")$name
ESP_sample_names_12S <- subset(ESP_CTD_12S_envt_metadata,CTD_or_ESP == "ESP")$name

setdiff(CTD_sample_names_12S, unique(colnames(ESP_CTD_12S_family_superheat_df_top10family)))

setdiff(CTD_sample_names_12S, unique(ESP_CTD_12S_family_superheat_df_top10family_gather$sample))
unique(ESP_CTD_12S_family_superheat_df_top10family_gather$sample)

### Plot in ggplot

#### Prep data

# Convert to long form
gathercols <- colnames(ESP_CTD_12S_family_superheat_df_top10family)
ESP_CTD_12S_family_superheat_df_top10family_v2 <- tibble::rownames_to_column(ESP_CTD_12S_family_superheat_df_top10family,"Family")

ESP_CTD_12S_family_superheat_df_top10family_gather <- gather(ESP_CTD_12S_family_superheat_df_top10family_v2, key = "sample", value = "reads", gathercols, factor_key=TRUE)

# Make family into a factor (for plotting)
ESP_CTD_12S_family_superheat_df_top10family_gather$Family <- as.factor(ESP_CTD_12S_family_superheat_df_top10family_gather$Family)

# Split into ESP and CTD DFs
CTD_12S_family_superheat_df_top10family_gather <- subset(ESP_CTD_12S_family_superheat_df_top10family_gather, sample %in% CTD_sample_names_12S)
ESP_12S_family_superheat_df_top10family_gather <- subset(ESP_CTD_12S_family_superheat_df_top10family_gather, sample %in% ESP_sample_names_12S)

# Reorder the CTD samples so that they're in the same order as the ESP samples
# rev(as.vector(unique(CTD_12S_family_superheat_df_top10family_gather$sample)))
CTD_12S_family_superheat_df_top10family_gather$sample <- factor(CTD_12S_family_superheat_df_top10family_gather$sample, 
                                                                  levels = c("CN18F1_CTD", "CN18F2_CTD", "CN18F3_CTD", "CN18F4_CTD", 
                                                                             "CN18F5_CTD", "CN18F6_CTD", "CN18F7_CTD", "CN18F8_CTD", 
                                                                             "CN18F9_CTD", "CN18F10_CTD", "CN18F11_CTD", "CN18F12_CTD", 
                                                                             "CN18F13_CTD", "CN18F14_CTD", "CN18F15_CTD", "CN18F16_CTD", 
                                                                             "CN18F17_CTD", "CN18F18_CTD", "CN18F19_CTD", "CN18F20_CTD", 
                                                                             "CN18F21_CTD", "CN18F22_CTD", "CN18S1_CTD", "CN18S2_CTD", 
                                                                             "CN18S3_CTD", "CN18S4_CTD", "CN18S5_CTD", "CN18S6_CTD", 
                                                                             "CN18S7_CTD", "CN18S8_CTD", "CN18S9_CTD", "CN18S10_CTD"))
# Change zeros to NAs for plotting purposes
CTD_12S_family_superheat_df_top10family_gather[CTD_12S_family_superheat_df_top10family_gather$reads == 0,]$reads <- NA
ESP_12S_family_superheat_df_top10family_gather[ESP_12S_family_superheat_df_top10family_gather$reads == 0,]$reads <- NA


#### Make the heatmap - two panel

CTD_12S_family_heatmap_gg <- ggplot(CTD_12S_family_superheat_df_top10family_gather, aes(x = sample, y = Family, fill= reads)) + 
  geom_tile(color = "white")+
  scale_fill_viridis(discrete = FALSE, breaks = c(0.001,0.25, 0.5, 0.75, 1), limits = c(0,1), na.value = "gray50")+
  scale_y_discrete(limits = rev(levels(CTD_12S_family_superheat_df_top10family_gather$Family)))+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 15),
        axis.title = element_blank(),
        axis.ticks.length = unit(0,"cm"),
        legend.position = "none",
        plot.title = element_text(size = 15, face = "bold",hjust = 0.5),
        panel.background = element_rect(fill = "white", color = "white"))

ESP_12S_family_heatmap_gg <- ggplot(ESP_12S_family_superheat_df_top10family_gather, aes(x = sample, y = Family, fill= reads)) + 
  geom_tile(color = "white")+
  scale_fill_viridis(discrete = FALSE, breaks = c(0.001,0.25, 0.5, 0.75, 1), limits = c(0,1), na.value = "gray50")+
  labs(fill = "Relative\nabundance")+
  scale_y_discrete(limits = rev(levels(CTD_12S_family_superheat_df_top10family_gather$Family)))+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_blank(),
        axis.ticks.length = unit(0,"cm"),
        legend.position = "none",
        plot.title = element_text(size = 15, face = "bold",hjust = 0.5),
        panel.background = element_rect(fill = "white", color = "white"),
        # plot.margin = margin(0, 1.5, 0, 0, "cm"))
        plot.margin = margin(1, 0, 0, 0, "cm"))

gA <- ggplotGrob(CTD_12S_family_heatmap_gg)
gB <- ggplotGrob(ESP_12S_family_heatmap_gg)
ESP_CTD_12S_gg_heatmap <- arrangeGrob(cbind(gA, gB))


ggsave(paste0(fig_dir, "/ESP_CTD_12S_gg_heatmap.png"), ESP_CTD_12S_gg_heatmap, width = 8.5, height = 2.5)


## ----Plot 16S ----------------------------------------------------------------

### Prep data

# Subset environmental samples for comparison
ESP_CTD_16S_envt_phyloseq <- prune_samples(sample_data(ESP_CTD_16S_phyloseq)$CTD_or_ESP %in% c("CTD","ESP"),ESP_CTD_16S_phyloseq)

# Merge at the phylum level
ESP_CTD_16S_envt_phyloseq_phylum <- tax_glom(ESP_CTD_16S_envt_phyloseq,taxrank=rank_names(ESP_CTD_16S_envt_phyloseq)[2], NArm = FALSE)

# Extract tax table
ESP_CTD_16S_phylum_tax_table <- as.data.frame(as(tax_table(ESP_CTD_16S_envt_phyloseq_phylum),"matrix"))
ESP_CTD_16S_phylum_tax_table %>% tibble::rownames_to_column("ASV") -> ESP_CTD_16S_phylum_tax_table
ESP_CTD_16S_phylum_tax_table <- ESP_CTD_16S_phylum_tax_table[,c("ASV","Kingdom","Phylum")]

# Extract asv table
ESP_CTD_16S_phylum_asv_table <- as.data.frame(as(otu_table(ESP_CTD_16S_envt_phyloseq_phylum),"matrix"))
ESP_CTD_16S_phylum_asv_table %>% tibble::rownames_to_column("ASV") -> ESP_CTD_16S_phylum_asv_table

# Convert all asv counts to relative abundance
ESP_CTD_16S_phylum_asv_table <- tibble::column_to_rownames(ESP_CTD_16S_phylum_asv_table, "ASV")
ESP_CTD_16S_phylum_asv_table_matrix <- as.matrix(ESP_CTD_16S_phylum_asv_table)
ESP_CTD_16S_phylum_sample_sums <- as.vector(colSums(ESP_CTD_16S_phylum_asv_table))
t(t(ESP_CTD_16S_phylum_asv_table_matrix)/ESP_CTD_16S_phylum_sample_sums) -> ESP_CTD_16S_phylum_normalized_asv_table
ESP_CTD_16S_phylum_normalized_asv_table <- as.data.frame(ESP_CTD_16S_phylum_normalized_asv_table)
ESP_CTD_16S_phylum_normalized_asv_table <- rownames_to_column(ESP_CTD_16S_phylum_normalized_asv_table,"ASV")


# Replace ASVs with the phylum it belongs to in asv table
ESP_CTD_16S_phylum_asv_table_combined <- left_join(ESP_CTD_16S_phylum_normalized_asv_table,ESP_CTD_16S_phylum_tax_table,by = "ASV")
ESP_CTD_16S_phylum_asv_table_combined <- dplyr::select(ESP_CTD_16S_phylum_asv_table_combined,-c(Kingdom,ASV))
# combine no_hit, unknown, and unassigned together
ESP_CTD_16S_phylum_asv_table_combined$Phylum <- as.character(ESP_CTD_16S_phylum_asv_table_combined$Phylum)
ESP_CTD_16S_phylum_asv_table_combined %>% mutate(., Phylum = ifelse(Phylum %in% c("no_hit","unassigned","unknown"),"unassigned",Phylum)) -> ESP_CTD_16S_phylum_asv_table_combined
# add the values for all "unassigned" of the same phyla together
ESP_CTD_16S_phylum_asv_table_combined %>% group_by(Phylum) %>% summarise_all(funs(sum)) -> ESP_CTD_16S_phylum_asv_table_combined
# remove index row names
rownames(ESP_CTD_16S_phylum_asv_table_combined) <- NULL
ESP_CTD_16S_phylum_asv_table_combined <- tibble::column_to_rownames(ESP_CTD_16S_phylum_asv_table_combined,"Phylum")

# Transpose
ESP_CTD_16S_phylum_asv_table_combined <- as.data.frame(t(ESP_CTD_16S_phylum_asv_table_combined))
ESP_CTD_16S_phylum_asv_table_combined <- tibble::rownames_to_column(ESP_CTD_16S_phylum_asv_table_combined,"sample_name")

# Determine which phyla are the most abundant
phyla_sums_16S <- as.data.frame(colSums(dplyr::select(ESP_CTD_16S_phylum_asv_table_combined,-sample_name)))
phyla_sums_16S <- tibble::rownames_to_column(phyla_sums_16S,"Phylum")
colnames(phyla_sums_16S)[2] <- "total_reads"
phyla_sums_16S <- phyla_sums_16S[order(-phyla_sums_16S$total_reads),]
top_10_16S_phyla <- subset(phyla_sums_16S, Phylum != "unassigned")$Phylum[1:10] 

# Extract metadata
ESP_CTD_16S_envt_metadata <- as.data.frame(as(sample_data(ESP_CTD_16S_envt_phyloseq_phylum),"matrix"))
ESP_CTD_16S_envt_metadata <- tibble::rownames_to_column(ESP_CTD_16S_envt_metadata,"sample_name")
ESP_CTD_16S_envt_metadata %>% mutate(.,name = paste0(matching_ID,"_",CTD_or_ESP)) -> ESP_CTD_16S_envt_metadata # add a new name that combines the matching ID and sample type
ESP_CTD_16S_envt_metadata_subset <- dplyr::select(ESP_CTD_16S_envt_metadata,c(sample_name,CTD_or_ESP,matching_ID,name)) # subset only the data that we want

# Join together
ESP_CTD_16S_phylum_superheat_df <- left_join(ESP_CTD_16S_phylum_asv_table_combined,ESP_CTD_16S_envt_metadata_subset,by = "sample_name")

ESP_CTD_16S_phylum_superheat_df <- tibble::column_to_rownames(ESP_CTD_16S_phylum_superheat_df,"name")
# Subset only the top 10 phyla

ESP_CTD_16S_phylum_superheat_df_top10phyla <- ESP_CTD_16S_phylum_superheat_df[,top_10_16S_phyla]
ESP_CTD_16S_phylum_superheat_df_top10phyla <- as.data.frame(t(ESP_CTD_16S_phylum_superheat_df_top10phyla))


# Determine the sample names of the ESP and CTD samples
CTD_sample_names_16S <- subset(ESP_CTD_16S_envt_metadata,CTD_or_ESP == "CTD")$name
ESP_sample_names_16S <- subset(ESP_CTD_16S_envt_metadata,CTD_or_ESP == "ESP")$name


### Plot in ggplot

#### Prep data

# Convert to long form
gathercols <- colnames(ESP_CTD_16S_phylum_superheat_df_top10phyla)
ESP_CTD_16S_phylum_superheat_df_top10phyla_v2 <- tibble::rownames_to_column(ESP_CTD_16S_phylum_superheat_df_top10phyla,"Phylum")

ESP_CTD_16S_phylum_superheat_df_top10phyla_gather <- gather(ESP_CTD_16S_phylum_superheat_df_top10phyla_v2, key = "sample", value = "reads", gathercols, factor_key=TRUE)

# Change two phyla for plotting
ESP_CTD_16S_phylum_superheat_df_top10phyla_gather %>%
  mutate(.,Phylum = ifelse(Phylum == " Marinimicrobia_(SAR406_clade)", "Marinimicrobia",
                         ifelse(Phylum == " SAR324_clade(Marine_group_B)", "SAR324_clade", 
                         Phylum))) -> ESP_CTD_16S_phylum_superheat_df_top10phyla_gather

# Make phyla into a factor (for plotting)
ESP_CTD_16S_phylum_superheat_df_top10phyla_gather$Phylum <- as.factor(ESP_CTD_16S_phylum_superheat_df_top10phyla_gather$Phylum)

# Split into ESP and CTD DFs
CTD_16S_phylum_superheat_df_top10phyla_gather <- subset(ESP_CTD_16S_phylum_superheat_df_top10phyla_gather, sample %in% CTD_sample_names_16S)
ESP_16S_phylum_superheat_df_top10phyla_gather <- subset(ESP_CTD_16S_phylum_superheat_df_top10phyla_gather, sample %in% ESP_sample_names_16S)

# Reorder the CTD samples so that they're in the same order as the ESP samples
# rev(as.vector(unique(CTD_16S_phylum_superheat_df_top10phyla_gather$sample)))
CTD_16S_phylum_superheat_df_top10phyla_gather$sample <- factor(CTD_16S_phylum_superheat_df_top10phyla_gather$sample, levels = c("CN18F1_CTD", "CN18F2_CTD", "CN18F3_CTD", "CN18F4_CTD", "CN18F5_CTD", "CN18F6_CTD", "CN18F7_CTD", "CN18F8_CTD", "CN18F9_CTD", "CN18F10_CTD", "CN18F11_CTD", "CN18F12_CTD", "CN18F13_CTD", "CN18F14_CTD", "CN18F15_CTD", "CN18F16_CTD", "CN18F17_CTD", "CN18F18_CTD", "CN18F19_CTD", "CN18F20_CTD", "CN18F21_CTD", "CN18F22_CTD", "CN18S1_CTD", "CN18S2_CTD", "CN18S3_CTD", "CN18S4_CTD", "CN18S5_CTD", "CN18S6_CTD", "CN18S7_CTD", "CN18S8_CTD", "CN18S9_CTD", "CN18S10_CTD"))

# Change zeros to NAs for plotting purposes
# CTD_16S_phylum_superheat_df_top10phyla_gather[CTD_16S_phylum_superheat_df_top10phyla_gather$reads == 0,]$reads <- NA
# No zeros in CTD data
ESP_16S_phylum_superheat_df_top10phyla_gather[ESP_16S_phylum_superheat_df_top10phyla_gather$reads == 0,]$reads <- NA

#### Make the heatmap - two panel

CTD_16S_phylum_heatmap_gg <- ggplot(CTD_16S_phylum_superheat_df_top10phyla_gather, aes(x = sample, y = Phylum, fill= reads)) + 
  geom_tile(color = "white")+
  scale_fill_viridis(discrete = FALSE, breaks = c(0.001,0.25, 0.5, 0.75, 1), limits = c(0,1), na.value = "gray50")+
  scale_y_discrete(limits = rev(levels(CTD_16S_phylum_superheat_df_top10phyla_gather$Phylum)))+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 15),
        axis.title = element_blank(),
        axis.ticks.length = unit(0,"cm"),
        legend.position = "none",
        plot.title = element_text(size = 15, face = "bold",hjust = 0.5),
        panel.background = element_rect(fill = "white", color = "white"))

ESP_16S_phylum_heatmap_gg <- ggplot(ESP_16S_phylum_superheat_df_top10phyla_gather, aes(x = sample, y = Phylum, fill= reads)) + 
  geom_tile(color = "white")+
  scale_fill_viridis(discrete = FALSE, breaks = c(0.001,0.25, 0.5, 0.75, 1), limits = c(0,1), na.value = "gray50")+
  labs(fill = "Relative\nabundance")+
  scale_y_discrete(limits = rev(levels(CTD_16S_phylum_superheat_df_top10phyla_gather$Phylum)))+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_blank(),
        axis.ticks.length = unit(0,"cm"),
        legend.position = "none",
        plot.title = element_text(size = 15, face = "bold",hjust = 0.5),
        panel.background = element_rect(fill = "white", color = "white"),
        # plot.margin = margin(0, 1.5, 0, 0, "cm"))
        plot.margin = margin(1, 0, 0, 0, "cm"))

gA <- ggplotGrob(CTD_16S_phylum_heatmap_gg)
gB <- ggplotGrob(ESP_16S_phylum_heatmap_gg)
ESP_CTD_16S_gg_heatmap <- arrangeGrob(cbind(gA, gB))


ggsave(paste0(fig_dir, "/ESP_CTD_16S_gg_heatmap.png"), ESP_CTD_16S_gg_heatmap, width = 8.5, height = 2.5)

## ----- Extract the legend ----------------------------------------------------
ESP_CTD_heatmap_for_legend_gg <- ggplot(ESP_16S_phylum_superheat_df_top10phyla_gather, aes(x = sample, y = Phylum, fill= reads)) + 
  geom_tile(color = "white")+
  scale_fill_viridis(discrete = FALSE, breaks = c(0.001,0.25, 0.5, 0.75, 1), limits = c(0,1), na.value = "gray50",   
                     guide = guide_colorbar(title.position = "top", direction = "horizontal",label.position = "top",
                                            barwidth = 40, barheight = 2.5, title.hjust = 0.5))+
  labs(fill = "Relative abundance")+
  scale_y_discrete(limits = rev(levels(CTD_16S_phylum_superheat_df_top10phyla_gather$Phylum)))+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_blank(),
        axis.ticks.length = unit(0,"cm"),
        plot.title = element_text(size = 15, face = "bold",hjust = 0.5),
        panel.background = element_rect(fill = "white", color = "white"),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 16),
        legend.position = c(0.575, 0.5))

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

heatmap_legend <- g_legend(ESP_CTD_heatmap_for_legend_gg)


## ----Arrange heatmaps for all markers with ggarrange--------------------------

combined_heatmap_gg <- annotate_figure(ggpubr::ggarrange(egg::ggarrange(plots = list(CTD_18S_phylum_heatmap_gg, ESP_18S_phylum_heatmap_gg,
                                                                   CTD_COI_phylum_heatmap_gg, ESP_COI_phylum_heatmap_gg,
                                                                   CTD_12S_family_heatmap_gg, ESP_12S_family_heatmap_gg,
                                                                   CTD_16S_phylum_heatmap_gg, ESP_16S_phylum_heatmap_gg), ncol = 2, nrow = 4,
                                                      # labels = c(" ", "(a)", " ", "(b)", " ", "(c)", " ", "(d)"), 
                                                      labels = c("(a) 18S", " ", "(b) COI", " ", "(c) 12S", " ", "(d) 16S", " "), 
                                                      # label.args = list(gp = grid::gpar(fontsize = 25, fontface = "bold"),
                                                      #                   x = 0.9, y = 0.95)),
                                                      label.args = list(gp = grid::gpar(fontsize = 20, fontface = "bold"),
                                                                        x = 0.01, y = 0.95)),
                                                 heatmap_legend, nrow = 2, ncol = 1, heights = c(10,1)),
                                       top = text_grob("Manual                                         Automated             ", 
                                      color = "black", face = "bold", size = 25, vjust = 2))


ggsave(combined_heatmap_gg, height = 16, width = 12, path = fig_dir, filename = "combined_heatmap_fig.png", device = "png")
