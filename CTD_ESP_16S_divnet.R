#load
library(phyloseq)
library(breakaway)
library(DivNet)
library(ggplot2)
library(ggfortify)
library(RColorBrewer)
library(data.table)
library(textshape)
library(magrittr)
library(vegan)
library(doParallel)
library(dplyr)
library(foreach)
library(doSNOW)

# import ASV table and ID-to-taxonomy table

microbe <- column_to_rownames(table_phyloseq, 'ID')
tax_table <- column_to_rownames(taxonomy_phyloseq, 'ID')
metadata <- column_to_rownames(eCruises2018_metadata_phyloseq_nocontrols, 'sample_name')

# convert all to matrices
microbe <- as.matrix(microbe)
tax_table <- as.matrix(tax_table)
metadata <- as.data.frame(metadata)

# combine into phyloseq object
ASV = otu_table(microbe, taxa_are_rows = TRUE)
TAX = tax_table(tax_table)
MET = sample_data(metadata)
physeq = phyloseq(ASV, TAX, MET)

# Filter out positive controls
physeq_noPCs <- subset_samples(physeq, control_type!="PC")

# Filter out all controls
physeq_exp <- subset_samples(physeq, control_type=="exp")

### RUN DIVNET ###

# Run DivNet at the ASV levels
# If you use X=variable, estimates will be made for each group rather than each sample

# unfiltered base asv
divnet_asv <- divnet(physeq_noPCs, base="10f86f692e9965c610fe1442ec51e246", ncores=4)

# filtered base asv
divnet_asv <- divnet(physeq_noPCs, base="ecccfdc6bc637875988540bab042e56c", ncores=6)

# filtered, experimental samples only
divnet_asv_exp <- divnet(physeq_exp, base="ecccfdc6bc637875988540bab042e56c", ncores=4)

################ ALPHA DIVERSITY #################

# Calculate alpha diversity metrics for each sample
shan = divnet_asv$shannon
simp = divnet_asv_exp$simpson

# To get diversity indices for sample groups rather than individual samples   
divnet_asv_exp <- physeq_exp %>%
  divnet(X = "sample_type", ncores = 2)

# To export values as table
shan_tib <- bind_rows(shan)
write.csv(shan_tib, file="DivNet_exp_shannon.csv")

################ BETA DIVERSITY #################

# Pull out Bray-Curtis distances as a square distance matrix
bc <- divnet_asv$'bray-curtis'
bc <- divnet_asv_exp$'bray-curtis'

# Base R PCA function from distance matrix
PCA <- prcomp(x=bc)
PCAi <- data.frame(PCA$x)

# Base R plotting for quick check
autoplot(PCA)

# copy data frame and add any metadata columns
bc_meta <- bc
# make sure the metadata is in the correct sample order
metadata <- metadata[ order(row.names(metadata)), ]

# For plot with NCs only
bc_meta <- cbind.data.frame(bc_meta, control = metadata$control_status)

# For plot without controls
bc_meta <- cbind.data.frame(bc_meta, type = metadata$sample_type)
bc_meta <- cbind.data.frame(bc_meta, season = metadata$season)
bc_meta <- cbind.data.frame(bc_meta, depth = metadata$depth_group)

##### PLOT WITH NCs #######
PCAi <- data.frame(PCA$x, group=bc_meta$control, shape=bc_meta$season)

# Shape guide for season, color guide for control/exp
season_shapes <- c("Fall" = 21, "Spring" = 22)
control_colors = c("control" = "#6a3d9a", "exp" = "#33a02c")

p <- ggplot(PCAi, aes(PC1, PC2, col=group, shape=shape)) + 
  geom_point(size = 5, stroke = 2) + 
  scale_color_manual(values = control_colors) +
  scale_shape_manual(values = season_shapes)
  
p + theme_classic() + theme(axis.ticks.x=element_blank(), axis.ticks.y=element_blank(),
          axis.text.x=element_blank(), axis.text.y=element_blank(), legend.position="none")+ 
  labs(x = "PC1 (77.5%)", y = "PC2 (17.4%)")

##### PLOT ONLY EXPERIMENTAL SAMPLES #######
PCAi <- data.frame(PCA$x, group=bc_meta$type, shape=bc_meta$depth, 
                   season=bc_meta$season)

# Shape guide for depths, color guide for CTD/ESP
depth_shapes <- c("0-25" = 21, "27-50" = 22, "195-200" = 24)
CTD_ESP_colors = c("CTD" = "#6a3d9a", "ESP" = "#33a02c")

p <- ggplot(PCAi, aes(PC1, PC2, color=group, shape=shape)) + 
  # Make fill for geom point be white if fall; if spring, then make it the color corresponding to ESP or CTD
  geom_point(size = 5, stroke = 2, fill = ifelse(PCAi$season == "Fall", "white", 
                                                 ifelse(PCAi$group == "ESP", "#33a02c", "#6a3d9a"))) +
  scale_color_manual(values = CTD_ESP_colors) +
  scale_shape_manual(values = depth_shapes)
   
p + theme_classic() + labs(x = "PC1 (82.2%)", y = "PC2 (16.0%)") +
  theme(axis.ticks.x=element_blank(), axis.ticks.y=element_blank(),
        axis.text.x=element_blank(), axis.text.y=element_blank(), legend.position="none")

## Permanova
library(vegan)

# convert bray-curtis matrix to "dist" class object
bc_dist <- as.dist(bc)

# PERMANOVA for control vs exp

permanova <- adonis2(bc_dist ~ depth_group, data = metadata, permutations=999)
permanova
