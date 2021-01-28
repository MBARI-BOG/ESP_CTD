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
metadata <- column_to_rownames(eCruises2018_metadata_phyloseq_noPCs, 'sample_name')

# convert all to matrices
microbe <- as.matrix(microbe)
tax_table <- as.matrix(tax_table)
metadata <- as.data.frame(metadata)

# combine into phyloseq object
ASV = otu_table(microbe, taxa_are_rows = TRUE)
TAX = tax_table(tax_table)
MET = sample_data(metadata)
physeq = phyloseq(ASV, TAX, MET)

# Filter out controls
physeq_exp <- subset_samples(physeq, control_type=="exp")

# Filter out positive controls
physeq_exp <- subset_samples(physeq, control_type!="PC")

### ALPHA DIVERSITY ###

## 1. Breakaway###

# Aggregate taxa at the order, family, or genus level
water <- physeq %>% tax_glom("Order")
water <- physeq %>% tax_glom("Family")
water <- physeq %>% tax_glom("Genus") 

ba <- breakaway(water)

plot(ba, water, color = "control_type")

### 2. DivNet ###

# Run DivNet at various taxonomic levels: phylum, genus, and ASV
# If you use X=variable, estimates will be made for each group rather than each sample

divnet_phylum <- divnet(tax_glom(physeq, taxrank="Phylum"), X="Location", ncores=6)
divnet_genus <- divnet(tax_glom(physeq, taxrank="Genus"), ncores=4) # X="control_type"

# unfiltered base asv
divnet_asv <- divnet(physeq, base="10f86f692e9965c610fe1442ec51e246", ncores=4)

#filtered base asv
divnet_asv <- divnet(physeq, base="ecccfdc6bc637875988540bab042e56c", ncores=4)

#filtered, experimental samples only
divnet_asv_exp <- divnet(physeq_exp, base="ecccfdc6bc637875988540bab042e56c", ncores=4)

divnet_asv %>% names
physeq %>% sample_data

################ ALPHA DIVERSITY #################

# Calculate alpha diversity metrics for each sample
shan = divnet_asv_exp$shannon
simp = divnet_asv_exp$simpson

# Plot each sample metric
plot(divnet_asv_exp$shannon, 
     physeq_exp, 
     col = "sample_type")

plot(divnet_asv_exp$simpson, 
     physeq_exp, 
     col = "sample_type")

# To get diversity indices for sample groups rather than individual samples   
divnet_exp <- physeq_exp %>%
  divnet(X = "sample_type", ncores = 2)

# To export values as table
shan_tib <- bind_rows(shan)
write.csv(shan_tib, file="DivNet_exp_shannon.csv")

testDiversity(divnet_asv, "shannon")

################ BETA DIVERSITY #################

simplifyBeta(divnet_asv, physeq_exp, "bray-curtis", "control_type")

simplifyBeta(divnet_genus, physeq_exp, "bray-curtis", "control_type") %>%
  ggplot(aes(x = interaction(Covar1, Covar2), 
             y = beta_est,
             col = interaction(Covar1, Covar2))) +
  geom_point() +
  geom_linerange(aes(ymin = lower, ymax = upper)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("") + ylab("Estimates of Bray-Curtis distance")

# Pull out Bray-Curtis distances as a square distance matrix
bc <- divnet_asv_exp$'bray-curtis'

# Base R PCA function from distance matrix
PCA <- prcomp(x=bc)
PCAi <- data.frame(PCA$x)

# Base R plotting for quick check
autoplot(PCA)
p <- ggplot(PCAi, aes(x=PC1, y=PC2)) + geom_point(size=2)

# Pretty plot in ggplot2
# copy data frame and add any metadata columns
bc_meta <- bc
# make sure the metadata is in the correct sample order
metadata <- metadata[ order(row.names(metadata)), ]

bc_meta <- cbind.data.frame(bc_meta, control = metadata$sample_type)
bc_meta <- cbind.data.frame(bc_meta, control = metadata$control_type)
bc_meta <- cbind.data.frame(bc_meta, sampler = metadata$sampling_device)
bc_meta <- cbind.data.frame(bc_meta, season = metadata$season)
bc_meta <- cbind.data.frame(bc_meta, depth = metadata$depth_group)
bc_meta <- cbind.data.frame(bc_meta, name = metadata$ESP_CTD_match)

# For plot without any controls
PCAi <- data.frame(PCA$x, group=bc_meta$depth, shape=bc_meta$season, 
                   size=bc_meta$sample, name=bc_meta$name)

# For plot with NCs
PCAi <- data.frame(PCA$x, group=bc_meta$control, shape=bc_meta$season, 
                   name=bc_meta$name)

# Make a custom color pallette
cbp1 <- c("#009E73","#F0E442", "#0072B2") # ", "#D55E00", "#CC79A7", "#999999", "#E69F00", "#56B4E9"

# For plot with NCs
p <- ggplot(PCAi, aes(PC1, PC2, col=group, shape=shape)) + geom_point(alpha=0.75, size=3) + 
  scale_size_manual(values=c(3, 4)) + scale_color_brewer(palette = "Set1") + 
  labs(col="Sample Type", shape="Season") +
  guides(color = guide_legend(override.aes = list(size = 3)), 
         shape = guide_legend(override.aes = list(size = 3)))
    
p + theme_classic() + labs(x = "PC1 (77.5%)", y = "PC2 (17.4%)")
  # + geom_text(aes(label=name),hjust=0, vjust=0)

## Rarefaction curves with phyloseq object

library(ranacapa)

sample_variables(physeq)
physeq %>% sample_names

r <- ggrare(physeq, step=1000, se=FALSE, color="control_type")
r <- ggrare(physeq, step=1000, se=FALSE, color="control_type", label=sample_names(physeq))
r <- ggrare(physeq, step=1000, se=FALSE, color="control_type", label="sampling_device")
  
r + scale_color_brewer(palette="Set1") + labs(col="Sample Type") +
  theme_classic() + scale_x_continuous(limits = c(0, 350000), expand = c(0, 0)) + 
  scale_y_continuous(limits = c(0,1200), expand = c(0, 0))

