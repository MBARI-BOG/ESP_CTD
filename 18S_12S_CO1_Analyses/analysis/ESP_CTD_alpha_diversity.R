## Alpha diversity (Shannon diversity) analyses and plots

library(tidyverse)
library(phyloseq)
library(ggpubr)
library(here)

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
asv_table_path_16S <- here::here("data", "16S", "table_phyloseq.csv")
tax_table_path_16S <- here::here("data", "16S", "taxonomy_phyloseq.csv")
metadata_path_16S <- here::here("data", "16S", "ESP_CTD_16S_metadata.csv")

ESP_CTD_16S_phyloseq <- merge_phyloseq(otu_table(read.csv(asv_table_path_16S, row.names = 1),taxa_are_rows = TRUE),
                                       tax_table(as.matrix(read.csv(tax_table_path_16S, row.names = 1))),
                                       sample_data(read.csv(metadata_path_16S, row.names = 1)))

# Remove sample that has no reads
ESP_CTD_16S_phyloseq <- subset_samples(ESP_CTD_16S_phyloseq, sample_names(ESP_CTD_16S_phyloseq) != "CN18Fc35_5_eDNA")


# Set figure directory
fig_dir <- here::here("figures")


##----Extract Metadata----------------------------------------------------------

# COI
ESP_CTD_COI_envt_phyloseq <- prune_samples(sample_data(ESP_CTD_COI_phyloseq)$CTD_or_ESP %in% c("CTD","ESP"),ESP_CTD_COI_phyloseq)
ESP_CTD_COI_envtsamples_metadata <- as.data.frame(as.matrix(sample_data(ESP_CTD_COI_envt_phyloseq)))
ESP_CTD_COI_envtsamples_metadata <- rownames_to_column(ESP_CTD_COI_envtsamples_metadata, "sample_name")

# 18S
ESP_CTD_18S_envt_phyloseq <- prune_samples(sample_data(ESP_CTD_18S_phyloseq)$CTD_or_ESP %in% c("CTD","ESP"),ESP_CTD_18S_phyloseq)
ESP_CTD_18S_envtsamples_metadata <- as.data.frame(as.matrix(sample_data(ESP_CTD_18S_envt_phyloseq)))
ESP_CTD_18S_envtsamples_metadata <- rownames_to_column(ESP_CTD_18S_envtsamples_metadata, "sample_name")

# 12S
ESP_CTD_12S_envt_phyloseq <- prune_samples(sample_data(ESP_CTD_12S_phyloseq)$CTD_or_ESP %in% c("CTD","ESP"),ESP_CTD_12S_phyloseq)
ESP_CTD_12S_envtsamples_metadata <- as.data.frame(as.matrix(sample_data(ESP_CTD_12S_envt_phyloseq)))
ESP_CTD_12S_envtsamples_metadata <- rownames_to_column(ESP_CTD_12S_envtsamples_metadata, "sample_name")

# 16S
ESP_CTD_16S_envt_phyloseq <- prune_samples(sample_data(ESP_CTD_16S_phyloseq)$CTD_or_ESP %in% c("CTD","ESP"),ESP_CTD_16S_phyloseq)
ESP_CTD_16S_envtsamples_metadata <- as.data.frame(as.matrix(sample_data(ESP_CTD_16S_envt_phyloseq)))
ESP_CTD_16S_envtsamples_metadata <- rownames_to_column(ESP_CTD_16S_envtsamples_metadata, "sample_name")






##----Calculate Shannon for COI-------------------------------------------------

ESP_CTD_COI_shannon <- estimate_richness(ESP_CTD_COI_phyloseq, split = TRUE, measures = c("Shannon"))

ESP_CTD_COI_shannon %>% tibble::rownames_to_column(.,"sample_name") -> ESP_CTD_COI_shannon

# Add metadata
ESP_CTD_COI_shannon <- left_join(ESP_CTD_COI_shannon,ESP_CTD_COI_envtsamples_metadata,by = "sample_name")

ESP_CTD_COI_shannon <- ESP_CTD_COI_shannon[,c("sample_name","Shannon","CTD_or_ESP","matching_ID")]

# subset ESP and CTD samples
ESP_CTD_COI_shannon <- subset(ESP_CTD_COI_shannon, CTD_or_ESP %in% c("ESP", "CTD"))

# Add a marker name for plotting
ESP_CTD_COI_shannon$marker <- "COI"


##----Calculate Shannon for 18S-------------------------------------------------

ESP_CTD_18S_shannon <- estimate_richness(ESP_CTD_18S_phyloseq, split = TRUE, measures = c("Shannon"))

ESP_CTD_18S_shannon %>% tibble::rownames_to_column(.,"sample_name") -> ESP_CTD_18S_shannon

# Add metadata
ESP_CTD_18S_shannon <- left_join(ESP_CTD_18S_shannon,ESP_CTD_18S_envtsamples_metadata,by = "sample_name")

ESP_CTD_18S_shannon <- ESP_CTD_18S_shannon[,c("sample_name","Shannon","CTD_or_ESP","matching_ID")]

# subset ESP and CTD samples
ESP_CTD_18S_shannon <- subset(ESP_CTD_18S_shannon, CTD_or_ESP %in% c("ESP", "CTD"))

# Add a marker name for plotting
ESP_CTD_18S_shannon$marker <- "18S"



##----Calculate Shannon for 12S-------------------------------------------------

ESP_CTD_12S_shannon <- estimate_richness(ESP_CTD_12S_phyloseq, split = TRUE, measures = c("Shannon"))

ESP_CTD_12S_shannon %>% tibble::rownames_to_column(.,"sample_name") -> ESP_CTD_12S_shannon

# Add metadata
ESP_CTD_12S_shannon <- left_join(ESP_CTD_12S_shannon,ESP_CTD_12S_envtsamples_metadata,by = "sample_name")

ESP_CTD_12S_shannon <- ESP_CTD_12S_shannon[,c("sample_name","Shannon","CTD_or_ESP","matching_ID")]

# subset ESP and CTD samples
ESP_CTD_12S_shannon <- subset(ESP_CTD_12S_shannon, CTD_or_ESP %in% c("ESP", "CTD"))

# Add a marker name for plotting
ESP_CTD_12S_shannon$marker <- "12S"



##----Calculate Shannon for 16S-------------------------------------------------

ESP_CTD_16S_shannon <- estimate_richness(ESP_CTD_16S_phyloseq, split = TRUE, measures = c("Shannon"))

ESP_CTD_16S_shannon %>% tibble::rownames_to_column(.,"sample_name") -> ESP_CTD_16S_shannon

# Add metadata
ESP_CTD_16S_shannon <- left_join(ESP_CTD_16S_shannon,ESP_CTD_16S_envtsamples_metadata,by = "sample_name")

ESP_CTD_16S_shannon <- ESP_CTD_16S_shannon[,c("sample_name","Shannon","CTD_or_ESP","matching_ID")]

# subset ESP and CTD samples
ESP_CTD_16S_shannon <- subset(ESP_CTD_16S_shannon, CTD_or_ESP %in% c("ESP", "CTD"))

# Add a marker name for plotting
ESP_CTD_16S_shannon$marker <- "16S"



##----Measure significance using t-test------------------------------------------
shannon_ttest_16S <- t.test(subset(ESP_CTD_16S_shannon, CTD_or_ESP == "ESP")$Shannon, 
                            subset(ESP_CTD_16S_shannon, CTD_or_ESP == "CTD")$Shannon)

shannon_ttest_18S <- t.test(subset(ESP_CTD_18S_shannon, CTD_or_ESP == "ESP")$Shannon, 
                            subset(ESP_CTD_18S_shannon, CTD_or_ESP == "CTD")$Shannon)

shannon_ttest_COI <- t.test(subset(ESP_CTD_COI_shannon, CTD_or_ESP == "ESP")$Shannon, 
                            subset(ESP_CTD_COI_shannon, CTD_or_ESP == "CTD")$Shannon)

shannon_ttest_12S <- t.test(subset(ESP_CTD_12S_shannon, CTD_or_ESP == "ESP")$Shannon, 
                            subset(ESP_CTD_12S_shannon, CTD_or_ESP == "CTD")$Shannon)

# All are > 0.05

##----Join all markers together-------------------------------------------------

ESP_CTD_shannon <- rbind(ESP_CTD_COI_shannon, ESP_CTD_18S_shannon, ESP_CTD_12S_shannon, ESP_CTD_16S_shannon)

# Change CTD to "Shipboard" and ESP to "Autonomous"
ESP_CTD_shannon %>%
  mutate(CTD_or_ESP = ifelse(CTD_or_ESP == "CTD", "Shipboard", "Autonomous")) -> ESP_CTD_shannon

##----Violin Plots--------------------------------------------------------------

CTD_ESP_colors = c("Shipboard" = "#6a3d9a", "Autonomous" = "#33a02c")

## Plot four markers in a for loop
markers <- c("16S", "18S", "COI", "12S")

for (i in 1:length(markers)){
  data = subset(ESP_CTD_shannon, marker == markers[i])
  plot <- ggplot(data, aes(CTD_or_ESP, Shannon, colour=CTD_or_ESP)) + 
    geom_point(size=1, position=position_dodge(width=1)) +
    geom_violin(alpha=0, position=position_dodge(width=1)) + 
    scale_color_manual(values = CTD_ESP_colors) + 
    ylim(0,7) +
    theme(legend.position = "none",
          panel.background = element_rect(fill = "white"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line.x = element_line(size = 0.5, color = "black"),
          axis.line.y = element_line(size = ifelse(i %in% c(2,4), 0, 0.5), color = "black"),
          axis.ticks.length.x = unit(0,"cm"),
          axis.ticks.length.y = unit(ifelse(i %in% c(2,4), 0, 0.2),"cm"),
          axis.text.y = element_text(size = ifelse(i %in% c(2,4), 0, 15)),
          axis.text.x = element_text(size = ifelse(i %in% c(1,2), 0, 15)),
          axis.title = element_blank())+
    annotate(geom = "text", label = paste0(markers[i], ": ", "p > 0.05"), x = 2.4, y = 6.8, size = 5, hjust = 1)
  assign(paste0("shannon_plot_",markers[i]), plot)  
}



## Arrange all plots

shannon_plot_combined <- annotate_figure(ggpubr::ggarrange(shannon_plot_16S, shannon_plot_18S,shannon_plot_COI,shannon_plot_12S, 
                                          ncol = 2, nrow = 2,labels = c("(a)", "(b)", "(c)", "(d)"), 
                                          label.x = 0.1, label.y = 0.985,
                                          font.label = list(size = 20, color = "black", face = "bold")),
                                         # annotate_figure
                                         left = text_grob("Shannon Index", color = "black", size = 15, rot = 90),
                                         bottom = text_grob("Sampling Method", size = 15))


ggsave(shannon_plot_combined, height = 8, width = 8, path = fig_dir, filename = "ESP_CTD_shannon_plot.png", device = "png")

