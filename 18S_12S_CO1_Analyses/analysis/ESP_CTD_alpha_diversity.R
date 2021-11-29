## Alpha diversity (Shannon diversity) analyses and plots

library(tidyverse)
library(phyloseq)
library(ggpubr)
library(here)
library(rstatix)
library(egg)
library(grid)

fig_dir <- here("resubmission")

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



##----Check all markers for normality------------------------------------------
qqplot(ESP_CTD_COI_shannon$Shannon, runif(1000))
qqplot(ESP_CTD_18S_shannon$Shannon, runif(1000))
qqplot(ESP_CTD_12S_shannon$Shannon, runif(1000))
qqplot(ESP_CTD_16S_shannon$Shannon, runif(1000))





##----Measure significance using t-test------------------------------------------

# UNPAIRED
shannon_ttest_16S <- t.test(subset(ESP_CTD_16S_shannon, CTD_or_ESP == "ESP")$Shannon, 
                            subset(ESP_CTD_16S_shannon, CTD_or_ESP == "CTD")$Shannon)

shannon_ttest_18S <- t.test(subset(ESP_CTD_18S_shannon, CTD_or_ESP == "ESP")$Shannon, 
                            subset(ESP_CTD_18S_shannon, CTD_or_ESP == "CTD")$Shannon)

shannon_ttest_COI <- t.test(subset(ESP_CTD_COI_shannon, CTD_or_ESP == "ESP")$Shannon, 
                            subset(ESP_CTD_COI_shannon, CTD_or_ESP == "CTD")$Shannon)

shannon_ttest_12S <- t.test(subset(ESP_CTD_12S_shannon, CTD_or_ESP == "ESP")$Shannon, 
                            subset(ESP_CTD_12S_shannon, CTD_or_ESP == "CTD")$Shannon)


# PAIRED

# Confirm that paired samples are in the same order
ESP_CTD_COI_shannon %>% 
  arrange(matching_ID) -> ESP_CTD_COI_shannon
subset(ESP_CTD_COI_shannon, CTD_or_ESP == "ESP")$matching_ID
subset(ESP_CTD_COI_shannon, CTD_or_ESP == "CTD")$matching_ID

ESP_CTD_18S_shannon %>% 
  arrange(matching_ID) -> ESP_CTD_18S_shannon
subset(ESP_CTD_18S_shannon, CTD_or_ESP == "ESP")$matching_ID
subset(ESP_CTD_18S_shannon, CTD_or_ESP == "CTD")$matching_ID

ESP_CTD_12S_shannon %>% 
  arrange(matching_ID) -> ESP_CTD_12S_shannon
subset(ESP_CTD_12S_shannon, CTD_or_ESP == "ESP")$matching_ID
subset(ESP_CTD_12S_shannon, CTD_or_ESP == "CTD")$matching_ID
# 12S has the issue where there are some missing CTD samples - remove the following:
setdiff(subset(ESP_CTD_12S_shannon, CTD_or_ESP == "ESP")$matching_ID, 
        subset(ESP_CTD_12S_shannon, CTD_or_ESP == "CTD")$matching_ID)
# CN18F1, CN18F2, CN18F21, CN18F3

ESP_CTD_16S_shannon %>% 
  arrange(matching_ID) -> ESP_CTD_16S_shannon
subset(ESP_CTD_16S_shannon, CTD_or_ESP == "ESP")$matching_ID
subset(ESP_CTD_16S_shannon, CTD_or_ESP == "CTD")$matching_ID

shannon_ttest_16S_paired <- t.test(subset(ESP_CTD_16S_shannon, CTD_or_ESP == "ESP")$Shannon, 
                            subset(ESP_CTD_16S_shannon, CTD_or_ESP == "CTD")$Shannon,
                            paired = TRUE)

shannon_ttest_18S_paired <- t.test(subset(ESP_CTD_18S_shannon, CTD_or_ESP == "ESP")$Shannon, 
                            subset(ESP_CTD_18S_shannon, CTD_or_ESP == "CTD")$Shannon,
                            paired = TRUE)

shannon_ttest_COI_paired <- t.test(subset(ESP_CTD_COI_shannon, CTD_or_ESP == "ESP")$Shannon, 
                            subset(ESP_CTD_COI_shannon, CTD_or_ESP == "CTD")$Shannon,
                            paired = TRUE)

shannon_ttest_12S_paired <- t.test(subset(ESP_CTD_12S_shannon, CTD_or_ESP == "ESP" & 
                                     !(matching_ID %in% c("CN18F1", "CN18F2", "CN18F21", "CN18F3")))$Shannon, 
                            subset(ESP_CTD_12S_shannon, CTD_or_ESP == "CTD" & 
                                     !(matching_ID %in% c("CN18F1", "CN18F2", "CN18F21", "CN18F3")))$Shannon,
                            paired = TRUE)

shannon_ttest_COI_paired
shannon_ttest_18S_paired
shannon_ttest_12S_paired
shannon_ttest_16S_paired

# 12S is significant (p = 0.006), all others are not
# 12S: mean of the differences = -0.2826637 (ESP vs. CTD)
mean(subset(ESP_CTD_12S_shannon, CTD_or_ESP == "ESP")$Shannon)
mean(subset(ESP_CTD_12S_shannon, CTD_or_ESP == "CTD")$Shannon)
# 12S: CTD has mean difference of 0.28 more Shannon diversity


##----Add metadata to shannon stats-------------------------------------------------
ESP_CTD_COI_shannon %>%
  left_join(., ESP_CTD_COI_envtsamples_metadata, by = "sample_name") -> ESP_CTD_COI_shannon
ESP_CTD_18S_shannon %>%
  left_join(., ESP_CTD_18S_envtsamples_metadata, by = "sample_name") -> ESP_CTD_18S_shannon
ESP_CTD_12S_shannon %>%
  left_join(., ESP_CTD_12S_envtsamples_metadata, by = "sample_name") -> ESP_CTD_12S_shannon
ESP_CTD_16S_shannon %>%
  left_join(., ESP_CTD_16S_envtsamples_metadata, by = "sample_name") -> ESP_CTD_16S_shannon


##----Join all markers together-------------------------------------------------

# ESP_CTD_shannon <- rbind(ESP_CTD_COI_shannon, ESP_CTD_18S_shannon, ESP_CTD_12S_shannon, ESP_CTD_16S_shannon)
ESP_CTD_COI_shannon %>% 
  bind_rows(., ESP_CTD_18S_shannon) %>% 
  bind_rows(., ESP_CTD_12S_shannon) %>% 
  bind_rows(., ESP_CTD_16S_shannon) -> ESP_CTD_shannon

# Change CTD to "Shipboard" and ESP to "Autonomous"
# Change depth to shallow_deep
ESP_CTD_shannon %>%
  mutate(CTD_or_ESP = ifelse(CTD_or_ESP.x == "CTD", "Shipboard", "Autonomous")) %>%
  # mutate(CTD_or_ESP = ifelse(CTD_or_ESP == "CTD", "Shipboard", "Autonomous")) %>% 
  mutate(shallow_deep = ifelse(depth < 100, "Shallow", "Deep")) %>% 
  mutate(shallow_deep = ifelse(is.na(shallow_deep), ifelse(Depth == "Deep_200m", "Deep", "Shallow"), shallow_deep)) -> ESP_CTD_shannon

##----Violin Plots--------------------------------------------------------------

cruise_colors <- c("CN18F" = "#ff7f00", "CN18S" = "#1f78b4")
depth_shapes <- c("Shallow" = 19, "Deep" = 1)

## Plot three markers in a for loop; 12S is separate because it's significant
# markers <- c("16S", "18S", "COI", "12S")
markers <- c("16S", "18S", "COI")

for (i in 1:length(markers)){
  data = subset(ESP_CTD_shannon, marker == markers[i])
  plot <- ggplot(data, aes(CTD_or_ESP, Shannon, color=SAMPLING_cruise, shape = shallow_deep)) + 
    geom_point(size=1.5, position=position_dodge(width=1)) +
    scale_shape_manual(values = depth_shapes) +
    geom_violin(alpha=0, position=position_dodge(width=1)) + 
    scale_color_manual(values = cruise_colors) + 
    # ylim(0,7) +
    theme(legend.position = "none",
          panel.background = element_rect(fill = "white"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line.x = element_line(size = 0.5, color = "black"),
          # axis.line.y = element_line(size = ifelse(i %in% c(2,4), 0, 0.5), color = "black"),
          axis.line.y = element_line(size = 0.5, color = "black"),
          axis.ticks.length.x = unit(0,"cm"),
          # axis.ticks.length.y = unit(ifelse(i %in% c(2,4), 0, 0.2),"cm"),
          axis.ticks.length.y = unit(0.2,"cm"),
          axis.text.x = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          axis.title = element_blank())+
    # annotate(geom = "text", label = paste0(markers[i], ": ", "p > 0.05"), x = 2.4, y = max(data$Shannon)+0.3, size = 5, hjust = 1)
    annotate(geom = "text", label = paste0("p > 0.05"), x = 2.4, y = max(data$Shannon)+0.3, size = 5, hjust = 1)
  assign(paste0("shannon_plot_",markers[i]), plot)  
}

data = subset(ESP_CTD_shannon, marker == "12S")
shannon_plot_12S <- ggplot(data, aes(CTD_or_ESP, Shannon, color=SAMPLING_cruise, shape = shallow_deep)) + 
  scale_shape_manual(values = depth_shapes) +
  geom_point(size=1.5, position=position_dodge(width=1)) +
  geom_violin(alpha=0, position=position_dodge(width=1)) + 
  scale_color_manual(values = cruise_colors) + 
  # ylim(0,7) +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(size = 0.5, color = "black"),
        # axis.line.y = element_line(size = ifelse(i %in% c(2,4), 0, 0.5), color = "black"),
        axis.line.y = element_line(size = 0.5, color = "black"),
        axis.ticks.length.x = unit(0,"cm"),
        # axis.ticks.length.y = unit(ifelse(i %in% c(2,4), 0, 0.2),"cm"),
        axis.ticks.length.y = unit(0.2,"cm"),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        # axis.text.y = element_text(size = ifelse(i %in% c(2,4), 0, 15)),
        # axis.text.x = element_text(size = ifelse(i %in% c(1,2), 0, 15)),
        axis.title = element_blank())+
  # annotate(geom = "text", label = paste0("12S", ": ", "p = 0.006"), x = 2.4, y = max(data$Shannon)+0.3, size = 5, hjust = 1)
  annotate(geom = "text", label = paste0("p = 0.006"), x = 2.4, y = max(data$Shannon)+0.3, size = 5, hjust = 1)


## Arrange all plots

shannon_plot_combined <- annotate_figure(ggpubr::ggarrange(shannon_plot_16S, shannon_plot_18S,shannon_plot_COI,shannon_plot_12S,
                                          ncol = 2, nrow = 2,labels = c("(a)", "(b)", "(c)", "(d)"), 
                                          label.x = 0.1, label.y = 0.985,
                                          font.label = list(size = 20, color = "black", face = "plain")),
                                         # annotate_figure
                                         left = text_grob("Shannon Index", color = "black", size = 15, rot = 90),
                                         bottom = text_grob("Sampling Method", size = 15))

shannon_plot_combined <- annotate_figure(egg::ggarrange(shannon_plot_16S, shannon_plot_18S,shannon_plot_COI,shannon_plot_12S, 
                                                        ncol = 2, nrow = 2,labels = c("(a) 16S", "(b) 18S", "(c) COI", "(d) 12S"),
                                                      label.args = list(gp=gpar(font=1, cex = 1.5), x=unit(1.5,"cm"), y=unit(8.75,"cm"), hjust = 0)),
                                         # annotate_figure
                                         left = text_grob("Shannon Index", color = "black", size = 15, rot = 90),
                                         bottom = text_grob("Sampling Method", size = 15))

# Make a legend manually
shannon_legend_violinplot <- ggplot(ESP_CTD_shannon) +
  # Season and depth
  annotate("text",label = "Season/Depth:", x = 1, y = 1.1,size = 6, adj = 0)+ # Title 
  annotate("point", x = 3.2, y = 1, shape = 21, colour = "#1f78b4", fill = "#1f78b4", size = 5, stroke = 3)+ # circle point
  annotate("text",label = "Spring (shallow)", x = 3.5, y = 1,size = 4,adj = 0)+
  annotate("point", x = 5.4, y = 1, shape = 1, colour = "#1f78b4",  size = 6, stroke = 2)+ # circle open point
  annotate("text",label = "Spring (deep)", x = 5.7, y = 1,size = 4,adj = 0)+
  annotate("point", x = 7.6, y = 1, shape = 21, colour = "#ff7f00", fill = "#ff7f00", size = 5, stroke = 3)+ # circle point
  annotate("text",label = "Fall (shallow)", x = 7.9, y = 1,size = 4,adj = 0)+
  
  xlim(1,9)+
  # ylim(0.2,3.2)+
  ylim(0,2)+
  theme(panel.background = element_rect(fill="white"),
        panel.border = element_rect(colour = "black", fill=NA, size=2),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none",
        plot.margin = unit(c(0.2,1,0,1),"cm"))


shannon_violin_plot <- ggpubr::ggarrange(shannon_plot_combined, shannon_legend_violinplot, ncol = 1,
                                                    heights = c(12,1))


# ggsave(shannon_plot_combined, height = 8, width = 8, path = fig_dir, filename = "ESP_CTD_shannon_plot.png", device = "png")

# Resave according to figure guidelines of eDNA journal
final_fig_dir <- here("resubmission")
ggsave(shannon_violin_plot, height = 8, width = 8, path = final_fig_dir, filename = "Fig4_alpha_diversity_v3.pdf", device = "pdf")



##### Run three-way ANOVA #####

##### COI #####

# Add metadata to all data frames of shannon div
ESP_CTD_COI_shannon %>% 
  # left_join(., dplyr::select(ESP_CTD_COI_envtsamples_metadata, c(sample_name, SAMPLING_cruise, depth)), by = "sample_name") %>%
  mutate(depth = as.numeric(depth)) %>% 
  mutate(shallow_deep = ifelse(depth > 50, "deep", "shallow")) -> ESP_CTD_COI_shannon_meta


# Test ANOVA assumptions
# Plot data distributions
ggplot(ESP_CTD_COI_shannon_meta, aes(x = Shannon)) +
  geom_histogram()

# Build the linear model
model  <- lm(Shannon ~ CTD_or_ESP.x, data = ESP_CTD_COI_shannon_meta)
# Create a QQ plot of residuals
ggqqplot(residuals(model))

# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))

# Not normally distributed

# Run the ANOVA anyway
res.aov.COI <- aov(Shannon ~ SAMPLING_cruise + shallow_deep + CTD_or_ESP.x,  data = ESP_CTD_COI_shannon_meta)
summary(res.aov.COI)


##### 18S #####

# Add metadata to all data frames of shannon div
ESP_CTD_18S_shannon %>% 
  # left_join(., dplyr::select(ESP_CTD_18S_envtsamples_metadata, c(sample_name, SAMPLING_cruise, depth)), by = "sample_name") %>% 
  mutate(depth = as.numeric(depth)) %>% 
  mutate(shallow_deep = ifelse(depth > 50, "deep", "shallow")) -> ESP_CTD_18S_shannon_meta


# Test ANOVA assumptions
# Plot data distributions
ggplot(ESP_CTD_18S_shannon_meta, aes(x = Shannon)) +
  geom_histogram()

# Build the linear model
model  <- lm(Shannon ~ CTD_or_ESP.x, data = ESP_CTD_18S_shannon_meta)
# Create a QQ plot of residuals
ggqqplot(residuals(model))

# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))

# Not normally distributed

# Run the ANOVA anyway
res.aov.18S <- aov(Shannon ~ SAMPLING_cruise + shallow_deep + CTD_or_ESP.x,  data = ESP_CTD_18S_shannon_meta)
summary(res.aov.18S)


##### 12S #####

# Add metadata to all data frames of shannon div
ESP_CTD_12S_shannon %>% 
  # left_join(., dplyr::select(ESP_CTD_12S_envtsamples_metadata, c(sample_name, SAMPLING_cruise, depth)), by = "sample_name") %>% 
  mutate(depth = as.numeric(depth)) %>% 
  mutate(shallow_deep = ifelse(depth > 50, "deep", "shallow")) -> ESP_CTD_12S_shannon_meta


# Test ANOVA assumptions
# Plot data distributions
ggplot(ESP_CTD_12S_shannon_meta, aes(x = Shannon)) +
  geom_histogram()

# Build the linear model
model  <- lm(Shannon ~ CTD_or_ESP.x, data = ESP_CTD_12S_shannon_meta)
# Create a QQ plot of residuals
ggqqplot(residuals(model))

# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))

# Not normally distributed

# Run the ANOVA anyway
res.aov.12S <- aov(Shannon ~ SAMPLING_cruise + shallow_deep + CTD_or_ESP.x,  data = ESP_CTD_12S_shannon_meta)
summary(res.aov.12S)


##### 16S #####

# Add metadata to all data frames of shannon div
ESP_CTD_16S_shannon %>% 
  # left_join(., dplyr::select(ESP_CTD_16S_envtsamples_metadata, c(sample_name, SAMPLING_cruise, Depth)), by = "sample_name") %>% 
  mutate(Depth = as.numeric(parse_number(Depth))) %>% 
  mutate(shallow_deep = ifelse(Depth > 50, "deep", "shallow")) -> ESP_CTD_16S_shannon_meta


# Test ANOVA assumptions
# Plot data distributions
ggplot(ESP_CTD_16S_shannon_meta, aes(x = Shannon)) +
  geom_histogram()

# Build the linear model
model  <- lm(Shannon ~ CTD_or_ESP.x, data = ESP_CTD_16S_shannon_meta)
# Create a QQ plot of residuals
ggqqplot(residuals(model))

# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))

# Not normally distributed

# Run the ANOVA anyway
res.aov.16S <- aov(Shannon ~ SAMPLING_cruise + shallow_deep + CTD_or_ESP.x,  data = ESP_CTD_16S_shannon_meta)
summary(res.aov.16S)


##### Examine all ANOVA results #####

summary(res.aov.COI)
summary(res.aov.18S)
summary(res.aov.12S)
summary(res.aov.16S)


##----Shannon paired plot----

# Reformat data function
reformat_shannon_data <- function(input_data){
  input_data %>% 
    dplyr::select(Shannon, CTD_or_ESP.x, matching_ID.x, SAMPLING_cruise, shallow_deep) %>% 
    dplyr::rename(CTD_or_ESP = CTD_or_ESP.x, matching_ID = matching_ID.x) %>% 
    pivot_longer(., cols = "Shannon", names_to = "Shannon") %>%
    dplyr::select(-Shannon) %>% 
    group_by(matching_ID) %>% 
    pivot_wider(names_from = "CTD_or_ESP", values_from = "value") -> output_data
  return(output_data)
}

# Plot data function
plot_shannon_paired <- function(input_data, marker){
  # Remove all rows with NA values
  input_data <- subset(input_data, !(is.na(ESP)) & !(is.na(CTD)))
  
  
  # Get R^2 value
  shannon_lm <- lm(formula = ESP~CTD, data = input_data)
  lm_summary <- summary(shannon_lm)
  r_squared <- lm_summary$r.squared
  r_squared <- round(r_squared, 4)
  
  
  # Generate plot
  cruise_colors <- c("CN18F" = "#ff7f00", "CN18S" = "#1f78b4")
  depth_shapes <- c("shallow" = 19, "deep" = 17)
  # Generate plot
  shannon_gg <- ggplot(input_data,aes(x = CTD, y = ESP, color = SAMPLING_cruise, shape = shallow_deep))+
    geom_point(size = 2.5)+
    scale_color_manual(values = cruise_colors)+
    scale_shape_manual(values = depth_shapes)+
    geom_abline(intercept = 0, slope = 1,lty = 4,color = "gold2", size = 1)+
    geom_smooth(aes(group = 1), method='lm', formula= y~x, color = "black", lty = 2,se = FALSE)+
    # annotate("text",x = 2.5, y =5.75, label = expression(paste("R"^"2"," = "), eval(parse(text = "r_squared"))),size = 6,hjust = 0)+ 
    annotate("text",x = round(min(input_data$CTD-0.5), 0)+0.5, y = round(max(input_data$ESP+0.5), 0)-0.75, label = as.expression(bquote(R^2 ~ "=" ~ .(r_squared))),size = 6,hjust = 0)+
    ylab("Autonomous\n")+
    xlab("\nShipboard")+
    scale_y_continuous(breaks = seq(round(min(input_data$ESP-0.5), 0), round(max(input_data$ESP+0.5), 0), 1), 
                       limits = c(round(min(input_data$ESP-0.5), 0), round(max(input_data$ESP+0.5), 0)), expand = c(0,0))+
    scale_x_continuous(breaks = seq(round(min(input_data$CTD-0.5), 0), round(max(input_data$CTD+0.5), 0), 1), 
                       limits = c(round(min(input_data$CTD-0.5), 0), round(max(input_data$CTD+0.5), 0)), expand = c(0,0))+
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          plot.title = element_text(face = "bold",hjust = 0.5,size = 18),
          axis.text = element_text(size = 12),
          axis.title = element_text(face = "bold", size = 15),
          axis.ticks.length = unit(0.3, "cm"),
          legend.position = "none")
    # ggtitle(paste0("Shannon Diversity: ", marker))
  
  return(shannon_gg)
}

# Run through all of the markers

ESP_CTD_COI_shannon_paired <- reformat_shannon_data(input_data = ESP_CTD_COI_shannon_meta)
ESP_CTD_18S_shannon_paired <- reformat_shannon_data(input_data = ESP_CTD_18S_shannon_meta) 
ESP_CTD_12S_shannon_paired <- reformat_shannon_data(input_data = ESP_CTD_12S_shannon_meta) 
ESP_CTD_16S_shannon_paired <- reformat_shannon_data(input_data = ESP_CTD_16S_shannon_meta) 

shannon_paired_COI <- plot_shannon_paired(input_data = ESP_CTD_COI_shannon_paired, marker = "COI")
shannon_paired_18S <- plot_shannon_paired(input_data = ESP_CTD_18S_shannon_paired, marker = "18S")
shannon_paired_12S <- plot_shannon_paired(input_data = ESP_CTD_12S_shannon_paired, marker = "12S")
shannon_paired_16S <- plot_shannon_paired(input_data = ESP_CTD_16S_shannon_paired, marker = "16S")

# Plot 16S differently because the limits are so different
# Remove all rows with NA values
input_data <- ESP_CTD_16S_shannon_paired
input_data <- subset(input_data, !(is.na(ESP)) & !(is.na(CTD)))


# Get R^2 value
shannon_lm <- lm(formula = ESP~CTD, data = input_data)
lm_summary <- summary(shannon_lm)
r_squared <- lm_summary$r.squared
r_squared <- round(r_squared, 4)


# Generate plot
cruise_colors <- c("CN18F" = "#ff7f00", "CN18S" = "#1f78b4")
depth_shapes <- c("shallow" = 19, "deep" = 17)
# Generate plot
shannon_paired_16S <- ggplot(input_data,aes(x = CTD, y = ESP, color = SAMPLING_cruise, shape = shallow_deep))+
  geom_point(size = 2.5)+
  scale_color_manual(values = cruise_colors)+
  scale_shape_manual(values = depth_shapes)+
  geom_abline(intercept = 0, slope = 1,lty = 4,color = "gold2", size = 1)+
  geom_smooth(aes(group = 1), method='lm', formula= y~x, color = "black", lty = 2,se = FALSE)+
  # annotate("text",x = 2.5, y =5.75, label = expression(paste("R"^"2"," = "), eval(parse(text = "r_squared"))),size = 6,hjust = 0)+ 
  annotate("text",x = 3.9, y = 4.9, label = as.expression(bquote(R^2 ~ "=" ~ .(r_squared))),size = 6,hjust = 0)+
  ylab("Autonomous\n")+
  xlab("\nShipboard")+
  # scale_y_continuous(breaks = seq(round(min(input_data$ESP-0.5), 0), round(max(input_data$ESP+0.5), 0), 1), 
  #                    limits = c(round(min(input_data$ESP-0.5), 0), round(max(input_data$ESP+0.5), 0)), expand = c(0,0))+
  # scale_x_continuous(breaks = seq(round(min(input_data$CTD-0.5), 0), round(max(input_data$CTD+0.5), 0), 1), 
  #                    limits = c(round(min(input_data$CTD-0.5), 0), round(max(input_data$CTD+0.5), 0)), expand = c(0,0))+
  scale_y_continuous(breaks = seq(4, 5, 0.5), limits = c(3.8, 5.2), expand = c(0,0))+
  scale_x_continuous(breaks = seq(4, 5, 0.5), limits = c(3.8, 5.2), expand = c(0,0))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(face = "bold",hjust = 0.5,size = 18),
        axis.text = element_text(size = 12),
        axis.title = element_text(face = "bold", size = 15),
        axis.ticks.length = unit(0.3, "cm"),
        legend.position = "none")
  # ggtitle(paste0("Shannon Diversity: 16S"))

# Plot 12S differently because the limits are so different
# Remove all rows with NA values
input_data <- ESP_CTD_12S_shannon_paired
input_data <- subset(input_data, !(is.na(ESP)) & !(is.na(CTD)))


# Get R^2 value
shannon_lm <- lm(formula = ESP~CTD, data = input_data)
lm_summary <- summary(shannon_lm)
r_squared <- lm_summary$r.squared
r_squared <- round(r_squared, 4)


# Generate plot
cruise_colors <- c("CN18F" = "#ff7f00", "CN18S" = "#1f78b4")
depth_shapes <- c("shallow" = 19, "deep" = 17)
# Generate plot
shannon_paired_12S <- ggplot(input_data,aes(x = CTD, y = ESP, color = SAMPLING_cruise, shape = shallow_deep))+
  geom_point(size = 2.5)+
  scale_color_manual(values = cruise_colors)+
  scale_shape_manual(values = depth_shapes)+
  geom_abline(intercept = 0, slope = 1,lty = 4,color = "gold2", size = 1)+
  geom_smooth(aes(group = 1), method='lm', formula= y~x, color = "black", lty = 2,se = FALSE)+
  # annotate("text",x = 2.5, y =5.75, label = expression(paste("R"^"2"," = "), eval(parse(text = "r_squared"))),size = 6,hjust = 0)+ 
  annotate("text",x = 0.38, y = 2.45, label = as.expression(bquote(R^2 ~ "=" ~ .(r_squared))),size = 6,hjust = 0)+
  ylab("Autonomous\n")+
  xlab("\nShipboard")+
  # scale_y_continuous(breaks = seq(round(min(input_data$ESP-0.5), 0), round(max(input_data$ESP+0.5), 0), 1), 
  #                    limits = c(round(min(input_data$ESP-0.5), 0), round(max(input_data$ESP+0.5), 0)), expand = c(0,0))+
  # scale_x_continuous(breaks = seq(round(min(input_data$CTD-0.5), 0), round(max(input_data$CTD+0.5), 0), 1), 
  #                    limits = c(round(min(input_data$CTD-0.5), 0), round(max(input_data$CTD+0.5), 0)), expand = c(0,0))+
  scale_y_continuous(breaks = seq(0,3,1), limits = c(0,3), expand = c(0,0))+
  scale_x_continuous(breaks = seq(0,3,1), limits = c(0,3), expand = c(0,0))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(face = "bold",hjust = 0.5,size = 18),
        axis.text = element_text(size = 12),
        axis.title = element_text(face = "bold", size = 15),
        axis.ticks.length = unit(0.3, "cm"),
        legend.position = "none")
# ggtitle(paste0("Shannon Diversity: 12S"))


shannon_paired_COI
shannon_paired_18S
shannon_paired_12S
shannon_paired_16S

# Join plots, export

# shannon_paired_combined <- egg::ggarrange(shannon_paired_16S, shannon_paired_18S, shannon_paired_COI, shannon_paired_12S,
#                                      ncol = 2, nrow = 2,
#                                      labels = c("(a) 16S", "(b) 18S", "(c) COI", "(d) 12S"), label.x = 0.15, label.y = 0.985,
#                                      font.label = list(size = 24, color = "black", face = "plain"), hjust = -0.12)

shannon_paired_combined <- egg::ggarrange(shannon_paired_16S, shannon_paired_18S, shannon_paired_COI, shannon_paired_12S, ncol = 2, nrow = 2,
                                          labels = c("(a) 16S", "(b) 18S", "(c) COI", "(d) 12S"),
                                          label.args = list(gp=gpar(font=1, cex = 2), x=unit(2.7,"cm"), y=unit(11,"cm"), hjust = 0))


# Create a legend manually
shannon_legend <- ggplot(ESP_CTD_COI_shannon_paired) +
  # Season
  annotate("text",label = "Season:", x = 1, y = 2,size = 7, adj = 0)+ # Title 
  annotate("point", x = 4, y = 2, shape = 21, colour = "#1f78b4", fill = "#1f78b4", size = 5, stroke = 3)+ # circle point
  annotate("text",label = "Spring", x = 4.5, y = 2,size = 5.5,adj = 0)+
  annotate("point", x = 7, y = 2, shape = 21, colour = "#ff7f00", fill = "#ff7f00", size = 5, stroke = 3)+ # circle point
  annotate("text",label = "Fall", x = 7.5, y = 2,size = 5.5,adj = 0)+
  
  # Depth
  annotate("text",label = "Depth:", x = 1, y = 1,size = 7, adj = 0)+ # Title 
  annotate("point", x = 4, y = 1, shape = 19, colour = "black", size = 5, stroke = 3)+ # circle point
  annotate("text",label = "Shallow", x = 4.5, y = 1,size = 5.5,adj = 0)+
  annotate("point", x = 7, y = 1, shape = 17, colour = "black", size = 5, stroke = 3)+ # circle point
  annotate("text",label = "Deep", x = 7.5, y = 1,size = 5.5,adj = 0)+
  
  xlim(1,9)+
  # ylim(0.2,3.2)+
  ylim(0.4,2.6)+
  theme(panel.background = element_rect(fill="white"),
        panel.border = element_rect(colour = "black", fill=NA, size=2),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none",
        plot.margin = unit(c(0.2,6,0.2,6),"cm"))
  # theme(panel.background = element_rect(fill="white"),
  #       panel.border = element_rect(colour = "black", fill=NA, size=2),
  #       axis.text = element_blank(),
  #       axis.title = element_blank(),
  #       axis.ticks = element_blank(), 
  #       legend.position = "none",
  #       plot.margin = unit(c(0,10,0,10),"cm"))

# shannon_legend

shannon_paired_combined_legend <- ggpubr::ggarrange(shannon_paired_combined, shannon_legend, ncol = 1,
                  heights = c(10,1.5))




ggsave(shannon_paired_combined_legend, height = 10, width = 10, path = final_fig_dir, filename = "shannon_paired_combined_legend_v2.pdf", device = "pdf")








