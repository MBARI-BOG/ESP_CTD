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

# Load data

## COI
asv_table_path_COI <- here("data", "COI", "ESP_CTD_COI_asv_table.csv")
tax_table_path_COI <- here("data", "COI", "ESP_CTD_COI_tax_table.csv")
metadata_path_COI <- here("data", "COI", "ESP_CTD_COI_metadata.csv")

ESP_CTD_COI_phyloseq <- merge_phyloseq(otu_table(read.csv(asv_table_path_COI, row.names = 1),taxa_are_rows = TRUE),
                                       tax_table(as.matrix(read.csv(tax_table_path_COI, row.names = 1))),
                                       sample_data(read.csv(metadata_path_COI, row.names = 1)))

## 18S
asv_table_path_18S <- here("data", "18S", "ESP_CTD_18S_asv_table.csv")
tax_table_path_18S <- here("data", "18S", "ESP_CTD_18S_tax_table.csv")
metadata_path_18S <- here("data", "18S", "ESP_CTD_18S_metadata.csv")

ESP_CTD_18S_phyloseq <- merge_phyloseq(otu_table(read.csv(asv_table_path_18S, row.names = 1),taxa_are_rows = TRUE),
                                       tax_table(as.matrix(read.csv(tax_table_path_18S, row.names = 1))),
                                       sample_data(read.csv(metadata_path_18S, row.names = 1)))

## 12S
asv_table_path_12S <- here("data", "12S", "ESP_CTD_12S_asv_table.csv")
tax_table_path_12S <- here("data", "12S", "ESP_CTD_12S_tax_table.csv")
metadata_path_12S <- here("data", "12S", "ESP_CTD_12S_metadata.csv")

ESP_CTD_12S_phyloseq <- merge_phyloseq(otu_table(read.csv(asv_table_path_12S, row.names = 1),taxa_are_rows = TRUE),
                                       tax_table(as.matrix(read.csv(tax_table_path_12S, row.names = 1))),
                                       sample_data(read.csv(metadata_path_12S, row.names = 1)))


