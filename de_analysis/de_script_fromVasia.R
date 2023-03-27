rm(list=ls(all=TRUE))

## Setups 
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(DESeq2)
library(pheatmap)
#install.packages(<path to downloaded file with DEGreport tar.gz>, repos = NULL)
library(DEGreport)
#BiocManager::install("vsn")
library(vsn)
library(plyr) 
library(dplyr)

## reading file:
df <- read_tsv("data/fromVasia/aipl6_gfp.tsv") 
sig_res <- filter(df, baseMean>5 & padj < 0.001 & log2FoldChange > 1)
