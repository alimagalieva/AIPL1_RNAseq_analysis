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
library(splitstackshape)

df <- read_tsv("../kallisto_res/abundance.tsv")

annot <- read.delim("../Homo_sapiens.GRCh38.107.gff3", header=F, comment.char="#") 
annotation <- cSplit(annot, 'V9', ';')
colnames(annotation)[1] ="chromosome"

