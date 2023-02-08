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

## reaading file:
df_CC_AIPLopt <- read.csv("data/de.Sirius_CC_AIPLopt.csv", sep = ";") 
df_CC_AIPLwt <- read.csv("data/de.Sirius_CC_AIPLwt.csv", sep = ";") 
df_GFP_AIPLopt <- read.csv("data/de.Sirius_GFP_AIPLopt.csv", sep = ";") 
df_GFP_AIPLwt <- read.csv("data/de.Sirius_GFP_AILPwt.csv", sep = ";") 
df_AIPLwt_AIPLopt <- read.csv("data/de.Sirius_AIPLwt_AIPLopt.csv", sep = ";") 

