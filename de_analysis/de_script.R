#install.packages(<path to downloaded file with DEGreport tar.gz>, repos = NULL)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.16")
#BiocManager::install("vsn")
#BiocManager::install("tximport")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("AnnotationDbi")
BiocManager::install("enrichplot")
BiocManager::install("DESeq2")
BiocManager::install("ComplexHeatmap")
packages <- c("readr", "ggplot2", "dplyr", "magrittr")
install.packages(packages, dependencies = TRUE)

source("https://bioconductor.org/biocLite.R")
biocLite("SummarizedExperiment", dependencies = TRUE)
biocLite("DESeq2", dependencies = TRUE)
biocLite("airway", dependencies = TRUE)

################################################################################
########### Gene-level differential expression analysis using DESeq2 ###########
################################################################################

rm(list=ls(all=TRUE))

## Setups 
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(DESeq2)
library(pheatmap)
library(plyr) 
library(dplyr)
library(stringr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(DOSE)
library(enrichplot)
library(ComplexHeatmap)
library(VennDiagram)
library(RColorBrewer)
#library(DEGreport)
#library(vsn)
#library(tximeta)

## reading files:
dfFC <- read.delim("../countsFromFeatureCounts.out", sep = "\t", header = T, skip = 1)
df <- read.table("data/counts.txt", sep = "\t", header = F)
df <- df[1:(nrow(df)-5),]
df <- merge(x = df, y = data.frame(dfFC$Geneid, dfFC$Length), by = 1, all.x = TRUE)
colnames(df) <- c('geneID', 'name', 'AIPL1_6_2', 'AIPL1_6_3', 'AIPL_3_1', 
                  'AIPL_3_3', 'CC1', 'GFP_1', 'GFP_3', 'AIPL1_6_1', 
                  'AIPL_3_2', 'CC2', 'CC3', 'GFP_2', 'length')
df[df[1]=='optAIPL1g', 'length'] <- '1158'
df[df[1]=='wtAIPL1g', 'length'] <- '1158'
df[df[1]=='GFPg', 'length'] <- '741'
df[df[1]=='GFP_CBg', 'length'] <- '741'
df[df[1]=='optAIPL1g', 'name'] <- 'optAIPL1'
df[df[1]=='wtAIPL1g', 'name'] <- 'wtAIPL1'
df[df[1]=='GFPg', 'name'] <- 'GFP'
df[df[1]=='GFP_CBg', 'name'] <- 'GFP_CB'

## Trying to reorder columns by replicates alphabetically:
colnS <- sort(c(colnames(df[3:14])))
colnT <- c(colnames(df[1:2]), colnames(df[15]))
coln = c(colnT, colnS)
df <- df[, coln]

#row.names(data) <- data[, 'name'] ### can't do because some of them are non-unique

rownames(df) <- df$geneID
data = df[, (names(df) %in% colnS)]
meta <- data.frame(colnames(data))
meta['sampletype'] <- c('AAV9_wtAIPL1', 'AAV9_wtAIPL1', 'AAV9_wtAIPL1', 
                        'AAV9_optAIPL1', 'AAV9_optAIPL1', 'AAV9_optAIPL1', 
                        'ctrl', 'ctrl', 'ctrl', 
                        'AAV9_GFP', 'AAV9_GFP', 'AAV9_GFP')
row.names(meta) <- meta[,1]
meta <- subset(meta, select = c('sampletype'))

all(colnames(data) %in% rownames(meta))
all(colnames(data) == rownames(meta))


################################# PRE-FILTERING ################################

## creating DESeq DataSet
#dds <- DESeqDataSetFromMatrix(countData = data, colData = meta, design = ~ sampletype)
## Run DESeq2 differential expression analysis
#dds <- DESeq(dds)

##  **Optional step** - Output normalized counts to save as a file to access outside RStudio
#normalized_counts <- counts(dds, normalized=TRUE) 
#write.table(normalized_counts, file="data/normalized_counts.txt", sep="\t", quote=F, col.names=NA) 

## Total number of raw counts per sample
#colSums(counts(dds))

## Total number of normalized counts per sample
#colSums(counts(dds, normalized=T))

## log2fold change calculation and MAplot: - what for?
#plotMA(dds, main = "dds")
#res <- results(dds)
#plotMA(res, alpha = 0.05, main = "res", ylim=c(-7,7))

################################################################################
###################################### QC ######################################
################################################################################ 

## Transform counts for data visualization
#rld <- rlog(dds)
#pca <- plotPCA(rld, intgroup="sampletype") + geom_text_repel(aes(label=name), size = 3) + geom_point(size = 1)
#ggsave(pca, file='plots/PCA.png', width = 7, height = 4, dpi = 600)
## Extract the rlog matrix from the object
#rld_mat <- assay(rld)
## Compute pairwise correlation values
#rld_cor <- cor(rld_mat)
## Plot heatmap
#heat.colors <- brewer.pal(6, "Blues")
#heatmap <- pheatmap(rld_cor, color = heat.colors, fontsize = 10, fontsize_row = 10, height=20)
#ggsave(heatmap, file='plots/heatmap.png', width = 7, height = 4, dpi = 600)

#remove CC1 from analysis, cause it doesn't cluster well and re-do dds object
data = data[, !(names(data) %in% c('CC1'))] 
meta <- data.frame(colnames(data))
meta['sampletype'] <- c('AAV9_wtAIPL1', 'AAV9_wtAIPL1', 'AAV9_wtAIPL1', 
                        'AAV9_optAIPL1', 'AAV9_optAIPL1', 'AAV9_optAIPL1', 
                        'ctrl', 'ctrl',
                        'AAV9_GFP', 'AAV9_GFP', 'AAV9_GFP')
row.names(meta) <- meta[,1]
meta <- subset(meta, select = c('sampletype'))

all(colnames(data) %in% rownames(meta))
all(colnames(data) == rownames(meta))

## creating DESeq DataSet
dds <- DESeqDataSetFromMatrix(countData = data, colData = meta, design = ~ sampletype)
## Run DESeq2 differential expression analysis
dds <- DESeq(dds)

df_genes <- df[, 2:3]

#### QC
#rld <- rlog(dds)
#pca <- plotPCA(rld, intgroup="sampletype") + geom_text_repel(aes(label=name), size = 3) + geom_point(size = 1)
#ggsave(pca, file='plots/PCA1.png', width = 7, height = 4, dpi = 600)
## Extract the rlog matrix from the object
#rld_mat <- assay(rld)
## Compute pairwise correlation values
#rld_cor <- cor(rld_mat)
## Plot heatmap
#heat.colors <- brewer.pal(6, "Blues")
#heatmap <- pheatmap(rld_cor, color = heat.colors, fontsize = 10, fontsize_row = 10, height=20)
#ggsave(heatmap, file='plots/heatmap1.png', width = 7, height = 4, dpi = 600)

################################################################################
############################## ///////////////// ###############################
############################## /// tentative /// ###############################
############################## ///////////////// ###############################
################################################################################

'''
res_unshrunken_wt_gfp <- results(dds, contrast = contrast_wt_gfp,  alpha = pval, lfcThreshold = lfc)
res_shrunken_wt_gfp <- lfcShrink(dds, contrast = contrast_wt_gfp, res=res_unshrunken_wt_gfp, type='norm')
res_df_wt_gfp <- data.frame(res_shrunken_wt_gfp)
sig_res_wt_gfp <- filter(res_df_wt_gfp, padj < pval & abs(log2FoldChange) > lfc)

hist(sig_res_wt_gfp$padj) ## looks norm

res_df_wt_gfp_up_TF <- res_df_wt_gfp %>% 
  mutate(dif_expr = padj < pval & log2FoldChange >= lfc)

de_up_wt_gfp <- res_df_wt_gfp_up_TF[which(res_df_wt_gfp_up_TF$dif_expr == 'TRUE'), ]

## Merge to see geneNames:
de_up_wt_gfp <- merge(x = de_up_wt_gfp, y = df_genes, by = 0, all.x = TRUE)
de_up_wt_gfp <- de_up_wt_gfp %>% arrange(desc(log2FoldChange))
'''

################################################################################
################################### Shrinkage ##################################
################################################################################

contrast_wt_gfp <- c("sampletype", "AAV9_wtAIPL1", "AAV9_GFP")
contrast_opt_gfp <- c("sampletype", "AAV9_optAIPL1", "AAV9_GFP")
contrast_opt_wt <- c("sampletype", "AAV9_optAIPL1", "AAV9_wtAIPL1")
contrast_gfp_ctrl <- c("sampletype", "AAV9_GFP", "ctrl")
contrast_wt_ctrl <- c("sampletype", "AAV9_wtAIPL1", "ctrl")
contrast_opt_ctrl <- c("sampletype", "AAV9_optAIPL1", "ctrl")

lfc <- 2
pval <- 0.001
samp <- 'opt_ctrl'

shrink <- function(samp, pval, lfc) {
  unshrunken <- results(dds, contrast = get(paste('contrast_', samp, sep='')))
  shrunken <- lfcShrink(dds, contrast = get(paste('contrast_', samp, sep='')), 
                        res=unshrunken, type='normal')
  res_df <<- data.frame(shrunken)
  res_df$ens <- gsub("\\..*","", rownames(res_df))
  res_df$symbol <- mapIds(org.Hs.eg.db, keys = res_df$ens, keytype = 'ENSEMBL', column = 'SYMBOL')
  res_df[rownames(res_df[is.na(res_df$symbol), ]), 'symbol'] <- rownames(res_df[is.na(res_df$symbol),])
  sigs <<- filter(res_df, padj < pval & abs(log2FoldChange) > lfc)
  res_df[rownames(res_df[is.na(res_df$symbol), ]), 'symbol'] <<- rownames(res_df[is.na(res_df$symbol),])
}

volcano <- function(samp, pval, lfc) {
  sigs_top5 <- subset(sigs %>% arrange(log2FoldChange) %>% head(n=5))
  sigs_low5 <- subset(sigs %>% arrange(desc(log2FoldChange)) %>% head(n=5))
  sigs_10 <- rbind(sigs_top5, sigs_low5)
  res_df <- merge(res_df, sigs_10, by = 'row.names', all.x = T)
  row.names(res_df) <- res_df$Row.names
  res_df <- res_df[c('baseMean.x', 'log2FoldChange.x', 'lfcSE.x', 'stat.x', 'pvalue.x', 'padj.x', 'symbol.y')]
  colnames(res_df) <- c('baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'padj', 'name10')
  
  res_df$super <- res_df$log2FoldChange >  lfc  & res_df$padj < pval
  res_df$sub   <- res_df$log2FoldChange < -lfc  & res_df$padj < pval
  res_df$threshold <- as.factor(abs(res_df$log2FoldChange) > lfc & res_df$padj < pval)
  
  volc <- ggplot(data=res_df, aes(x=log2FoldChange, y=-log10(padj))) +
    geom_point(data=res_df, size=1, colour="gray") +
    geom_point(data=res_df[res_df$super==TRUE, ], size=1.5, colour="#CC0000") +
    geom_point(data=res_df[res_df$sub  ==TRUE, ], size=1.5, colour="#000099") +
    geom_text_repel(aes(label = name10)) +
    xlab("log2 fold change") +
    ylab("-log10 p-value adjusted") +
    ggtitle(paste('Differential expression ', samp, '\nl2fc=', lfc, ', padj=', pval)) +
    scale_x_continuous() +
    scale_y_continuous() +
    theme_bw() +
    theme(axis.title.y = element_text(size=16),
          axis.title.x = element_text(size=16, colour="black"),
          axis.text = element_text(size=12),
          legend.title =element_blank() ,
          legend.text = element_text(size = 12)) +
    theme(plot.title = element_text(face="bold", size = 17, hjust = 0.5)) +
    geom_hline(yintercept = -log10(pval)) + 
    geom_vline(xintercept = c(lfc, -lfc))
  
  ggsave(volc, file=paste('plots/volcano_', samp, '_lfc=', lfc,'_padj=', pval, '.png', sep=""), 
         width = 7, height = 7, dpi = 600)
}

heat <- function(samp, pval, lfc) {
  matr <- counts(dds, normalized=T)[rownames(sigs),]
  matr.z <- t(apply(matr, 1, scale))
  colnames(matr.z) <- colnames(matr)
  fig <- Heatmap(matr.z, cluster_rows = T, cluster_columns = T, column_labels = colnames(matr.z), 
                 name = 'Z.score', row_labels = sigs[rownames(matr.z), 'symbol'], 
                 column_title = paste('Heatmap of ', as.character(samp),', lfc=',as.character(lfc),', padj=',as.character(pval), sep=''))
  png(paste('plots/heat_', samp, '_lfc=', lfc, '_padj=', pval, '.png', sep=''), 
      res = 250, width = 1500, height = 2000)
  print(fig)
  dev.off()
}

go_bar <- function(samp, pval, lfc) {
  for (dir in c('up', 'down')) {
    genes_to_test <- get(paste('sigs_', dir, '_', samp, sep=''))[,1]
    GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont="BP")
    as.data.frame(GO_results)
    if (nrow(GO_results) != 0) {
      fig <- barplot(GO_results, showCategory = 20, main = paste('GO of ', dir, '-regulated ', samp, ', \npadj=', pval, 'l2fc=', lfc)) + 
        ggtitle(paste('GO of ', samp, dir, '-regulated', ', padj=', pval, ', l2fc=', lfc, sep='')) +
        theme(plot.title = element_text(face = 'bold', hjust = 0.5), 
              axis.title.x = element_text(size = 16), 
              axis.text.x = element_text(size = 12), 
              axis.text.y = element_text(size = 10))
      ggsave(fig, 
             file=paste('plots/GO_', samp, '_', dir, '_lfc=', lfc, '_padj=', pval, '.png', sep=''), 
             width = 14, height = 7, dpi = 600)
    }
  }
}


for (samp in c('wt_ctrl','opt_ctrl', 'opt_wt', 'wt_gfp', 'opt_gfp')) {
  shrink(samp, pval, lfc)
  
  assign(paste('sigs_', samp, sep=''), sigs[,c('symbol', 'baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'padj')])
#  write.table(get(paste('sigs_', samp, sep=''))[c('symbol', 'padj', 'log2FoldChange')], 
#              file=paste('results/sigs_', samp, '_l2fc=', lfc, '_padj=', pval, '.txt', sep=''), sep="\t", quote=F, col.names=NA)
  assign(paste('sigs_down_', samp, sep=''), filter(get(paste('sigs_', samp, sep='')), padj < pval & log2FoldChange < lfc))
#  write.table(get(paste('sigs_down_', samp, sep=''))[c('symbol', 'padj', 'log2FoldChange')], 
#              file=paste('results/sigs_', samp, '_down_l2fc=', lfc, '_padj=', pval, '.txt', sep=''), sep="\t", quote=F, col.names=NA)
  assign(paste('sigs_up_', samp, sep=''), filter(get(paste('sigs_', samp, sep='')), padj < pval & log2FoldChange > lfc))
#  write.table(get(paste('sigs_up_', samp, sep=''))[c('symbol', 'padj', 'log2FoldChange')], 
#              file=paste('results/sigs_', samp, '_up_l2fc=', lfc, '_padj=', pval, '.txt', sep=''), sep="\t", quote=F, col.names=NA)
  nrow(get(paste('sigs_', samp, sep=''))) == nrow(get(paste('sigs_up_', samp, sep='')))+nrow(get(paste('sigs_down_', samp, sep='')))
  
  #volcano(samp, pval, lfc)
  
  #heat(samp, pval, lfc)
  
  #go_bar(samp, pval, lfc)
} ### never delete middle strings!


i <- 1
for (samp in list(c('wt_ctrl', 'wt_gfp'), c('opt_ctrl', 'opt_gfp'), c('wt_gfp', 'opt_gfp'), c('wt_ctrl', 'opt_ctrl')) ) {
  #print(samp[2])
  for (dir in c('up', 'down')) {
    venn.diagram(
      x = list(row.names(get(paste('sigs_', dir, '_', samp[1], sep=''))), row.names(get(paste('sigs_', dir, '_', samp[2], sep='')))),
      category.names = c(samp[1], samp[2]),
      filename = paste('plots/venn_', as.character(i), '_', dir, '.png', sep=''),
      output=TRUE, 
      fill = c(2,4), 
      main = paste('Venn diagram of ', samp[1], ' and ', samp[2], ', ', dir, '-regulated', sep='')
    )
  }
  i <- i + 1
}

### mega-venn
for (dir in c('up', 'down')) {
  venn.diagram(
    x = list(row.names(get(paste('sigs_', dir, '_wt_ctrl', sep=''))), row.names(get(paste('sigs_', dir, '_wt_gfp', sep=''))), row.names(get(paste('sigs_', dir, '_opt_ctrl', sep=''))), row.names(get(paste('sigs_', dir, '_opt_gfp', sep='')))),
    category.names = c("wt_ctrl", "wt_gfp", 'opt_ctrl', 'opt_gfp'),
    filename = paste('plots/venn_', dir, '.png', sep=''),
    output=TRUE, 
    fill = c(1,2,3,4), 
    main = paste('Venn diagram, ', dir, '-regulated', sep='')
  )
}

intersect1_wt_down <- intersect(sigs_down_wt_gfp$symbol, sigs_down_wt_ctrl$symbol)
intersect1_wt_up <- intersect(sigs_up_wt_gfp$symbol, sigs_up_wt_ctrl$symbol)

intersect2_opt_down <- intersect(sigs_down_opt_gfp$symbol, sigs_down_opt_ctrl$symbol)
intersect2_opt_up <- intersect(sigs_up_opt_gfp$symbol, sigs_up_opt_ctrl$symbol)

intersect3_gfp_down <- intersect(sigs_down_opt_gfp$symbol, sigs_down_wt_gfp$symbol)
intersect3_gfp_up <- intersect(sigs_up_opt_gfp$symbol, sigs_up_wt_gfp$symbol)

intersect4_ctrl_down <- intersect(sigs_down_opt_ctrl$symbol, sigs_down_wt_ctrl$symbol)
intersect4_ctrl_up <- intersect(sigs_up_opt_ctrl$symbol, sigs_up_wt_ctrl$symbol)

intersect_all_down <- intersect(intersect3_gfp_down, intersect4_ctrl_down)
intersect_all_up <- intersect(intersect3_gfp_up, intersect4_ctrl_up)



deg <- data.frame(matrix(ncol = 4, nrow = 0))
x <- c("sample", "padj", "l2fc", 'N of DEG')
colnames(deg) <- x
N_of_DEG <- function(samp, pval, lfc) {
  assign(paste('sigs_', samp, sep=''), filter(get(paste('res_df_', samp, sep='')), padj < pval & abs(log2FoldChange) > lfc))
  assign(paste('sigs_', samp, sep=''), merge(get(paste('sigs_', samp, sep='')), y = df_genes, by = 0, all.x = T))
  genes_to_test <- get(paste('sigs_', samp, sep=''))[,7]
  new <- c(samp,pval,lfc,length(genes_to_test))
  deg[nrow(deg) + 1, ] <<- new 
}

for (pval in c(0.05, 0.01, 0.005, 0.001)) {
  for (lfc in c(0.5, 1, 2)) {
    for (samp in c('wt_ctrl','opt_ctrl', 'opt_wt', 'wt_gfp', 'opt_gfp')) {
      N_of_DEG(samp, pval, lfc)
    }
  }
}

pval <- 0.005 
lfc <- 1
samp <- 'opt_ctrl'
## trying to remove 'virus'
for (dir in c('up', 'down')) {
  genes_to_test <- get(paste('sigs_', dir, '_', samp, sep=''))[,1]
  GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont="BP")
  GO_results <- as.data.frame(GO_results)
  
}



################################################################################
########################### CHOOSING GENES TO TEST #############################
################################################################################

hist(sig_res_gfp_ctrl$log2FoldChange, breaks = 100)


sig_res_wt_gfp <- filter(res_df_wt_gfp, padj < pval & abs(log2FoldChange) > lfc)
sig_res_opt_gfp <- filter(res_df_opt_gfp, padj < pval & abs(log2FoldChange) > lfc)
sig_res_opt_wt <- filter(res_df_opt_wt, padj < pval & abs(log2FoldChange) > lfc)

sig_res_wt_gfp <- subset(sig_res_wt_gfp %>% arrange(desc(abs(log2FoldChange))) %>% head(n=20))
sig_res_opt_gfp <- subset(sig_res_opt_gfp %>% arrange(desc(abs(log2FoldChange))) %>% head(n=20))
sig_res_opt_wt <- subset(sig_res_opt_wt %>% arrange(desc(abs(log2FoldChange))) %>% head(n=20))

sig_res_wt_gfp <- merge(x = sig_res_wt_gfp, y = df_genes, by = 0, all.x = TRUE)
sig_res_opt_gfp <- merge(x = sig_res_opt_gfp, y = df_genes, by = 0, all.x = TRUE)
sig_res_opt_wt <- merge(x = sig_res_opt_wt, y = df_genes, by = 0, all.x = TRUE)

write.table(sig_res_wt_gfp[c('Row.names','name', 'length', 'padj', 'log2FoldChange')], file=str_replace_all("results/top20_wt_gfp=lfcTr_padj=padjTr.txt", c(lfcTr=as.character(lfc), padjTr=as.character(pval))), sep="\t", quote=F, col.names=NA)
write.table(sig_res_opt_gfp[c('Row.names','name', 'length', 'padj', 'log2FoldChange')], file=str_replace_all("results/top20_opt_gfp=lfcTr_padj=padjTr.txt", c(lfcTr=as.character(lfc), padjTr=as.character(pval))), sep="\t", quote=F, col.names=NA)
write.table(sig_res_opt_wt[c('Row.names','name', 'length', 'padj', 'log2FoldChange')], file=str_replace_all("results/top20_opt_wt=lfcTr_padj=padjTr.txt", c(lfcTr=as.character(lfc), padjTr=as.character(pval))), sep="\t", quote=F, col.names=NA)

genes_to_test <- get(paste('sig_res_opt_gfp', sep=''))[,8]
GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", pvalueCutoff = 0.1, ont="BP")
fig <- barplot(GO_results, showCategory = 20, main = paste("GO top20 of opt_gfp",', \npadj=',as.character(pval),'l2fc=',as.character(lfc))) + 
  ggtitle(paste('GO of top20 opt_gfp, padj=',as.character(pval),', l2fc=',as.character(lfc), sep='')) + 
  theme(plot.title = element_text(face = 'bold', hjust = 0.5),
        axis.title.x = element_text(size = 16), 
        axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 10))
ggsave(fig, 
         file=paste('plots/GO_top20_opt_gfp',as.character(samp),'_lfc=',as.character(lfc),'_padj=',as.character(pval),'.png', sep=''), 
         width = 14, height = 7, dpi = 600)
}


################################################################################
###############  EXTRACTING GENES PRESENTED IN AILP1 MINUS GFP #################
################################################################################

de_down_opt_minus_gfp <- filter(de_down_opt, !Row.names %in% de_down_gfp_ctrl$Row.names)
de_down_wt_minus_gfp <- filter(de_down_wt, !Row.names %in% de_down_gfp_ctrl$Row.names)
de_up_opt_minus_gfp <- filter(de_up_opt, !Row.names %in% de_up_gfp_ctrl$Row.names)
de_up_wt_minus_gfp <- filter(de_up_wt, !Row.names %in% de_up_gfp_ctrl$Row.names)

#write.table(de_down_opt_minus_gfp[c('Row.names','Symbol', 'Gene.name', 'padj', 'log2FoldChange')], file="results/de_down_opt_minus_gfp.txt", sep="\t", quote=F, col.names=NA) 
#write.table(de_down_wt_minus_gfp[c('Row.names','Symbol', 'Gene.name', 'padj', 'log2FoldChange')], file="results/de_down_wt_minus_gfp.txt", sep="\t", quote=F, col.names=NA) 
#write.table(de_up_opt_minus_gfp[c('Row.names','Symbol', 'Gene.name', 'padj', 'log2FoldChange')], file="results/de_up_opt_minus_gfp.txt", sep="\t", quote=F, col.names=NA) 
#write.table(de_up_wt_minus_gfp[c('Row.names','Symbol', 'Gene.name', 'padj', 'log2FoldChange')], file="results/de_up_wt_minus_gfp.txt", sep="\t", quote=F, col.names=NA) 

################################################################################
#################### Dot-plot expression of a single gene ######################
################################################################################
### ENSG00000129221.16 - AIPL1 gene 
### now makes no sense
# Save plotcounts to a data frame object
d1 <- plotCounts(dds, gene="ENSG00000129221.16", intgroup="sampletype", returnData=TRUE)
# Plotting AIPL1 normalized counts, using the samplenames (rownames of d as labels)
ggplot(d1, aes(x = sampletype, y = count, color = sampletype)) + 
  geom_point(position=position_jitter(w = 0.1,h = 0)) +
  geom_text_repel(aes(label = rownames(d1))) + 
  theme_bw() +
  ggtitle("AIPL1 gene expression",) +
  theme(plot.title = element_text(hjust = 0.5, face='bold'))

### barplot:
## calculating sd:
d2 <- ddply(d1, .(sampletype), summarize, mean = round(mean(count), 2), 
            sd=round(sd(count), 2)) 

ggplot(d2, aes(x = sampletype, y = mean, color = sampletype)) + 
  geom_bar(position=position_dodge(), stat="identity") + 
  ylab("normalized Counts") + xlab("sample") + ggtitle("AIPL1 gene expression")+ 
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, 
                position=position_dodge(.9))

################################################################################
################################################################################
################################################################################

### !!! Difficult to opperate due to (and, mb, fuckuped numeration of my rows)
### Using ggplot2 to plot multiple genes (e.g. top 20):
## Order results by padj values
top20_sig_genes <- sig_res_opt %>% 
  arrange(padj) %>% 	#Arrange rows by padj values
  head(n=20) 		#Extract the first 20 genes
#top20_sig_genes$gene <- rownames(top20_sig_genes)
## Then, we can extract the normalized count values for these top 20 genes:
## normalized counts for top 20 significant genes
top20_sig_norm <- normalized_counts %>% filter(gene %in% top20_sig_genes)
# Gathering the columns to have normalized counts to a single column
gathered_top20_sigOE <- top20_sigOE_norm %>%
  gather(colnames(top20_sigOE_norm)[2:9], key = "samplename", value = "normalized_counts")
## check the column header in the "gathered" data frame
View(gathered_top20_sigOE)
gathered_top20_sigOE <- inner_join(mov10_meta, gathered_top20_sigOE)
## plot using ggplot2
ggplot(gathered_top20_sigOE) +
  geom_point(aes(x = gene, y = normalized_counts, color = sampletype)) +
  scale_y_log10() +
  xlab("Genes") +
  ylab("log10 Normalized Counts") +
  ggtitle("Top 20 Significant DE Genes") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5))
