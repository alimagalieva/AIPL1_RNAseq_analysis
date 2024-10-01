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
library(AnnotationDbi)
library(org.Hs.eg.db)
library(DOSE)
library(enrichplot)
library(ComplexHeatmap)
library(VennDiagram)
library(RColorBrewer)
library(igraph)
library(STRINGdb)
library(tximport)

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
df[df[1]=='optAIPL1g', 'name'] <- 'AIPL1co'
df[df[1]=='wtAIPL1g', 'name'] <- 'AIPL1wt'
df[df[1]=='GFPg', 'name'] <- 'GFP'
df[df[1]=='GFP_CBg', 'name'] <- 'GFP_CB'

## Trying to reorder columns by replicates alphabetically:
colnS <- sort(c(colnames(df[3:14])))
colnT <- c(colnames(df[1:2]), colnames(df[15]))
coln = c(colnT, colnS)
df <- df[, coln]

df[df[2]=='AIPL1co', 'geneID'] <- 'AIPL1co'
df[df[2]=='AIPL1wt', 'geneID'] <- 'AIPL1wt'
df[df[2]=='GFP', 'geneID'] <- 'GFP'
df[df[2]=='GFP_CB', 'geneID'] <- 'GFP_CB'
rownames(df) <- df$geneID
data = df[, (names(df) %in% colnS)]
meta <- data.frame(row.names=colnames(data), sampletype= c('AAV9_AIPL1wt', 'AAV9_AIPL1wt', 'AAV9_AIPL1wt', 
                        'AAV9_AIPL1co', 'AAV9_AIPL1co', 'AAV9_AIPL1co', 
                        'non_transd_ctrl', 'non_transd_ctrl', 'non_transd_ctrl',
                        'AAV9_GFP', 'AAV9_GFP', 'AAV9_GFP'), repl=rep(c(1,2,3), 4))
all(colnames(data) %in% rownames(meta))
all(colnames(data) == rownames(meta))

hist(rowSums(data), breaks = 1000)


################################################################################
###################################### QC ######################################
################################################################################ 

#remove CC1 from analysis, cause it doesn't cluster well and re-do dds object
data = data[, !(names(data) %in% c('CC1'))] 
meta <- data.frame(colnames(data))
meta['sampletype'] <- c('AAV9_AIPL1wt', 'AAV9_AIPL1wt', 'AAV9_AIPL1wt', 
                        'AAV9_AIPL1co', 'AAV9_AIPL1co', 'AAV9_AIPL1co', 
                        'non_transd_ctrl', 'non_transd_ctrl', 
                        'AAV9_GFP', 'AAV9_GFP', 'AAV9_GFP')
row.names(meta) <- meta[,1]
meta <- subset(meta, select = c('sampletype'))
meta$repl <- c(1,2,3,1,2,3,2,3,1,2,3)

all(colnames(data) %in% rownames(meta))
all(colnames(data) == rownames(meta))

summary(data)

## creating DESeq DataSet
dds <- DESeqDataSetFromMatrix(countData = data, colData = meta, design = ~ sampletype)
## Run DESeq2 differential expression analysis
dds <- DESeq(dds)

sum(apply(counts(dds, normalized = TRUE), 1, function(row) all(row ==0)))


################################################################################
################################### Shrinkage ##################################
################################################################################

contr_wt_gfp <- c("sampletype", "AAV9_AIPL1wt", "AAV9_GFP")
contr_opt_gfp <- c("sampletype", "AAV9_AIPL1co", "AAV9_GFP")
contr_opt_wt <- c("sampletype", "AAV9_AIPL1co", "AAV9_AIPL1wt")
contr_gfp_ctrl <- c("sampletype", "AAV9_GFP", "non_transd_ctrl")
contr_wt_ctrl <- c("sampletype", "AAV9_AIPL1wt", "non_transd_ctrl")
contr_opt_ctrl <- c("sampletype", "AAV9_AIPL1co", "non_transd_ctrl")

lfc <- 1.58
pval <- 0.005
samp <- 'opt_wt'

shrink <- function(samp, pval, lfc) {
  unshrunken <- results(dds, contrast = get(paste('contr_', samp, sep='')))
  shrunken <- lfcShrink(dds, contrast = get(paste('contr_', samp, sep='')), 
                        res=unshrunken, type='normal')
  res_df <<- data.frame(shrunken)
  #print(samp, res_df)
}

filtr <- function(samp, pval, lfc){
  res_df$ens <- gsub("\\..*","", rownames(res_df))
  res_df$symbol <- mapIds(org.Hs.eg.db, keys = res_df$ens, keytype = 'ENSEMBL', column = 'SYMBOL')
  res_df$map <- mapIds(org.Hs.eg.db, keys = res_df$ens, keytype = 'ENSEMBL', column = 'MAP')
  res_df[rownames(res_df[is.na(res_df$symbol), ]), 'symbol'] <- rownames(res_df[is.na(res_df$symbol),])
  sigs <<- filter(res_df, padj < pval & abs(log2FoldChange) > lfc)
  res_df[rownames(res_df[is.na(res_df$symbol), ]), 'symbol'] <<- rownames(res_df[is.na(res_df$symbol),])
  res_df <<- res_df
}

volcano_new <- function(samp, pval, lfc) {
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
    geom_point(data=res_df[!is.na(res_df$name10), ], size=2, shape=1) +
    geom_text_repel(aes(label = name10), colour='#00000000') +
    xlab("log2 fold change") +
    ylab("-log10 p-value adjusted") +
    ggtitle(paste('Differential expression of ', samp, '\nl2fc = ', lfc, ', padj = ', pval, sep='')) +
    scale_x_continuous() +
    scale_y_continuous() +
    theme_bw() +
    theme(axis.title.y = element_text(size=16),
          axis.title.x = element_text(size=16, colour="black"),
          axis.text = element_text(size=12),
          legend.title =element_blank() ,
          legend.text = element_text(size = 12), 
          plot.title = element_text(face="bold", size = 17, hjust = 0.5)) +
    geom_hline(yintercept = -log10(pval), linetype="dotted") + 
    geom_vline(xintercept = c(lfc, -lfc), linetype="dotted")
  
  #ggsave(volc, file=paste('plots/volcano_', samp, '_lfc=', lfc,'_padj=', pval, '.png', sep=""), 
  #       width = 7, height = 7, dpi = 600)
  ggsave(volc, file=paste('plots/volcano_pl_', samp, '_lfc=', lfc,'_padj=', pval, '.svg', sep=""), width=10, height=10)
}

heat <- function(samp, pval, lfc) {
  genesToMap <- c(row.names(sigs))
  matr <- counts(dds, normalized=T)[genesToMap, ]
  matr.z <- t(apply(matr, 1, scale))
  matr.z <- na.omit(matr.z)
  colnames(matr.z) <- paste(meta[rownames(meta)==colnames(matr),'sampletype'], meta[rownames(meta)==colnames(matr),'repl'], sep='-')
  fig <- Heatmap(matr.z, cluster_rows = T, cluster_columns = T, 
                 column_labels = colnames(matr.z), column_names_rot = 45, column_names_gp = gpar(fontsize = 8),
                 name = 'Z.score', row_labels = sigs[as.factor(rownames(matr.z)), 'symbol'],
                 row_names_gp = gpar(fontsize = 3), 
                 column_title = paste('Heatmap of ', as.character(samp),' DEGs, l2fc = ',as.character(lfc),', padj = ',as.character(pval), sep=''))
  png(paste('plots/heat_', samp, '_lfc=', lfc, '_padj=', pval, '.png', sep=''), 
      res = 2500, width = 15000, height = 20000)
  print(fig)
  dev.off()
}

go_bar <- function(samp, pval, lfc) {
  for (dir in c('up', 'down')) {
    genes_to_test <- get(paste('sigs_', dir, '_', samp, sep=''))[,1]
    GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont="BP")
    as.data.frame(GO_results)
    if (nrow(GO_results) != 0) {
      fig <- barplot(GO_results, showCategory = 20) + 
        ggtitle(paste('GO of ', samp, ', ', dir, '-regulated', ', l2fc = ', lfc, ', padj = ', pval, sep='')) +
        theme(plot.title = element_text(face = 'bold', hjust = 0.5), 
              axis.title.x = element_text(size = 16), 
              axis.text.x = element_text(size = 12), 
              axis.text.y = element_text(size = 9))
      ggsave(fig, 
             file=paste('plots/GO_', samp, '_', dir, '_lfc=', lfc, '_padj=', pval, '.png', sep=''), 
             width = 14, height = 7, dpi = 600)
    }
  }
}

for (samp in c('wt_ctrl','opt_ctrl', 'opt_wt', 'wt_gfp', 'opt_gfp', 'gfp_ctrl')) {
  
  shrink(samp, pval, lfc)
  
  filtr(samp, pval, lfc)
  
  #assign(paste('res_', samp, sep=''), res_df[,c(-7)])
  
  assign(paste('sigs_', samp, sep=''), sigs[,c('symbol', 'map', 'baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'padj')])
  #write.table(get(paste('sigs_', samp, sep=''))[c('symbol', 'map', 'padj', 'log2FoldChange')], 
  #            file=paste('results/sigs_', samp, '_l2fc=', lfc, '_padj=', pval, '.txt', sep=''), sep="\t", quote=F, col.names=NA)
  assign(paste('sigs_down_', samp, sep=''), filter(get(paste('sigs_', samp, sep='')), padj < pval & log2FoldChange < lfc))
  #write.table(get(paste('sigs_down_', samp, sep=''))[c('symbol', 'map', 'padj', 'log2FoldChange')], 
  #            file=paste('results/sigs_', samp, '_down_l2fc=', lfc, '_padj=', pval, '.txt', sep=''), sep="\t", quote=F, col.names=NA)
  assign(paste('sigs_up_', samp, sep=''), filter(get(paste('sigs_', samp, sep='')), padj < pval & log2FoldChange > lfc))
  #write.table(get(paste('sigs_up_', samp, sep=''))[c('symbol', 'map', 'padj', 'log2FoldChange')], 
  #            file=paste('results/sigs_', samp, '_up_l2fc=', lfc, '_padj=', pval, '.txt', sep=''), sep="\t", quote=F, col.names=NA)
  #print(nrow(get(paste('sigs_', samp, sep=''))) == nrow(get(paste('sigs_up_', samp, sep=''))) + nrow(get(paste('sigs_down_', samp, sep=''))))
  
  #volcano(samp, pval, lfc)
  
  #heat(samp, pval, lfc)
  
  #go_bar(samp, pval, lfc)
} ### never delete middle strings!


png(file='plots/plotCounts_AIPL1co.png', width=7, height=4, units="in", res=1200)
plotCounts(dds, gene='AIPL1co', intgroup='sampletype', 
                  normalized = TRUE, main = 'AIPL1co', xlab = '')
dev.off()
png(file='plots/plotCounts_AIPL1wt.png', width=7, height=4, units="in", res=1200)
plotCounts(dds, gene='ENSG00000129221.16', intgroup='sampletype', 
           normalized = TRUE, main = 'AIPL1wt', xlab = '')
dev.off()

png(file='plots/plotCounts_IL1RN.png', width=7, height=4, units="in", res=1200)
plotCounts(dds, gene='ENSG00000136689.20', intgroup='sampletype', 
           normalized = TRUE, main = 'IL1RN', xlab = '')
dev.off()

png(file='plots/plotCounts_H2AZ2.png', width=7, height=4, units="in", res=1200)
plotCounts(dds, gene='ENSG00000105968.19', intgroup='sampletype', 
           normalized = TRUE, main = 'H2AZ2', xlab = '')
dev.off()
png(file='plots/plotCounts_H2AZ1.png', width=7, height=4, units="in", res=1200)
plotCounts(dds, gene='ENSG00000164032.12', intgroup='sampletype', 
           normalized = TRUE, main = 'H2AZ1', xlab = '')
dev.off()

d1 <- plotCounts(dds, gene="ENSG00000129221.16", intgroup="sampletype", returnData=TRUE)
d2 <- ddply(d1, .(sampletype), summarize, mean = round(mean(count), 2), 
            sd=round(sd(count), 2)) 
d1 <- plotCounts(dds, gene="AIPL1co", intgroup="sampletype", returnData=TRUE)
d3 <- ddply(d1, .(sampletype), summarize, mean = round(mean(count), 2), 
            sd=round(sd(count), 2)) 
d <- rbind(d2, d3)
d[,'AIPL1_variant'] <- c(factor(rep(c('AIPL1wt', 'AIPL1co'), each=4)))
counts_AIPL1 <- ggplot(data=d, aes(x = sampletype, y = mean, fill = AIPL1_variant)) + 
  geom_bar(stat="identity", position=position_dodge(), colour='black', linewidth=0.1) + 
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, 
                position=position_dodge(.9)) +
  ylab("normalized Counts") + 
  xlab("") + 
  ggtitle("AIPL1 expression") + 
  theme(axis.text=element_text(colour="black"), axis.title.y=element_text(size=18), 
        axis.text.x=element_text(size=16), legend.text=element_text(size=16), 
        legend.title=element_text(size=18), 
        plot.title = element_text(hjust = 0.5, size=20)) +
  scale_fill_brewer(palette="Blues") 
ggsave(counts_AIPL1, file='plots/counts_AIPL1.png', width = 12, height = 7, dpi = 1200)


d1 <- plotCounts(dds, gene="ENSG00000105968.19", intgroup="sampletype", returnData=TRUE)
d2 <- plotCounts(dds, gene="ENSG00000164032.12", intgroup="sampletype", returnData=TRUE)
d <- rbind(d2, d1)
d[,'H2AZ'] <- c(factor(rep(c('H2AZ2', 'H2AZ1'), each=11)))
ggplot(data=d, aes(x = sampletype, y = count, fill = H2AZ)) + 
  geom_boxplot()

d1 <- plotCounts(dds, gene="AIPL1co", intgroup="sampletype", returnData=TRUE)
d2 <- plotCounts(dds, gene="ENSG00000129221.16", intgroup="sampletype", returnData=TRUE)
d <- rbind(d2, d1)
d[,'AIPL1_variant'] <- c(factor(rep(c('AIPL1co', 'AIPL1wt'), each=11)))
ggplot(data=d, aes(x = sampletype, y = count, fill=AIPL1_variant)) 
  

################################################################################
############################# compare with genes list ##########################
################################################################################

for (samp in c('wt_ctrl','opt_ctrl', 'opt_wt', 'wt_gfp', 'opt_gfp')) {
  res <- get(as.name(paste('res_', samp, sep='')))
  
  print(samp)
  print('up')
  for (g in c('MAF', 'ABCB5', 'GAS1', 'WHRN', 'UAP1L1', 'ANGPTL2', 'UCN')) {
    print(res[(res$symbol == g), c('symbol', 'padj', 'log2FoldChange')] )
  }
  print('down')
  for (g in c('RRM2')) {
    print(res[(res$symbol == g), c('symbol', 'padj', 'log2FoldChange')] )
  }
}

Andreys_genes <- c('MAF', 'ABCB5', 'GAS1', 'WHRN', 'UAP1L1', 'ANGPTL2', 'UCN', 'RRM2')

normCounts <- as.data.frame(counts(dds, normalized=T))
normCounts$ens <- gsub("\\..*","", rownames(normCounts))
normCounts$symbol <- mapIds(org.Hs.eg.db, keys = normCounts$ens, keytype = 'ENSEMBL', column = 'SYMBOL')
normCounts$map <- mapIds(org.Hs.eg.db, keys = normCounts$ens, keytype = 'ENSEMBL', column = 'MAP')
normCounts[rownames(normCounts[is.na(normCounts$symbol), ]), 'symbol'] <- rownames(normCounts[is.na(normCounts$symbol),])
normCounts = subset(normCounts, select = -c(ens))

#write.table(((res_df))[c('symbol', 'map', 'padj', 'log2FoldChange')], file=paste('results/res_df.txt'), sep="\t", quote=F, col.names=NA)
h1 <- strsplit(normCounts[grep("^H1", normCounts$symbol),]$symbol, split=' ')
h2a <- strsplit(normCounts[grep("^H2A", normCounts$symbol),]$symbol, split=' ')
h2b <- strsplit(normCounts[grep("^H2B", normCounts$symbol),]$symbol, split=' ')
h3 <- strsplit(normCounts[grep("^H3", normCounts$symbol),]$symbol, split=' ')
h4 <- strsplit(normCounts[grep("^H4", normCounts$symbol),]$symbol, split=' ')
histones <- c(h1, h2a, h2b, h3, h4)

ifis <- strsplit(normCounts[grep("^IFI", normCounts$symbol),]$symbol, split=' ')
irfs <- strsplit(normCounts[grep("^IRF", normCounts$symbol),]$symbol, split=' ')
ifs <- c(ifis, irfs)

ASs <- strsplit(normCounts[grep("-AS", normCounts$symbol),]$symbol, split=' ')

linc <- strsplit(normCounts[grep("^LINC", normCounts$symbol),]$symbol, split=' ')

hsp <- strsplit(normCounts[grep("^HSP", normCounts$symbol),]$symbol, split=' ')
dnaj <- strsplit(normCounts[grep("^DNAJ", normCounts$symbol),]$symbol, split=' ')
fkbp <- strsplit(normCounts[grep("^FKBP", normCounts$symbol),]$symbol, split=' ')
hsp_dnaj_fkbp <- c('AIPL1', hsp, dnaj, fkbp)

IL1R <- strsplit(normCounts[grep("^IL1R", normCounts$symbol),]$symbol, split=' ')
IL1 <- strsplit(normCounts[grep("^IL1", normCounts$symbol),]$symbol, split=' ')
ILs <- strsplit(normCounts[grep("^IL", normCounts$symbol),]$symbol, split=' ')
IL_sel <- c('IL1R1', 'IL1RL2', 'IL1RL1', 'IL1A', 'IL1B', 'IL1RN', 'IL36RN', 'IL36B', 'IL1F10', 'IL1RAPL1', 'IL6')

POT1 <- c('ACD', 'ACOT7', 'ACTB', 'ACTN4', 'ACY1', 'AFAP1L2', 'AHCY', 'AHNAK', 'AIPL1', 'ALDH1A1', 'ALDH3A1', 'AMPD2', 'ANKMY2', 'ANXA2', 'ANXA4', 'APPL2', 'ARHGDIA', 'ARID3B', 'ARRB1', 'BAG3', 'BCAS2', 'BIN2', 'C15orf57', 'C2orf74', 'CALD1', 'CAMK1D', 'CCDC9', 'CCM2', 'CFL1', 'CFL2', 'CKB', 'CLIC3', 'CNST', 'CORO1A', 'COX6A2', 'CPNE3', 'CPPED1', 'CRK', 'CRYGS', 'CSNK2B', 'CYP4F11', 'DBN1', 'DBNL', 'DCX', 'DDX19B', 'DNPH1', 'DOK2', 'DPP3', 'DPYSL3', 'ECI1', 'EEF1D', 'EIF3G', 'EIF4B', 'ENO2', 'ENSA', 'EPB41L1', 'EVL', 'FAM131B', 'FBP1', 'FES', 'GAMT', 'GAPDH', 'GAS2L1', 'GFPT2', 'GNMT', 'GPA33', 'GPR52', 'GRN', 'HAAO', 'HLCS', 'HMOX1', 'HNMT', 'HOXA3', 'HSP90AB1', 'HSPA1A', 'IFRD2', 'IL1RN', 'ISYNA1', 'IVL', 'KHDRBS1', 'KIAA1191', 'KRT18', 'LAMC3', 'LASP1', 'LDHA', 'LDHB', 'MADD', 'MAGEA4', 'MAP4', 'MAP4K2', 'MAP7', 'MDM2', 'MICA', 'MT1X', 'MVK', 'MVP', 'MYO5C', 'NAP1L1', 'NCDN', 'NOL3', 'NUDC', 'NUDCD2', 'NXNL1', 'PACSIN1', 'PACSIN2', 'PAGE2', 'PAGE5', 'PAK4', 'PALM', 'PCP4', 'PDE1B', 'PDLIM2', 'PEX5', 'PFKP', 'PGLS', 'PGM1', 'PGM2', 'PHYHD1', 'PHYKPL', 'PIPOX', 'PRMT7', 'PROSER2', 'RBKS', 'RECQL4', 'RGS14', 'RIF1', 'RPAP1', 'RPSA', 'RTN4', 'SARS', 'SBDS', 'SERTAD1', 'SH3BP1', 'SNCG', 'STIP1', 'STUB1', 'SULT1B1', 'SULT1C2', 'SULT4A1', 'SYAP1', 'TAGLN', 'TBCD', 'TERF1', 'TERF2', 'TINF2', 'TMSB10', 'TMSB4Y', 'TNKS', 'TOMM34', 'TPI1', 'TPP1', 'TRIM16', 'TRIP10', 'TUBB2A', 'TUBB4B', 'TWF2', 'WIBG', 'WIPI2', 'XAGE2', 'YWHAE', 'YWHAG', 'ZBED2', 'ZBTB49', 'ZFP36L1', 'ZNF32', 'ZNF790')
TINF2 <- c('ACD', 'ACOT7', 'ACTB', 'ADA', 'AFAP1L2', 'AIPL1', 'AK1', 'ANKMY2', 'ANXA4', 'ANXA5', 'APOBEC3F', 'APPL2', 'ARHGDIA', 'ARID3B', 'BAG3', 'BIN2', 'CALD1', 'CCDC43', 'CCDC9', 'CKB', 'CLK3', 'CPNE3', 'CRYGS', 'CTTN', 'DBN1', 'DCX', 'DHFRL1', 'DPP3', 'EEF1D', 'EIF4B', 'ENO2', 'FAM131B', 'FERMT3', 'FKBP6', 'GAPDH', 'GNMT', 'HNMT', 'HOXA3', 'HSP90AB1', 'IPO5', 'KCTD17', 'LASP1', 'LDHA', 'LGALSL', 'LTA4H', 'MAP2K3', 'NCDN', 'NOL3', 'NUDC', 'NUDCD2', 'OR2H1', 'PAGE2', 'PAK1IP1', 'PEX5', 'PFKP', 'PGM1', 'PGM2', 'PHPT1', 'PKM', 'POT1', 'PPP1R2', 'PPP6R3', 'PRKCB', 'REM2', 'RGS14', 'RPSA', 'SARS', 'SCRN2', 'SIAH2', 'STUB1', 'SULT1C2', 'TAGLN', 'TBL1X', 'TERF1', 'TERF2', 'TOMM34', 'TPP1', 'TRIM15', 'TRIM16', 'TRIP10', 'TUBB', 'TUBB4B', 'TXNDC17', 'UCHL1', 'VIPR1', 'YWHAG', 'ZFP36L1', 'ZNF790')
ACD <- c('ACTB', 'ADPRH', 'AFAP1L2', 'AIPL1', 'ANKMY2', 'APOBEC3F', 'BAG3', 'CKB', 'CTTN', 'DBN1', 'DBNL', 'DCX', 'DPYSL3', 'EIF3G', 'EIF4B', 'ENO2', 'ENSA', 'FAM131B', 'FKBP6', 'GAPDH', 'HSP90AB1', 'IPO5', 'IVL', 'LASP1', 'LDHA', 'LGALSL', 'LLGL1', 'LRRC25', 'NCDN', 'NUDC', 'NUDCD2', 'PAGE2', 'PDLIM2', 'PEX5', 'PFKP', 'PGM2', 'POT1', 'RBKS', 'RECQL4', 'RGS3', 'RPSA', 'SBDS', 'STIP1', 'STUB1', 'SULT1C2', 'TAGLN', 'TBC1D10A', 'TBCD', 'TINF2', 'TOMM34', 'TRIM15', 'TRIM16', 'TUBB2A', 'TUBB4B', 'USP7', 'XRCC6', 'YWHAE', 'ZFP36L1', 'ZNF790')
RAP1A <- c('HNRNPD', 'PPP4R2', 'SYVN1', 'UQCRC2', 'CLPP', 'CHMP4C', 'MARCKS', 'GANAB', 'TRIM67', 'ECHS1', 'KIFAP3', 'TSC2', 'RADIL', 'RAPGEF1', 'SMARCA2', 'RABIF', 'DGKI', 'FAS', 'PCDHA8', 'CDC34', 'CALCOCO2', 'RAB7A', 'RGL4', 'HNRNPL', 'CCDC85A', 'HSPA5', 'RAPGEF5', 'RAB8B', 'PARL', 'RGS14', 'BRAF', 'ARHGDIA', 'RUNDC3A', 'PPP2R1A', 'RASGRP4', 'TMEM31', 'HYPM', 'TCEB3', 'SOD1', 'LIMA1', 'CDC42', 'NTRK1', 'RAB5B', 'RAB5C', 'KRIT1', 'PRIM1', 'DUSP22', 'KLRC4', 'IFT22', 'PFN1', 'OR2T10', 'GPR113', 'CPLX2', 'TRIM25', 'RAPGEF4', 'ILF3', 'CENPM', 'MLLT4', 'RAP1GAP', 'MEF2BNB', 'EXOC6', 'DUSP9', 'DPF2', 'ESR1', 'FADD', 'SAAL1', 'CACUL1', 'RAB8A', 'TBRG4', 'KIF14', 'IQGAP1', 'DLST', 'CAV2', 'GNB1', 'PDE6D', 'POLE3', 'INSC', 'KBTBD4', 'RAF1', 'PON2', 'PARK2', 'KLK10', 'RBM39', 'RGL1', 'HNRNPM', 'PCTP', 'AKAP1', 'FEN1', 'TNFRSF10C', 'TACR3', 'OSGEP', 'HDAC1', 'RAPGEF6', 'RAP1GDS1', 'LGALS9', 'BMX', 'RAPGEF2', 'BIN1', 'DNAJC8', 'FGFR1', 'ARHGEF1', 'MMGT1', 'ELAVL1', 'HSPA4', 'XRCC6', 'C10ORF91', 'GNG12', 'RHOB', 'RYBP', 'HNRNPDL', 'CUL4A', 'LYN', 'DUSP19', 'RHEB', 'CCDC53', 'IMMP2L', 'CD53', 'ESR2', 'MRPL34', 'RAB1A', 'RALGDS', 'RAB35', 'MGST3', 'CHST12', 'RAB5A', 'GABARAPL2', 'OR6N1', 'TAX1BP1', 'MTOR', 'KIAA1429', 'GPC3', 'FAF1', 'UNK', 'AR', 'MAGEA3', 'PMAIP1', 'STX7', 'BSG', 'RAPGEF3', 'RPL35A', 'MTNR1A', 'SIRT1', 'LAMP2', 'KRAS', 'CDK11B', 'PIGO')
POT1_TINF2_ACD_RAP1A <- unique(c(POT1, TINF2, ACD, RAP1A))

P53 <- c('MDM2', 'MDM4', 'RFWD2', 'TP53', 'TP63', 'TP73')
P53_panther <- c('AKT1', 'AKT2', 'AKT3', 'ATM', 'ATR', 'CCNA1', 'CCNA2', 'CCNE1', 'CCNG1', 'CDK2', 'CDKN1A', 'CTNNB1', 'HRAS', 'KRAS', 'MAPK11', 'MAPK12', 'MAPK13', 'MAPK14', 'MYC', 'NRAS', 'PDPK1', 'PIK3C2A', 'PIK3C2B', 'PIK3C2G', 'PIK3C3', 'PIK3CA', 'PIK3CB', 'PIK3CD', 'PIK3CG', 'PIK3R1', 'PIK3R2', 'PIK3R3', 'PIK3R5', 'PPP2CA', 'PPP2CB', 'PTEN', 'RB1', 'RBL1', 'SIAH1', 'TP53', 'TP63', 'TP73')
YWHA <- strsplit(normCounts[grep("^YWHA", normCounts$symbol),]$symbol, split=' ')
YWHA <- c("YWHAE","YWHAH","YWHAQ","YWHAZ","YWHAB","YWHAG")
TP53 <- c(P53, P53_panther, YWHA)

H2AZ <- strsplit(normCounts[grep("^H2AZ", normCounts$symbol),]$symbol, split=' ')

AIPL1_pittEdu <- c('ACD', 'HSP90AA1', 'HSPA8', 'NUB1', 'POT1', 'TINF2')
AIPL1_BioGrid <- c('optAIPL1g', 'NUB1', 'PDE5A', 'BLM', 'EDRF1', 'HSPA8', 'A2ML1', 'ACD', 'ACPP', 'ALOX12B', 'ANXA8', 'ARHGAP1	', 'C3', 'C4A', 'CALML5', 'CAPNS2', 'CASP14', 'CBR1', 'CPA4', 'CXORF57	', 'DTX2', 'EVPL', 'FGB', 'FGFR1OP', 'FGG', 'FKBP15', 'FLG', 'HAL', 'HERC2', 'HIST2H3PS2', 'HMGCS1', 'HMOX1', 'HPGD', 'HSP90AA1', 'HSP90AA1', 'HSP90AB1', 'HSPA4', 'HSPA8', 'IGH', 'IGHG1', 'IGHG2', 'IGHG3', 'IGHM', 'IGLC6', 'IL37', 'IVL', 'KLK7', 'KPNA5', 'KPRP', 'LBX1', 'LCN2', 'LCP1', 'LOR', 'MIOS', 'MTMR3', 'NADSYN1', 'NCCRP1', 'NUB1', 'ORM1', 'PATZ1', 'PKP1', 'PLSCR3', 'POF1B', 'POT1', 'PPL', 'RAB14', 'RAB5A', 'RAI14', 'RNASE7', 'S100A7A', 'SCPEP1', 'SDR9C7', 'SEC16A', 'SEPT2', 'SERPINA3', 'SERPINB4', 'SMAP1', 'SMTNL2', 'SPRR1A', 'SPRR1B', 'STAT3', 'SUCLG2', 'SULT2B1', 'SUPT5H', 'TF', 'TGM3', 'TINF2', 'TP53', 'TREX2', 'TXNL4B', 'TYMP', 'UBA6', 'UBD', 'VIM', 'WASF1', 'WDR24')

NFKB1 <- c('NFKB1', 'PIK3CA', 'CTNNB1', 'EP300', 'HDAC2', 'MYC', 'NFKBIZ', 'KAT5', 'PML', 'NCOA1', 'PTGS2', 'RXRA', 'NR4A1', 'AURKA', 'CHEK1', 'CUL4A', 'TP63', 'CUL3', 'CUL1', 'HDAC1', 'ETS1', 'CEBPB', 'NOTCH1', 'RELA', 'IL10RA', 'STAT6', 'STAT3', 'JUN', 'GSK3B', 'AURKB', 'JAK1', 'CUL4B', 'PARP1')

OAS <- c('OAS1', 'OAS2', 'OAS3', 'OASL', 'PKR', 'MDA5', 'MAVS')

all <- strsplit(normCounts$symbol, split=' ')

heatList <- function(sigs, NameOfList) {
  sigs <- normCounts[normCounts$symbol %in% NameOfList, ]
  #keep <- rowSums(sigs >= 5) >= 4
  #sigs <- sigs[keep,]
  
  #sigs <- sigs[order(sigs$map),]
  matr.z <- t(apply(sigs[,1:11], 1, scale))
  matr.z <- na.omit(matr.z)
  colnames(matr.z) <- paste(meta[rownames(meta)==colnames(sigs),'sampletype'], meta[rownames(meta)==colnames(sigs),'repl'], sep='-')
  #colnames(matr.z) <- colnames(sigs[,1:11])
  fig <- Heatmap(matr.z, cluster_rows = T, cluster_columns = T, column_labels = colnames(matr.z), 
                 name = 'Z.score', column_names_rot = 45, 
                 row_labels = paste(sigs[rownames(matr.z), 'symbol'], sigs[rownames(matr.z), 'map'], sep=' - '),
                 row_names_gp = gpar(fontsize = 12), 
                 column_title = paste('Heatmap of ', deparse(substitute(NameOfList)), sep=''))
  png(paste('plots/heat_', deparse(substitute(NameOfList)), '.png', sep=''), res = 1200, width = 7000, height = 10000)
  print(fig)
  dev.off()
}
heatList(sigs=normCounts, NameOfList=ILs)

##for patent:
heatList <- function(sigs, NameOfList) {
  sigs <- normCounts[normCounts$symbol %in% NameOfList, ]
  #keep <- rowSums(sigs >= 5) >= 4
  #sigs <- sigs[keep,]
  
  #sigs <- sigs[order(sigs$map),]
  matr.z <- t(apply(sigs[,1:11], 1, scale))
  matr.z <- na.omit(matr.z)
  colnames(matr.z) <- paste(meta[rownames(meta)==colnames(sigs),'sampletype'], meta[rownames(meta)==colnames(sigs),'repl'], sep='-')
  #colnames(matr.z) <- colnames(sigs[,1:11])
  fig <- Heatmap(matr.z, cluster_rows = T, cluster_columns = T, column_labels = colnames(matr.z), 
                 name = 'Z.score', column_names_rot = 45, 
                 row_labels = paste(sigs[rownames(matr.z), 'symbol'], sigs[rownames(matr.z), 'map'], sep=' - '),
                 row_names_gp = gpar(fontsize = 7), 
                 column_title = paste('Heatmap of ', deparse(substitute(NameOfList)), sep=''), 
                 col = colorRamp2(breaks=c(-1, 1), colors = c("#F5f5f5", "black"), transparency = 0, space = "LAB", hcl_palette = NULL, reverse = FALSE))
  png(paste('plots/heat_', deparse(substitute(NameOfList)), '.png', sep=''), res = 1200, width = 7000, height = 10000)
  print(fig)
  dev.off()
}
heatList(sigs=normCounts, NameOfList=ILs)

hist(rowSums(sigs[,1:11]), breaks = 1000)

intersect_AIPL1_POT1 <- intersect(POT1, AIPL1_BioGrid)
intersect_AIPL1_shelterin <- intersect(AIPL1_BioGrid, c('TERF1', 'TERF2', 'TERF2IP', 'TINF2', 'ACD', 'POT1'))

genes_to_test <- AIPL1_BioGrid
GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont="BP")
as.data.frame(GO_results)
if (nrow(GO_results) != 0) {
  fig <- barplot(GO_results, showCategory = 20) + ggtitle('GO of AIPL1_BioGrid') + theme(plot.title = element_text(face = 'bold', hjust = 0.5), axis.title.x = element_text(size = 16), axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 10))
  ggsave(fig, file='plots/GO_AIPL1_BioGrid.png', width = 14, height = 7, dpi = 1200)
}


################################################################################
################################### N of DEGs ##################################
################################################################################

deg <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(deg) <- c("sample", "padj", "l2fc", 'N')
N_of_DEG <- function(samp, pval, lfc) {
  new <- c(samp,pval,lfc, nrow(sigs))
  deg[nrow(deg) + 1, ] <<- new 
}

for (pval in c(0.05, 0.01, 0.005, 0.001)) {
  for (lfc in c(0.58, 1, 1.58, 2, 2.32)) {
    for (samp in c('wt_ctrl','opt_ctrl', 'wt_gfp', 'opt_gfp', 'opt_wt', 'gfp_ctrl')) {
      shrink(samp, pval, lfc)
      filtr(samp, pval, lfc)
      N_of_DEG(samp, pval, lfc)
    }
  }
}
#write.table(deg, file=paste('results/deg.txt'), sep="\t", quote=F, col.names=T)

### plot:
for (samp in c('wt_ctrl','opt_ctrl', 'wt_gfp', 'opt_gfp', 'opt_wt', 'gfp_ctrl')) {
  deg_sub <- filter(deg, sample==samp)
  N <- c(deg$N)
  plt <- ggplot(deg_sub, aes(x = padj, y = l2fc)) + ggtitle(samp) +
    theme(plot.title = element_text(hjust = 0.5)) +
    geom_point(shape=1, size = 3*log10(as.integer(deg_sub$N))) + 
    geom_text(aes(label = N), vjust = -2.5) 
  ggsave(plt, file=paste('plots/degs_', samp, '.png', sep=""), 
         width = 7, height = 7, dpi = 600)
}

ggplot(deg, aes(x = padj, y = l2fc, shape=deg$samp, color=deg$samp)) + 
  geom_point(size = 3*log10(as.integer(deg$N))) + 
  #geom_point(shape=1, size = 3*log10(as.integer(deg$N))) + 
  geom_text(aes(label = N), vjust = deg$samp) +
  geom_vline(xintercept = deg$padj) +
  geom_hline(yintercept = deg$l2fc)

ggplot(NULL) +
  geom_vline(xintercept = c(0.05, 0.01, 0.005, 0.001)) +
  geom_hline(yintercept = c(0.58, 1, 1.58, 2, 2.32)) +
  xlab("p-value adjusted") +
  ylab("log2 fold change")

#wt/ctrl
ggplot(res_df, aes(x=abs(log2FoldChange), y=padj)) + 
  geom_point(size=1, shape=23)

################################################################################
##################################### Venns ####################################
################################################################################

i <- 1
#list(c('wt_ctrl', 'wt_gfp'), c('opt_ctrl', 'opt_gfp'), c('wt_gfp', 'opt_gfp'), c('wt_ctrl', 'opt_ctrl'))
for (samp in list(c('wt_ctrl', 'opt_ctrl')) ) {
  #print(samp[2])
  for (dir in c('up', 'down')) {
    venn.diagram(
      x = list(row.names(get(paste('sigs_', dir, '_', samp[1], sep=''))), row.names(get(paste('sigs_', dir, '_', samp[2], sep='')))),
      category.names = c(samp[1], samp[2]),
      filename = paste('plots/venn_', as.character(i), '_', dir, '.png', sep=''),
      output=TRUE, 
      fill = c(2,4), 
      main = paste('Venn diagram of ', samp[1], ' and ', samp[2], ', ', dir, '-regulated, p-adj = ', pval, ', l2fc = ', lfc, sep='')
    )
  }
  i <- i + 1
}

### mega-venns
for (dir in c('up', 'down')) {
  venn.diagram(
    x = list(row.names(get(paste('sigs_', dir, '_wt_ctrl', sep=''))), row.names(get(paste('sigs_', dir, '_wt_gfp', sep=''))), row.names(get(paste('sigs_', dir, '_opt_ctrl', sep=''))), row.names(get(paste('sigs_', dir, '_opt_gfp', sep='')))),
    category.names = c("AAV9-AIPL1wt/\nnon-transd.ctrl", "AAV9-AIPL1wt/\nAAV9-GFP", 'AAV9-AIPL1co/\nnon-transd.ctrl', 'AAV9-AIPL1co/\nAAV9-GFP'),
    filename = paste('plots/venn_', dir, '.png', sep=''),
    output=TRUE, 
    fill = c(1,2,3,4), force.unique=T, cat.cex=0.7,
    main = paste('Venn diagram, ', dir, '-regulated', sep=''), 
    height = 3000, width = 3000, resolution = 600
  )
}

intersect1_wt_down <- intersect(sigs_down_wt_gfp$symbol, sigs_down_wt_ctrl$symbol)
intersect1_wt_up <- intersect(sigs_up_wt_gfp$symbol, sigs_up_wt_ctrl$symbol)

intersect2_opt_down <- intersect(sigs_down_opt_gfp$symbol, sigs_down_opt_ctrl$symbol)
intersect2_opt_up <- intersect(sigs_up_opt_gfp$symbol, sigs_up_opt_ctrl$symbol)

intersect3_gfp_down <- intersect(sigs_down_opt_gfp$symbol, sigs_down_wt_gfp$symbol)
intersect3_gfp_up <- intersect(sigs_up_opt_gfp$symbol, sigs_up_wt_gfp$symbol)

intersect4_ctrl_down <- str_split(intersect(sigs_down_opt_ctrl$symbol, sigs_down_wt_ctrl$symbol), ' ')
intersect4_ctrl_up <- str_split(intersect(sigs_up_opt_ctrl$symbol, sigs_up_wt_ctrl$symbol), ' ')
#exclusion <- list(list(str_split(sigs_down_opt_ctrl$symbol, ' ')) - list(str_split(sigs_down_wt_ctrl$symbol, ' ')))

intersect_all_down <- intersect(intersect3_gfp_down, intersect4_ctrl_down)
intersect_all_up <- intersect(intersect3_gfp_up, intersect4_ctrl_up)


sigs_opt_up_ctrl_minus_gfp <- filter(sigs_up_opt_ctrl[!rownames(sigs_up_opt_ctrl) %in% rownames(sigs_up_opt_gfp), ])
write.table(sigs_opt_up_ctrl_minus_gfp, 
            file=paste('results/sigs_opt_up_ctrl_minus_gfp.txt'), sep="\t", quote=F, col.names=NA)

matr <- counts(dds, normalized=T)[rownames(sigs_opt_up_ctrl_minus_gfp),]
matr.z <- t(apply(matr, 1, scale))
colnames(matr.z) <- colnames(matr)
fig <- Heatmap(matr.z, cluster_rows = T, cluster_columns = T, column_labels = colnames(matr.z), 
               name = 'Z.score', row_labels = sigs_opt_up_ctrl_minus_gfp[rownames(matr.z), 'symbol'], 
               column_title = paste('Heatmap of ', as.character(samp),', lfc=',as.character(lfc),', padj=',as.character(pval), sep=''))
#png(paste('plots/heat_', samp, '_lfc=', lfc, '_padj=', pval, '.png', sep=''), 
#    res = 250, width = 1500, height = 2000)
#print(fig)
#dev.off()
matr.zz <- matr.z[rowSums(matr.z[] > 1.5) > 2, ]
sigs_opt_up_ctrl_minus_gfp[rownames(matr.zz), 'symbol']

genes_to_test <- intersect4_ctrl_down
GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont="BP")
as.data.frame(GO_results)
if (nrow(GO_results) != 0) {
  fig <- barplot(GO_results, showCategory = 20) + 
    ggtitle(paste('down')) +
    theme(plot.title = element_text(face = 'bold', hjust = 0.5), 
          axis.title.x = element_text(size = 16), 
          axis.text.x = element_text(size = 12), 
          axis.text.y = element_text(size = 10))
  ggsave('plots/GO_down.png', plot=fig, width = 14, height = 7, dpi = 600) 
}


################################################################################
################################## article figs ################################
################################################################################

fig1 <- data.frame(genes=rep(c('OPN3', 'KLF4', 'TUBB3'), 3), 
                   sample = factor(rep(c('control','FEN', 'bFGF'), each=3), levels=c('control','FEN', 'bFGF')), 
                   RQ = c(1, 1, 1.02, 0.46, 1.45, 3.26, 1.51, 1.79, 3.18),
                   sd = c(0.03, 0.09, 0.27, 0.06, 0.11, 0.05, 0.21, 0.22, 0.14))
fig1_bar <- ggplot(fig1, aes(x=sample, y=RQ, fill=genes)) + 
  geom_bar(stat="identity", position=position_dodge(), colour='black', linewidth=0.1) +
  geom_errorbar(aes(ymin=RQ-sd, ymax=RQ+sd), width=.2,
                position=position_dodge(.9)) + 
  ylab('Relative gene expression') + xlab('') +
  theme(axis.text=element_text(colour="black"), axis.title.y=element_text(size=20), axis.text.x=element_text(size=16)) +
  scale_fill_brewer(palette="Blues") 
ggsave(fig1_bar, file='plots/fig1_bar.png', width = 12, height = 7, dpi = 1200)

fig3 <- data.frame(genes=rep(c('NUB1', 'AHR', 'LCA5', 'RDH12', 'ND4', 'OPN3', 'RPE65', 'KLF4', 'AIPLI1', 'TUBB3', 'PDE6B'), each = 2),
                   sample=factor(rep(c('non-\ndifferetiated', 'differetiated\nwith fenretinide'), length.out = 22)),
                   RQ=c(1.00, 0.21, 1.00, 0.23, 1.00, 0.42, 1.00, 1.23, 1.00, 4.15, 1.00, 0.02, 1.02, 113.76, 1.01, 8.06, 1.01, 9.84, 1.00, 130.66, 1.00, 0.78),
                   sd=c(0.11, 0.06, 0.00, 0.04, 0.05, 0.02, 0.08, 0.13, 0.06, 0.11, 0.00, 0.00, 0.30, 21.25, 0.18, 0.44, 0.17, 1.41, 0.11, 29.40, 0.11, 0.03))
fig3 <- fig3[fig3$sample=='differetiated\nwith fenretinide',]
fig3_bar <- ggplot(fig3, aes(x=reorder(genes,RQ), y=RQ)) + 
  geom_bar(stat="identity", position=position_dodge(), fill="#2171B5") + 
  ylab('Relative gene expression') + xlab('') +
  theme(axis.text=element_text(colour="black"), axis.title.y=element_text(size=16), 
        axis.text.x=element_text(size=12)) +
  geom_errorbar(aes(ymin=RQ-sd, ymax=RQ+sd), width=.2,
                position=position_dodge(.9)) + scale_y_continuous(trans = "log10")
ggsave(fig3_bar, file='plots/fig3_bar.png', width = 10, height = 7, dpi = 1200)

################################################################################
########################### Lists of genes for article #########################
################################################################################

genes_intro <- c('AIPL1 ', 'DNAJA2', 'PDE6', 'FKBP', 'PDEα', 'FKBP51', 
                    'FKBP52', 'HSP90', 'SUB1', 'NUB1', 'EB1', 'NEDD8', 'UBD', 
                    'HSP70', 'GNAT1', 'GNAT2', 'GRK1')
genes_string <- c('NUB1', 'GUCY2D', 'AHR', 'CRB1', 'TULP1', 'RDH12', 
                    'AIP', 'LCA5', 'CRX', 'RPE65', 'TULP1', 'LBX1', 'TXNL4B', 
                    'PLSCR3', 'SMAP1')
genes_OMIM <- c('PATZ1', 'POT1', 'TXNL4B', 'WDR24', 'HSPA8', 'HSP90AA1', 'UBA6', 
                'MIOS', 'TINF2', 'FGFR1OP', 'WASF1', 'SEC16A', 'SUPT5H', 'SMAP1', 
                'TP53', 'UBD', 'LBX1', 'FKBP15', 'HSPA4', 'BLM', 'ACD', 'PDE5A')

intersect <- intersect(genes_OMIM, sigs_wt_gfp$symbol)


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