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
#install.packages(<path to downloaded file with DEGreport tar.gz>, repos = NULL)
library(DEGreport)
#BiocManager::install("vsn")
library(vsn)
library(plyr) 
library(dplyr)

## reading file:
df <- read.csv("data/counts.Sirius.csv",sep = ";") 

## Trying to reorder columns by replicates alphabetically:
coln1 <- sort(c(colnames(df[5:16]))) 
coln2 <- c(colnames(df[1:4]))
coln = c(coln2, coln1)
coln
data <- df[, coln]

## create column with gene names/ID:
data_geneNames <- data %>%  mutate(geneLabel = ifelse(is.na(Symbol), 
                                                      name, Symbol)) 

## some preparations of data dataframe (try to remove duplicates of geneNames):
#drop <- c("gene", "Refseq","Gene.name", "Symbol")
#data = data[,!(names(data) %in% drop)]
#data <- data[complete.cases(data['name']),]

#data <- data %>%  mutate(dupl = duplicated(data['gene'])) ## find duplicates of Refseq number 

#data[data$dupl == 'TRUE', ] == #... need to finish or delete

row.names(data) <- data[, 1]

drop <- c('CC1', coln2, 'gene', 'dupl') #remove CC1 from analysis, cause it doesn't cluster well
data = data[,!(names(data) %in% drop)]


meta <- data.frame(colnames(data[, 1:11]))
meta['sampletype'] <- c('AAV9_AIPL1wt', 'AAV9_AIPL1wt', 'AAV9_AIPL1wt', 
                        'AAV9_AIPL1opt', 'AAV9_AIPL1opt', 'AAV9_AIPL1opt', 
                        'control', 'control', 
                        'AAV9_GFP', 'AAV9_GFP', 'AAV9_GFP')
#meta
row.names(meta) <- meta[,1]
meta <- subset(meta, select = c('sampletype'))

all(colnames(data) %in% rownames(meta))
all(colnames(data) == rownames(meta))


################################# PRE-FILTERING ################################

## creating DESeq DataSet
dds <- DESeqDataSetFromMatrix(countData = data, colData = meta, design = ~ sampletype)
## Run DESeq2 differential expression analysis
dds <- DESeq(dds)

##  **Optional step** - Output normalized counts to save as a file to access outside RStudio
normalized_counts <- counts(dds, normalized=TRUE) 
write.table(normalized_counts, file="data/normalized_counts.txt", sep="\t", quote=F, col.names=NA) 

## Total number of raw counts per sample
colSums(counts(dds))

## Total number of normalized counts per sample
colSums(counts(dds, normalized=T))

## log2fold change calculation and MAplot:
plotMA(dds, main = "dds")
res <- results(dds)
plotMA(res, alpha = 0.05, main = "res", ylim=c(-7,7))

####################################### QC #####################################

# Transform counts for data visualization
rld <- rlog(dds, blind=TRUE)
plotPCA(rld, intgroup="sampletype")

# Extract the rlog matrix from the object
rld_mat <- assay(rld)
# Compute pairwise correlation values
rld_cor <- cor(rld_mat)
# Plot heatmap
heat.colors <- brewer.pal(6, "Blues")
pheatmap(rld_cor, color = heat.colors, fontsize = 10, fontsize_row = 10, height=20)


#################################### Shrinkage #################################

contrast_wt <- c("sampletype", "AAV9_AIPL1wt", "control")
contrast_opt <- c("sampletype", "AAV9_AIPL1opt", "control")
contrast_wt_opt <- c("sampletype", "AAV9_AIPL1wt", "AAV9_AIPL1opt")
res_unshrunken_wt <- results(dds, contrast = contrast_wt, alpha = 0.05, pAdjustMethod = 'bonferroni')
res_unshrunken_opt <- results(dds, contrast = contrast_opt, alpha = 0.05, pAdjustMethod = 'bonferroni')
res_unshrunken_wt_opt <- results(dds, contrast = contrast_wt_opt, alpha = 0.05, pAdjustMethod = 'bonferroni')
res_shrunken_wt <- lfcShrink(dds, contrast = contrast_wt, res=res_unshrunken_wt, type='normal')
res_shrunken_opt <- lfcShrink(dds, contrast = contrast_opt, res=res_unshrunken_opt, type='normal')
res_shrunken_wt_opt <- lfcShrink(dds, contrast = contrast_wt_opt, res=res_unshrunken_wt_opt, type='normal')

#res_shrunken_wt %>% data.frame() %>% View()
summary(res_shrunken_wt_opt)
sum(res_shrunken_wt_opt$padj < 1.580078e-06, na.rm=TRUE)

0.05/31644

## How many genes are differentially expressed compared to control:
padj.cutoff <- 0.05/31644
lfc.cutoff <- 0.58
# Turn the results object into a data frame
res_df_wt <- data.frame(res_shrunken_wt)
res_df_opt <- data.frame(res_shrunken_opt)
res_df_wt_opt <- data.frame(res_shrunken_wt_opt)

# Subset the significant results
sig_res_wt <- filter(res_df_wt, padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff)
sig_res_opt <- filter(res_df_opt, padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff)
sig_res_wt_opt <- filter(res_df_wt_opt, padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff)

#################### Dot-plot expression of a single gene ######################
# ENSG00000129221 - AIPL1 gene 
# Save plotcounts to a data frame object
d1 <- plotCounts(dds, gene="ENSG00000129221", intgroup="sampletype", returnData=TRUE)
# Plotting AIPL1 normalized counts, using the samplenames (rownames of d as labels)
ggplot(d, aes(x = sampletype, y = count, color = sampletype)) + 
  geom_point(position=position_jitter(w = 0.1,h = 0)) +
  geom_text_repel(aes(label = rownames(d))) + 
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






## extract df with differentially expressed up-regulated genes:
res_df_wt_up_TF <- res_df_wt %>% 
  mutate(dif_expr = padj < 1.580078e-06 & log2FoldChange >= 0.58)
res_df_opt_up_TF <- res_df_opt %>% 
  mutate(dif_expr = padj < 1.580078e-06 & log2FoldChange >= 0.58)
res_df_wt_opt_up_TF <- res_df_opt %>% 
  mutate(dif_expr = padj < 1.580078e-06 & log2FoldChange >= 0.58)

de_up_wt <- res_df_wt_up_TF[which(res_df_wt_up_TF$dif_expr == 'TRUE'), ]
de_up_opt <- res_df_opt_up_TF[which(res_df_opt_up_TF$dif_expr == 'TRUE'), ]
de_up_wt_opt <- res_df_wt_opt_up_TF[which(res_df_wt_opt_up_TF$dif_expr == 'TRUE'), ]

df_genes <- df[, 1:4]
row.names(df_genes) <- df_genes[,1]
drop <- c('name') 
df_genes = df_genes[,!(names(df_genes) %in% drop)]
## Merge:
de_up_wt <- merge(x = de_up_wt, y = df_genes, by = 0, all.x = TRUE)
de_up_opt <- merge(x = de_up_opt, y = df_genes, by = 0, all.x = TRUE)
de_up_wt_opt <- merge(x = de_up_wt_opt, y = df_genes, by = 0, all.x = TRUE)

de_up_wt <- de_up_wt %>% arrange(desc(log2FoldChange))
de_up_opt <- de_up_opt %>% arrange(desc(log2FoldChange))
de_up_wt_opt <- de_up_opt %>% arrange(desc(log2FoldChange))

write.table(de_up_wt[c('Row.names','Symbol', 'Gene.name', 'padj', 'log2FoldChange')], file="results/de_up_wt_vs_ctrl.txt", sep="\t", quote=F, col.names=NA) 
write.table(de_up_opt[c('Row.names','Symbol', 'Gene.name', 'padj', 'log2FoldChange')], file="results/de_up_opt_vs_ctrl.txt", sep="\t", quote=F, col.names=NA) 
write.table(de_up_wt_opt[c('Row.names','Symbol', 'Gene.name', 'padj', 'log2FoldChange')], file="results/de_up_wt_vs_opt.txt", sep="\t", quote=F, col.names=NA) 


## extract df with differentially expressed down-regulated genes:
res_df_wt_down_TF <- res_df_wt %>% 
  mutate(dif_expr = padj < 1.580078e-06 & log2FoldChange <= -0.58)
res_df_opt_down_TF <- res_df_opt %>% 
  mutate(dif_expr = padj < 1.580078e-06 & log2FoldChange <= -0.58)
res_df_wt_opt_down_TF <- res_df_wt_opt %>% 
  mutate(dif_expr = padj < 1.580078e-06 & log2FoldChange <= -0.58)

de_down_wt <- res_df_wt_down_TF[which(res_df_wt_down_TF$dif_expr == 'TRUE'), ]
de_down_opt <- res_df_opt_down_TF[which(res_df_opt_down_TF$dif_expr == 'TRUE'), ]
de_down_wt_opt <- res_df_wt_opt_down_TF[which(res_df_wt_opt_down_TF$dif_expr == 'TRUE'), ]

df_genes <- df[, 1:4]
row.names(df_genes) <- df_genes[,1]
drop <- c('name') 
df_genes = df_genes[,!(names(df_genes) %in% drop)]
## Merge:
de_down_wt <- merge(x = de_down_wt, y = df_genes, by = 0, all.x = TRUE)
de_down_opt <- merge(x = de_down_opt, y = df_genes, by = 0, all.x = TRUE)
de_down_wt_opt <- merge(x = de_down_wt_opt, y = df_genes, by = 0, all.x = TRUE)

de_down_wt <- de_down_wt %>% arrange(log2FoldChange)
de_down_opt <- de_down_opt %>% arrange(log2FoldChange)
de_down_wt_opt <- de_down_wt_opt %>% arrange(log2FoldChange)

write.table(de_down_wt[c('Row.names','Symbol', 'Gene.name', 'padj', 'log2FoldChange')], file="results/de_down_wt_vs_ctrl.txt", sep="\t", quote=F, col.names=NA) 
write.table(de_down_opt[c('Row.names','Symbol', 'Gene.name', 'padj', 'log2FoldChange')], file="results/de_down_opt_vs_ctrl.txt", sep="\t", quote=F, col.names=NA) 
write.table(de_down_wt_opt[c('Row.names','Symbol', 'Gene.name', 'padj', 'log2FoldChange')], file="results/de_down_wt_vs_opt.txt", sep="\t", quote=F, col.names=NA) 


################################## Volcano plot ################################
################################# for opt vs ctrl ##############################
## Obtain logical vector where TRUE values denote padj values < 0.05 and fold change > 1.5 in either direction
res_df_opt_TF <- res_df_opt %>% 
  mutate(dif_expr = padj < 1.580078e-06 & abs(log2FoldChange) >= 0.58)

## Create a column to indicate which genes to label 
res_df_opt_TF$geneLabels <- data_geneNames$geneLabel

res_df_opt_TF_top20 <- subset(res_df_opt_TF %>% arrange(log2FoldChange) %>% head(n=20))
res_df_opt_TF_low20 <- subset(res_df_opt_TF %>% arrange(desc(log2FoldChange)) %>% head(n=20))

res_df_opt_TF <- res_df_opt_TF[,!names(res_df_opt_TF) %in% c("geneLabels")]

res_df_opt_TF_40 <- rbind(res_df_opt_TF_top20, res_df_opt_TF_low20)
res_df_opt_TF_40 <- res_df_opt_TF_40[,!names(res_df_opt_TF_top20) %in% 
                                       c("baseMean", 'lfcSE', 'stat', 'pvalue', 'padj', 'dif_expr')]
res_df_opt_TF <- merge(res_df_opt_TF, res_df_opt_TF_40, by = 'row.names', all.x = T)
row.names(res_df_opt_TF) <- res_df_opt_TF$Row.names
drop = c('Row.names', 'log2FoldChange.y')
res_df_opt_TF = res_df_opt_TF[,!(names(res_df_opt_TF) %in% drop)]
colnames(res_df_opt_TF)[2] ="log2FoldChange"

res_df_opt_TF$super <- res_df_opt_TF$log2FoldChange >  0.58  & res_df_opt_TF$padj < 1.580078e-06
res_df_opt_TF$sub   <- res_df_opt_TF$log2FoldChange < -0.58  & res_df_opt_TF$padj < 1.580078e-06
res_df_opt_TF$threshold <- as.factor(abs(res_df_opt_TF$log2FoldChange) > 0.58 & res_df_opt_TF$padj < 1.580078e-06)

#PLOT!
volcano_plot <- ggplot(data=res_df_opt_TF, aes(x=log2FoldChange, y=-log10(padj))) +
  geom_point(data=res_df_opt_TF, size=1, colour="gray") +
  geom_point(data=res_df_opt_TF[res_df_opt_TF$super==TRUE, ], size=2, colour="#CC0000") +
  geom_point(data=res_df_opt_TF[res_df_opt_TF$sub  ==TRUE, ], size=2, colour="#000099") +
  geom_text_repel(aes(label = geneLabels)) +
  xlab("log2 fold change") +
  ylab("-log10 p-value adjusted") +
  ggtitle("Differential expression optAIPL1 vs CC") +
  scale_x_continuous() +
  scale_y_continuous() +
  theme_bw() +
  theme(axis.title.y = element_text(size=16),
        axis.title.x = element_text(size=16, colour="black"),
        axis.text = element_text(size=12),
        legend.title =element_blank() ,
        legend.text = element_text(size = 12)) +
  theme(plot.title = element_text(face="bold", size = 17, hjust = 0.5))
volcano_plot

################################# for wt vs ctrl ###############################
res_df_wt_TF <- res_df_wt %>% 
  mutate(dif_expr = padj < 1.580078e-06 & abs(log2FoldChange) >= 0.58)

## Create a column to indicate which genes to label 
res_df_wt_TF$geneLabels <- data_geneNames$geneLabel

res_df_wt_TF_top20 <- subset(res_df_wt_TF %>% arrange(log2FoldChange) %>% head(n=20))
res_df_wt_TF_low20 <- subset(res_df_wt_TF %>% arrange(desc(log2FoldChange)) %>% head(n=20))

res_df_wt_TF <- res_df_wt_TF[,!names(res_df_wt_TF) %in% c("geneLabels")]

res_df_wt_TF_40 <- rbind(res_df_wt_TF_top20, res_df_wt_TF_low20)
res_df_wt_TF_40 <- res_df_wt_TF_40[,!names(res_df_wt_TF_top20) %in% 
                                       c("baseMean", 'lfcSE', 'stat', 'pvalue', 'padj', 'dif_expr')]
res_df_wt_TF <- merge(res_df_wt_TF, res_df_wt_TF_40, by = 'row.names', all.x = T)
row.names(res_df_wt_TF) <- res_df_wt_TF$Row.names
drop = c('Row.names', 'log2FoldChange.y')
res_df_wt_TF = res_df_wt_TF[,!(names(res_df_wt_TF) %in% drop)]
colnames(res_df_wt_TF)[2] ="log2FoldChange"

res_df_wt_TF$super <- res_df_wt_TF$log2FoldChange >  0.58  & res_df_wt_TF$padj < 1.580078e-06
res_df_wt_TF$sub   <- res_df_wt_TF$log2FoldChange < -0.58  & res_df_wt_TF$padj < 1.580078e-06
res_df_wt_TF$threshold <- as.factor(abs(res_df_wt_TF$log2FoldChange) > 0.58 & res_df_wt_TF$padj < 1.580078e-06)

volcano_plot <- ggplot(data=res_df_wt_TF, aes(x=log2FoldChange, y=-log10(padj))) +
  geom_point(data=res_df_wt_TF, size=1, colour="gray") +
  geom_point(data=res_df_wt_TF[res_df_wt_TF$super==TRUE, ], size=2, colour="#CC0000") +
  geom_point(data=res_df_wt_TF[res_df_wt_TF$sub  ==TRUE, ], size=2, colour="#000099") +
  geom_text_repel(aes(label = geneLabels)) +
  xlab("log2 fold change") +
  ylab("-log10 p-value adjusted") +
  ggtitle("Differential expression wtAIPL1 vs CC")+
  scale_x_continuous() +
  scale_y_continuous() +
  theme_bw() +
  theme(axis.title.y = element_text(size=16),
        axis.title.x = element_text(size=16, colour="black"),
        axis.text = element_text(size=12),
        legend.title =element_blank() ,
        legend.text = element_text(size = 12)) +
  theme(plot.title = element_text(face="bold", size = 17, hjust = 0.5))
volcano_plot




################################## for wt vs opt ###############################
## Obtain logical vector where TRUE values denote padj values < 0.05 and fold change > 1.5 in either direction
res_df_opt_TF <- res_df_opt %>% 
  mutate(dif_expr = padj < 1.580078e-06 & abs(log2FoldChange) >= 0.58)

## Create a column to indicate which genes to label 
res_df_opt_TF$geneLabels <- data_geneNames$geneLabel

res_df_opt_TF_top20 <- subset(res_df_opt_TF %>% arrange(log2FoldChange) %>% head(n=20))
res_df_opt_TF_low20 <- subset(res_df_opt_TF %>% arrange(desc(log2FoldChange)) %>% head(n=20))

res_df_opt_TF <- res_df_opt_TF[,!names(res_df_opt_TF) %in% c("geneLabels")]

res_df_opt_TF_40 <- rbind(res_df_opt_TF_top20, res_df_opt_TF_low20)
res_df_opt_TF_40 <- res_df_opt_TF_40[,!names(res_df_opt_TF_top20) %in% 
                                       c("baseMean", 'lfcSE', 'stat', 'pvalue', 'padj', 'dif_expr')]
res_df_opt_TF <- merge(res_df_opt_TF, res_df_opt_TF_40, by = 'row.names', all.x = T)
row.names(res_df_opt_TF) <- res_df_opt_TF$Row.names
drop = c('Row.names', 'log2FoldChange.y')
res_df_opt_TF = res_df_opt_TF[,!(names(res_df_opt_TF) %in% drop)]
colnames(res_df_opt_TF)[2] ="log2FoldChange"

res_df_opt_TF$super <- res_df_opt_TF$log2FoldChange >  0.58  & res_df_opt_TF$padj < 1.580078e-06
res_df_opt_TF$sub   <- res_df_opt_TF$log2FoldChange < -0.58  & res_df_opt_TF$padj < 1.580078e-06
res_df_opt_TF$threshold <- as.factor(abs(res_df_opt_TF$log2FoldChange) > 0.58 & res_df_opt_TF$padj < 1.580078e-06)

#PLOT!
volcano_plot <- ggplot(data=res_df_opt_TF, aes(x=log2FoldChange, y=-log10(padj))) +
  geom_point(data=res_df_opt_TF, size=1, colour="gray") +
  geom_point(data=res_df_opt_TF[res_df_opt_TF$super==TRUE, ], size=2, colour="#CC0000") +
  geom_point(data=res_df_opt_TF[res_df_opt_TF$sub  ==TRUE, ], size=2, colour="#000099") +
  geom_text_repel(aes(label = geneLabels)) +
  xlab("log2 fold change") +
  ylab("-log10 p-value adjusted") +
  ggtitle("Differential expression optAIPL1 vs CC") +
  scale_x_continuous() +
  scale_y_continuous() +
  theme_bw() +
  theme(axis.title.y = element_text(size=16),
        axis.title.x = element_text(size=16, colour="black"),
        axis.text = element_text(size=12),
        legend.title =element_blank() ,
        legend.text = element_text(size = 12)) +
  theme(plot.title = element_text(face="bold", size = 17, hjust = 0.5))
volcano_plot





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








#reset par
par(mfrow=c(1,1))
# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot"))
# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))


