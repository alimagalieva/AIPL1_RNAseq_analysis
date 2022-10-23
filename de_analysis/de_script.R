## Gene-level differential expression analysis using DESeq2 

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
df <- read.csv("data/counts.Sirius.csv",sep = ";") #, row.names = 2) 

## Trying to reorder columns by replicates alphabetically:
coln1 <- c(colnames(df[1:4]))
coln2 <- c(colnames(df[5:16]))
coln2 <- sort(coln2)
coln = c(coln1, coln2)
coln
data <- df[, coln]

### Check classes of the data we just brought in
class(data)

pdf("Raw_expression_counts_plots.pdf", width = 10, height = )
ggplot(data) +
  geom_histogram(aes(x = AIPL1_6_1), stat = "bin", bins = 200) +
  xlab("Raw expression counts") +
  ylab("Number of genes")
#jpeg('rplot.jpg')

ggplot(data) +
  geom_histogram(aes(x = AIPL1_6_1), stat = "bin", bins = 200) + 
  xlim(-5, 500)  +
  xlab("Raw expression counts") +
  ylab("Number of genes")
dev.off() 

## So
pdf("mean_counts_plots.pdf")
mean_counts_AIPL3 <- apply(data[, 5:7], 1, mean)
variance_counts <- apply(data[, 5:7], 1, var)
df <- data.frame(mean_counts_AIPL3, variance_counts)

ggplot(df) +
  geom_point(aes(x=mean_counts_AIPL3, y=variance_counts)) + 
  geom_line(aes(x=mean_counts_AIPL3, y=mean_counts_AIPL3, color="red")) +
  scale_y_log10() +
  scale_x_log10()


mean_counts_AIPL16 <- apply(data[, 8:10], 1, mean)
variance_counts <- apply(data[, 8:10], 1, var)
df <- data.frame(mean_counts_AIPL16, variance_counts)

ggplot(df) +
  geom_point(aes(x=mean_counts_AIPL16, y=variance_counts)) + 
  geom_line(aes(x=mean_counts_AIPL16, y=mean_counts_AIPL16, color="red")) +
  scale_y_log10() +
  scale_x_log10()


mean_counts_CC <- apply(data[, 11:13], 1, mean)
variance_counts <- apply(data[, 11:13], 1, var)
df <- data.frame(mean_counts_CC, variance_counts)

ggplot(df) +
  geom_point(aes(x=mean_counts_CC, y=variance_counts)) + 
  geom_line(aes(x=mean_counts_CC, y=mean_counts_CC, color="red")) +
  scale_y_log10() +
  scale_x_log10()


mean_counts_GFP <- apply(data[, 14:16], 1, mean)
variance_counts <- apply(data[, 14:16], 1, var)
df <- data.frame(mean_counts_GFP, variance_counts)

ggplot(df) +
  geom_point(aes(x=mean_counts_GFP, y=variance_counts)) + 
  geom_line(aes(x=mean_counts_GFP, y=mean_counts_GFP, color="red")) +
  scale_y_log10() +
  scale_x_log10()
dev.off() ## Close writing in pdf-file


## some preparations of data dataframe:
drop <- c("Refseq","Gene.name", "Symbol")
data = data[,!(names(data) %in% drop)]
data <- data[complete.cases(data['name']),]
row.names(data) <- data[,1]
drop <- c("name", 'CC1') #remove CC1 from analysis, cause it doesn't cluster well
data = data[,!(names(data) %in% drop)]

meta <- data.frame(colnames(data))
meta['sampletype'] <- c('wt', 'wt', 'wt', 'opt', 'opt', 'opt', 
                        'control', 'control', 'GFP', 'GFP', 'GFP')
row.names(meta) <- meta[,1]
meta <- subset(meta, select = c('sampletype'))

all(colnames(data) %in% rownames(meta))
all(colnames(data) == rownames(meta))


## creating DESeq DataSet
dds <- DESeqDataSetFromMatrix(countData = data, colData = meta, design = ~ sampletype)
## Run DESeq2 differential expression analysis
dds <- DESeq(dds)

##  **Optional step** - Output normalized counts to save as a file to access outside RStudio
normalized_counts <- counts(dds, normalized=TRUE) 
write.table(normalized_counts, file="data/normalized_counts.txt", sep="\t", quote=F, col.names=NA) 

## QC
# Transform counts for data visualization
rld <- rlog(dds, blind=TRUE)
# Plot PCA 
plotPCA(rld, intgroup="sampletype")
# Extract the rlog matrix from the object
rld_mat <- assay(rld)
# Compute pairwise correlation values
rld_cor <- cor(rld_mat)
# Plot heatmap
heat.colors <- brewer.pal(6, "Blues")
pheatmap(rld_cor, color = heat.colors, fontsize = 10, fontsize_row = 10, height=20)

## Shrinkage 
contrast <- c("sampletype", "wt", "control")
contrast <- c("sampletype", "opt", "control")
res_unshrunken_opt <- results(dds, contrast = contrast)
res_unshrunken_wt <- results(dds, contrast = contrast)
res_shrunken_opt <- lfcShrink(dds, contrast = contrast, res=res_unshrunken_opt, type='normal')
res_shrunken_wt <- lfcShrink(dds, contrast = contrast, res=res_unshrunken_wt, type='normal')

## How many genes are differentially expressed compared to control:
padj.cutoff <- 0.05
lfc.cutoff <- 0.58
# Turn the results object into a data frame
res_df_wt <- data.frame(res_shrunken_wt)
res_df_opt <- data.frame(res_shrunken_opt)
# Subset the significant results
sig_res_wt <- filter(res_df_wt, padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff)
sig_res_opt <- filter(res_df_opt, padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff)

#### Plot expression of a single gene 
# ENSG00000129221 - AIPL1 gene 
# Save plotcounts to a data frame object
d <- plotCounts(dds, gene="ENSG00000129221", intgroup="sampletype", returnData=TRUE)
# Plotting AIPL1 normalized counts, using the samplenames (rownames of d as labels)
ggplot(d, aes(x = sampletype, y = count, color = sampletype)) + 
  geom_point(position=position_jitter(w = 0.1,h = 0)) +
  geom_text_repel(aes(label = rownames(d))) + 
  theme_bw() +
  ggtitle("AIPL1") +
  theme(plot.title = element_text(hjust = 0.5))


!!! 
### Using ggplot2 to plot multiple genes (e.g. top 20):
## Order results by padj values
top20_sig_genes <- res %>% 
  arrange(padj) %>% 	#Arrange rows by padj values
  pull(gene) %>% 		#Extract character vector of ordered genes
  head(n=20) 		#Extract the first 20 genes

## Then, we can extract the normalized count values for these top 20 genes:
## normalized counts for top 20 significant genes
top20_sig_norm <- normalized_counts %>%
  filter(gene %in% top20_sig_genes)

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






res <- results(dds)
head(results(dds, tidy=TRUE)) #let's look at the results table
summary(res)
res <- res[order(res$padj),]
head(res)


#reset par
par(mfrow=c(1,1))
# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot"))
# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))


## QC:
pdf("QC_plots.pdf", width = 10, height = )

## PCA:
### Transform counts for data visualization
rld <- rlog(dds, blind=TRUE)
plotPCA(rld, intgroup="sampletype")

## Heatmap:
### Extract the rlog matrix from the object
rld_mat <- assay(rld)    ## assay() is function from the "SummarizedExperiment" package that was loaded when you loaded DESeq2
### Compute pairwise correlation values
rld_cor <- cor(rld_mat)    ## cor() is a base R function
head(rld_cor)   ## check the output of cor(), make note of the rownames and colnames
### Plot heatmap
pheatmap(rld_cor)

dev.off() ## Close writing in pdf-file



## Total number of raw counts per sample
colSums(counts(dds))

## Total number of normalized counts per sample
colSums(counts(dds, normalized=T))

## Plot dispersion estimates
plotDispEsts(dds)

## log2fold change calculation and MAplot:
dds <- DESeq(dds)
plotMA(dds, main = "dds")
res <- results(dds)
plotMA(res, alpha = 0.05, main = "res", ylim=c(-7,7))