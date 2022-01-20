rm(list=ls())

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

library(DESeq2)
library(pheatmap)
library("genefilter")
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library("cowplot")


##read the resulting featureCounts file created in the 'Count.sh' step, change the row names with the genes name and remove the column containing the genes name

cts <- read.table("C:/Users/Maëlle/Documents/Master/RNA-seq/Project/cut_read_count", h=T, sep='\t')
row.names(cts) <- cts$Geneid
cts <- as.matrix(cts[,-1])


## Create a dataframe with the condition (Normal or HER2-positive) and the type (paired-end) for each sample

condition <- as.factor(c(rep("HER2",3),rep("Normal",3)))
type <- as.factor(c(rep("Paired_end",6)))
coldata <-  data.frame(row.names = c("HER21","HER22","HER23", "Normal1","Normal2","Normal3"),condition, type)


## Create a DESeqDataSet object from the featureCounts arranged matrix and the condition/type dataframe. Separate (design) the results by condition (HER2 or Normal)
## and do a differential expression analysis with DESeq()

dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design=~condition)
dds <- DESeq(dds)


## Do the mean of each row from dds with normalized counts and take the 20 with the bigger means.
## Estimate the dispersion trend and plot it as an heatmap for the genes selected above in 'select' by samples. 
## This allows to check quality by seeing if the samples of same group have the same pattern.

select <- order(rowMeans(counts(dds,normalized=TRUE)), decreasing=TRUE)[1:20]
vst <- vst(dds, blind = T) # blind the transformation to the experimental design
pheatmap(assay(vst)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=coldata)


## Create a Principal component analysis plot with the dispersion data (vst)

plotPCA(vst)


## extract the results from the differential analysis with a significance level of 5% and show the summary of the results (how many up/down regulated)

res <- results(dds, contrast=c("condition","HER2","Normal"),alpha=0.05)
summary(res)


## create a MA plot with the results 

plotMA(res)


## extract the 20 more differentially expressed genes from the analysis and plot their dispersion on an heatmap
## This also permitted to find the most differentially expressed genes for the next step

topVarGenes <- head(order(rowVars(assay(vst)), decreasing = TRUE), 20)
mat  <- assay(vst)[ topVarGenes, ]
mat  <- mat - rowMeans(mat)
pheatmap(mat, annotation_col = coldata)


## plot 3 graphs with the counts per sample of 3 genes of interest separated by the condition 

par(mfrow=c(1,3))
plotCounts(dds, gene='ENSG00000108342', intgroup="condition",main = 'CSF3')
plotCounts(dds, gene='ENSG00000141736', intgroup="condition",main = 'ERBB2')
plotCounts(dds, gene='ENSG00000135446', intgroup="condition",main = 'CDK4')


## Create a vector with the name of all the genes with a significant padj and omit all the the NA
## Do the same with the log2FoldChange values
## Create a dataframe with the two vectors
## Create a vector with all the genes name

gene_DE <- rownames(res)[res$padj < 0.05]
gene_DE <- na.omit(gene_DE)
log2 <- res$log2FoldChange[res$padj < 0.05]
log2 <- na.omit(log2)
df <- data.frame(gene_DE,log2)
gene <- rownames(res)


## Create two vectors with the name of the genes with a significant padj(gene_DE) that have a log2FoldChange greater/lower than 0 to keep only the up/down regulated genes

up_regulated <- df$gene_DE[log2 > 0]
down_regulated <- df$gene_DE[log2 < 0]


## Do 3 overexpression analysis with first all the genes differentially expressed, then the differentially expressed genes up-regulated and finally the differentially expressed genes down-regulated.

ego_BP <- enrichGO(gene       = gene_DE, #gene differentially expressed
                universe      = gene, #all the genes in our analysis
                OrgDb         = org.Hs.eg.db, #genome wide annotation 
                qvalueCutoff  = 0.05, 
                ont           = 'BP', #filter results for the Bioprocesses ontology 
                keyType       = 'ENSEMBL') #the key of the genes are from ensembl (ENSG...)
ego_up <- enrichGO(gene       = up_regulated,
                   universe      = gene,
                   OrgDb         = org.Hs.eg.db,
                   qvalueCutoff  = 0.05,
                   ont           = 'BP',
                   keyType       = 'ENSEMBL')

ego_down <- enrichGO(gene       = down_regulated,
                   universe      = gene,
                   OrgDb         = org.Hs.eg.db,
                   qvalueCutoff  = 0.05,
                   ont           = 'BP',
                   keyType       = 'ENSEMBL')


## Create a new dataframe with the 10 first result of enrichGO on all genes ordered by counts (decreasing)
## Create one dotplot with the 10 most frequent categorys (derived from most_enriched dataframe) and their counts from the overexpression analysis of all the differentially expressed genes
## Create 2 barplot with the 10 most frequent ontology and their counts for the overexpression analysis of the up/down regulated genes 

most_enriched = head(ego_BP@result[order(ego_BP@result$Count, decreasing = T),],10)
A <- dotplot(ego_BP, showCategory=most_enriched$Description, label_format=20)
B <- barplot(ego_up, showCategory=10,label_format=15, order='T') 
C <- barplot(ego_down, showCategory=10,label_format=15, order='T') 


## Draw the 3 plots from above in one plot.

ggdraw() +
  draw_plot(A, x = 0, y = 0, width = 0.5, height = 1) +
  draw_plot(B, x = .5, y = .5, width = .5, height = .5) +
  draw_plot(C, x = .5, y = 0, width = .5, height = .5) +
  draw_plot_label(label = c("A", "B", "C"), size = 25,
                  x = c(0, 0.5, 0.5), y = c(1, 1, 0.5))+
  draw_plot_label(label = c('Up-regulated','Down-regulated'),size = 20,x = c(0.55,0.52),y=c(1.005,0.505))
