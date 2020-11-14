
####Salmon/Tximport####

library(Leonardi4.0)

library(GenomicFeatures)
library(tximport)
library(readr)
library(rjson)
library(dplyr)
library(AnnotationDbi)

library(tximportData)
dir <- system.file("extdata", package = "Leonardi4.0", mustWork=TRUE)
list.files(dir)
list.files(file.path(dir, "quants"))

csvfile <- file.path(dir, "sample_table_Mouse_Leonardi_F.csv")
coldata <- read.csv(csvfile, row.names=1, stringsAsFactors=FALSE)
coldata

coldata <- coldata[1:16,]
coldata$names <- coldata$Run
coldata$files <- file.path(dir, "quants", coldata$names, "quant.sf")
file.exists(coldata$files)


library(tximeta)
se <- tximeta(coldata)
dim(se)
head(rownames(se))
gse <- summarizeToGene(se)
dim(gse)

library(SummarizedExperiment)
data(gse)

gse
assayNames(gse)
head(assay(gse), 3)
colSums(assay(gse))
rowRanges(gse)
seqinfo(rowRanges(gse))
colData(gse)

#files <- file.path(dir, "quants", coldata$Run, "quant.sf")
#names(files) <- paste0("sample", 1:32)
#all(file.exists(files))

#library(EnsDb.Mmusculus.v79)
#TxDb <- EnsDb.Mmusculus.v79
#k <- keys(TxDb, keytype = "TXNAME")
#tx2gene <- select(TxDb, k, "GENEID", "TXNAME")
#head(tx2gene)

#txi <- tximport(file, design = ~ genotype + diet + genotype:diet)es, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
#names(txi)

#head(txi$counts)
#txi.tx <- tximport(files, type = "salmon", txOut = TRUE, ignoreTxVersion = TRUE)
#txi.sum <- summarizeToGene(txi.tx, tx2gene, ignoreTxVersion = TRUE)
#all.equal(txi$counts, txi.sum$counts)

####Make Counts Table####

#assay(gse)
#resOrderedDF <- as.data.frame(assay(gse))
#write.csv(resOrderedDF, file = "/home/john/Documents/Bioinformatics_Seq_Files/Fastq/Leonardi_mRNA/Output/Counts_Table_Male_mRNA_gene.csv")

####DESeq2####

#library(magrittr)
#levels(gse$genotype)
#gse$genotype <- relevel(gse$genotype, "WT")
#gse$genotype %<>% relevel("WT")
#gse$diet

round( colSums(assay(gse)) / 1e6, 1 )


library(DESeq2)

#Use group for general comparison(1st), use interaction term for interactions (2nd), pairwise comparisions (3rd) 

#dds <- DESeqDataSet(gse, design = ~ genotype + diet + sex)
dds <- DESeqDataSet(gse, design = ~ genotype + diet + genotype:diet)
#dds <- DESeqDataSet(gse, design = ~ group)

nrow(dds)

keep <- rowSums(counts(dds) >= 10) >= 4
dds <- dds[keep,]
nrow(dds)

dds$genotype <- relevel(dds$genotype, ref = "WT")
dds$genotype
levels(dds$genotype)

#rownames(coldata) <- colnames(txi$counts)

#dds <- DESeqDataSetFromTximport(txi, samples, ~dex)

#nrow(dds)
#keep <- rowSums(counts(dds)) > 1
#dds <- dds[keep,]
#nrow(dds)

#keep <- rowSums(counts(dds) >= 5) >= 3


##Data Normalization Assessment##
vsd <- vst(dds, blind = FALSE)
head(assay(vsd), 3)

colData(vsd)

rld <- rlog(dds, blind = FALSE)
head(assay(rld), 3)

library("dplyr")
library("ggplot2")

dds <- estimateSizeFactors(dds)

df <- bind_rows(
  as_tibble(log2(counts(dds, normalized=TRUE)[, 1:2]+1)) %>%
    mutate(transformation = "log2(x + 1)"),
  as_tibble(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"),
  as_tibble(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"))

#df <- bind_rows(
#  as_data_frame(log2(counts(dds, normalized=TRUE)[, 1:2]+1)) %>%
#    mutate(transformation = "log2(x + 1)"),
#  as_data_frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"),
#  as_data_frame(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"))

colnames(df)[1:2] <- c("x", "y")  

lvls <- c("log2(x + 1)", "vst", "rlog")
df$transformation <- factor(df$transformation, levels=lvls)

ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid( . ~ transformation)


##Sample Distances##
sampleDists <- dist(t(assay(vsd)))
sampleDists

library("pheatmap")
library("RColorBrewer")
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( vsd$genotype, vsd$diet, sep = " - " )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)


library("PoiClaClu")
poisd <- PoissonDistance(t(counts(dds)))

samplePoisDistMatrix <- as.matrix( poisd$dd )
rownames(samplePoisDistMatrix) <- paste( dds$genotype, dds$diet, sep=" - " )
colnames(samplePoisDistMatrix) <- NULL
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows = poisd$dd,
         clustering_distance_cols = poisd$dd,
         col = colors)


##PCA Plot##
#genotype and diet#
plotPCA(vsd, intgroup = c("genotype", "diet"))
pcaData <- plotPCA(vsd, intgroup = c("genotype", "diet"), returnData = TRUE)
pcaData
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2, color = genotype, shape = diet)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  ggtitle("PCA with VST data")

#genotype, diet, and sex#
plotPCA(vsd, intgroup = c("genotype", "diet", "sex"))
pcaData <- plotPCA(vsd, intgroup = c("genotype", "diet", "sex"), returnData = TRUE)
pcaData
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2, color = genotype, shape = diet, size = sex)) + geom_point() +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  ggtitle("PCA with VST data")


library("glmpca")
gpca <- glmpca(counts(dds), L=2)
gpca.dat <- gpca$factors
gpca.dat$genotype <- dds$genotype
gpca.dat$diet <- dds$diet

ggplot(gpca.dat, aes(x = dim1, y = dim2, color = genotype, shape = diet)) +
  geom_point(size =3) + coord_fixed() + ggtitle("glmpca - Generalized PCA")


##MDS Plot##
mds <- as.data.frame(colData(vsd))  %>%
  cbind(cmdscale(sampleDistMatrix))
ggplot(mds, aes(x = `1`, y = `2`, color = genotype, shape = diet)) +
  geom_point(size = 3) + coord_fixed() + ggtitle("MDS with VST data")

mdsPois <- as.data.frame(colData(dds)) %>%
  cbind(cmdscale(samplePoisDistMatrix))
ggplot(mdsPois, aes(x = `1`, y = `2`, color = genotype, shape = diet)) +
  geom_point(size = 3) + coord_fixed() + ggtitle("MDS with PoissonDistances")


##Differential Expression Analysis##
#library(ashr)
dds <- DESeq(dds)
resultsNames(dds)


#Middle group is "treatment", end group is "control"
##General##
res01 <- results(dds, contrast=c("genotype","KO","WT"), alpha = 0.05)
mcols(res01, use.names = TRUE)
summary(res01)
table(res01$padj < 0.05)

res02 <- results(dds, contrast=c("diet","WD","CD"), alpha = 0.05)
mcols(res02, use.names = TRUE)
summary(res02)
table(res02$padj < 0.05)

res03 <- results(dds, contrast=c("sex","Female","Male"), alpha = 0.05)
mcols(res03, use.names = TRUE)
summary(res03)
table(res03$padj < 0.05)

##More Specific##
#Interaction Sets#
#Do individually for male and female
resWT <- results(dds, contrast=c("diet", "WD", "CD"), alpha = 0.05)
resKO <- results(dds, list(c("diet_WD_vs_CD", "genotypeKO.dietWD")), alpha = 0.05)
resI <- results(dds, name="genotypeKO.dietWD", alpha = 0.05)

#resWT2 <- lfcShrink(dds, res = resWT, type="ashr")
#resKO2 <- lfcShrink(dds, res = resKO, type="ashr")
#resI2 <- lfcShrink(dds, res = resI, type="ashr")

ixWT = which.min(resWT$padj) # most significaggplot(pcaData, aes(x = PC1, y = PC2, color = genotype, shape = diet)) +
geom_point(size=3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  nt
barplot(assay(dds)[ixWT,],las=2,names=coldata$Run,cex.names=0.9, main=rownames(dds)[ ixWT  ]  )

ixKO = which.min(resKO$padj) # most significant
barplot(assay(dds)[ixKO,],las=2,names=coldata$Run,cex.names=0.9, main=rownames(dds)[ ixKO  ]  )

ixI = which.min(resI$padj) # most significant
barplot(assay(dds)[ixI,],las=2,names=coldata$Run,cex.names=0.9, main=rownames(dds)[ ixI  ]  )
#barplot(assay(dds)[ixI,],las=2,names=coldata$Run,cex.names=0.9,beside=TRUE,col=rep(c("black","white"),each=8), main=rownames(dds)[ ixI  ]  )

##Most Specific##
res1 <- results(dds, contrast=c("group","KOFWD","WTFWD"), alpha = 0.05)
mcols(res1, use.names = TRUE)
summary(res1)
table(res1$padj < 0.05)

res2 <- results(dds, contrast=c("group","KOMWD","WTMWD"), alpha = 0.05)
mcols(res2, use.names = TRUE)
summary(res2)
table(res2$padj < 0.05)

res3 <- results(dds, contrast=c("group","KOFCD","WTFCD"), alpha = 0.05)
mcols(res3, use.names = TRUE)
summary(res3)
table(res3$padj < 0.05)

res4 <- results(dds, contrast=c("group","KOMCD","WTMCD"), alpha = 0.05)
mcols(res4, use.names = TRUE)
summary(res4)
table(res4$padj < 0.05)

#res5 <- results(dds, contrast=c("group","KOFWD","KOMWD"), alpha = 0.05)
#mcols(res5, use.names = TRUE)
#summary(res5)
#table(res5$padj < 0.05)

#res6 <- results(dds, contrast=c("group","WTFWD","WTMWD"), alpha = 0.05)
#mcols(res6, use.names = TRUE)
#summary(res6)
#table(res6$padj < 0.05)

#res7 <- results(dds, contrast=c("group","KOFCD","KOMCD"), alpha = 0.05)
#mcols(res7, use.names = TRUE)
#summary(res7)
#table(res7$padj < 0.05)

#res8 <- results(dds, contrast=c("group","WTFCD","WTMCD"), alpha = 0.05)
#mcols(res8, use.names = TRUE)
#summary(res8)
#table(res8$padj < 0.05)



##VennDiagram##
library(limma)

res1 <- results(dds, contrast=c("group","KOFWD","WTFWD"), alpha = 0.05)
res1 <- res1[which(res1$padj < 0.05),]
res1.genes <- row.names(res1)

res2 <- results(dds, contrast=c("group","KOMWD","WTMWD"), alpha = 0.05)
res2 <- res2[which(res2$padj < 0.05),]
res2.genes <- row.names(res2)

res3 <- results(dds, contrast=c("group","KOFCD","WTFCD"), alpha = 0.05)
res3 <- res3[which(res3$padj < 0.05),]
res3.genes <- row.names(res3)

res4 <- results(dds, contrast=c("group","KOMCD","WTMCD"), alpha = 0.05)
res4 <- res4[which(res4$padj < 0.05),]
res4.genes <- row.names(res4)

#res5 <- results(dds, contrast=c("group","KOFWD","KOMWD"), alpha = 0.05)
#res5 <- res5[which(res5$padj < 0.05),]
#res5.genes <- row.names(res5)

#res6 <- results(dds, contrast=c("group","WTFWD","WTMWD"), alpha = 0.05)
#res6 <- res6[which(res6$padj < 0.05),]
#res6.genes <- row.names(res6)

#res7 <- results(dds, contrast=c("group","KOFCD","KOMCD"), alpha = 0.05)
#res7 <- res7[which(res7$padj < 0.05),]
#res7.genes <- row.names(res7)

#res8 <- results(dds, contrast=c("group","WTFCD","WTMCD"), alpha = 0.05)
#res8 <- res8[which(res8$padj < 0.05),]
#res8.genes <- row.names(res8)

#Unique <- sort(unique(c(res1.genes, res2.genes)))
#Unique <- sort(unique(c(res1.genes, res2.genes, res3.genes)))
Unique <- sort(unique(c(res1.genes, res2.genes, res3.genes, res4.genes)))
#Unique <- sort(unique(c(res5.genes, res6.genes, res7.genes, res8.genes)))
#Unique <- sort(unique(c(res1.genes, res2.genes, res3.genes, res4.genes, res5.genes, res6.genes, res7.genes, res8.genes)))


res1.genes.2 <- Unique %in% res1.genes
res2.genes.2 <- Unique %in% res2.genes
res3.genes.2 <- Unique %in% res3.genes
res4.genes.2 <- Unique %in% res4.genes
#res5.genes.2 <- Unique %in% res5.genes
#res6.genes.2 <- Unique %in% res6.genes
#res7.genes.2 <- Unique %in% res7.genes
#res8.genes.2 <- Unique %in% res8.genes


#counts.1 <- cbind(res1.genes.2,res2.genes.2)
#counts.1 <- cbind(res1.genes.2,res2.genes.2,res3.genes.2)
counts.1 <- cbind(res1.genes.2,res2.genes.2,res3.genes.2,res4.genes.2)
#counts.1 <- cbind(res5.genes.2,res6.genes.2,res7.genes.2,res8.genes.2)
#counts.1 <- cbind(res1.genes.2,res2.genes.2,res3.genes.2,res4.genes.2,res5.genes.2,res6.genes.2,res7.genes.2,res8.genes.2)

results.1 <- vennCounts(counts.1, include="both")

#vennDiagram(results.1, include="both", cex = 1, names = c("",""), circle.col = c("blue", "red"))
#vennDiagram(results.1, include="both", cex = 1, names = c("","",""), circle.col = c("blue", "red", "green"))
vennDiagram(results.1, include="both", cex = 1, names = c("","","",""), circle.col = c("blue", "red", "green", "black"))

#library(UpSetR)
#UpSetHisto <- as.data.frame(array(as.numeric(unlist(results.1)), dim=c(256, 9)))
#upset(UpSetHisto, sets = c("V1", "V2", "V3", "V4", "V5", 
#                           "V6", "V7", "V8"), order.by = "freq", nsets = 8)
#upset(UpSetHisto)


##Other Visualizations##
library(vidger)

vsBoxPlot(dds, d.factor = "group", type = "deseq")
vsDEGMatrix(dds, d.factor = "group", padj = 0.05, type = "deseq")
#vsFourWay("trt2", "trt", "untrt2", d.factor = "group", dds, type = "deseq")
vsMAMatrix(dds, d.factor = "group", type = "deseq", y.lim = c(-10,10))
vsScatterMatrix(dds, d.factor = "group", type = "deseq")
#vsScatterPlot("untrt2", "trt2", dds, d.factor = "group", type = "deseq")
#vsVolcano("untrt2", "trt2", dds, d.factor = "group", type = "deseq", x.lim = c(-10,10), padj = 0.05)
vsVolcanoMatrix(dds, d.factor = "group", type = "deseq", lfc = 2, padj = 0.05, x.lim = c(-8,8),
                title = FALSE, legend = TRUE, grid = TRUE, counts = FALSE, facet.title.size = 10)


##Volcano Plots##
library(EnhancedVolcano)

#write.csv(Unique, file = "/home/john/Documents/Bioinformatics_Seq_Files/Fastq/Salik_mRNA/Salmon_1.1.0_Output/Salik_mRNA_unique_1_trt6.csv")
#write.csv(counts.1, file = "/home/john/Documents/Bioinformatics_Seq_Files/Fastq/Salik_mRNA/Salmon_1.1.0_Output/Salik_mRNA_unique_2_trt6.csv")

#resUNIQUE <- read.csv("/home/john/Documents/Bioinformatics_Seq_Files/Fastq/Salik_mRNA/Salmon_1.1.0_Output/Salik_mRNA_unique_3_trt6.csv", row.names = 1)

#EnhancedVolcano(resUNIQUE, lab = rownames(resUNIQUE), x = "log2FoldChange", y = "padj", pCutoff = 0.05,
#                xlim = c(-10, 8),  ylim = c(0, -log10(10e-11)))

EnhancedVolcano(res01, lab = rownames(res01), x = "log2FoldChange", y = "padj", pCutoff = 0.05,
                xlim = c(-7.5, 7.5))

EnhancedVolcano(res02, lab = rownames(res02), x = "log2FoldChange", y = "padj", pCutoff = 0.05,
                xlim = c(-7.5, 7.5))

EnhancedVolcano(res3, lab = rownames(res3), x = "log2FoldChange", y = "padj", pCutoff = 0.05,
                xlim = c(-7.5, 7.5))

EnhancedVolcano(res4, lab = rownames(res4), x = "log2FoldChange", y = "padj", pCutoff = 0.05,
                xlim = c(-7.5, 7.5))




##Plotting Results##
topGene <- rownames(res01)[which.min(res01$padj)]
plotCounts(dds, gene = topGene, intgroup=c("genotype"))

#GeneofChoice <- "ENSMUSG00000092341"
#plotCounts(dds, gene = GeneofChoice, intgroup=c("genotype"))


library("ggbeeswarm")
geneCounts <- plotCounts(dds, gene = topGene, intgroup = c("genotype","diet"),
                         returnData = TRUE)
ggplot(geneCounts, aes(x = genotype, y = count, color = diet)) +
  scale_y_log10() +  geom_beeswarm(cex = 3)

ggplot(geneCounts, aes(x = genotype, y = count, color = diet, group = diet)) +
  scale_y_log10() + geom_point(size = 3) + geom_line()


##MA-Plot##
library("apeglm")
resultsNames(dds)

#res <- lfcShrink(dds, coef="genotype_WT_vs_KO", type="apeglm")
#plotMA(res01, ylim = c(-5, 5))
#res.noshr <- results(dds, name="genotype_WT_vs_KO")
#plotMA(res.noshr, ylim = c(-5, 5))

plotMA(res01, ylim = c(-5,5))
topGene <- rownames(res01)[which.min(res01$padj)]
with(res01[topGene, ], {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, topGene, pos=2, col="dodgerblue")
})


plotMA(res02, ylim = c(-5,5))
topGene <- rownames(res02)[which.min(res02$padj)]
with(res02[topGene, ], {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, topGene, pos=2, col="dodgerblue")
})


plotMA(res03, ylim = c(-10,10))
topGene <- rownames(res03)[which.min(res03$padj)]
with(res03[topGene, ], {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, topGene, pos=2, col="dodgerblue")
})

##Histogram of p vaules for genes with mean normalized count larger than 2##
hist(res1$pvalue[res$baseMean > 1], breaks = 0:20/20,
     col = "grey50", border = "white")

##Gene Clustering##

library("genefilter")
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 1000)
#or#
#topGenes <- head(order(res01$padj),decreasing = TRUE, 50)

mat  <- assay(vsd)[ topVarGenes, ]
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vsd)[, c("diet","genotype","sex")])
pheatmap(mat, annotation_col = anno)

#Annotating and Exporting
library("AnnotationDbi")
library("org.Mm.eg.db")

columns(org.Mm.eg.db)

ens.str01 <- rownames(res01)
res01$symbol <- mapIds(org.Mm.eg.db,
                      keys=ens.str01,
                      column="SYMBOL",
                      keytype="ENSEMBL",
                      multiVals="first")
res01$entrez <- mapIds(org.Mm.eg.db,
                      keys=ens.str01,
                      column="ENTREZID",
                      keytype="ENSEMBL",
                      multiVals="first")

resOrdered01 <- res01[order(res01$padj),]
head(resOrdered01)

ens.str02 <- rownames(res02)
res02$symbol <- mapIds(org.Mm.eg.db,
                      keys=ens.str02,
                      column="SYMBOL",
                      keytype="ENSEMBL",
                      multiVals="first")
res02$entrez <- mapIds(org.Mm.eg.db,
                      keys=ens.str02,
                      column="ENTREZID",
                      keytype="ENSEMBL",
                      multiVals="first")

resOrdered02 <- res02[order(res02$padj),]
head(resOrdered02)

ens.str03 <- rownames(res03)
res03$symbol <- mapIds(org.Mm.eg.db,
                      keys=ens.str03,
                      column="SYMBOL",
                      keytype="ENSEMBL",
                      multiVals="first")
res03$entrez <- mapIds(org.Mm.eg.db,
                      keys=ens.str03,
                      column="ENTREZID",
                      keytype="ENSEMBL",
                      multiVals="first")

resOrdered03 <- res03[order(res03$padj),]
head(resOrdered03)

ens.strWT <- rownames(resWT)
resWT$symbol <- mapIds(org.Mm.eg.db,
                       keys=ens.strWT,
                       column="SYMBOL",
                       keytype="ENSEMBL",
                       multiVals="first")
resWT$entrez <- mapIds(org.Mm.eg.db,
                       keys=ens.strWT,
                       column="ENTREZID",
                       keytype="ENSEMBL",
                       multiVals="first")

resOrderedWT <- resWT[order(resWT$padj),]
head(resOrderedWT)

ens.strKO <- rownames(resKO)
resKO$symbol <- mapIds(org.Mm.eg.db,
                       keys=ens.strKO,
                       column="SYMBOL",
                       keytype="ENSEMBL",
                       multiVals="first")
resKO$entrez <- mapIds(org.Mm.eg.db,
                       keys=ens.strKO,
                       column="ENTREZID",
                       keytype="ENSEMBL",
                       multiVals="first")

resOrderedKO <- resKO[order(resKO$padj),]
head(resOrderedKO)

ens.strI <- rownames(resI)
resI$symbol <- mapIds(org.Mm.eg.db,
                      keys=ens.strI,
                      column="SYMBOL",
                      keytype="ENSEMBL",
                      multiVals="first")
resI$entrez <- mapIds(org.Mm.eg.db,
                      keys=ens.strI,
                      column="ENTREZID",
                      keytype="ENSEMBL",
                      multiVals="first")

resOrderedI <- resI[order(resI$padj),]
head(resOrderedI)


ens.str1 <- rownames(res1)
res1$symbol <- mapIds(org.Mm.eg.db,
                     keys=ens.str1,
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
res1$entrez <- mapIds(org.Mm.eg.db,
                     keys=ens.str1,
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")

resOrdered1 <- res1[order(res1$padj),]
head(resOrdered1)

ens.str2 <- rownames(res2)
res2$symbol <- mapIds(org.Mm.eg.db,
                      keys=ens.str2,
                      column="SYMBOL",
                      keytype="ENSEMBL",
                      multiVals="first")
res2$entrez <- mapIds(org.Mm.eg.db,
                      keys=ens.str2,
                      column="ENTREZID",
                      keytype="ENSEMBL",
                      multiVals="first")

resOrdered2 <- res2[order(res2$padj),]
head(resOrdered2)

ens.str3 <- rownames(res3)
res3$symbol <- mapIds(org.Mm.eg.db,
                      keys=ens.str3,
                      column="SYMBOL",
                      keytype="ENSEMBL",
                      multiVals="first")
res3$entrez <- mapIds(org.Mm.eg.db,
                      keys=ens.str3,
                      column="ENTREZID",
                      keytype="ENSEMBL",
                      multiVals="first")

resOrdered3 <- res3[order(res3$padj),]
head(resOrdered3)

ens.str4 <- rownames(res4)
res4$symbol <- mapIds(org.Mm.eg.db,
                      keys=ens.str4,
                      column="SYMBOL",
                      keytype="ENSEMBL",
                      multiVals="first")
res4$entrez <- mapIds(org.Mm.eg.db,
                      keys=ens.str4,
                      column="ENTREZID",
                      keytype="ENSEMBL",
                      multiVals="first")

resOrdered4 <- res4[order(res4$padj),]
head(resOrdered4)

ens.str5 <- rownames(res5)
res5$symbol <- mapIds(org.Mm.eg.db,
                      keys=ens.str5,
                      column="SYMBOL",
                      keytype="ENSEMBL",
                      multiVals="first")
res5$entrez <- mapIds(org.Mm.eg.db,
                      keys=ens.str5,
                      column="ENTREZID",
                      keytype="ENSEMBL",
                      multiVals="first")

resOrdered5 <- res5[order(res5$padj),]
head(resOrdered5)

ens.str6 <- rownames(res6)
res6$symbol <- mapIds(org.Mm.eg.db,
                      keys=ens.str6,
                      column="SYMBOL",
                      keytype="ENSEMBL",
                      multiVals="first")
res6$entrez <- mapIds(org.Mm.eg.db,
                      keys=ens.str6,
                      column="ENTREZID",
                      keytype="ENSEMBL",
                      multiVals="first")

resOrdered6 <- res6[order(res6$padj),]
head(resOrdered6)

ens.str7 <- rownames(res7)
res7$symbol <- mapIds(org.Mm.eg.db,
                      keys=ens.str7,
                      column="SYMBOL",
                      keytype="ENSEMBL",
                      multiVals="first")
res7$entrez <- mapIds(org.Mm.eg.db,
                      keys=ens.str7,
                      column="ENTREZID",
                      keytype="ENSEMBL",
                      multiVals="first")

resOrdered7 <- res7[order(res7$padj),]
head(resOrdered7)

ens.str8 <- rownames(res8)
res8$symbol <- mapIds(org.Mm.eg.db,
                      keys=ens.str8,
                      column="SYMBOL",
                      keytype="ENSEMBL",
                      multiVals="first")
res8$entrez <- mapIds(org.Mm.eg.db,
                      keys=ens.str8,
                      column="ENTREZID",
                      keytype="ENSEMBL",
                      multiVals="first")

resOrdered8 <- res8[order(res8$padj),]
head(resOrdered8)

##Exporting Results##

resOrderedDF <- as.data.frame(resOrdered01)[1:15000, ]
write.csv(resOrderedDF, file = "/home/john/Documents/Bioinformatics_Seq_Files/Fastq/Leonardi_mRNA/Output/Final/Female_Genotype.csv")

resOrderedDF <- as.data.frame(resOrdered02)[1:15000, ]
write.csv(resOrderedDF, file = "/home/john/Documents/Bioinformatics_Seq_Files/Fastq/Leonardi_mRNA/Output/Final/Female_Diet.csv")

resOrderedDF <- as.data.frame(resOrdered03)[1:15000, ]
write.csv(resOrderedDF, file = "/home/john/Documents/Bioinformatics_Seq_Files/Fastq/Leonardi_mRNA/Output/Final/Sex.csv")

resOrderedDF <- as.data.frame(resOrderedWT)[1:15000, ]
write.csv(resOrderedDF, file = "/home/john/Documents/Bioinformatics_Seq_Files/Fastq/Leonardi_mRNA/Output/Final/Diet_Effect_on_WT_M.csv")

resOrderedDF <- as.data.frame(resOrderedKO)[1:15000, ]
write.csv(resOrderedDF, file = "/home/john/Documents/Bioinformatics_Seq_Files/Fastq/Leonardi_mRNA/Output/Final/Diet_Effect_on_KO_M.csv")

resOrderedDF <- as.data.frame(resOrderedI)[1:15000, ]
write.csv(resOrderedDF, file = "/home/john/Documents/Bioinformatics_Seq_Files/Fastq/Leonardi_mRNA/Output/Final/Interaction_M.csv")

resOrderedDF <- as.data.frame(resOrdered1)[1:15000, ]
write.csv(resOrderedDF, file = "/home/john/Documents/Bioinformatics_Seq_Files/Fastq/Leonardi_mRNA/Output/Final/KOFWDvsWTFWD.csv")

resOrderedDF <- as.data.frame(resOrdered2)[1:15000, ]
write.csv(resOrderedDF, file = "/home/john/Documents/Bioinformatics_Seq_Files/Fastq/Leonardi_mRNA/Output/Final/KOMWDvsWTMWD.csv")

resOrderedDF <- as.data.frame(resOrdered3)[1:15000, ]
write.csv(resOrderedDF, file = "/home/john/Documents/Bioinformatics_Seq_Files/Fastq/Leonardi_mRNA/Output/Final/KOFCDvsWTFCD.csv")

resOrderedDF <- as.data.frame(resOrdered4)[1:15000, ]
write.csv(resOrderedDF, file = "/home/john/Documents/Bioinformatics_Seq_Files/Fastq/Leonardi_mRNA/Output/Final/KOMCDvsWTMCD.csv")

#resOrderedDF <- as.data.frame(resOrdered5)[1:15000, ]
#write.csv(resOrderedDF, file = "/home/john/Documents/Bioinformatics_Seq_Files/Fastq/Leonardi_mRNA/Output/KOFWDvsKOMWD.csv")

#resOrderedDF <- as.data.frame(resOrdered6)[1:15000, ]
#write.csv(resOrderedDF, file = "/home/john/Documents/Bioinformatics_Seq_Files/Fastq/Leonardi_mRNA/Output/WTFWDvsWTMWD.csv")

#resOrderedDF <- as.data.frame(resOrdered7)[1:15000, ]
#write.csv(resOrderedDF, file = "/home/john/Documents/Bioinformatics_Seq_Files/Fastq/Leonardi_mRNA/Output/KOFCDvsKOMCD.csv")

#resOrderedDF <- as.data.frame(resOrdered7LfcS)[1:15000, ]
#write.csv(resOrderedDF, file = "/home/john/Documents/Bioinformatics_Seq_Files/Fastq/Leonardi_mRNA/Output/KOFCDvsKOMCD_Shrink.csv")

#resOrderedDF <- as.data.frame(resOrdered8)[1:15000, ]
#write.csv(resOrderedDF, file = "/home/john/Documents/Bioinformatics_Seq_Files/Fastq/Leonardi_mRNA/Output/WTFCDvsWTMCD.csv")

##Genomic Space##
resGR <- lfcShrink(dds, coef="genotype_WT_vs_KO", type="apeglm", format="GRanges")
resGR

ens.str <- rownames(res1)
resGR$symbol <- mapIds(org.Mm.eg.db, ens.str, "SYMBOL", "ENSEMBL")

library("Gviz")

window <- resGR[topGene] + 1e6
strand(window) <- "*"
resGRsub <- resGR[resGR %over% window]
naOrDup <- is.na(resGRsub$symbol) | duplicated(resGRsub$symbol)
resGRsub$group <- ifelse(naOrDup, names(resGRsub), resGRsub$symbol)

status <- factor(ifelse(resGRsub$padj < 0.05 & !is.na(resGRsub$padj),
                        "sig", "notsig"))

options(ucscChromosomeNames = FALSE)
g <- GenomeAxisTrack()
a <- AnnotationTrack(resGRsub, name = "gene ranges", feature = status)
d <- DataTrack(resGRsub, data = "log2FoldChange", baseline = 0,
               type = "h", name = "log2 fold change", strand = "+")
plotTracks(list(g, d, a), groupAnnotation = "group",
           notsig = "grey", sig = "hotpink")



##Gene Ontology##
library(topGO)
library(KEGG.db)

#resUNIQUE2 <- read.csv("/home/john/Documents/Bioinformatics_Seq_Files/Fastq/Salik_mRNA/Salmon_1.1.0_Output/Salik_mRNA_unique_4_trt6.csv", row.names = 1)

rowsum.threshold <- 1 # user chosen
fdr.threshold <- 0.05 # user chosen
rs <- rowSums(counts(dds))
dds <- dds[ rs > rowsum.threshold ,]
dds <- DESeq(dds)
res <- results(dds, contrast=c("group","KOFWD","KOFCD"), independentFiltering=FALSE, alpha = 0.05) # use count threshold instead of IF
assayed.genes <- rownames(res)
de.genes <- rownames(res)[ which(res$padj < fdr.threshold) ]
gene.vector=as.integer(assayed.genes%in%de.genes)
names(gene.vector)=assayed.genes
head(gene.vector)

library(goseq)
pwf=nullp(gene.vector,"mm9","ensGene")
head(pwf)

GO.wall=goseq(pwf,"mm9","ensGene")
head(GO.wall)

#Random Re-sampling
#GO.samp=goseq(pwf,"mm9","ensGene",method="Sampling",repcnt=1000)
#head(GO.samp)

write.csv(GO.wall, file = "/home/john/Documents/Bioinformatics_Seq_Files/Fastq/Salik_mRNA/Salmon_1.1.0_Output/Salik_GO_Unique_CB-O3_Day4.csv")

KEGG=goseq(pwf,'mm9','ensGene',test.cats="KEGG")
head(KEGG)
