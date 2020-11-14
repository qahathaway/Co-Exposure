
####Salmon/Tximport####

library(Co-Exposure)

library(GenomicFeatures)
library(tximport)
library(readr)
library(rjson)
library(dplyr)
library(AnnotationDbi)

library(tximportData)
dir <- system.file("extdata", package = "Co-Exposure", mustWork=TRUE)
list.files(dir)
list.files(file.path(dir, "quants"))

csvfile <- file.path(dir, "sample_table.csv")
coldata <- read.csv(csvfile, row.names=1, stringsAsFactors=FALSE)
coldata

coldata <- coldata[1:28,]
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


####Make Counts Table####

assay(gse)
resOrderedDF <- as.data.frame(assay(gse))
write.csv(resOrderedDF, file = ".../counts_table.csv")


####DESeq2####

round( colSums(assay(gse)) / 1e6, 1 )

library(DESeq2)
dds <- DESeqDataSet(gse, design = ~ treatment + duration)

nrow(dds)

keep <- rowSums(counts(dds) >= 3) >= 1
dds <- dds[keep,]
nrow(dds)


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
rownames(sampleDistMatrix) <- paste( vsd$treatment, vsd$duration, sep = " - " )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)

library("PoiClaClu")
poisd <- PoissonDistance(t(counts(dds)))

samplePoisDistMatrix <- as.matrix( poisd$dd )
rownames(samplePoisDistMatrix) <- paste( dds$treatment, dds$duration, sep=" - " )
colnames(samplePoisDistMatrix) <- NULL
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows = poisd$dd,
         clustering_distance_cols = poisd$dd,
         col = colors)


##PCA Plot##
#treatment and duration#
plotPCA(vsd, intgroup = c("treatment", "duration"))
pcaData <- plotPCA(vsd, intgroup = c("treatment", "duration"), returnData = TRUE)
pcaData
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2, color = treatment, shape = duration)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  ggtitle("PCA with VST data")

library("glmpca")
gpca <- glmpca(counts(dds), L=2)
gpca.dat <- gpca$factors
gpca.dat$treatment <- dds$treatment
gpca.dat$duration <- dds$duration

ggplot(gpca.dat, aes(x = dim1, y = dim2, color = treatment, shape = duration)) +
  geom_point(size =3) + coord_fixed() + ggtitle("glmpca - Generalized PCA")


##MDS Plot##
mds <- as.data.frame(colData(vsd))  %>%
  cbind(cmdscale(sampleDistMatrix))
ggplot(mds, aes(x = `1`, y = `2`, color = treatment, shape = duration)) +
  geom_point(size = 3) + coord_fixed() + ggtitle("MDS with VST data")

mdsPois <- as.data.frame(colData(dds)) %>%
  cbind(cmdscale(samplePoisDistMatrix))
ggplot(mdsPois, aes(x = `1`, y = `2`, color = treatment, shape = duration)) +
  geom_point(size = 3) + coord_fixed() + ggtitle("MDS with PoissonDistances")


##Differential Expression Analysis##
#library(ashr)
dds <- DESeq(dds)
resultsNames(dds)

###untrt = Sham-1, untrt2 = Sham-4, trt = CB-1, trt2 = CB-4, trt3 = O3-1, trt4 = O3-4, trt5 = CB-O3-1, trt6 = CB-O3-4###

res1 <- results(dds, contrast=c("treatment","trt","untrt"), alpha = 0.05)
mcols(res1, use.names = TRUE)
summary(res1)
table(res1$padj < 0.05)

res2 <- results(dds, contrast=c("treatment","trt3","untrt"), alpha = 0.05)
mcols(res2, use.names = TRUE)
summary(res2)
table(res2$padj < 0.05)

res3 <- results(dds, contrast=c("treatment","trt5","untrt"), alpha = 0.05)
mcols(res3, use.names = TRUE)
summary(res3)
table(res3$padj < 0.05)

res4 <- results(dds, contrast=c("treatment","trt2","untrt2"), alpha = 0.05)
mcols(res4, use.names = TRUE)
summary(res4)
table(res4$padj < 0.05)

res5 <- results(dds, contrast=c("treatment","trt4","untrt2"), alpha = 0.05)
mcols(res5, use.names = TRUE)
summary(res5)
table(res5$padj < 0.05)

res6 <- results(dds, contrast=c("treatment","trt6","untrt2"), alpha = 0.05)
mcols(res6, use.names = TRUE)
summary(res6)
table(res6$padj < 0.05)


##VennDiagram##
library(limma)

res2 <- results(dds, contrast=c("treatment","trt3","untrt"), alpha = 0.05)
res2 <- res2[which(res2$padj < 0.05),]
res2.genes <- row.names(res2)

res5 <- results(dds, contrast=c("treatment","trt4","untrt2"), alpha = 0.05)
res5 <- res5[which(res5$padj < 0.05),]
res5.genes <- row.names(res5)

res3 <- results(dds, contrast=c("group","KOFCD","WTFCD"), alpha = 0.05)
res3 <- res3[which(res3$padj < 0.05),]
res3.genes <- row.names(res3)

res6 <- results(dds, contrast=c("treatment","trt6","untrt2"), alpha = 0.05)
res6 <- res6[which(res6$padj < 0.05),]
res6.genes <- row.names(res6)

Unique <- sort(unique(c(res2.genes, res5.genes, res3.genes, res6.genes)))

res2.genes.2 <- Unique %in% res2.genes
res5.genes.2 <- Unique %in% res5.genes
res3.genes.2 <- Unique %in% res3.genes
res6.genes.2 <- Unique %in% res6.genes

counts.1 <- cbind(res2.genes.2,res5.genes.2,res3.genes.2,res6.genes.2)
results.1 <- vennCounts(counts.1, include="both")
vennDiagram(results.1, include="both", cex = 1, names = c("","","",""), circle.col = c("blue", "red", "green", "black"))


##Other Visualizations##
library(vidger)

vsBoxPlot(dds, d.factor = "treamtment", type = "deseq")
vsDEGMatrix(dds, d.factor = "treamtment", padj = 0.05, type = "deseq")
vsFourWay("trt2", "trt", "untrt2", d.factor = "group", dds, type = "deseq")
vsMAMatrix(dds, d.factor = "treamtment", type = "deseq", y.lim = c(-10,10))
vsScatterMatrix(dds, d.factor = "treamtment", type = "deseq")
vsScatterPlot("untrt2", "trt2", dds, d.factor = "treamtment", type = "deseq")
vsVolcano("untrt2", "trt2", dds, d.factor = "treamtment", type = "deseq", x.lim = c(-10,10), padj = 0.05)
vsVolcanoMatrix(dds, d.factor = "treamtment", type = "deseq", lfc = 2, padj = 0.05, x.lim = c(-8,8),
                title = FALSE, legend = TRUE, grid = TRUE, counts = FALSE, facet.title.size = 10)


##Volcano Plots##
library(EnhancedVolcano)
EnhancedVolcano(res1, lab = rownames(res1), x = "log2FoldChange", y = "padj", pCutoff = 0.05,
                xlim = c(-7.5, 7.5))


##Plotting Results##
topGene <- rownames(res1)[which.min(res1$padj)]
plotCounts(dds, gene = topGene, intgroup=c("treatment"))

library("ggbeeswarm")
geneCounts <- plotCounts(dds, gene = topGene, intgroup = c("treatment","duration"),
                         returnData = TRUE)

ggplot(geneCounts, aes(x = treatment, y = count, color = duration)) +
  scale_y_log10() +  geom_beeswarm(cex = 3)

ggplot(geneCounts, aes(x = treatment, y = count, color = duration, group = duration)) +
  scale_y_log10() + geom_point(size = 3) + geom_line()


##MA-Plot##
library("apeglm")
resultsNames(dds)

res1 <- lfcShrink(dds, coef="genotype_trt_vs_untrt", type="apeglm")
plotMA(res01, ylim = c(-5, 5))

res2 <- lfcShrink(dds, coef="genotype_trt3_vs_untrt", type="apeglm")
plotMA(res01, ylim = c(-5, 5))

res3 <- lfcShrink(dds, coef="genotype_trt5_vs_untrt", type="apeglm")
plotMA(res01, ylim = c(-5, 5))

res4 <- lfcShrink(dds, coef="genotype_trt2_vs_untrt2", type="apeglm")
plotMA(res01, ylim = c(-5, 5))

res5 <- lfcShrink(dds, coef="genotype_trt4_vs_untrt2", type="apeglm")
plotMA(res01, ylim = c(-5, 5))

res6 <- lfcShrink(dds, coef="genotype_trt6_vs_untrt2", type="apeglm")
plotMA(res01, ylim = c(-5, 5))


##Histogram of p vaules for genes with mean normalized count larger than 2##

hist(res1$pvalue[res$baseMean > 2], breaks = 0:20/20,
     col = "grey50", border = "white")


##Gene Clustering##

library("genefilter")
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 1000)
#or#
topGenes <- head(order(res01$padj),decreasing = TRUE, 50)

mat  <- assay(vsd)[ topVarGenes, ]
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vsd)[, c("treatment","duration")])
pheatmap(mat, annotation_col = anno)


##Annotating and Exporting##
library("AnnotationDbi")
library("org.Mm.eg.db")

columns(org.Mm.eg.db)

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


##Exporting Results##

resOrderedDF <- as.data.frame(resOrdered1)[1:25000, ]
write.csv(resOrderedDF, file = ".../res1.csv")

resOrderedDF <- as.data.frame(resOrdered2)[1:25000, ]
write.csv(resOrderedDF, file = ".../res2.csv")

resOrderedDF <- as.data.frame(resOrdered3)[1:25000, ]
write.csv(resOrderedDF, file = ".../res3.csv")

resOrderedDF <- as.data.frame(resOrdered4)[1:25000, ]
write.csv(resOrderedDF, file = ".../res4.csv")

resOrderedDF <- as.data.frame(resOrdered5)[1:25000, ]
write.csv(resOrderedDF, file = ".../res5.csv")

resOrderedDF <- as.data.frame(resOrdered6)[1:25000, ]
write.csv(resOrderedDF, file = ".../res6.csv")

resOrderedDF <- as.data.frame(resOrdered7)[1:25000, ]
write.csv(resOrderedDF, file = ".../res7.csv")

resOrderedDF <- as.data.frame(resOrdered8)[1:25000, ]
write.csv(resOrderedDF, file = ".../res8.csv")


##Genomic Space##
resGR <- lfcShrink(dds, coef="genotype_trt_vs_untrt", type="apeglm", format="GRanges")
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

rowsum.threshold <- 1 # user chosen
fdr.threshold <- 0.05 # user chosen
rs <- rowSums(counts(dds))
dds <- dds[ rs > rowsum.threshold ,]
dds <- DESeq(dds)
res <- results(dds, contrast=c("treatment","trt","untrt"), independentFiltering=FALSE, alpha = 0.05) # use count threshold instead of IF
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

write.csv(GO.wall, file = ".../GO.csv")

KEGG=goseq(pwf,'mm9','ensGene',test.cats="KEGG")
head(KEGG)
