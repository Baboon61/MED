#!/usr/bin/env Rscript
#args = commandArgs(trailingOnly=TRUE)

# Start DESeq2
library("DESeq2")
library("apeglm")
library("vsn")
library("hexbin")
library("pheatmap")
library("pcaExplorer")
library("RColorBrewer")
library("stringr")

df <- read.table("/media/bastien/e3dfe106-a690-43a4-a5a0-b387d2d04861/Data/MED/479_featurecounts.txt", header = TRUE, sep = "\t", row.names= 1)
df$Start <- str_split_fixed(df$Start, ";",2)[,1]
df$End <- str_split_fixed(df$End, ";",2)[,1]
df$Chr <- str_split_fixed(df$Chr, ";",2)[,1]
df$Strand <- str_split_fixed(df$Strand, ";",2)[,1]
igh_genes <- rownames(df[(df$Chr == "12" & df$Start >= 113225000 & df$End <= 116010000),])


non_informative_genes <- rownames(df[(nchar(df$Chr) > 2 | df$Chr == "X" | df$Chr == "Y" | df$Chr == "MT"),])


df <- df[,-c(1,2,3,4,5)]
df <- df[,c(5,6,7,1,2,3,4)]

# Create design
df_infos <- data.frame(c('untreated', 'untreated','untreated', 'treated', 'treated', 'treated', 'treated'), c('paired-end', 'paired-end', 'paired-end','paired-end', 'paired-end', 'paired-end', 'paired-end'))
row.names(df_infos) <- c('WT7','WT8','WT9','s3145','s3146','s3182','s3184')
colnames(df_infos) <- c("condition", "type")

# Check of we have the same colnames in df & design
all(rownames(df_infos) %in% colnames(df))

dds <- DESeqDataSetFromMatrix(countData = df,
                              colData = df_infos,
                              design = ~ condition)

# Pre-filtering removing rows that have only 0 or 1 read

dds$condition <- relevel(dds$condition, ref = "untreated")

dds <- DESeq(dds)

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

#PRINT FIGURE TO GET SOME FILTER
png("/media/bastien/e3dfe106-a690-43a4-a5a0-b387d2d04861/Data/MED/sizeFactor.png")
plot(sizeFactors(dds), colSums(counts(dds)))
abline(lm(colSums(counts(dds)) ~ sizeFactors(dds) + 0))
dev.off()

res <- results(dds, name="condition_treated_vs_untreated")

png("/media/bastien/e3dfe106-a690-43a4-a5a0-b387d2d04861/Data/MED/MA.png")
plotMA(res, ylim=c(-2,2))
dev.off()

resLFC <- lfcShrink(dds, coef="condition_treated_vs_untreated", type="apeglm")

png("/media/bastien/e3dfe106-a690-43a4-a5a0-b387d2d04861/Data/MED/MA_normalized.png")
plotMA(resLFC, ylim=c(-2,2))
dev.off()

rld <- rlog(dds, blind=FALSE)

png("/media/bastien/e3dfe106-a690-43a4-a5a0-b387d2d04861/Data/MED/MeanSDplot.png")
meanSdPlot(assay(rld), ranks = FALSE)
dev.off()

png("/media/bastien/e3dfe106-a690-43a4-a5a0-b387d2d04861/Data/MED/violin.png")
distro_expr(rld, plot_type = "violin")
dev.off()

png("/media/bastien/e3dfe106-a690-43a4-a5a0-b387d2d04861/Data/MED/PCA.png")
plotPCA(rld, intgroup=c("condition", "type"))
dev.off()

#SHRINK DDS WITH RESULT OF FILTER
keep <- rowSums(counts(dds)) >= 100
dds <- dds[keep,]

res <- results(dds, name="condition_treated_vs_untreated")

rld <- rlog(dds, blind=FALSE)

res <- subset(res, !(rownames(res) %in% non_informative_genes))

res05 <- subset(res, res$padj < 0.05)

#BEST 50
topVarGenes <- head(order(-abs(res05$log2FoldChange)),50)

mat <- assay(dds)[ c(rownames(res05[topVarGenes,])), ]

df <- as.data.frame(colData(rld)[,c("condition")])
rownames(df) <- colnames(mat)
colnames(df) <- "condition"
df$condition <- c("WT", "WT", "WT", "KO", "KO", "KO", "KO")
pheatmap(log2(mat+1), annotation_col=df, filename="/media/bastien/e3dfe106-a690-43a4-a5a0-b387d2d04861/Data/MED/heatmap_50.png", border_color=NA, treeheight_row = 150, treeheight_col = 50, fontsize = 5, display_numbers=TRUE)

#BEST 300
topVarGenes <- head(order(-abs(res05$log2FoldChange)),300)

mat <- assay(dds)[ c(rownames(res05[topVarGenes,])), ]

df <- as.data.frame(colData(rld)[,c("condition")])
rownames(df) <- colnames(mat)
colnames(df) <- "condition"
df$condition <- c("WT", "WT", "WT", "KO", "KO", "KO", "KO")
pheatmap(log2(mat+1), annotation_col=df, filename="/media/bastien/e3dfe106-a690-43a4-a5a0-b387d2d04861/Data/MED/heatmap_300.png", width = 15, height = 20, border_color=NA, cellwidth=40, treeheight_row = 150, treeheight_col = 50, fontsize = 3, display_numbers=TRUE)

#BEST 50 IGH
res_igh <- res[rownames(res) %in% igh_genes ,]
topVarGenes_igh <- head(order(-abs(res_igh$log2FoldChange)),50)

mat_igh <- assay(dds)[ c(rownames(res_igh[topVarGenes_igh,])), ]

df <- as.data.frame(colData(rld)[,c("condition")])
rownames(df) <- colnames(mat_igh)
colnames(df) <- "condition"
df$condition <- c("WT", "WT", "WT", "KO", "KO", "KO", "KO")

pheatmap(log2(mat_igh+1), annotation_col= df, filename="/media/bastien/e3dfe106-a690-43a4-a5a0-b387d2d04861/Data/MED/heatmap_igh_50.png", border_color=NA, treeheight_row = 150, treeheight_col = 50, fontsize = 5, display_numbers=TRUE)

#ALL IGH
res_igh <- res[rownames(res) %in% igh_genes ,]
topVarGenes_igh <- head(order(-abs(res_igh$log2FoldChange)),200)

mat_igh <- assay(dds)[ c(rownames(res_igh)), ]

df <- as.data.frame(colData(rld)[,c("condition")])
rownames(df) <- colnames(mat_igh)
colnames(df) <- "condition"
df$condition <- c("WT", "WT", "WT", "KO", "KO", "KO", "KO")
pheatmap(log2(mat_igh+1), annotation_col= df, filename="/media/bastien/e3dfe106-a690-43a4-a5a0-b387d2d04861/Data/MED/heatmap_igh_all.png", width = 15, height = 20, border_color=NA, cellwidth=40, treeheight_row = 150, treeheight_col = 50, fontsize = 3, display_numbers=TRUE)

#TEST
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- colnames(rld)
colnames(sampleDistMatrix) <- colnames(rld)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
png("/media/bastien/e3dfe106-a690-43a4-a5a0-b387d2d04861/Data/MED/sampleDists.png")
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.off()