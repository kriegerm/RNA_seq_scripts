#Tutorial:
#http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

install.packages("fdrtool")
install.packages("pheatmap")
library("fdrtool")
library("DESeq2")
library("Biobase")

getwd()
setwd("/Users/kriegema/Box/RNA-seq/Hua_RNAseq_11.2021/DESeq2_rRNAdepreads/")

#######CHOOSE WHICH FILE YOU WANT TO LOAD FOR FUTURE ANALYSIS
countData = read.csv( "/Users/kriegema/Box/RNA-seq/Hua_RNAseq_11.2021/DESeq2_rRNAdepreads/countsandconditions/WT_1375-CDM.csv", header=TRUE, row.names = 1, check.names=FALSE, sep="," )
colData = read.table("countsandconditions/conditions_WT_1375_CDM.txt", header=TRUE, check.names=FALSE, row.names=1 )

countData = read.csv( "countsandconditions/WT_1377-CDM.csv", header=TRUE, row.names = 1, check.names=FALSE, sep="," )
colData = read.table( "countsandconditions/conditions_WT_1377_CDM.txt", header=TRUE, check.names=FALSE, row.names=1 )

countData = read.csv( "countsandconditions/WT_1375-THYE.csv", header=TRUE, row.names = 1, check.names=FALSE, sep="," )
colData = read.table( "countsandconditions/conditions_WT_1375_THYE.txt", header=TRUE, check.names=FALSE, row.names=1 )

countData = read.csv( "countsandconditions/WT_1377-THYE.csv", header=TRUE, row.names = 1, check.names=FALSE, sep="," )
colData = read.table( "countsandconditions/conditions_WT_1377_THYE.txt", header=TRUE, check.names=FALSE, row.names=1 )

countData = read.csv( "totalcounts_norRNAnotRNA_CDM.csv", header=TRUE, row.names = 1, check.names=FALSE, sep="," )
colData = read.table( "conditions_CDM.txt", header=TRUE, check.names=FALSE, row.names=1 )

countData = read.csv( "totalcounts_norRNAnotRNA_THYE.csv", header=TRUE, row.names = 1, check.names=FALSE, sep="," )
colData = read.table( "conditions_THYE.txt", header=TRUE, check.names=FALSE, row.names=1 )


head(countData)
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ condition)
dds <- rowSums(counts(dds)) >= 100 #Keep only rows with more than 100 gene counts
dds <- dds[keep,]

#######CHOOSE LINE TO RUN BASED ON THE INPUT FILES
#dds$condition <- factor(dds$condition, levels=c("WT_CDM", "1375_CDM", "1377_CDM"))

dds$condition <- factor(dds$condition, levels=c("WT_CDM", "1375_CDM"))
dds$condition <- factor(dds$condition, levels=c("WT_CDM", "1377_CDM"))

#dds$condition <- factor(dds$condition, levels=c("WT_THYE", "1375_THYE", "1377_THYE"))
dds$condition <- factor(dds$condition, levels=c("WT_THYE", "1375_THYE"))
dds$condition <- factor(dds$condition, levels=c("WT_THYE", "1377_THYE"))

dds$condition <- factor(dds$condition, levels=c("WT_CDM", "WT_THYE"))
dds$condition <- factor(dds$condition, levels=c("1375_CDM", "1375_THYE"))
dds$condition <- factor(dds$condition, levels=c("1377_CDM", "1377_THYE"))

#dds$condition <- factor(dds$condition, levels=c("WT_CDM", "1375_CDM", "1377_CDM","WT_THYE", "1375_THYE", "1377_THYE"))

dds <- DESeq(dds)


##Data Transformations and PCAs
vsd <- vst(dds, blind=FALSE)
plotPCA(vsd, intgroup="condition")

rld <- rlog(dds, blind=FALSE)
plotPCA(rld, intgroup="condition")


##Outliers - boxplot
par(mar=c(8,5,2,2))
boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)


##Heat Map
library("RColorBrewer")
library("pheatmap")
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$type, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)


#######CHOOSE LINE TO RUN BASED ON THE INPUT FILES
res_WT_1375_CDM <- results(dds, contrast=c("condition","WT_CDM", "1375_CDM"))
res_WT_1377_CDM <- results(dds, contrast=c("condition","WT_CDM", "1377_CDM"))
res_WT_1375_THYE <- results(dds, contrast=c("condition","WT_THYE", "1375_THYE"))
res_WT_1377_THYE <- results(dds, contrast=c("condition","WT_THYE", "1377_THYE"))

res_WT_CDM_THYE <- results(dds, contrast=c("condition","WT_THYE", "WT_CDM"))
res_1375_CDM_THYE <- results(dds, contrast=c("condition","1375_THYE", "1375_CDM"))
res_1377_CDM_THYE <- results(dds, contrast=c("condition","1377_THYE", "1377_CDM"))

#######CHOOSE LINE TO RUN BASED ON THE INPUT FILES
#P-value hist
hist(res_WT_1375_CDM$pvalue, col = "lavender",  xlab = "p-values")
hist(res_WT_1377_CDM$pvalue, col = "lavender",  xlab = "p-values")
hist(res_WT_1375_THYE$pvalue, col = "lavender",  xlab = "p-values")
hist(res_WT_1377_THYE$pvalue, col = "lavender",  xlab = "p-values")

hist(res_WT_CDM_THYE$pvalue, col = "lavender",  xlab = "p-values")
hist(res_1375_CDM_THYE$pvalue, col = "lavender",  xlab = "p-values")
hist(res_1377_CDM_THYE$pvalue, col = "lavender",  xlab = "p-values")

#######CHOOSE LINE TO RUN BASED ON THE INPUT FILES
#Write the results if they look ok based on p-vals
write.csv(res_WT_1375_CDM, file="DEGs_WT_1375_CDM.csv")
write.csv(res_WT_1377_CDM, file="DEGs_WT_1377_CDM.csv")


write.csv(res_WT_1375_THYE, file="DEGs_WT_1375_THYE.csv")
write.csv(res_WT_1377_THYE, file="DEGs_WT_1377_THYE.csv")
write.csv(res_WT_CDM_THYE, file="DEGs_WT_CDM_THYE.csv")
write.csv(res_1375_CDM_THYE, file="DEGs_1375_CDM_THYE.csv")
write.csv(res_1377_CDM_THYE, file="DEGs_1377_CDM_THYE.csv")


##if the graphs of your p-values look wonky....
###############################################################################################################
###############################################################################################################
###FDR STUFF from Shawn




FDR.ddsRes <- fdrtool(ddsRes$stat, statistic= "normal", plot = T)
ddsRes[,"padj"]  <- p.adjust(FDR.ddsRes$pval, method = "BH")
table(ddsRes[,"padj"] < 0.1) #27 TRUE 20136 FALSE

#######EDIT THIS "res" VARIABLE BASED ON THE INPUT FILES
res <- res_WT_1377_THYE
  
res <- res[ !is.na(res$padj), ]
res <- res[ !is.na(res$pvalue), ]
res <- res[, -which(names(res) == "pvalue")] #remove pvalue
res <- res[, -which(names(res) == "padj")] #remove adjusted pvalue
FDR.res <- fdrtool(res$stat, statistic= "normal", plot = T) #calculate new pvalues WTF IS GOING ON HERE
FDR.res$param[1, "sd"]
res[,"pvalue"] <- FDR.res$pval #add back in new pvalue
res[,"padj"]  <- p.adjust(FDR.res$pval, method = "BH") #add back in new adjusted pvalues
hist(res$pvalue, col = "lavender",  xlab = "p-values", main = "Histogram of FDR corrected p-vals")

table(res$padj<0.05)
res <- res[order(res$padj), ]
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "Gene"
head(resdata)
write.csv(resdata, file="WT_1377_THYE_FDRtool_DEGs.csv") ############CHANGE THE NAME OF YOUR FILE BASED ON INPUT

par(mfrow=c(1,1))
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))
with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
