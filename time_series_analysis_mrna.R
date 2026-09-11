#install packages
BiocManager::install("edgeR")
BiocManager::install("apeglm")
#loading libraries
library(BiocManager)
library(DESeq2)
library(edgeR)
#library(clusterProfiler)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(apeglm)

#importing files
counts = read.csv("gene_counts_matrix.csv", header = TRUE, row.names = 1, check.names = FALSE)
metadata = read.csv("metadata.csv",header = TRUE, row.names = 1, check.names = FALSE)

#convert counts to matrix
counts = as.matrix(counts)

#making sure the count matrix contains numeric values
mode(counts) = "numeric"

#check sample names
colnames(counts)
rownames(metadata)
identical(colnames(counts), rownames(metadata))

setdiff(colnames(counts), rownames(metadata))
setdiff(rownames(metadata), colnames(counts))

#matching the rownames of metadata with colnames of matrix
metadata <- metadata[match(colnames(counts), rownames(metadata)), , drop = FALSE]
identical(colnames(counts), rownames(metadata))


#create deseq2 object
dds = DESeqDataSetFromMatrix(countData = counts,
                             colData = metadata,
                             design = ~individual + time_point)
dds

#prefiltering
keep = rowSums(counts(dds)) >= 10
dds = dds[keep, ]

dds$condition = factor(dds$time_point, levels = c("0", "6", "24", "72"))

#important command :)
dds = DESeq(dds)
res = results(dds)
res

head(counts(dds))
colSums(counts(dds))
#colSums(counts(dds)) %>% barplot

counts(dds) %>% str
assay(dds) %>% str

#normalizing seq depth
gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x>0]), na.rm=na.rm)/length(x)) }
pseudo_refs = counts(dds) %>% apply(., 1, gm_mean)
pseudo_ref_ratios = counts(dds) %>% apply(., 2, function(cts){ cts/pseudo_refs})

counts(dds)[1,]/pseudo_refs[1]
apply(pseudo_ref_ratios, 2, median)

#calculate and apply the size factor
dds = estimateSizeFactors(dds)
plot(sizeFactors(dds), colSums(counts(dds)), 
     ylab = "library sizes", xlab = "size factors", cex = .6)

par(mfrow=c(1,2))
counts.sf_normalized = counts(dds, normalized=TRUE)
boxplot(counts.sf_normalized, main = "SF normalized", cex = .6)
boxplot(counts(dds), main = "read counts only", cex = .6)
