# DESeq2 analysis for lung toxoplasma RNA-seq data
# Following workflow from course materials

# load packages
library(DESeq2)
library(pheatmap)
library(ggplot2)

# check version
packageVersion("DESeq2")

# set working directory - change this to wherever you have the counts file
setwd("C:/Users/yi/Desktop/rnaseq_content")

# read in counts
counts <- read.delim("counts_matrix_clean.tsv", row.names = 1, check.names = FALSE)

# the column names have the full path, need to clean them up
colnames(counts) <- basename(colnames(counts))
colnames(counts) <- sub("\\.sorted\\.bam$", "", colnames(counts))

dim(counts)
colnames(counts)

# create sample metadata based on README file
# from /data/courses/rnaseq_course/toxoplasma_de/reads_Lung/README
# 
# WT infected: SRR7821918, 919, 920, 921, 922
# DKO infected: SRR7821923, 924, 925, 927
# WT uninfected: SRR7821937, 938, 939
# DKO uninfected: SRR7821940, 941, 942

sample_info <- data.frame(
  sample = c("SRR7821918", "SRR7821919", "SRR7821920", "SRR7821921", "SRR7821922",
             "SRR7821923", "SRR7821924", "SRR7821925", "SRR7821927",
             "SRR7821937", "SRR7821938", "SRR7821939",
             "SRR7821940", "SRR7821941", "SRR7821942"),
  genotype = c(rep("WT", 5), rep("DKO", 4), rep("WT", 3), rep("DKO", 3)),
  infection = c(rep("Infected", 5), rep("Infected", 4), rep("Uninfected", 3), rep("Uninfected", 3))
)

# combined group for plotting
sample_info$group <- paste(sample_info$genotype, sample_info$infection, sep = "_")
rownames(sample_info) <- sample_info$sample

# reorder to match counts
sample_info <- sample_info[colnames(counts), ]

# check that everything lines up
all(rownames(sample_info) == colnames(counts))

print(sample_info)


### Step 5: Exploratory analysis ###

# create DESeq object
dds <- DESeqDataSetFromMatrix(countData = counts, 
                               colData = sample_info, 
                               design = ~ genotype + infection)

# filter low counts
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep,]
nrow(dds)

# variance stabilizing transform for visualization
vsd <- vst(dds, blind = TRUE)

# PCA plot with 4 groups
pca_data <- plotPCA(vsd, intgroup = "group", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

ggplot(pca_data, aes(PC1, PC2, color = group)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw() +
  ggtitle("PCA of lung samples")

ggsave("PCA_4groups.png", width = 8, height = 6)

# also make one with shape for genotype
pca_data2 <- plotPCA(vsd, intgroup = c("infection", "genotype"), returnData = TRUE)

ggplot(pca_data2, aes(PC1, PC2, color = infection, shape = genotype)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw()

ggsave("PCA_infection_genotype.png", width = 8, height = 6)

# sample distance heatmap
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)

annotation <- data.frame(infection = sample_info$infection, 
                         genotype = sample_info$genotype,
                         row.names = rownames(sample_info))

pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         annotation_col = annotation,
         main = "Sample-to-sample distances (VST)")


### Step 6: Differential expression ###

# run DESeq
dds <- DESeq(dds)

resultsNames(dds)

# main comparison: infected vs uninfected
res <- results(dds, contrast = c("infection", "Infected", "Uninfected"))
res <- res[order(res$padj),]

summary(res)

# how many significant?
sum(res$padj < 0.05, na.rm = TRUE)

# up vs down
sum(res$padj < 0.05 & res$log2FoldChange > 0, na.rm = TRUE)  # up
sum(res$padj < 0.05 & res$log2FoldChange < 0, na.rm = TRUE)  # down

# top genes
head(res, 20)

# save results
write.csv(as.data.frame(res), "DE_results_infection.csv")

# MA plot
plotMA(res, main = "Infected vs Uninfected")

# volcano plot
res_df <- as.data.frame(res)
res_df$significant <- res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1

ggplot(res_df, aes(log2FoldChange, -log10(pvalue), color = significant)) +
  geom_point(alpha = 0.5) +
  scale_color_manual(values = c("grey", "red")) +
  theme_bw() +
  ggtitle("Volcano plot: Infected vs Uninfected")

ggsave("volcano_plot.png", width = 8, height = 6)


### Step 7: GO enrichment ###

library(clusterProfiler)
library(org.Mm.eg.db)

# get significant genes
sig_genes <- rownames(res)[which(res$padj < 0.05)]
length(sig_genes)

# all genes as background
all_genes <- rownames(res)[!is.na(res$padj)]

# GO biological process
ego <- enrichGO(gene = sig_genes,
                universe = all_genes,
                OrgDb = org.Mm.eg.db,
                keyType = "ENSEMBL",
                ont = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05)

head(ego, 15)

# save
write.csv(as.data.frame(ego), "GO_enrichment.csv")

# plots
dotplot(ego, showCategory = 15) + ggtitle("GO Biological Process")
ggsave("GO_dotplot.png", width = 10, height = 8)

barplot(ego, showCategory = 15)
ggsave("GO_barplot.png", width = 10, height = 8)

# session info for reproducibility
sessionInfo()
