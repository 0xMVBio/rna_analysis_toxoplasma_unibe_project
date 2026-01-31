# step 7 - GO enrichment analysis
# run after step 6, needs res_infection object

library(clusterProfiler)
library(org.Mm.eg.db)

# all genes tested (universe)
all_genes <- rownames(res_infection)[!is.na(res_infection$padj)]

# all DE genes (up + down combined)
de_genes_all <- rownames(res_infection)[which(res_infection$padj < 0.05)]
length(de_genes_all)

# upregulated genes
de_genes_up <- rownames(res_infection)[which(res_infection$padj < 0.05 & res_infection$log2FoldChange > 0)]
length(de_genes_up)

# downregulated genes
de_genes_down <- rownames(res_infection)[which(res_infection$padj < 0.05 & res_infection$log2FoldChange < 0)]
length(de_genes_down)


## GO enrichment - upregulated genes ##

go_up <- enrichGO(
  gene = de_genes_up,
  universe = all_genes,
  OrgDb = org.Mm.eg.db,
  keyType = "ENSEMBL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05
)

head(go_up, 15)
dotplot(go_up, showCategory = 20, font.size = 5, title = "Upregulated genes - GO BP")

write.csv(as.data.frame(go_up), "GO_upregulated_BP.csv")


## GO enrichment - downregulated genes ##

go_down <- enrichGO(
  gene = de_genes_down,
  universe = all_genes,
  OrgDb = org.Mm.eg.db,
  keyType = "ENSEMBL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05
)

head(go_down, 15)
dotplot(go_down, showCategory = 20, font.size = 5, title = "Downregulated genes - GO BP")

write.csv(as.data.frame(go_down), "GO_downregulated_BP.csv")


## GO enrichment - all DE genes (combined) ##

go_all <- enrichGO(
  gene = de_genes_all,
  universe = all_genes,
  OrgDb = org.Mm.eg.db,
  keyType = "ENSEMBL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05
)

head(go_all, 15)
dotplot(go_all, showCategory = 20, font.size = 5, title = "All DE genes - GO BP")

write.csv(as.data.frame(go_all), "GO_all_DE_BP.csv")
