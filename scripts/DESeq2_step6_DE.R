# step 6 - differential expression
# run after step 5, needs dds object

library(DESeq2)

# run the DE analysis
dds <- DESeq(dds)

resultsNames(dds)


## comparison 1: infected vs uninfected (main effect) ##

res_infection <- results(dds, contrast = c("infection", "Infected", "Uninfected"))
res_infection <- res_infection[order(res_infection$padj), ]

summary(res_infection)

# how many DE
sum(res_infection$padj < 0.05, na.rm = TRUE)

# up and down
sum(res_infection$padj < 0.05 & res_infection$log2FoldChange > 0, na.rm = TRUE)
sum(res_infection$padj < 0.05 & res_infection$log2FoldChange < 0, na.rm = TRUE)

# top genes
head(res_infection, 20)

write.csv(as.data.frame(res_infection), "DE_infected_vs_uninfected.csv")


## comparison 2: genotype effect (DKO vs WT) ##

res_genotype <- results(dds, contrast = c("genotype", "DKO", "WT"))
summary(res_genotype)

write.csv(as.data.frame(res_genotype), "DE_DKO_vs_WT.csv")


## interaction: does infection effect differ between genotypes? ##
# need to refit with interaction term

dds_int <- dds
design(dds_int) <- ~ genotype + infection + genotype:infection
dds_int <- DESeq(dds_int)

resultsNames(dds_int)

# interaction effect - name from resultsNames(dds_int)
res_interaction <- results(dds_int, name = "genotypeWT.infectionUninfected")
summary(res_interaction)

# genes where infection response differs between WT and DKO
sum(res_interaction$padj < 0.05, na.rm = TRUE)

write.csv(as.data.frame(res_interaction), "DE_interaction_genotype_infection.csv")


## look at genes of interest from singhania paper ##
# interferon response genes - check a few

# normalised counts
norm_counts <- counts(dds, normalized = TRUE)

# plot counts for top DE gene
plotCounts(dds, gene = rownames(res_infection)[1], intgroup = c("infection", "genotype"))
