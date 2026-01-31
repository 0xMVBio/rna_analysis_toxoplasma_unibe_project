# step 6 answering questions part - genes of interest from singhania paper
# run after step 6 main script, needs dds object

# these genes are highlighted in singhania et al 2019 as key responding genes to toxoplasma infection via IFN-gamma pathway

# Create explicit group variable for clearer plotting
dds$group <- factor(paste(dds$genotype, dds$infection, sep = "_"))

# Gbp2 - guanylate binding protein 2
# one of the main IFN-gamma induced genes, helps kill intracellular parasites
plotCounts(dds, gene = "ENSMUSG00000028270", intgroup = "group")

# Stat1 - signal transducer and activator of transcription 1
# central transcription factor for both type I and II interferon signaling
plotCounts(dds, gene = "ENSMUSG00000026104", intgroup = "group")

# Irgm1 - immunity related GTPase family M member 1
# important for autophagy and killing of intracellular pathogens
plotCounts(dds, gene = "ENSMUSG00000020009", intgroup = "group")


# check if these genes are DE in our results
res_infection["ENSMUSG00000028270", ]  # Gbp2
res_infection["ENSMUSG00000026104", ]  # Stat1
res_infection["ENSMUSG00000020009", ]  # Irgm1
