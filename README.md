# rna-seq course project

analyzing toxoplasma infection in mouse lung samples

## data

from Singhania et al 2019, GEO accession GSE119855

15 lung samples:
- wildtype + infected (5)
- wildtype + uninfected (3)  
- DKO (double knockout) + infected (4)
- DKO + uninfected (3)

## workflow

1. fastqc for quality check // optional multiQC for compilation
2. hisat2 for mapping (mouse genome GRCm39)
3. featurecounts to get gene counts
4. deseq2 for differential expression
5. clusterprofiler for GO enrichment

## files

- `scripts/` - bash scripts for cluster jobs + R scripts for analysis
- `analysis_rnaseq/01_fastqc/` - fastqc reports
- `analysis_rnaseq/02_hisat/counts/` - count matrix from featurecounts

## running it

cluster jobs submitted with sbatch, see scripts folder

R analysis done locally, needs DESeq2 and clusterProfiler packages
