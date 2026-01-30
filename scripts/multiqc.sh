#!/bin/bash
# usage: sbatch multiqc.sh (paths hardcoded)

#SBATCH -J multiqc
#SBATCH -c 1
#SBATCH --time=01:00:00        # Time limit hh:mm:ss 
#SBATCH --partition=pibu_el8    # our compute node 
#SBATCH --ntasks=1             # n of tasks
#SBATCH --cpus-per-task=1 
#SBATCH --mem=4G

#overwrite if output alrdy exists
# -B flag makes "outside'  data visible on the bound mount
apptainer exec -B /data:/data -B /home:/home /containers/apptainer/multiqc-1.19.sif \
multiqc /data/users/mvaldivia/rna_seq_course/analysis_rnaseq/01_fastqc \
-o /data/users/mvaldivia/rna_seq_course/analysis_rnaseq/01_fastqc/multiqc --force
