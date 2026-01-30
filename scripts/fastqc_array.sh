#!/usr/bin/env bash
# usage: sbatch fastqc_array.sh <samplelist.tsv> <outdir>

#SBATCH -J fastqc
#SBATCH --array=1-15
#SBATCH -c 1
#SBATCH --time=02:00:00        # would be too short if not launched as array i think 
#SBATCH --partition=pibu_el8    # compute partition 
#SBATCH --ntasks=1             
#SBATCH --mem=4G

# below are my err and output files with correct naming
#SBATCH -o /data/users/mvaldivia/rna_seq_course/analysis_rnaseq/01_fastqc/logs/%x_%A_%a.out
#SBATCH -e /data/users/mvaldivia/rna_seq_course/analysis_rnaseq/01_fastqc/logs/%x_%A_%a.err


# trick to get the whole pipe to fail if anything fails so we don't get residues
set -euo pipefail

# paths, uses fields for list and output but I like hardcoding my logs
SAMPLELIST="$1"
OUTDIR="$2"
LOGDIR="/data/users/mvaldivia/rna_seq_course/analysis_rnaseq/01_fastqc/logs"
IMG="/containers/apptainer/fastqc-0.12.1.sif"

mkdir -p "${OUTDIR}" "${LOGDIR}"

# get the line for this array task from samplelist; had some help here
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLELIST")
SAMPLE=$(echo "$LINE" | cut -f1)
R1=$(echo "$LINE" | cut -f2)
R2=$(echo "$LINE" | cut -f3)


#run fastqc and print whats happening to output; -B flag so apptainer sees wht's outside
echo "[$(date)] FastQC ${SAMPLE}"
apptainer exec -B /data:/data -B /home:/home "${IMG}" \
  fastqc --outdir "${OUTDIR}" "${R1}" "${R2}"
echo "[$(date)] Done ${SAMPLE}"
