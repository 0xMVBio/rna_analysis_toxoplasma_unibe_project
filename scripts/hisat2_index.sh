#!/bin/bash
# usage: sbatch hisat2_index.sh <genome.fa> <index_prefix>

#SBATCH -J hisat2_index
#SBATCH -p pibu_el8
#SBATCH -c 1
#SBATCH --mem=8000
#SBATCH -t 03:00:00
#SBATCH -o hisat2_index_%j.out
#SBATCH -e hisat2_index_%j.err

# build index from reference genome; fields for fasta path and prefix; it's GRcM39 but we stay agnostic, its prettier that way

FASTA="$1"
PREFIX="$2"
CONTAINER="/containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif"

echo "Indexing: $FASTA"
echo "Prefix: $PREFIX"

#runs the script; using apptainer for reproducibility
apptainer exec --bind /data:/data "$CONTAINER" \
  hisat2-build -f "$FASTA" "$PREFIX"

echo "Done. Created index files: ${PREFIX}.*.ht2"
