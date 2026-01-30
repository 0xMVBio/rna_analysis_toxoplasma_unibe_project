#!/bin/bash
# usage: sbatch hisat2_map.sh <index_prefix> <reads_dir> <output_dir>

#SBATCH -J hisat2_map
#SBATCH -p pibu_el8
#SBATCH -c 4                   # more cores necessary this time;mapping is compute intensive
#SBATCH --mem=8000
#SBATCH -t 12:00:00
#SBATCH -o hisat2_map_%j.out
#SBATCH -e hisat2_map_%j.err

# map paired-end reads to reference genome with hisat2 / they're all paired end - two files per

# fail if anything goes wrong
set -euo pipefail

INDEX_PREFIX="$1"
READS_DIR="$2"
OUT_DIR="$3"
CONTAINER="/containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif"

mkdir -p "$OUT_DIR"

# loop through all R1 files
for R1 in "$READS_DIR"/*_1.fastq.gz; do
  # get sample name and find R2
  base=$(basename "$R1")
  sample="${base%%_1*}"
  R2="${R1/_1.fastq.gz/_2.fastq.gz}"
  
  echo "Mapping $sample"
  
  # run hisat2 (RF = reverse stranded) and pipe to samtools to make bam
  apptainer exec --bind /data:/data "$CONTAINER" \
    bash -c "hisat2 -p 4 --rna-strandness RF -x $INDEX_PREFIX -1 $R1 -2 $R2 \
      | samtools view -b -o $OUT_DIR/${sample}.bam -"
done

echo "Done. BAMs in: $OUT_DIR"
