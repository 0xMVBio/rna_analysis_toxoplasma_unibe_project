#!/bin/bash
# usage: sbatch sort_and_index_BAM.sh <input_dir> <output_dir>

#SBATCH -J samtools_sort_index
#SBATCH -p pibu_el8
#SBATCH -c 4
#SBATCH --mem=25000
#SBATCH -t 04:00:00
#SBATCH -o sort_index_%j.out
#SBATCH -e sort_index_%j.err

# sort and index bam files

set -euo pipefail

IN_DIR="$1"
OUT_DIR="$2"
CONTAINER="/containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif"

mkdir -p "$OUT_DIR"

# loop through bam files
for BAM in "$IN_DIR"/*.bam; do
  sample=$(basename "$BAM" .bam)
  echo "Sorting and indexing $sample"
  
  apptainer exec --bind /data:/data "$CONTAINER" \
    bash -c "samtools sort -@ 4 -o $OUT_DIR/${sample}.sorted.bam $BAM && \
             samtools index -@ 4 $OUT_DIR/${sample}.sorted.bam"
done

echo "Done. Sorted BAMs in: $OUT_DIR"
