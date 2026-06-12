#!/bin/bash

BASE_DIR="/g/data/rg47/mw9045/BLAST"
QUERY_DIR="/g/data/rg47/mw9045/BLAST/UKBombus_Megahit-contigs_MW"

# Create the sample list first
echo "Creating sample list..."
find "$QUERY_DIR" -name "*_renamed_contigs.fa" -type f | sort > "$BASE_DIR/all_samples.txt"

TOTAL_SAMPLES=$(wc -l < "$BASE_DIR/all_samples.txt")
echo "Total samples to process: $TOTAL_SAMPLES"

if [ $TOTAL_SAMPLES -eq 0 ]; then
    echo "Error: No samples found in $QUERY_DIR"
    exit 1
fi

echo "Sample list created at: $BASE_DIR/all_samples.txt"
echo ""
echo "Submitting individual jobs..."

# Submit individual jobs (one per sample, each with unique job ID)
for i in $(seq 1 $TOTAL_SAMPLES); do
    SAMPLE_FILE=$(sed -n "${i}p" "$BASE_DIR/all_samples.txt")
    SAMPLE=$(basename $(dirname "$SAMPLE_FILE"))
    
    JOB_ID=$(qsub -r y -v JOB_NUM=$i 1.2_blastn_taxa_parallel_MW_Version2.pbs)
    echo "  [$JOB_ID] Sample $i: $SAMPLE"
done

echo ""
echo "Submitted $TOTAL_SAMPLES individual jobs to express queue"
echo "Each job has its own unique job ID"