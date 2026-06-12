#!/bin/bash

# ============================================================================
# 2.1 - Submit Diamond BLASTX array job
# ----------------------------------------------------------------------------
# Generates sample list from unknown contigs, then submits one array task
# per sample.
# ============================================================================

BASE_DIR="/g/data/rg47/mw9045/BLAST"
SAMPLE_LIST="${BASE_DIR}/diamond_samples.txt"

# Build sample list from BLASTn unknown contigs
find "${BASE_DIR}/blast_results" -name "*_unknown_contigs.fa" -type f -size +0c | sort > "${SAMPLE_LIST}"

TOTAL_SAMPLES=$(wc -l < "${SAMPLE_LIST}")

if [ "$TOTAL_SAMPLES" -eq 0 ]; then
    echo "ERROR: No unknown contig files found."
    exit 1
fi

echo "Total samples to process: $TOTAL_SAMPLES"
echo "Sample list: ${SAMPLE_LIST}"

# Submit array job
qsub -J 1-${TOTAL_SAMPLES} 2.1_blastx_diamond_parallel_MW_Version2.pbs

echo "Submitted array job for $TOTAL_SAMPLES samples"
