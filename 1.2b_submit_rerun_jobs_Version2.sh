#!/bin/bash

BASE_DIR="/g/data/rg47/mw9045/BLAST"

BAD_SAMPLES=(
    "K10F"
    "K11F"
    "K17M"
    "K1M"
    "K23M"
    "K24F"
    "K24M"
    "K26M"
    "K2F"
    "K6M"
    "K6F"
)

echo "=========================================="
echo "Submitting BLAST re-run jobs for corrupted samples"
echo "Total samples to process: ${#BAD_SAMPLES[@]}"
echo "=========================================="
echo ""

for SAMPLE in "${BAD_SAMPLES[@]}"; do
    echo "Submitting job for $SAMPLE..."
    
    PBS_FILE="${BASE_DIR}/temp_rerun_${SAMPLE}.pbs"
    
    cat > "$PBS_FILE" << 'EOFPBS'
#!/bin/bash
#PBS -P rg47
#PBS -l ncpus=4
#PBS -l mem=32GB
#PBS -l walltime=48:00:00
#PBS -l wd
#PBS -q normal
#PBS -l storage=gdata/rg47+scratch/rg47+gdata/if89
#PBS -N rerun_SAMPLE_PLACEHOLDER

source /g/data/rg47/mw9045/miniconda3/etc/profile.d/conda.sh
conda activate /g/data/rg47/mw9045/miniconda3/envs/BLAST

DB_PATH="/g/data/if89/data_library/blast_db/10082025"
DB_NAME="nt"
QUERY_DIR="/g/data/rg47/mw9045/BLAST/UKBombus_Megahit-contigs_MW"
BASE_DIR="/g/data/rg47/mw9045/BLAST"
CHUNK_DIR="${BASE_DIR}/query_chunks"
RESULTS_DIR="${BASE_DIR}/blast_results"
CHUNK_SIZE=1000

# Use scratch for ALL temporary operations
SCRATCH_WORK="/scratch/rg47/mw9045/blast_rerun_SAMPLE_PLACEHOLDER_${PBS_JOBID}"
mkdir -p "${SCRATCH_WORK}"

export BLASTDB_LMDB_DISABLE=1
export BLASTDB="${DB_PATH}"
export TMPDIR="${SCRATCH_WORK}"
export TEMP="${SCRATCH_WORK}"
export TMP="${SCRATCH_WORK}"

SAMPLE="SAMPLE_PLACEHOLDER"

echo "[$(date)] =========================================="
echo "[$(date)] Re-running BLAST for:  ${SAMPLE}"
echo "[$(date)] Using scratch: ${SCRATCH_WORK}"
echo "[$(date)] Running SEQUENTIALLY (no parallelism)"
echo "[$(date)] =========================================="

QUERY_FILE="${QUERY_DIR}/${SAMPLE}/${SAMPLE}_renamed_contigs.fa"

if [ !  -f "${QUERY_FILE}" ]; then
    echo "[$(date)] ERROR: Query file not found"
    exit 1
fi

SAMPLE_CHUNK_DIR="${CHUNK_DIR}/${SAMPLE}_rerun"
SAMPLE_RESULT_DIR="${RESULTS_DIR}/${SAMPLE}"

echo "[$(date)] Backing up corrupted results..."
[ -f "${SAMPLE_RESULT_DIR}/${SAMPLE}_blast.txt" ] && mv "${SAMPLE_RESULT_DIR}/${SAMPLE}_blast.txt" "${SAMPLE_RESULT_DIR}/${SAMPLE}_blast.txt.corrupted"
[ -f "${SAMPLE_RESULT_DIR}/${SAMPLE}_sorted_blast.txt" ] && mv "${SAMPLE_RESULT_DIR}/${SAMPLE}_sorted_blast.txt" "${SAMPLE_RESULT_DIR}/${SAMPLE}_sorted_blast.txt.corrupted"
[ -f "${SAMPLE_RESULT_DIR}/${SAMPLE}_possible_viral_hits.txt" ] && mv "${SAMPLE_RESULT_DIR}/${SAMPLE}_possible_viral_hits.txt" "${SAMPLE_RESULT_DIR}/${SAMPLE}_possible_viral_hits.txt.corrupted"

[ -d "${CHUNK_DIR}/${SAMPLE}" ] && rm -rf "${CHUNK_DIR}/${SAMPLE}"
rm -rf "${SAMPLE_CHUNK_DIR}"

mkdir -p "${SAMPLE_CHUNK_DIR}" "${SAMPLE_RESULT_DIR}"

echo "[$(date)] Splitting FASTA..."
seqkit split2 "${QUERY_FILE}" -s ${CHUNK_SIZE} -O "${SAMPLE_CHUNK_DIR}" --force

find "${SAMPLE_CHUNK_DIR}" -name "*.fa" | sort > "${SAMPLE_RESULT_DIR}/${SAMPLE}_chunks_rerun.txt"
TOTAL_CHUNKS=$(wc -l < "${SAMPLE_RESULT_DIR}/${SAMPLE}_chunks_rerun.txt")
echo "[$(date)] Created ${TOTAL_CHUNKS} chunks"

if [ "${TOTAL_CHUNKS}" -eq 0 ]; then
    echo "[$(date)] ERROR: No chunks created"
    rm -rf "${SCRATCH_WORK}"
    exit 1
fi

find "${SAMPLE_RESULT_DIR}" -name "*_blast.txt" -type f !  -name "${SAMPLE}_blast.txt" -delete

echo "[$(date)] Running BLAST sequentially on ${TOTAL_CHUNKS} chunks..."
CHUNK_NUM=0
FAILED_CHUNKS=0

while IFS= read -r CHUNK; do
    CHUNK_NUM=$((CHUNK_NUM + 1))
    CHUNK_NAME=$(basename "${CHUNK}" .fa)
    BLAST_OUT="${SAMPLE_RESULT_DIR}/${CHUNK_NAME}_blast.txt"
    
    echo "[$(date)] [$CHUNK_NUM/$TOTAL_CHUNKS] Processing:  ${CHUNK_NAME}"
    
    blastn \
        -db "${DB_NAME}" \
        -query "${CHUNK}" \
        -out "${BLAST_OUT}" \
        -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sskingdoms" \
        -evalue 1e-5 \
        -num_threads 4 \
        2>&1
    
    EXIT_CODE=$?
    
    if [ ${EXIT_CODE} -ne 0 ]; then
        echo "[$(date)] ERROR: BLAST failed for ${CHUNK_NAME} with exit code ${EXIT_CODE}"
        FAILED_CHUNKS=$((FAILED_CHUNKS + 1))
        
        # Allow some failures but not too many
        if [ ${FAILED_CHUNKS} -gt 5 ]; then
            echo "[$(date)] ERROR: Too many failed chunks (${FAILED_CHUNKS}), aborting"
            rm -rf "${SCRATCH_WORK}"
            exit 1
        fi
        continue
    fi
    
    if [ ! -f "${BLAST_OUT}" ]; then
        echo "[$(date)] ERROR: Output file not created for ${CHUNK_NAME}"
        FAILED_CHUNKS=$((FAILED_CHUNKS + 1))
        continue
    fi
    
    # Show progress every 10 chunks
    if [ $((CHUNK_NUM % 10)) -eq 0 ]; then
        echo "[$(date)] Progress: ${CHUNK_NUM}/${TOTAL_CHUNKS} chunks completed"
    fi
    
done < "${SAMPLE_RESULT_DIR}/${SAMPLE}_chunks_rerun.txt"

echo "[$(date)] BLAST completed:  ${CHUNK_NUM} total, ${FAILED_CHUNKS} failed"

if [ ${FAILED_CHUNKS} -gt 0 ]; then
    echo "[$(date)] WARNING: ${FAILED_CHUNKS} chunks failed - results may be incomplete"
fi

echo "[$(date)] Merging results..."
find "${SAMPLE_RESULT_DIR}" -name "*_blast.txt" -type f !  -name "${SAMPLE}_blast.txt" -exec cat {} + > "${SAMPLE_RESULT_DIR}/${SAMPLE}_blast.txt"

if [ ! -s "${SAMPLE_RESULT_DIR}/${SAMPLE}_blast.txt" ]; then
    echo "[$(date)] ERROR: Merged file empty"
    rm -rf "${SCRATCH_WORK}"
    exit 1
fi

TOTAL_LINES=$(wc -l < "${SAMPLE_RESULT_DIR}/${SAMPLE}_blast.txt")
echo "[$(date)] Total hits: ${TOTAL_LINES}"

FIRST_LINE=$(head -n 1 "${SAMPLE_RESULT_DIR}/${SAMPLE}_blast.txt")
NUM_FIELDS=$(echo "${FIRST_LINE}" | awk -F'\t' '{print NF}')
echo "[$(date)] Fields: ${NUM_FIELDS}"

if [ "${NUM_FIELDS}" -ne 13 ]; then
    echo "[$(date)] ERROR: Malformed results"
    echo "[$(date)] First line: ${FIRST_LINE: 0:200}"
    rm -rf "${SCRATCH_WORK}"
    exit 1
fi

echo "[$(date)] Processing results..."

echo "[$(date)] Extracting blast hits..."
cut -f1 "${SAMPLE_RESULT_DIR}/${SAMPLE}_blast.txt" | sort -T "${SCRATCH_WORK}" | uniq > "${SAMPLE_RESULT_DIR}/${SAMPLE}_blast_hits.txt"

echo "[$(date)] Extracting unknown contigs..."
awk 'BEGIN {
    while (getline < "'"${SAMPLE_RESULT_DIR}"'/'"${SAMPLE}"'_blast_hits.txt") hits[$1] = 1
}
/^>/ {
    header = substr($1, 2)
    print_header = !(header in hits)
}
{
    if (print_header) print $0
}' "${QUERY_FILE}" > "${SAMPLE_RESULT_DIR}/${SAMPLE}_unknown_contigs.fa"

echo "[$(date)] Sorting BLAST results..."
sort -t $'\t' -k1,1 -k11,11g -T "${SCRATCH_WORK}" "${SAMPLE_RESULT_DIR}/${SAMPLE}_blast.txt" > "${SAMPLE_RESULT_DIR}/${SAMPLE}_sorted_blast.txt"

echo "[$(date)] Extracting viral hits..."
awk -F'\t' '! seen[$1]++ && ($13 == "Viruses" || $13 == "Unclassified" || $13 == "N/A" || $13 == "Other" || $13 == "")' "${SAMPLE_RESULT_DIR}/${SAMPLE}_sorted_blast.txt" > "${SAMPLE_RESULT_DIR}/${SAMPLE}_possible_viral_hits.txt"

echo "[$(date)] Extracting viral query IDs..."
cut -f1 "${SAMPLE_RESULT_DIR}/${SAMPLE}_possible_viral_hits.txt" | sort -T "${SCRATCH_WORK}" | uniq > "${SAMPLE_RESULT_DIR}/${SAMPLE}_possible_viral_query_ids.txt"

echo "[$(date)] Extracting viral contigs..."
awk 'BEGIN {
    while (getline < "'"${SAMPLE_RESULT_DIR}"'/'"${SAMPLE}"'_possible_viral_query_ids.txt") viral[$1] = 1
}
/^>/ {
    header = substr($1, 2)
    print_header = (header in viral)
}
{
    if (print_header) print $0
}' "${QUERY_FILE}" > "${SAMPLE_RESULT_DIR}/${SAMPLE}_possible_viral_contigs.fasta"

NUM_TOTAL_HITS=$(wc -l < "${SAMPLE_RESULT_DIR}/${SAMPLE}_blast_hits.txt")
NUM_VIRAL_HITS=$(wc -l < "${SAMPLE_RESULT_DIR}/${SAMPLE}_possible_viral_query_ids.txt")

echo "[$(date)] =========================================="
echo "[$(date)] Summary for ${SAMPLE}:"
echo "[$(date)]   Total hits: ${NUM_TOTAL_HITS}"
echo "[$(date)]   Viral hits: ${NUM_VIRAL_HITS}"
echo "[$(date)] =========================================="

# Cleanup
rm -rf "${SAMPLE_CHUNK_DIR}"
rm -rf "${SCRATCH_WORK}"
find "${SAMPLE_RESULT_DIR}" -name "*_blast.txt" -type f ! -name "${SAMPLE}_blast.txt" -delete

echo "[$(date)] ✓ Completed ${SAMPLE}"
EOFPBS

    # Replace placeholder with actual sample name
    sed -i "s/SAMPLE_PLACEHOLDER/${SAMPLE}/g" "$PBS_FILE"
    
    if [ -f "$PBS_FILE" ]; then
        JOB_ID=$(qsub "$PBS_FILE" 2>&1)
        echo "  Job:  $JOB_ID"
    else
        echo "  ERROR: PBS file not created"
    fi
    
    sleep 1
done

echo ""
echo "=========================================="
echo "All submitted!"
echo "=========================================="