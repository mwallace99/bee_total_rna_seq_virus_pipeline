#!/bin/bash

BASE_DIR="/g/data/${PROJECT}/${USER}/BLAST"

# Only samples that failed
FAILED_SAMPLES=(
    "K10F"
    "K11F"
    "K24F"
    "K24M"
)

echo "=========================================="
echo "MEMORY-SAFE PARALLEL BLAST RESUME"
echo "Total samples to resume: ${#FAILED_SAMPLES[@]}"
echo "=========================================="
echo ""

for SAMPLE in "${FAILED_SAMPLES[@]}"; do
    echo "Submitting SAFE parallel job for $SAMPLE..."
    
    PBS_FILE="${BASE_DIR}/temp_resume_${SAMPLE}.pbs"
    
    cat > "$PBS_FILE" << 'EOFPBS'
#!/bin/bash
#PBS -l ncpus=48
#PBS -l mem=190GB
#PBS -l walltime=24:00:00
#PBS -l wd
#PBS -q normal
#PBS -N resume_SAMPLE_PLACEHOLDER

# LOAD MODULES
source /g/data/${PROJECT}/${USER}/miniconda3/etc/profile.d/conda.sh
conda activate /g/data/${PROJECT}/${USER}/miniconda3/envs/BLAST

# CONFIGURATION
DB_PATH="/g/data/if89/data_library/blast_db/10082025"
DB_NAME="nt"
QUERY_DIR="/g/data/${PROJECT}/${USER}/BLAST/Bombus_Megahit-contigs"
BASE_DIR="/g/data/${PROJECT}/${USER}/BLAST"
CHUNK_DIR="${BASE_DIR}/query_chunks"
RESULTS_DIR="${BASE_DIR}/blast_results"
CHUNK_SIZE=1000

# --- SAFETY CONFIGURATION ---
# We limit to 6 concurrent jobs to ensure ~31GB RAM per job
# We increase threads to 8 to utilize the CPUs freed up by running fewer jobs.
# Calculation: 6 jobs * 8 threads = 48 CPUs (Full Node Usage)
THREADS_PER_JOB=8
MAX_JOBS=6 

# Use scratch for temporary operations
SCRATCH_WORK="/scratch/${PROJECT}/${USER}/blast_resume_SAMPLE_PLACEHOLDER_${PBS_JOBID}"
mkdir -p "${SCRATCH_WORK}"

export BLASTDB_LMDB_DISABLE=1
export BLASTDB="${DB_PATH}"

SAMPLE="SAMPLE_PLACEHOLDER"

echo "[$(date)] =========================================="
echo "[$(date)] RESUMING BLAST for:  ${SAMPLE}"
echo "[$(date)] Mode: MEMORY-SAFE PARALLEL"
echo "[$(date)] Concurrency: ${MAX_JOBS} jobs x ${THREADS_PER_JOB} threads"
echo "[$(date)] Memory per job: ~31GB"
echo "[$(date)] =========================================="

QUERY_FILE="${QUERY_DIR}/${SAMPLE}/${SAMPLE}_renamed_contigs.fa"

if [ ! -f "${QUERY_FILE}" ]; then
    echo "[$(date)] ERROR: Query file not found"
    rm -rf "${SCRATCH_WORK}"
    exit 1
fi

SAMPLE_CHUNK_DIR="${CHUNK_DIR}/${SAMPLE}_rerun"
SAMPLE_RESULT_DIR="${RESULTS_DIR}/${SAMPLE}"
SPLIT_MARKER="${SAMPLE_CHUNK_DIR}/splitting_complete.marker"

mkdir -p "${SAMPLE_CHUNK_DIR}" "${SAMPLE_RESULT_DIR}"

# --- STEP 1: PREPARE CHUNKS ---

if [ -d "${SAMPLE_CHUNK_DIR}" ] && [ -f "${SPLIT_MARKER}" ]; then
    echo "[$(date)] Found existing chunks and marker, resuming..."
else
    echo "[$(date)] Splitting FASTA..."
    rm -rf "${SAMPLE_CHUNK_DIR}"
    mkdir -p "${SAMPLE_CHUNK_DIR}"
    
    seqkit split2 "${QUERY_FILE}" -s ${CHUNK_SIZE} -O "${SAMPLE_CHUNK_DIR}" --force
    
    if [ $? -eq 0 ]; then
        touch "${SPLIT_MARKER}"
    else
        echo "[$(date)] ERROR: Seqkit split failed."
        rm -rf "${SCRATCH_WORK}"
        exit 1
    fi
fi

find "${SAMPLE_CHUNK_DIR}" -name "*.fa" | sort > "${SAMPLE_RESULT_DIR}/${SAMPLE}_chunks_rerun.txt"
TOTAL_CHUNKS=$(wc -l < "${SAMPLE_RESULT_DIR}/${SAMPLE}_chunks_rerun.txt")

echo "[$(date)] Total chunks to process: ${TOTAL_CHUNKS}"

# --- STEP 2: PARALLEL PROCESSING LOOP ---

echo "[$(date)] Starting parallel execution..."

while IFS= read -r CHUNK; do
    CHUNK_NAME=$(basename "${CHUNK}" .fa)
    BLAST_OUT="${SAMPLE_RESULT_DIR}/${CHUNK_NAME}_blast.txt"
    
    # Check if done (Atomic check)
    if [ -f "${BLAST_OUT}" ] && [ -s "${BLAST_OUT}" ]; then
        continue
    fi
    
    # WAIT if we have reached max concurrency
    # This loop pauses the script until a slot opens up
    while [ $(jobs -r | wc -l) -ge ${MAX_JOBS} ]; do
        sleep 2
    done

    echo "[$(date)] Launching: ${CHUNK_NAME}"
    
    # LAUNCH BACKGROUND JOB (Atomic Write Logic)
    (
        TEMP_OUT="${BLAST_OUT}.tmp"
        
        blastn \
            -db "${DB_NAME}" \
            -query "${CHUNK}" \
            -out "${TEMP_OUT}" \
            -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sskingdoms" \
            -evalue 1e-5 \
            -num_threads ${THREADS_PER_JOB} \
            2>/dev/null
        
        EXIT_VAL=$?
        
        if [ $EXIT_VAL -eq 0 ] && [ -f "${TEMP_OUT}" ]; then
            mv "${TEMP_OUT}" "${BLAST_OUT}"
        else
            # If it failed, we silently remove the temp file
            # The main script will catch the missing file at the end
            rm -f "${TEMP_OUT}"
        fi
    ) & 
    
    # Stagger launch slightly
    sleep 1

done < "${SAMPLE_RESULT_DIR}/${SAMPLE}_chunks_rerun.txt"

# WAIT for all remaining background jobs to finish
echo "[$(date)] All jobs submitted. Waiting for completion..."
wait

# --- STEP 3: VERIFICATION AND MERGE ---

COMPLETED_COUNT=$(find "${SAMPLE_RESULT_DIR}" -name "*_blast.txt" -type f ! -name "${SAMPLE}_blast.txt" | wc -l)
echo "[$(date)] Final count: ${COMPLETED_COUNT}/${TOTAL_CHUNKS} chunks completed."

if [ "${COMPLETED_COUNT}" -lt "${TOTAL_CHUNKS}" ]; then
    echo "[$(date)] WARNING: Not all chunks completed."
    echo "[$(date)] Re-submit this job to retry failed chunks."
    rm -rf "${SCRATCH_WORK}"
    exit 0
fi

echo "[$(date)] Merging results..."
rm -f "${SAMPLE_RESULT_DIR}/${SAMPLE}_blast.txt"
find "${SAMPLE_RESULT_DIR}" -name "*_blast.txt" -type f ! -name "${SAMPLE}_blast.txt" -exec cat {} + > "${SAMPLE_RESULT_DIR}/${SAMPLE}_blast.txt"

if [ ! -s "${SAMPLE_RESULT_DIR}/${SAMPLE}_blast.txt" ]; then
    echo "[$(date)] ERROR: Merged file empty."
    rm -rf "${SCRATCH_WORK}"
    exit 1
fi

# --- STEP 4: POST-PROCESSING ---

echo "[$(date)] Processing hits..."
cut -f1 "${SAMPLE_RESULT_DIR}/${SAMPLE}_blast.txt" | sort -T "${SCRATCH_WORK}" | uniq > "${SAMPLE_RESULT_DIR}/${SAMPLE}_blast_hits.txt"

echo "[$(date)] Extracting unknown..."
awk 'BEGIN {while (getline < "'"${SAMPLE_RESULT_DIR}"'/'"${SAMPLE}"'_blast_hits.txt") hits[$1] = 1} /^>/ {header = substr($1, 2); print_header = !(header in hits)} {if (print_header) print $0}' "${QUERY_FILE}" > "${SAMPLE_RESULT_DIR}/${SAMPLE}_unknown_contigs.fa"

sort -t $'\t' -k1,1 -k11,11g -T "${SCRATCH_WORK}" "${SAMPLE_RESULT_DIR}/${SAMPLE}_blast.txt" > "${SAMPLE_RESULT_DIR}/${SAMPLE}_sorted_blast.txt"

echo "[$(date)] Extracting viral..."
awk -F'\t' '! seen[$1]++ && ($13 == "Viruses" || $13 == "Unclassified" || $13 == "N/A" || $13 == "Other" || $13 == "")' "${SAMPLE_RESULT_DIR}/${SAMPLE}_sorted_blast.txt" > "${SAMPLE_RESULT_DIR}/${SAMPLE}_possible_viral_hits.txt"

cut -f1 "${SAMPLE_RESULT_DIR}/${SAMPLE}_possible_viral_hits.txt" | sort -T "${SCRATCH_WORK}" | uniq > "${SAMPLE_RESULT_DIR}/${SAMPLE}_possible_viral_query_ids.txt"

awk 'BEGIN {while (getline < "'"${SAMPLE_RESULT_DIR}"'/'"${SAMPLE}"'_possible_viral_query_ids.txt") viral[$1] = 1} /^>/ {header = substr($1, 2); print_header = (header in viral)} {if (print_header) print $0}' "${QUERY_FILE}" > "${SAMPLE_RESULT_DIR}/${SAMPLE}_possible_viral_contigs.fasta"

# CLEANUP
rm -rf "${SAMPLE_CHUNK_DIR}"
rm -rf "${SCRATCH_WORK}"
find "${SAMPLE_RESULT_DIR}" -name "*_blast.txt" -type f ! -name "${SAMPLE}_blast.txt" -delete

echo "[$(date)] ✓ SUCCESS: ${SAMPLE}"

EOFPBS

    sed -i "s/SAMPLE_PLACEHOLDER/${SAMPLE}/g" "$PBS_FILE"
    
    if [ -f "$PBS_FILE" ]; then
        # Project defaults to $PROJECT; storage is injected here (the #PBS line
        # in the generated script cannot expand shell variables).
        JOB_ID=$(qsub -l "storage=gdata/${PROJECT}+scratch/${PROJECT}+gdata/if89" "$PBS_FILE" 2>&1)
        echo "  Job submitted: $JOB_ID"
    fi
    sleep 1
done