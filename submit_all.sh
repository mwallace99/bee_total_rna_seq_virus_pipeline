#!/bin/bash
# ============================================================================
# submit_all.sh - submit the whole pipeline to PBS as a dependency chain
# ----------------------------------------------------------------------------
# Each stage starts only after the previous one finishes successfully
# (qsub -W depend=afterok). Run from the repo root:  ./submit_all.sh
#
# Stage 13 (blastn) is submitted as a PBS job ARRAY, one subjob per sample
# (qsub -J 1-N), so all samples blast in parallel across nodes instead of
# serializing into a single 48h job. The next chained stage waits for the
# whole array via depend=afterokarray. Each subjob is independently resumable:
# 13_blastn.pbs writes per-chunk outputs atomically and skips chunks already
# done, so re-submitting a failed/timed-out array picks up where it left off.
#
# Prefer running 00_setup once by hand first (it builds the RdRp-scan DBs and
# the results tree); it is included here too and is safe to re-run.
#
# To run a single stage instead, just:  qsub pbs/06_rrna_bowtie2.pbs
# To run stage 13 by hand as an array:   qsub -J 1-N pbs/13_blastn.pbs
# To preview a stage without running:    DRYRUN=1 bash pbs/06_rrna_bowtie2.pbs
# ============================================================================
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

# Sample count sizes the stage-13 blastn array (one subjob per sample).
source ./config.sh
source ./lib/parse_samples.sh
NSAMPLES=$(list_samples | grep -c .)
[ "$NSAMPLES" -ge 1 ] || { echo "ERROR: no samples found in $SAMPLE_SHEET" >&2; exit 1; }

# Inject project/storage from config so the #PBS directives in each stage (which
# cannot expand shell variables) are overridden with real values. PBS_ACCOUNT
# and PBS_STORAGE default to $PROJECT, which Gadi sets automatically.
QSUB_OPTS=(-P "$PBS_ACCOUNT" -l "storage=$PBS_STORAGE")

# Stage 13 (blastn) runs on hugemem as a GROUPED array: BLAST_NGROUPS subjobs,
# each warming the ~706 GB nt cache once and handling a slice of the samples
# (see config.sh). Cap the group count at the sample count. The stage that
# follows an array in the chain must depend on it with afterokarray (wait for
# ALL subjobs), not the plain afterok used between single jobs.
NGROUPS=${BLAST_NGROUPS:-5}
[ "$NGROUPS" -le "$NSAMPLES" ] || NGROUPS=$NSAMPLES
ARRAY_STAGES=" 13_blastn.pbs "

STAGES=(
    00_setup.pbs
    01_raw_fastqc.pbs
    02_trim_galore.pbs
    03_trimmed_fastqc.pbs
    04_concat_lanes.pbs
    05_concat_fastqc.pbs
    06_rrna_bowtie2.pbs
    07_rrna_fastqc.pbs
    08_star_build.pbs
    09_star_host.pbs
    10_post_star_rename.pbs
    11_star_fastqc.pbs
    12_megahit.pbs
    13_blastn.pbs
    # 14_blastn_add_names.pbs is intentionally NOT chained here: stage 13 already
    # emits scientific names, so 14 only back-fills names onto blastn tables that
    # lack them (e.g. older runs). Run it by hand when needed.
    15_blastx_diamond_nr.pbs
    16_add_taxid.pbs
    17_rdrpscan.pbs
    18_compare_virus_search.pbs
)

# Optionally restart partway through the chain (e.g. after qdel'ing a stage):
#   FROM=13_blastn ./submit_all.sh
# Submits from the first stage whose filename contains $FROM onward, with NO
# upstream dependency on the (assumed already-complete) earlier stages.
if [ -n "${FROM:-}" ]; then
    start=-1
    for i in "${!STAGES[@]}"; do
        [[ "${STAGES[$i]}" == *"$FROM"* ]] && { start=$i; break; }
    done
    [ "$start" -ge 0 ] || { echo "ERROR: FROM='$FROM' matches no stage in the chain" >&2; exit 1; }
    STAGES=("${STAGES[@]:$start}")
    echo "Starting from ${STAGES[0]} (FROM='$FROM'); earlier stages skipped."
fi

# START_DEP: make the FIRST submitted stage depend on an already-queued job,
# e.g. to chain 15-18 onto a stage-13 array submitted by an earlier run:
#   FROM=15_blastx START_DEP='afterokarray:171587169[].gadi-pbs' ./submit_all.sh
# Pass the full PBS depend value (type:jobid). Use the form your PBS accepts.
prev=""
prev_is_array=0
first=1
for stage in "${STAGES[@]}"; do
    # Dependency on the PREVIOUS job: afterokarray if it was an array, else afterok.
    dep=()
    if [ -n "$prev" ]; then
        if [ "$prev_is_array" -eq 1 ]; then
            dep=(-W "depend=afterokarray:${prev}")
        else
            dep=(-W "depend=afterok:${prev}")
        fi
    elif [ "$first" -eq 1 ] && [ -n "${START_DEP:-}" ]; then
        dep=(-W "depend=${START_DEP}")
    fi
    first=0

    # Submit THIS stage as a grouped array if it is listed above.
    #  -r y       : PBS requires array jobs to be rerunable (Gadi defaults -r n).
    #  -J 1-G     : G subjobs.
    #  -v ...     : pass NGROUPS (stride) into the job, and keep PROJECT/USER
    #               exported (our -v would otherwise replace Gadi's default
    #               `-v PROJECT`, leaving ${PROJECT} unset inside the job).
    arr=()
    this_is_array=0
    label=""
    if [[ "$ARRAY_STAGES" == *" $stage "* ]]; then
        arr=(-r y -J "1-${NGROUPS}" -v "NGROUPS=${NGROUPS},PROJECT,USER")
        this_is_array=1
        label=" (grouped array 1-${NGROUPS} over ${NSAMPLES} samples)"
    fi

    jid=$(qsub "${QSUB_OPTS[@]}" "${dep[@]}" "${arr[@]}" "pbs/${stage}")
    echo "Submitted ${stage}: ${jid}${label}"
    prev="$jid"
    prev_is_array=$this_is_array
done

echo "All stages submitted as a dependency chain."
