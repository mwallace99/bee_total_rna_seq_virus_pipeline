#!/bin/bash
# ============================================================================
# submit_all.sh - submit the whole pipeline to PBS as a dependency chain
# ----------------------------------------------------------------------------
# Each stage starts only after the previous one finishes successfully
# (qsub -W depend=afterok). Run from the repo root:  ./submit_all.sh
#
# The heavy per-sample stages (13 blastn, 15 diamond, 17 rdrpscan) are submitted
# as GROUPED PBS arrays: G subjobs, each striding over a slice of the samples, so
# they run in parallel across nodes instead of serializing into one 48h job. Gadi
# caps an array at 10 elements, so G = min(BLAST_NGROUPS, samples, 10). Arrays go
# out with -r y (rerunnable) and each sample is skipped once its <sdir>/.done
# marker exists, so a timed-out or re-submitted array resumes where it left off.
# The next chained stage waits for the whole array via depend=afterok on the
# array id (jobid[]) - Gadi's PBS has no afterokarray type.
#
# Prefer running 00_setup once by hand first (it builds the RdRp-scan DBs and
# the results tree); it is included here too and is safe to re-run.
#
# To run a single stage instead, just:  qsub pbs/06_rrna_bowtie2.pbs
# To run an array stage by hand (G <= 10):
#   qsub -r y -J 1-G -v NGROUPS=G,PROJECT,USER pbs/13_blastn.pbs
# To preview a stage without running:    DRYRUN=1 bash pbs/06_rrna_bowtie2.pbs
# ============================================================================
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

# Sample count sizes the grouped arrays (stages 13/15/17).
source ./config.sh
source ./lib/parse_samples.sh
NSAMPLES=$(list_samples | grep -c .)
[ "$NSAMPLES" -ge 1 ] || { echo "ERROR: no samples found in $SAMPLE_SHEET" >&2; exit 1; }

# Inject project/storage from config so the #PBS directives in each stage (which
# cannot expand shell variables) are overridden with real values. PBS_ACCOUNT
# and PBS_STORAGE default to $PROJECT, which Gadi sets automatically.
QSUB_OPTS=(-P "$PBS_ACCOUNT" -l "storage=$PBS_STORAGE")

# The heavy per-sample stages (13 blastn, 15 diamond, 17 rdrpscan) run as GROUPED
# arrays: G subjobs, each striding over a slice of the samples. Array width G
# starts from BLAST_NGROUPS (see config.sh), then is capped at the sample count
# and at Gadi's 10-element array limit. The stage that follows an array depends on
# it with afterok on the array id (jobid[]), which on Gadi waits for ALL subjobs
# to finish OK.
NGROUPS=${BLAST_NGROUPS:-5}
[ "$NGROUPS" -le "$NSAMPLES" ] || NGROUPS=$NSAMPLES
[ "$NGROUPS" -le 10 ] || NGROUPS=10          # Gadi caps a job array at 10 elements
ARRAY_STAGES=" 13_blastn.pbs 15_blastx_diamond_nr.pbs 17_rdrpscan.pbs "

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
#   FROM=15_blastx START_DEP='afterok:171587169[].gadi-pbs' ./submit_all.sh
# Pass the full PBS depend value (type:jobid). Use the form your PBS accepts.
prev=""
first=1
for stage in "${STAGES[@]}"; do
    # Dependency on the PREVIOUS job. Gadi's PBS has no afterokarray type;
    # instead afterok on an array job id (jobid[]) waits for ALL its subjobs to
    # finish OK, so the same afterok form works for both arrays and plain jobs
    # (prev already carries the [] for an array).
    dep=()
    if [ -n "$prev" ]; then
        dep=(-W "depend=afterok:${prev}")
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
done

echo "All stages submitted as a dependency chain."
