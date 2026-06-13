#!/bin/bash
# ============================================================================
# submit_all.sh - submit the whole pipeline to PBS as a dependency chain
# ----------------------------------------------------------------------------
# Each stage starts only after the previous one finishes successfully
# (qsub -W depend=afterok). Run from the repo root:  ./submit_all.sh
#
# Prefer running 00_setup once by hand first (it builds the RdRp-scan DBs and
# the results tree); it is included here too and is safe to re-run.
#
# To run a single stage instead, just:  qsub pbs/06_rrna_bowtie2.pbs
# To preview a stage without running:    DRYRUN=1 bash pbs/06_rrna_bowtie2.pbs
# ============================================================================
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

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

prev=""
for stage in "${STAGES[@]}"; do
    if [ -z "$prev" ]; then
        jid=$(qsub "pbs/${stage}")
    else
        jid=$(qsub -W "depend=afterok:${prev}" "pbs/${stage}")
    fi
    echo "Submitted ${stage}: ${jid}"
    prev="$jid"
done

echo "All stages submitted as a dependency chain."
