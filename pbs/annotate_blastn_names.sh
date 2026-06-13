#!/bin/bash
#PBS -P rg47
#PBS -l ncpus=1
#PBS -l mem=8gb
#PBS -l walltime=04:00:00
#PBS -l wd
#PBS -q normal
#PBS -l storage=gdata/rg47+gdata/if89
# ============================================================================
# annotate_blastn_names.sh - add scientific names to EXISTING blastn output
# ----------------------------------------------------------------------------
# For blastn tables that carry the subject accession in column 2 (sseqid) but no
# readable name (e.g. results produced before stage 13 added sscinames/stitle),
# this maps each subject accession to its scientific name + title with
# `blastdbcmd` (NO re-BLAST) and writes <input>_named.txt, appending two columns:
#       ... <sscinames>\t<stitle>
#
# Usage (run from the repo dir; uses NT_DB_PATH/NT_DB_NAME + the BLAST env from
# config.sh). Pass one or more blast TSV files:
#
#   bash pbs/annotate_blastn_names.sh FILE [FILE ...]
#   # e.g. every sample's possible-viral hits from an earlier run:
#   bash pbs/annotate_blastn_names.sh \
#        /g/data/rg47/mw9045/BLAST/blast_results/*/*_possible_viral_hits.txt
#
# It prints, per file, how many subjects resolved to a name so you can tell at a
# glance whether the accession format matched.
# ============================================================================
set -euo pipefail

_d="${PBS_O_WORKDIR:-$PWD}"
while [ "$_d" != "/" ] && [ ! -f "$_d/config.sh" ]; do _d="$(dirname "$_d")"; done
[ -f "$_d/config.sh" ] || { echo "ERROR: cannot find config.sh; run from the repo dir"; exit 1; }
source "$_d/config.sh"
cd "$PROJECT_DIR"

source "$CONDA_SH"
conda activate "$CONDA_ENV_BLAST"
export BLASTDB="$NT_DB_PATH"

# Files from positional args, or from `qsub -v FILES="f1 f2 ..."`.
if [ "$#" -gt 0 ]; then
    INPUTS=( "$@" )
elif [ -n "${FILES:-}" ]; then
    read -r -a INPUTS <<< "$FILES"
else
    echo "Usage: $0 <blast_tsv> [<blast_tsv> ...]   (subject accession must be column 2)" >&2
    exit 1
fi

tmp="$(mktemp -d "${SCRATCH_BASE:-/tmp}/annot.XXXXXX")"
trap 'rm -rf "$tmp"' EXIT

for f in "${INPUTS[@]}"; do
    if [ ! -s "$f" ]; then echo "skip (empty/missing): $f"; continue; fi
    echo "Annotating ${f}"

    # unique subject accessions (column 2)
    cut -f2 "$f" | sort -u > "$tmp/ids.txt"

    # accession -> scientific name <tab> title  (-target_only = one row per acc)
    blastdbcmd -db "$NT_DB_NAME" -entry_batch "$tmp/ids.txt" -target_only \
        -outfmt $'%a\t%S\t%t' 2>"$tmp/err" > "$tmp/map.tsv" || true

    out="${f%.txt}_named.txt"
    awk -F'\t' '
        # build lookup keyed by accession.version AND by bare accession
        NR==FNR {
            name[$1] = $2 "\t" $3
            base = $1; sub(/\.[0-9]+$/, "", base)
            if (!(base in name)) name[base] = $2 "\t" $3
            next
        }
        {
            key = $2
            if (!(key in name)) {                 # strip gi|..|ref|ACC| decorations
                n = split($2, a, "|"); for (i=n; i>=1; i--) if (a[i] != "") { key = a[i]; break }
            }
            if (!(key in name)) { sub(/\.[0-9]+$/, "", key) }   # try without version
            print $0 "\t" ((key in name) ? name[key] : "NA\tNA")
        }' "$tmp/map.tsv" "$f" > "$out"

    total=$(wc -l < "$f")
    na=$(grep -c $'\tNA\tNA$' "$out" || true)
    echo "  -> ${out}   (named $((total - na))/${total} rows; ${na} unresolved)"
done

echo "Done."
