#!/bin/bash
# ============================================================================
# lib/parse_samples.sh - sample-sheet helpers (source this after config.sh)
# ----------------------------------------------------------------------------
# Reads $SAMPLE_SHEET (one fastq filename per line, R1 and R2 on separate
# lines; blank lines and #-comments ignored) and exposes functions the stage
# scripts use instead of any hard-coded sample list:
#
#   list_samples              -> unique sample IDs, sorted
#   r1_files_for  <sample>    -> R1 fastqs for that sample, sorted (one per line)
#   r2_files_for  <sample>    -> R2 fastqs for that sample, sorted (one per line)
#   read_of       <file>      -> "R1" | "R2" | ""  (which read a filename is)
#   assert_pairs              -> exit 1 if any sample has #R1 != #R2
#   resolve_fastq <file>      -> absolute path (prepends $RAW_DIR for bare names)
#
# Requires config.sh to have been sourced first (for $SAMPLE_SHEET, $RAW_DIR,
# and sample_id_of()).
# ============================================================================

if [ -z "${SAMPLE_SHEET:-}" ]; then
    echo "ERROR: parse_samples.sh sourced before config.sh (SAMPLE_SHEET unset)" >&2
    return 1 2>/dev/null || exit 1
fi
if [ ! -f "${SAMPLE_SHEET}" ]; then
    echo "ERROR: sample sheet not found: ${SAMPLE_SHEET}" >&2
    return 1 2>/dev/null || exit 1
fi

# Resolve a sheet entry to a usable path: if it already has a directory
# component or is absolute, use it as-is; otherwise look for it in $RAW_DIR.
resolve_fastq() {
    local f="$1"
    case "$f" in
        /*|*/*) printf '%s\n' "$f" ;;
        *)      printf '%s\n' "${RAW_DIR}/$f" ;;
    esac
}

# Determine whether a filename is the R1 or R2 mate. Matches _R1 / _R2 followed
# by . _ or end-of-string, so it tolerates _R1.fastq.gz and _R1_001.fastq.gz.
read_of() {
    local b
    b="$(basename "$1")"
    if   [[ "$b" =~ _R1([._]|$) ]]; then printf 'R1\n'
    elif [[ "$b" =~ _R2([._]|$) ]]; then printf 'R2\n'
    else printf '\n'
    fi
}

# Emit the cleaned sheet (no comments / blank lines), one entry per line.
_sheet_entries() {
    grep -vE '^[[:space:]]*(#|$)' "${SAMPLE_SHEET}"
}

list_samples() {
    local line
    _sheet_entries | while IFS= read -r line; do
        sample_id_of "$line"
    done | sort -u
}

# Internal: emit fastqs for a given sample ($1) and read mate ($2 = R1|R2),
# resolved to paths and sorted (so lanes concatenate in a stable order).
_files_for() {
    local sample="$1" want="$2" line sid
    _sheet_entries | while IFS= read -r line; do
        sid="$(sample_id_of "$line")"
        [ "$sid" = "$sample" ] || continue
        [ "$(read_of "$line")" = "$want" ] || continue
        resolve_fastq "$line"
    done | sort
}

r1_files_for() { _files_for "$1" R1; }
r2_files_for() { _files_for "$1" R2; }

# Fail fast if the sheet is malformed (unequal R1/R2 counts, or unclassifiable
# entries). Call this at the top of stages that pair reads.
assert_pairs() {
    local rc=0 s n1 n2 line
    # warn about entries we could not classify as R1 or R2
    while IFS= read -r line; do
        if [ -z "$(read_of "$line")" ]; then
            echo "WARNING: cannot tell R1/R2 from: $line" >&2
            rc=1
        fi
    done < <(_sheet_entries)

    for s in $(list_samples); do
        n1="$(r1_files_for "$s" | grep -c .)"
        n2="$(r2_files_for "$s" | grep -c .)"
        if [ "$n1" -ne "$n2" ]; then
            echo "ERROR: sample ${s} has ${n1} R1 file(s) but ${n2} R2 file(s)" >&2
            rc=1
        elif [ "$n1" -eq 0 ]; then
            echo "ERROR: sample ${s} has no fastq files" >&2
            rc=1
        fi
    done
    return $rc
}
