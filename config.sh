#!/bin/bash
# ============================================================================
# config.sh - single source of truth for the bee total-RNAseq virus pipeline
# ----------------------------------------------------------------------------
# Sourced by every pbs/NN_*.pbs stage and by lib/parse_samples.sh.
# Edit this file (and your ${STUDY}_rna_samples.txt sample sheet) - you should
# not need to edit any of the stage scripts to run a new project.
#
# Paths are kept RELATIVE to $PROJECT_DIR wherever possible. Large shared
# resources (reference genomes, nt / nr databases) cannot be relative, so they
# are set here as absolute paths in one place.
# ============================================================================

# ---- Project identity -------------------------------------------------------
# STUDY is the prefix of your sample sheet: ${STUDY}_rna_samples.txt
# NB: do NOT call this PROJECT - on NCI Gadi $PROJECT is a reserved variable that
# holds your allocation project (e.g. rg47) and is read by nqstat / nci_account.
STUDY="NZApis2026"

# Root of this checkout. Defaults to the directory containing config.sh, so the
# scripts work from wherever the repo is cloned. Override if you run the stages
# from a copy that lives elsewhere.
PROJECT_DIR="${PROJECT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"

# Directory holding the raw paired-end fastq files listed in the sample sheet.
# NZApis2026: AGRF delivery sits on /g/data, not under the repo, so point at it
# directly rather than symlinking into ${PROJECT_DIR}/raw.
RAW_DIR="/g/data/rg47/mw9045/NZApis_virus/AGRF_NXGSQCAGRF26050366-1_23MN3YLT3"

# The sample sheet: one fastq filename per line, R1 and R2 on separate lines.
# Blank lines and lines starting with # are ignored. See example_rna_samples.txt.
SAMPLE_SHEET="${PROJECT_DIR}/${STUDY}_rna_samples.txt"

# All pipeline outputs live under here, one sub-dir per stage.
# NB: the repo lives on /home (small quota), so results are sent to /g/data
# alongside the raw data instead - a 20-sample total-RNAseq run far exceeds the
# /home quota. Override per-run as needed.
RESULTS_DIR="/g/data/rg47/mw9045/NZApis_virus/results"

# ---- PBS / scheduler (NCI Gadi) ---------------------------------------------
# These are referenced in docs / submit_all.sh. The #PBS directives inside each
# stage script cannot read shell variables, so if you change project/storage you
# must also update the literal #PBS lines (they are identical in every stage).
PBS_ACCOUNT="rg47"
PBS_STORAGE="gdata/rg47+scratch/rg47+gdata/if89"
PBS_QUEUE="normal"

# Scratch area for sort/tmp during BLAST (large, fast, purgeable).
SCRATCH_BASE="/scratch/rg47/mw9045"

# ---- Conda environments -----------------------------------------------------
CONDA_SH="/g/data/rg47/mw9045/miniconda3/etc/profile.d/conda.sh"
CONDA_ENV="/g/data/rg47/mw9045/miniconda3/envs/Beeviromics"   # qc/trim/align/assemble
CONDA_ENV_BLAST="/g/data/rg47/mw9045/miniconda3/envs/BLAST"   # blast/diamond/seqkit
CONDA_ENV_RDRP="${CONDA_ENV_BLAST}"                            # getorf(EMBOSS)/hmmer/diamond

# ---- Host-depletion references ----------------------------------------------
# Apis mellifera (NZApis2026). All host refs live under host_rrna/apis.
# !! CONFIRM the exact filenames/prefix with:  ls /g/data/rg47/mw9045/host_rrna/apis
# and adjust the four paths below to match.
APIS_HOST_DIR="/g/data/rg47/mw9045/host_rrna/apis"

# rRNA bowtie2 index PREFIX (the path you pass to bowtie2 -x; the .bt2 files
# sit alongside it, e.g. Apis_rRNA.1.bt2).
RRNA_BT2_INDEX="${APIS_HOST_DIR}/Apis_rRNA"

# STAR host genome. If STAR_GENOME_DIR is empty/not yet built, stage 08 builds it
# from STAR_FASTA + STAR_GTF (Amel_HAv3.1 / GCF_003254395.2).
STAR_GENOME_DIR="${APIS_HOST_DIR}/Star-output_Apis"
STAR_FASTA="${APIS_HOST_DIR}/GCF_003254395.2_Amel_HAv3.1_genomic.fa"
STAR_GTF="${APIS_HOST_DIR}/GCF_003254395.2_Amel_HAv3.1_genomic.gtf"
STAR_SA_NBASES=13                # --genomeSAindexNbases (lower for small genomes)
STAR_LIMIT_BAM_SORT_RAM=4249358420

# ---- Sequence databases -----------------------------------------------------
# blastn nt database: directory holding the nt.* files, and its base name.
NT_DB_PATH="/g/data/if89/data_library/blast_db/10082025"
NT_DB_NAME="nt"

# DIAMOND protein nr database (.dmnd) for blastx of unknown contigs.
NR_DMND="/g/data/rg47/mw9045/BLAST/NR_db/nr.dmnd"

# taxid -> scientific-name table (tab separated: taxid <TAB> name) used by the
# blastx annotation R script.
NAMES_TXT="${PROJECT_DIR}/names.txt"

# ---- RdRp-scan (novel RNA virus discovery) ----------------------------------
# 00_setup.pbs clones the RdRp-scan repo into RDRPSCAN_DIR and builds the two
# canonical databases below into RDRPSCAN_DIR/db (locating the source files
# wherever they sit in the repo). The scan stage (16) uses these.
RDRPSCAN_REPO_URL="https://github.com/JustineCharon/RdRp-scan"
RDRPSCAN_DIR="/g/data/rg47/mw9045/RdRp-scan"
# DIAMOND db built by 00_setup from the RdRp-scan core protein fasta
# (RdRp-scan_0.90.fasta), so it matches your local DIAMOND version.
RDRPSCAN_DMND="${RDRPSCAN_DIR}/db/RdRp-scan.dmnd"
# RdRp HMM profile db (HMMER3). The repo ships it pre-pressed as
# Profile_db_and_alignments/RdRp_HMM_profile_CLUSTALO.db.h3{f,i,m,p}; this is the
# prefix hmmscan reads (the .h3* files sit alongside it). 00_setup verifies it.
RDRPSCAN_HMM="${RDRPSCAN_DIR}/Profile_db_and_alignments/RdRp_HMM_profile_CLUSTALO.db"

# ---- Tuning -----------------------------------------------------------------
# Threads used inside per-sample tool calls. Defaults to the PBS allocation if
# present, else 8.
THREADS="${PBS_NCPUS:-8}"

CHUNK_SIZE=1000                  # contigs per blastn chunk (memory control)
EVALUE_BLASTN="1e-5"
EVALUE_DIAMOND="1e-5"
EVALUE_HMM="1e-6"
GETORF_MINLEN=200                # aa, RdRp-scan ORF length cutoff
GETORF_TABLES="1 3 4 5 6 11 16"  # genetic codes used in viruses (RdRp-scan)

# Set DRYRUN=1 in the environment to make stages print the heavy tool commands
# they would run instead of executing them. e.g.:
#   DRYRUN=1 bash pbs/06_rrna_bowtie2.pbs
DRYRUN="${DRYRUN:-0}"

# run <cmd> [args...] : execute a command, or just print it when DRYRUN=1.
# Redirections written on the same line (e.g. `run bowtie2 ... 2> log`) bind to
# this function call, so the child still writes to the logfile when it runs.
run() {
    if [ "${DRYRUN:-0}" = "1" ]; then
        printf '+ %s\n' "$*"
    else
        "$@"
    fi
}

# ============================================================================
# sample_id_of <fastq filename or path>
# ----------------------------------------------------------------------------
# Maps a raw fastq filename to its sample ID. Default: everything before the
# first underscore, which reproduces the K01-F / K01-M / NegC convention. If
# your filenames use a different scheme, override this function here only.
# ============================================================================
sample_id_of() {
    local b
    b="$(basename "$1")"
    printf '%s\n' "${b%%_*}"
}
