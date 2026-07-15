# bee_total_rna_seq_virus_pipeline

Finds viral contigs in **total RNA-seq** data from bees (or any host):
QC → trim → rRNA + host depletion → de novo assembly → a two-pronged virus
search (**DIAMOND blastx vs nr** + **RdRp-scan** for divergent/novel RNA
viruses). Driven by one sample sheet and one `config.sh`; no stage script needs
editing. Written as PBS (`qsub`) jobs for NCI Gadi.

- [Quick start](#quick-start) · [Sample sheet](#sample-sheet) · [Inputs](#inputs) ·
  [Conda environments](#conda-environments) · [Databases & taxonomy](#databases--taxonomy) ·
  [Stages](#stages) · [Running](#running) · [Outputs](#outputs) ·
  [Testing](#testing-without-the-cluster) · [Reruns](#reruns--robustness) ·
  [Layout](#repository-layout)

## Quick start

1. Clone the repo to your run location (e.g. `/g/data/$PROJECT/$USER/<run>`).
2. Put raw fastqs in `raw/` (or set `RAW_DIR` in `config.sh`).
3. Create `${STUDY}_rna_samples.txt` (copy `example_rna_samples.txt`): one fastq
   filename per line, R1 and R2 on separate lines.
4. Edit `config.sh`: project name, reference indexes, DB paths, conda env paths,
   threads.
5. Build the [conda envs](#conda-environments) and supply your own
   [`nt`/`nr` databases](#databases--taxonomy) (not shipped).
6. Run:
   ```bash
   qsub pbs/00_setup.pbs   # results tree + RdRp-scan DBs + names.txt
   ./submit_all.sh         # stages 01–17 as a dependency chain
   ```
   Single stage, or dry run (prints commands without running them):
   ```bash
   qsub pbs/06_rrna_bowtie2.pbs
   DRYRUN=1 bash pbs/06_rrna_bowtie2.pbs
   ```

## Sample sheet

`${STUDY}_rna_samples.txt`: one raw fastq filename per line, R1 and R2 on
separate lines:

```
K01-F_237V23LT3_GCGTGGATGG-TTGACCAATG_L003_R1.fastq.gz
K01-F_237V23LT3_GCGTGGATGG-TTGACCAATG_L003_R2.fastq.gz
K01-F_237V23LT3_GCGTGGATGG-TTGACCAATG_L004_R1.fastq.gz
K01-F_237V23LT3_GCGTGGATGG-TTGACCAATG_L004_R2.fastq.gz
K01-M_237V23LT3_TCTGGTATCC-CGTTGCTTAC_L003_R1.fastq.gz
NegC_237V23LT3_TCATAGATTG-GATTAAGGTG_L004_R1.fastq.gz
```

- **Sample ID** = text before the first `_`. Files sharing an ID are treated as
  **lanes** and concatenated into one R1/R2 pair after trimming (stage 04).
  Override `sample_id_of()` in `config.sh` for a different scheme.
- **R1/R2** detected from `_R1`/`_R2` (`_R1_001` tolerated). Blank lines and `#`
  comments ignored. Bare names resolve under `$RAW_DIR`; entries with a `/` are
  used as-is.

Sample handling lives in `lib/parse_samples.sh`; `assert_pairs` fails fast on
unequal R1/R2 counts or an unmatchable mate.

## Inputs

| Input | `config.sh` | Notes |
|-------|-------------|-------|
| Raw fastqs + sample sheet | `RAW_DIR`, `SAMPLE_SHEET` | total RNA-seq |
| rRNA Bowtie2 index | `RRNA_BT2_INDEX` | index **prefix** |
| Host STAR genome | `STAR_GENOME_DIR` | or `STAR_FASTA` + `STAR_GTF` to build (stage 08) |
| `nt` BLAST database | `NT_DB_PATH`, `NT_DB_NAME` | see [Databases](#databases--taxonomy) |
| `nr` DIAMOND database | `NR_DMND` | build with taxonomy, see [Databases](#databases--taxonomy) |
| taxid → name table | `NAMES_TXT` | stage 16 only; auto-built by stage 00 |
| RdRp-scan databases | `RDRPSCAN_DIR` | built automatically by stage 00 |

## Conda environments

Two environments, committed under [`envs/`](envs). Build them anywhere:

```bash
conda env create -f envs/Beeviromics.yml   # stages 01–12: fastqc trim-galore bowtie2 star megahit
conda env create -f envs/BLAST.yml         # stages 00, 13–17: blast diamond hmmer emboss seqkit + R
```

Point `config.sh` at them (`CONDA_ENV_RDRP` defaults to `CONDA_ENV_BLAST`):

```bash
CONDA_SH=$(conda info --base)/etc/profile.d/conda.sh
CONDA_ENV=$(conda info --base)/envs/Beeviromics
CONDA_ENV_BLAST=$(conda info --base)/envs/BLAST
```

- `*.yml` are lean, portable specs (re-solve on any OS). `*.lock.yml` are pinned
  `linux-64` exports for byte-for-byte reproduction.
- **On Gadi:** install Miniconda under `/g/data` (not `/home`) and run
  `conda env create` on a **login node**; compute nodes have no internet. The
  `config.sh` defaults already assume `/g/data/$PROJECT/$USER/miniconda3/envs/…`.
- `fastqc`/`bowtie2` come from `module load` on Gadi; on sites without those
  modules add them to `Beeviromics` and delete the `module load` lines (stages
  01, 03, 05, 06, 07, 11). `git` (stage 00 clone) uses system/`module load git`.

## Databases & taxonomy

**Bring your own `nt` and `nr`**: they are not shipped.

| DB | `config.sh` | Size | Build/get |
|----|-------------|-----:|-----------|
| `nt` | `NT_DB_PATH` + `NT_DB_NAME` | ~700 GB | `update_blastdb.pl --decompress nt` |
| `nr` | `NR_DMND` | ~200 GB `.dmnd` | `diamond makedb` with taxonomy (below) |

On Gadi the `if89` project hosts a maintained `nt` (the default `NT_DB_PATH`), so
you may only need to build `nr`. Stage 13 needs `nt` to fit a hugemem node's page
cache (~700 GB → `hugemem` queue).

**Taxonomy metadata must sit beside the databases**, or the pipeline runs but
silently mislabels hits:

- **`nt`**: put `taxdb.btd` + `taxdb.bti` in `NT_DB_PATH`. Without them `blastn`
  returns `N/A` for `sscinames`/`sskingdoms`, so the stage-13 viral filter keeps
  *every* hit and `*_possible_viral_contigs.fasta` is meaningless.
  ```bash
  cd "$NT_DB_PATH" && wget https://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz && tar -xzf taxdb.tar.gz
  ```
- **`nr`**: build with the taxonomy maps so stage 15 gets `staxids`:
  ```bash
  diamond makedb --in nr.gz -d nr \
      --taxonmap prot.accession2taxid.FULL.gz --taxonnodes nodes.dmp --taxonnames names.dmp
  ```

**`names.txt`** (taxid→name, stage 16 only) is **auto-built by stage 00** from
the NCBI taxdump; nothing to do. It is the only place `names.txt` is used;
blastn names come from `taxdb.*` directly. To drop it entirely, build `nr.dmnd`
with `--taxonnames` and add `sscinames` to stage 15's `--outfmt`, making stage 16
redundant.

## Stages

| Stage | Script | What it does | Env / module |
|------:|--------|--------------|--------------|
| 00 | `00_setup.pbs` | results tree; clone RdRp-scan + build its DBs; build `names.txt` | BLAST · **copyq** (internet) |
| 01 | `01_raw_fastqc.pbs` | FastQC on raw reads | module fastqc |
| 02 | `02_trim_galore.pbs` | Trim Galore trimming (per lane) | Beeviromics |
| 03 | `03_trimmed_fastqc.pbs` | FastQC on trimmed reads | module fastqc |
| 04 | `04_concat_lanes.pbs` | Concatenate lanes → one R1/R2 pair per sample | - |
| 05 | `05_concat_fastqc.pbs` | FastQC on concatenated reads | module fastqc |
| 06 | `06_rrna_bowtie2.pbs` | rRNA depletion (Bowtie2) | Beeviromics + module bowtie2 |
| 07 | `07_rrna_fastqc.pbs` | FastQC on rRNA-depleted reads | module fastqc |
| 08 | `08_star_build.pbs` | Build host STAR index (auto-skips if present) | Beeviromics |
| 09 | `09_star_host.pbs` | Host depletion (STAR); keep unmapped | Beeviromics |
| 10 | `10_post_star_rename.pbs` | Add `.fq` extension to STAR unmapped mates | - |
| 11 | `11_star_fastqc.pbs` | FastQC on host-depleted reads | module fastqc |
| 12 | `12_megahit.pbs` | MEGAHIT assembly → `<sample>_renamed_contigs.fa` | Beeviromics |
| 13 | `13_blastn.pbs` | blastn vs `nt`; split known / **unknown** / possible-viral | BLAST |
| 14 | `14_blastn_add_names.pbs` | *(optional)* back-fill names onto older blastn tables; not auto-chained | BLAST |
| 15 | `15_blastx_diamond_nr.pbs` | DIAMOND blastx of **unknown** contigs vs `nr` | BLAST |
| 16 | `16_add_taxid.pbs` | Annotate nr hits with names (`blastx_add_taxid.r`) | BLAST |
| 17 | `17_rdrpscan.pbs` | **RdRp-scan** on the same unknown contigs | BLAST |
| 18 | `18_compare_virus_search.pbs` | Compare nr (15) vs RdRp-scan (17) | - |

```
raw fastqs
  └─02 trim ─04 concat ─06 rRNA-deplete ─09 host-deplete ─12 assemble
        └─ <sample>_renamed_contigs.fa
              └─13 blastn vs nt  (+14 optional: name the hits)
                    ├─ <sample>_possible_viral_contigs.fasta   (viral-kingdom nt hits)
                    └─ <sample>_unknown_contigs.fa  ─┬─15 DIAMOND blastx vs nr ─16 add names
                                                     └─17 RdRp-scan ─18 compare 15 vs 17
```

## Running

`submit_all.sh` submits all stages chained with `qsub -W depend=afterok:<prev>`,
each held until the previous succeeds; a failure leaves the rest held.

Gadi's 48 h walltime is **per job**, so the whole chain can exceed it. The heavy
stages **13/15/17** are submitted as **grouped PBS arrays**: `G` subjobs, each
striding over a slice of the samples. Gadi caps an array at **10 elements**, so
`submit_all.sh` uses `G = min(BLAST_NGROUPS, samples, 10)` (default
`BLAST_NGROUPS=5`). To submit one by hand (`G ≤ 10`):

```bash
qsub -r y -J 1-5 -v NGROUPS=5,PROJECT,USER pbs/13_blastn.pbs
```

Arrays go out `-r y` (rerunnable), and each stage writes a `.done` marker per
sample and skips samples already marked, so a timed-out or re-submitted array
resumes where it left off (`rm <sample>/.done` or set `FORCE=1` to redo one).
Without `-J`, a stage loops over all samples in one job.

## Outputs

Everything lands under `results/`, one sub-dir per stage:

- `12_megahit/<sample>/<sample>_renamed_contigs.fa`: contigs (headers prefixed
  with sample ID, unique across samples).
- `13_blastn/<sample>/`: `<sample>_unknown_contigs.fa` (no `nt` hit; input to
  15/17), `<sample>_possible_viral_contigs.fasta`, `<sample>_blast.txt`.
- `15_diamond_nr/<sample>/<sample>_diamond_nr.tsv`: nr blastx hits.
- `16_taxid/<sample>/<sample>_diamond_nr_named.tsv`: nr hits + names.
- `17_rdrpscan/<sample>/`: `_rdrpscan_diamond.tsv`, `_rdrpscan_hmm.tbl`,
  `_rdrp_candidate_contigs.fasta`.
- `18_comparison/`: `summary.tsv` (per sample: unknown_contigs, nr_hits,
  rdrp_diamond, rdrp_hmm, rdrp_total, both, nr_only, rdrp_only) and, per sample,
  `_both.txt` / `_nr_only.txt` / `_rdrp_only.txt` plus
  `_rdrp_only_contigs.fasta`.

`rdrp_only` is the payoff: RdRp-scan recovers divergent RdRps (~10 % identity to
known viruses) that a plain nr search misses.

## Testing without the cluster

- **Dry run** a stage (prints commands, runs no tools):
  `DRYRUN=1 bash pbs/06_rrna_bowtie2.pbs`
- **Syntax check:**
  `for f in config.sh lib/parse_samples.sh submit_all.sh pbs/*.pbs; do bash -n "$f"; done`

The orchestration is verified end-to-end on synthetic data with the tools
stubbed out; only real tool/DB behaviour needs the cluster.

## Reruns & robustness

- **Re-submittable:** every stage `mkdir -p`s its output dir. Stage 08 skips an
  existing STAR index; 12 clears a sample's dir before re-assembling.
- **Resumable arrays:** stages 13/15/17 write a `.done` marker per finished
  sample and skip it on re-run, so a timed-out or re-submitted array picks up
  where it left off (`rm <sample>/.done` or `FORCE=1` to redo). Stage 13 also
  skips individual chunks whose output already exists.
- **One bad sample won't kill the batch:** 13/15/17 log a `WARN` and continue.
- **Fail-fast inputs:** read-pairing stages call `assert_pairs` first.
- **Scratch:** BLAST sort/temp files go under `SCRATCH_BASE` and are cleaned per
  sample (stages 13/15/17 list `scratch/$PROJECT` in their `-l storage`).

## Repository layout

```
config.sh                          # all paths / refs / DBs / tuning (edit this)
example_rna_samples.txt            # sample-sheet template
envs/*.yml, *.lock.yml             # conda env specs (conda env create -f …)
names.txt                          # taxid<TAB>name (auto-built by stage 00; stage 16 only)
lib/parse_samples.sh               # sample-sheet helpers
blastx_add_taxid.r                 # nr-hit taxonomy annotation
pbs/00_setup.pbs … 18_*.pbs        # pipeline stages
submit_all.sh                      # qsub dependency chain
results/                           # all outputs (created on run)
```

[![DOI](https://zenodo.org/badge/1266813709.svg)](https://doi.org/10.5281/zenodo.21366317)

## Credits

Thanks to James Damayo ([@joimes-d](https://github.com/joimes-d)) who provided
the basis for this pipeline.

RdRp-scan: Charon *et al.* (2022), *Virus Evolution* 8(2):veac082,
<https://github.com/JustineCharon/RdRp-scan>. Built on FastQC, Trim Galore,
Bowtie2, STAR, MEGAHIT, BLAST+, DIAMOND, seqkit, EMBOSS, HMMER3 and R.
