# bee_total_rna_seq_virus_pipeline

A configurable pipeline that processes raw paired-end **total RNA-seq** data to
discover viral contigs in bees (or any host). It runs from raw fastq files
through QC, trimming, rRNA + host depletion, de novo assembly, and a two-pronged
virus search: standard **DIAMOND blastx vs nr** plus **RdRp-scan** to recover
divergent / novel RNA viruses that sequence-similarity search misses.

The pipeline is driven by a **sample sheet** and a single **`config.sh`** — you
should not need to edit any stage script to run a new project. It is written as
PBS (`qsub`) jobs for NCI Gadi and uses relative paths under the repo directory
wherever possible.

---

## Contents

- [Quick start](#quick-start)
- [The sample sheet](#the-sample-sheet-project_rna_samplestxt)
- [Inputs you provide](#inputs-you-provide)
- [Software & conda environments](#software--conda-environments)
- [Stages](#stages)
- [Running it](#running-it-submission--walltime)
- [Outputs](#outputs)
- [Comparing the two virus searches](#comparing-rdrp-scan-vs-diamond-vs-nr-stage-17)
- [Testing without the cluster](#testing-without-the-cluster)
- [Reruns & robustness](#reruns--robustness)
- [Repository layout](#repository-layout)
- [Legacy scripts](#legacy-scripts)

---

## Quick start

1. **Clone / copy the repo** to where you will run it (e.g. `/g/data/rg47/<you>/<project>`).
2. **Put your raw fastqs** in `raw/` (or set `RAW_DIR` in `config.sh`).
3. **Create the sample sheet** `${STUDY}_rna_samples.txt` — one fastq filename
   per line, R1 and R2 on separate lines. Copy `example_rna_samples.txt` as a
   template. With the default `STUDY=UKBombus` the file is
   `UKBombus_rna_samples.txt`.
4. **Edit `config.sh`** — project name, reference indexes, database paths, conda
   env paths, threads. Everything project-specific lives here.
5. **Add `names.txt`** — tab-separated `taxid<TAB>scientific_name`, used to
   annotate the nr hits.
6. **Build the conda envs** (see [Software](#software--conda-environments)) and
   point `config.sh` at them.
7. **Run setup once, then submit the pipeline:**
   ```bash
   qsub pbs/00_setup.pbs          # build results tree + RdRp-scan databases
   ./submit_all.sh                # submit stages 01..17 as a dependency chain
   ```
   …or run / inspect a single stage:
   ```bash
   qsub pbs/06_rrna_bowtie2.pbs
   DRYRUN=1 bash pbs/06_rrna_bowtie2.pbs   # print the commands without running them
   ```

---

## The sample sheet (`${STUDY}_rna_samples.txt`)

One **raw fastq filename per line, R1 and R2 on separate lines**:

```
K01-F_237V23LT3_GCGTGGATGG-TTGACCAATG_L003_R1.fastq.gz
K01-F_237V23LT3_GCGTGGATGG-TTGACCAATG_L003_R2.fastq.gz
K01-F_237V23LT3_GCGTGGATGG-TTGACCAATG_L004_R1.fastq.gz
K01-F_237V23LT3_GCGTGGATGG-TTGACCAATG_L004_R2.fastq.gz
K01-M_237V23LT3_TCTGGTATCC-CGTTGCTTAC_L003_R1.fastq.gz
K01-M_237V23LT3_TCTGGTATCC-CGTTGCTTAC_L003_R2.fastq.gz
NegC_237V23LT3_TCATAGATTG-GATTAAGGTG_L004_R1.fastq.gz
NegC_237V23LT3_TCATAGATTG-GATTAAGGTG_L004_R2.fastq.gz
```

- **Sample ID** = text before the first underscore (`K01-F`, `K01-M`, `NegC`).
  Files that share a sample ID are treated as **lanes** and concatenated into one
  R1/R2 pair after trimming (stage 04). Any number of lanes per sample is fine.
  Override `sample_id_of()` in `config.sh` for a different naming convention.
- **R1 / R2** are detected from the `_R1` / `_R2` token; `_R1_001` style is also
  tolerated.
- Blank lines and `#` comments are ignored. **Bare filenames** are looked up in
  `$RAW_DIR`; lines containing a `/` (or absolute paths) are used as-is.

All sample handling lives in `lib/parse_samples.sh` (`list_samples`,
`r1_files_for`, `r2_files_for`, `assert_pairs`) — **no sample names are
hard-coded in any stage**. `assert_pairs` fails fast if a sample has an unequal
number of R1 and R2 files or an entry whose mate can't be determined.

---

## Inputs you provide

| Input | `config.sh` variable | Notes |
|-------|----------------------|-------|
| Raw paired-end fastqs + sample sheet | `RAW_DIR`, `SAMPLE_SHEET` | total RNA-seq |
| rRNA Bowtie2 index | `RRNA_BT2_INDEX` | index **prefix** (the `-x` value) |
| Host STAR genome | `STAR_GENOME_DIR` | or `STAR_FASTA` + `STAR_GTF` to build with stage 08 |
| `nt` BLAST database | `NT_DB_PATH`, `NT_DB_NAME` | dir holding `nt.*` + base name |
| `nr` DIAMOND database | `NR_DMND` | `.dmnd`; default `/g/data/rg47/mw9045/BLAST/NR_db/nr.dmnd` |
| taxid → name table | `NAMES_TXT` | tab-separated `taxid<TAB>name` |
| RdRp-scan databases | `RDRPSCAN_DIR` | **built automatically** by stage 00 |

---

## Software & conda environments

The stages get their tools from **two conda environments** plus a couple of NCI
Gadi **environment modules**. Set the env paths in `config.sh`:

```bash
CONDA_SH=/g/data/.../miniconda3/etc/profile.d/conda.sh
CONDA_ENV=/g/data/.../envs/Beeviromics      # QC / trim / align / assemble  (stages 01–12)
CONDA_ENV_BLAST=/g/data/.../envs/BLAST       # blast / diamond / R          (stages 13–15)
CONDA_ENV_RDRP=$CONDA_ENV_BLAST              # getorf / hmmer / diamond      (stages 00, 16)
```

### Environment 1 — `Beeviromics` (stages 01–12)

Read QC, trimming, host/rRNA alignment and assembly:

```bash
mamba create -n Beeviromics -c bioconda -c conda-forge \
    fastqc trim-galore cutadapt bowtie2 star megahit pigz
```

### Environment 2 — `BLAST` (stages 00, 13–16)

Homology search, ORF calling, profile search, taxonomy annotation, and the
RdRp-scan DB build. `CONDA_ENV_RDRP` defaults to this env, so it must contain the
EMBOSS / HMMER tools too:

```bash
mamba create -n BLAST -c bioconda -c conda-forge \
    blast diamond hmmer emboss seqkit git \
    r-base r-dplyr r-tidyr r-readr r-argparse
```

| Tool | Used by | Purpose |
|------|---------|---------|
| `blastn` (BLAST+) | 13 | contigs vs `nt` |
| `seqkit` | 13, 16 | split contigs into chunks; dedup ORFs |
| `diamond` | 00, 14, 16 | blastx vs `nr` and vs the RdRp-scan db; `makedb` |
| `Rscript` + `dplyr/tidyr/readr/argparse` | 15 | join nr hits to scientific names |
| `getorf` (EMBOSS) | 16 | translate ORFs under viral genetic codes |
| `hmmsearch` / `hmmpress` (HMMER3) | 00, 16 | RdRp HMM-profile search / press |
| `git` | 00 | clone the RdRp-scan repository |

### Gadi environment modules

FastQC and Bowtie2 are loaded with `module load fastqc` / `module load bowtie2`
inside the relevant stages (NCI Gadi provides these). **If your site has no such
modules**, add `fastqc` and `bowtie2` to the `Beeviromics` env and delete the
`module load …` lines in the affected stages (01, 03, 05, 06, 07, 11).

> The default paths in `config.sh` point at an existing Gadi install under
> `/g/data/rg47/mw9045/…`. Change them to your own envs/DBs before running.

---

## Stages

| Stage | Script | What it does | Env / module |
|------:|--------|--------------|--------------|
| 00 | `pbs/00_setup.pbs` | results tree; clone RdRp-scan; build its DIAMOND + HMM databases | BLAST |
| 01 | `pbs/01_raw_fastqc.pbs` | FastQC on raw reads | module fastqc |
| 02 | `pbs/02_trim_galore.pbs` | Trim Galore adapter/quality trimming (per lane) | Beeviromics |
| 03 | `pbs/03_trimmed_fastqc.pbs` | FastQC on trimmed reads | module fastqc |
| 04 | `pbs/04_concat_lanes.pbs` | Concatenate lanes → one R1/R2 pair per sample | — |
| 05 | `pbs/05_concat_fastqc.pbs` | FastQC on concatenated reads | module fastqc |
| 06 | `pbs/06_rrna_bowtie2.pbs` | rRNA depletion (Bowtie2); keep non-rRNA pairs | Beeviromics + module bowtie2 |
| 07 | `pbs/07_rrna_fastqc.pbs` | FastQC on rRNA-depleted reads | module fastqc |
| 08 | `pbs/08_star_build.pbs` | Build host STAR index (once; auto-skips if present) | Beeviromics |
| 09 | `pbs/09_star_host.pbs` | Host depletion (STAR); keep unmapped reads | Beeviromics |
| 10 | `pbs/10_post_star_rename.pbs` | Add `.fq` extension to STAR unmapped mates | — |
| 11 | `pbs/11_star_fastqc.pbs` | FastQC on host-depleted reads | module fastqc |
| 12 | `pbs/12_megahit.pbs` | MEGAHIT assembly + `<sample>_renamed_contigs.fa` | Beeviromics |
| 13 | `pbs/13_blastn.pbs` | blastn vs `nt`; split into known / **unknown** / possible-viral | BLAST |
| 14 | `pbs/14_blastx_diamond_nr.pbs` | DIAMOND blastx of **unknown** contigs vs `nr` | BLAST |
| 15 | `pbs/15_add_taxid.pbs` | Annotate nr hits with scientific names (`blastx_add_taxid.r`) | BLAST |
| 16 | `pbs/16_rdrpscan.pbs` | **RdRp-scan** on the same unknown contigs (getorf + DIAMOND + hmmsearch) | BLAST |
| 17 | `pbs/17_compare_virus_search.pbs` | Compare nr (14) vs RdRp-scan (16); tabulate novel candidates | — |

### Read flow (what feeds what)

```
raw fastqs
  └─02 trim ─04 concat ─06 rRNA-deplete ─09 host-deplete ─12 assemble
        └─ <sample>_renamed_contigs.fa
              └─13 blastn vs nt
                    ├─ <sample>_possible_viral_contigs.fasta   (viral-kingdom nt hits)
                    └─ <sample>_unknown_contigs.fa  ─┬─14 DIAMOND blastx vs nr ─15 add names
                                                     └─16 RdRp-scan (getorf+diamond+hmmsearch)
                                                          └─17 compare 14 vs 16
```

---

## Running it (submission & walltime)

**`submit_all.sh`** submits **all** stages at once and chains them with
`qsub -W depend=afterok:<prev>`, so each job is held in the queue until the
previous one finishes successfully. A failure leaves the rest **held** rather
than running on bad input.

Gadi's **48 h walltime is per job**, not for the whole chain — each stage is a
separate job, so the full run can span well beyond 48 h. The one caveat: stages
**13/14/16 process all samples in a single job**, which for a large cohort can
approach the per-job limit. Those three therefore also accept **PBS array**
submission, one sample per task (each with its own 48 h):

```bash
# how many samples? (array range = number of unique sample IDs)
source ./config.sh && source lib/parse_samples.sh && list_samples | wc -l
# then submit one task per sample:
qsub -J 1-<N> pbs/13_blastn.pbs
```

Without `-J` the same scripts simply loop over every sample in one job.

---

## Outputs

Everything lands under `results/` (one sub-dir per stage). Highlights:

- `results/12_megahit/<sample>/<sample>_renamed_contigs.fa` — assembled contigs
  (headers prefixed with the sample ID so they are unique across samples).
- `results/13_blastn/<sample>/`
  - `<sample>_unknown_contigs.fa` — contigs with **no** `nt` hit (input to 14 & 16)
  - `<sample>_possible_viral_contigs.fasta`, `<sample>_possible_viral_hits.txt`
  - `<sample>_blast.txt`, `<sample>_sorted_blast.txt`, `<sample>_blast_hits.txt`
- `results/14_diamond_nr/<sample>/<sample>_diamond_nr.tsv` — nr blastx hits.
- `results/15_taxid/<sample>/<sample>_diamond_nr_named.tsv` — nr hits + names.
- `results/16_rdrpscan/<sample>/`
  - `<sample>_rdrpscan_diamond.tsv` — DIAMOND hits vs the RdRp-scan protein db
  - `<sample>_rdrpscan_hmm.tbl` — hmmsearch hits vs the RdRp HMM profiles
  - `<sample>_rdrp_candidate_contigs.fasta` — contigs flagged as candidate RdRps
- `results/17_comparison/` — see below.

### Comparing RdRp-scan vs DIAMOND-vs-nr (stage 17)

Stages 14 and 16 take the **same** blastn-unknown contigs, so stage 17 partitions
them and quantifies what each method recovers. It writes `results/17_comparison/`:

- `summary.tsv` — one row per sample:
  `sample, unknown_contigs, nr_hits, rdrp_diamond, rdrp_hmm, rdrp_total, both, nr_only, rdrp_only`
- `<sample>/<sample>_both.txt`, `_nr_only.txt`, `_rdrp_only.txt` — contig-ID lists
- `<sample>/<sample>_rdrp_only_contigs.fasta` — **candidate novel-virus
  sequences**: flagged by RdRp-scan but with no nr hit.

The `rdrp_only` set is the payoff: RdRp-scan recovers divergent RdRps with as
little as ~10 % identity to known viruses — sequences a plain nr search misses.

---

## Testing without the cluster

You don't need Gadi / the real databases to sanity-check wiring:

- **Per-stage dry run** — print the exact commands a stage would run, for every
  sample, without executing the heavy tools:
  ```bash
  DRYRUN=1 bash pbs/06_rrna_bowtie2.pbs
  ```
  This is the recommended pre-flight check: eyeball one sample's paths against
  your real references before launching the chain.

- **Syntax check everything:**
  ```bash
  for f in config.sh lib/parse_samples.sh submit_all.sh pbs/*.pbs; do bash -n "$f"; done
  ```

The orchestration (sample parsing, lane grouping, the blastn known/unknown/viral
split, and the stage-17 comparison) has been verified end-to-end on synthetic
data with the external tools stubbed out; only the actual tool/DB behaviour needs
the real cluster environment.

---

## Reruns & robustness

- **Idempotent-ish:** every stage `mkdir -p`s its own output dir and can be
  re-submitted. Stage 08 **auto-skips** if the STAR index already exists; stage
  12 clears a sample's MEGAHIT dir before re-assembling; stage 13 skips chunks
  whose blast output already exists.
- **One bad sample won't kill the batch:** in the looping (non-array) form,
  stages 13/14/16 log `WARN: <sample> failed; continuing with next sample` and
  move on, rather than aborting the whole job.
- **Fail-fast inputs:** stages that pair reads call `assert_pairs` first, so a
  malformed sample sheet stops the run immediately with a clear message.
- **Scratch usage:** blastn sorts/temp files go under `SCRATCH_BASE` (config) and
  are cleaned up per sample.

---

## Repository layout

```
config.sh                          # all paths / refs / DBs / tuning (edit this)
example_rna_samples.txt            # sample-sheet template
names.txt                          # taxid<TAB>name table (you provide)
lib/parse_samples.sh               # sample-sheet helpers (sourced by every stage)
blastx_add_taxid.r                 # nr-hit taxonomy annotation (parameterised by --names)
pbs/00_setup.pbs … pbs/17_*.pbs    # pipeline stages
submit_all.sh                      # qsub dependency chain
results/                           # all outputs (created on run)
```

---

## Legacy scripts

`1.2b_submit_rerun_jobs_Version2.sh` and `1.2c_idempotent_rerun_jobs.sh` are the
original hard-coded BLAST re-run helpers (absolute paths, a fixed sample list).
They are kept for reference; the canonical, config-driven blastn is now
`pbs/13_blastn.pbs`.

---

## Credits

RdRp-scan: Charon *et al.* (2022), *Virus Evolution* 8(2):veac082 —
<https://github.com/JustineCharon/RdRp-scan>. Built on FastQC, Trim Galore,
Bowtie2, STAR, MEGAHIT, BLAST+, DIAMOND, seqkit, EMBOSS, HMMER3 and R.
