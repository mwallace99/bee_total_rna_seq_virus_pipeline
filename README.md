# bee_total_rna_seq_virus_pipeline

This repository contains a pipeline that processes raw paired-end total RNAseq data to generate potential viral contigs for further metagenome assembly.

## What it does

The pipeline takes raw RNA sequencing reads, performs quality control and trimming, and removes host reads using aligners like Bowtie2 and STAR. The remaining non-host reads are assembled into contigs using MEGAHIT. Finally, the contigs are queried against sequence databases using BLASTN and Diamond BLASTX to identify potential viruses, and their taxonomic information is added.

## Inputs needed

- **Raw RNAseq data**: Paired-end fastq files containing total RNAseq reads.
- **Host reference genomes/indices**: For Bowtie2 and STAR to perform host depletion.
- **`names.txt`**: A required tab-separated text file containing NCBI tax ids and corresponding scientific species names. This is strictly required by the annotation scripts to map BLAST hit taxonomy IDs to readable species names.

## What will be produced

- **Quality-controlled, host-depleted reads**: Clean fastq files ready for assembly.
- **Assembled contigs**: Fasta files containing the de novo assembled contigs (potential viral genomes).
- **Taxonomic annotations**: Tabular output files from BLAST searches, annotated with scientific species names derived from the `names.txt` file, allowing for easy interpretation of the viral taxa present in the samples.
