# DSA4262 Assignment 1 - Task 5: Long-read RNA-Seq Nextflow pipeline

A Nextflow pipeline for processing SG-NEx long-read (Oxford Nanopore) RNA-Seq data, built by extending the example workflow provided in the course (`workflow_longReadRNASeq.nf`, [GoekeLab/sg-nex-data](https://github.com/GoekeLab/sg-nex-data)).

## Pipeline steps

1. **MINIMAP2_ALIGN** - aligns reads to the reference genome with [minimap2](https://github.com/lh3/minimap2), using protocol-specific flags:
   - direct RNA: `-ax splice -uf -k14`
   - cDNA: `-ax splice` (default k-mer size)
2. **SAM_TO_BAM** - sorts and indexes the alignment with [samtools](http://www.htslib.org/)
3. **QC** - computes each sample's mapping rate (`samtools view -c -F 4` / `samtools view -c`, plus a full `samtools flagstat`) and reports PASS/FAIL against a configurable threshold (default 50%). **This is diagnostic only**: `QC` is a side-branch off the sorted BAMs and does not filter or block any sample from reaching `BAMBU` - all samples in the samplesheet are always included in the joint Bambu run. Note also that `samtools view -c [-F 4]` counts alignment records (including secondary/supplementary), not distinct reads; `samtools flagstat`'s own "primary mapped" line is the closer per-read metric and is included in each `*.qc_report.txt`.
4. **BAMBU** - joint transcript discovery and quantification across all samples with [Bambu](https://github.com/GoekeLab/bambu); genome annotations (`--refGtf`) are optional - omitting them runs Bambu in fully de-novo discovery mode

1. `samples.csv` → `MINIMAP2_ALIGN` (per sample, protocol-aware)
2. → `SAM_TO_BAM` (sort + index)
3. From each sorted BAM, two independent things happen:
   - → `QC` (per sample; diagnostic only - does **not** gate anything downstream)
   - → collected across all samples (regardless of QC result) → `BAMBU` (joint discovery + quantification, all samples always included)

## Usage

Samplesheet format (`samples.csv`):
```
sample_id,fastq_path,protocol
A549_directRNA,/path/to/A549_directRNA.fastq.gz,directRNA
K562_cDNA,/path/to/K562_cDNA.fastq.gz,cDNA
```
(`protocol` must be `directRNA` or `cDNA`)

Run with genome annotations (transcript discovery is annotation-guided):
```bash
nextflow run workflow_longReadRNASeq_modified.nf \
  -with-report report.html -resume \
  --samplesheet /path/to/samples.csv \
  --refFa /path/to/genome.fa \
  --refGtf /path/to/annotations.gtf \
  --outdir results/
```

Run without genome annotations (fully de-novo transcript discovery) - simply omit `--refGtf`:
```bash
nextflow run workflow_longReadRNASeq_modified.nf \
  -with-report report.html -resume \
  --samplesheet /path/to/samples.csv \
  --refFa /path/to/genome.fa \
  --outdir results/
```

## Requirements

- [Nextflow](https://www.nextflow.io/) (tested with v26.04.6)
- [minimap2](https://github.com/lh3/minimap2) (tested with v2.26)
- [samtools](http://www.htslib.org/) (tested with v1.19.2)
- [Bambu](https://github.com/GoekeLab/bambu) (tested with v3.12.1, installed via bioconda in a `bambu` conda environment)

## Files in this repository

- `workflow_longReadRNASeq_modified.nf` - the Nextflow pipeline
- `samples.csv` - example samplesheet (4 SG-NEx samples used in this assignment)
- `report_scenario1.html` / `report_scenario2.html` - Nextflow execution reports from the two runs described in the assignment report (with vs. without annotations) - these are the primary evidence for reported runtimes and cache behaviour
- `verification_scripts/` - small R/Python scripts used to independently verify the numbers reported in Tasks 4 and 5:
  - `analyze_bambu.R`, `pick_novel.R`, `verify_task4.R` - re-derive novel transcript counts, top expressed genes, transcripts-per-gene statistics, and the 3 selected novel transcripts' expression, directly from `bambu_output/se.rds`
  - `bed12_junction_analysis.py` - decodes the submitted `reads.bb` track and checks read-level support for Bambu's predicted splice junctions at the 3 novel-transcript loci (Task 4.4)

## Known limitations / portability notes

- The `BAMBU` process's script uses a hardcoded interpreter path (`/home/ubuntu/miniforge3/envs/bambu/bin/Rscript`) matching the exact server this was run on; to run elsewhere, point this at wherever Bambu is installed (or add a `conda`/`container` directive).
- `samples.csv` (as submitted) contains this run's absolute server paths (`/home/ubuntu/workshop/fastq/...`); replace these with paths valid on your own system.
- `refGtfPath` is passed into `BAMBU` as a plain Nextflow `val`, not a staged `path`. This worked correctly for the runs reported here, but Nextflow's caching does not track *content changes* to a `val`-typed file path the way it does for a real `path` input - a future version should stage the GTF as a `path` for more robust caching.
- The BAM/index lists collected across samples are not currently sorted by sample name before being passed to `BAMBU`; the two scenario runs happened to produce the same 4 samples in a different column order in the embedded R script. This does not affect the results shown (samples are matched by name in the output tables, not column position), but a future version should sort the collected channel by sample ID for deterministic output-column order.

## Context

This pipeline was built as part of DSA4262 (NUS) Assignment 1, using long-read Nanopore RNA-Seq data from the [SG-NEx project](https://github.com/GoekeLab/sg-nex-data). See the accompanying assignment report for full methods, results, and discussion.
