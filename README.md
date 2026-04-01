```

     ██████╗ ███╗   ██╗██╗   ██╗██╗  ██╗
    ██╔═══██╗████╗  ██║╚██╗ ██╔╝╚██╗██╔╝
    ██║   ██║██╔██╗ ██║ ╚████╔╝  ╚███╔╝
    ██║   ██║██║╚██╗██║  ╚██╔╝   ██╔██╗
    ╚██████╔╝██║ ╚████║   ██║   ██╔╝ ██╗
     ╚═════╝ ╚═╝  ╚═══╝   ╚═╝   ╚═╝  ╚═╝
     
```

ONYX is an alignment-free method for inferring biological sex from sequencing data using unique k-mers derived from sex chromosomes.

The method identifies sex-specific k-mers from reference genomes and measures their presence in sequencing data without read alignment.

ONYX supports FASTQ, FASTA, and BAM inputs, and works with both XY and ZW sex determination systems.

---

# Features

- Alignment-free biological sex inference
- Fast and memory-efficient method
- Supports FASTQ / multi-FASTA / BAM
- Works with XY and ZW systems
- Preset reference database support

---

# Requirements

ONYX requires the following external tools:

- [KMC](https://github.com/refresh-bio/KMC)
- [seqkit](https://bioinf.shenwei.me/seqkit/)
- [samtools](https://github.com/samtools/samtools/)

Install them via bioconda:

```bash
conda install -c bioconda kmc seqkit samtools
```

Python dependency:

```
pip install tqdm requests
```

---

# Installation

## Installation via Bioconda

```bash
#Conda
#Create a new conda environment & install onyx
conda create -n onyx_env -c conda-forge -c bioconda onyx
conda activate onyx_env
onyx -h
```

```bash
#Mamba
#Create a new conda environment & install onyx
mamba create -n onyx_env -c conda-forge -c bioconda onyx
mamba activate onyx_env
onyx -h
```

## Installation via Anaconda

```bash
#Conda
#Create a new conda environment & install onyx
conda create -n onyx_env -c omics-tools onyx
conda activate onyx_env
onyx -h
```

```bash
#Mamba
#Create a new conda environment & install onyx
mamba create -n onyx_env -c omics-tools onyx
conda activate onyx_env
onyx -h
```
---

# Quick Start

## 0. Get the list of available preset databases
```bash
onyx download-db --list
```

Then, you can get the list.
```bash
Available ONYX preset databases:

human      Homo sapiens (hg38)
chicken    Gallus gallus (bGalGal1)
```

## 1. Download preset database
```bash
#Download the human preset db at the current directory. 
onyx download-db human --out-dir ./
```

The preset human database includes:

```
human_hg38_k33_v1/
 ├ build_info.json
 └ kmc/
```

---

## 2. Classify sequencing data

Example using paired-end FASTQ files:

```bash
wget https://github.com/omics-tools/onyx/releases/download/example/example.tar.gz
tar zxvf example.tar.gz
onyx classify \
  --seqs ./example/human.HG00138.male.n10000.R1.fq.gz ./example/human.HG00138.male.n10000.R2.fq.gz \
  --db human_hg38_k33_v1 \
  --system XY \
  --out human.HG00138.result.tsv
```

---

# Bootstrap confidence estimation

Bootstrap sampling can be used to estimate classification stability.

Example:

```bash
onyx classify \
  --seqs sample_R1.fastq.gz sample_R2.fastq.gz \
  --db human_hg38_k33_v1 \
  --system XY \
  --bootstrap 20 \
  --bootstrap-fraction 0.7 \
  --out result.tsv
```

Parameters:

| Option | Description |
|------|-------------|
| `--bootstrap` | number of bootstrap replicates |
| `--bootstrap-fraction` | subsampling fraction |
| `--bootstrap-seed` | base random seed |

---

# Input formats

| Format | Description | Extensions |
|------|-------------|-------------|
| FASTQ | sequencing reads (single or paired reads) | `.fq`, `.fastq`, `.fq.gz`, `.fastq.gz` |
| BAM | aligned reads (single or paired reads) |`.bam` (with bam.bai) |
| FASTA | multi-FASTA sequences | `.fa`, `.fasta`, `.fa.gz`, `.fasta.gz`  |

ONYX supports several common sequencing file formats as input. All files provided to `--seqs` must have the same format.

#### Example:
```bash
#Single reads
--seqs single.reads.fq.gz
#Paired reads
--seqs paired_R1.fq.gz paired_R2.fq.gz
#FASTA
--seqs input.fa
#BAM
--seqs input.bam
```

---

# Output format (.tsv)

### Standard output columns

| Column | Description |
|------|-------------|
| `sample` | Sample identifier. If `--sample-id` is provided, that value is used; otherwise the first input filename is used. |
| `inputs` | Comma-separated list of input sequencing files used for the analysis. |
| `k` | k-mer size used for the analysis. |
| `KR_hom` | Normalized ratio of homologous sex chromosome k-mer hits. |
| `KR_het` | Normalized ratio of heterologous sex chromosome k-mer hits. |
| `class` | Sex class inferred from `KR_het` relative to the threshold (`HET` or `HOM`). |
| `sex` | Inferred biological sex (e.g. `XX`, `XY`, `ZZ`, `ZW`). |
| `sex_system` | Sex determination system used for classification (`XY` or `ZW`). |

### Additional columns (bootstrap enabled)

If the `--bootstrap` option is used, additional statistics are reported:

| Column | Description |
|------|-------------|
| `ci_low` | Lower bound of the 95% bootstrap confidence interval for `KR_het`. |
| `ci_high` | Upper bound of the 95% bootstrap confidence interval for `KR_het`. |
| `confidence` | Fraction of bootstrap replicates that produced the same classification as the original estimate. |
| `bootstrap_n` | Number of bootstrap iterations performed. |
| `bootstrap_fraction` | Fraction of reads sampled in each bootstrap replicate. |
| `bootstrap_seed` | Base random seed used for bootstrap sampling. |
| `threshold` | Threshold value used for sex classification. |

---

## Build a custom k-mer database

Create a database of unique k-mers from a reference genome. The reference genome must include both autosomes and sex chromosomes. It is recommended that contigs and ALT sequences with unknown chromosomal locations be excluded from the reference sequence.

Please specify the FASTA header ID for --sex_hom (homologous sex chromosome, e.g., chrX or chrZ) and --sex_het (heterologous sex chromosome, e.g., chrY or chrW), respectively.

Example using a reference sequence:

```bash
onyx build \
  --ref ref.fa \
  --sex_hom chrX \
  --sex_het chrY \
  --out custom_db
```

This generates:

```
custom_db/
 ├ build_info.json
 └ kmc/
```

If you want to set thresholds for a custom database, please create it as a preset database.

```bash
onyx build \
  --ref ref.fa \
  --sex_hom chrX \
  --sex_het chrY \
  --out custom_db_preset
  --preset
  --threshold your_KR_het_value
```

---

# Citation

If you use ONYX in your research, please cite the associated publication (to be released).

---

# License

MIT License

---

# Developer
Koji Ishiya  
Sapiens-LEM, Kanazawa University
