# Freshwater Viral Conservation - Comprehensive Metagenome Pipeline

A robust, scalable Snakemake pipeline for viral identification, clustering, and abundance quantification from freshwater metagenomic samples. This pipeline integrates multiple state-of-the-art viral discovery tools with advanced clustering algorithms to provide comprehensive viral community analysis.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.16737615.svg)](https://doi.org/10.5281/zenodo.16737615)
[![Snakemake](https://img.shields.io/badge/snakemake-≥7.0-brightgreen.svg)](https://snakemake.readthedocs.io)
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

## Table of Contents

- [Overview](#overview)
- [Pipeline Workflow](#pipeline-workflow)
- [Detailed Pipeline Steps](#detailed-pipeline-steps)
- [Key Features](#key-features)
- [Requirements](#requirements)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Configuration Guide](#configuration-guide)
- [Output Structure](#output-structure)
- [Advanced Usage](#advanced-usage)
- [Troubleshooting](#troubleshooting)
- [Performance Optimization](#performance-optimization)
- [Citations](#citations)
- [Version History](#version-history)

## Overview

This pipeline performs comprehensive viral metagenome analysis from raw sequencing reads through to final abundance quantification and taxonomic annotation. The workflow is designed for freshwater environmental samples but can be adapted for other sample types.

### Key Capabilities

- **Multi-tool viral identification**: Integrates GeNomad, VirSorter2, and DeepVirFinder with union-based approach for comprehensive viral detection
- **ANI-based clustering**: Uses CheckV's clustering method for accurate viral population delineation
- **Dual-level analysis**: Performs both per-sample and cross-sample clustering
- **Spike-in normalization**: Maps reads to both T7 internal control and viral contigs for normalization
- **Comprehensive abundance metrics**: Calculates TPM, mean coverage, and read counts
- **SLURM-optimized**: Ready for high-performance computing environments
- **Fully reproducible**: Conda/Mamba-managed dependencies

## Pipeline Workflow

![Pipeline Workflow](pipeline_rulegraph.png)

The complete pipeline consists of 21 integrated steps organized into five major phases:

### Phase 1: Quality Control (Steps 1-2)
- Adapter trimming and quality filtering with fastp
- T7 phage spike-in control removal for assembly with Bowtie2

### Phase 2: Assembly (Steps 3-4)
- Metagenomic assembly with SPAdes metaSPAdes mode
- Contig filtering (≥5kb) and systematic renaming

### Phase 3: Viral Identification (Steps 5-9)
- Parallel identification using three complementary tools
- Consensus-based merging of predictions
- Quality assessment with CheckV

### Phase 4: Clustering (Steps 10-17)
- ANI-based clustering at 95% identity, 85% coverage
- Representative sequence selection
- Combined T7+viral index building

### Phase 5: Quantification (Steps 18-21)
- Read alignment to combined references
- Multi-metric abundance calculation
- Taxonomic annotation with GeNomad

## Detailed Pipeline Steps

### Phase 1: Quality Control and Preprocessing

#### Step 1: Quality Control with fastp
**Tool**: fastp v0.23.4
**Purpose**: Remove adapters, trim low-quality bases, and filter short reads

```yaml
Parameters:
  - qualified_quality_phred: 20 (Phred score threshold)
  - length_required: 20 (Minimum read length)
  - trim_poly_g: true (Remove poly-G tails from NovaSeq data)
  - trim_poly_x: true (Remove poly-X tails)
```

**Outputs**:
- Paired-end clean reads: `{sample}_1P.fq.gz`, `{sample}_2P.fq.gz`
- Unpaired reads: `{sample}_U1.fq.gz`, `{sample}_U2.fq.gz`
- QC reports: HTML and JSON formats

#### Step 2: T7 Phage Spike-in Control Removal (for Assembly)
**Tool**: Bowtie2 v2.5.0
**Purpose**: Remove reads mapping to T7 phage spike-in control before assembly

**Background**: T7 phage is added as an **internal spike-in control** during wet lab sample processing to:
- Monitor DNA extraction efficiency
- Assess sequencing depth consistency
- Provide normalization standard across samples
- Quality control for the experimental workflow

**Note**: T7 reads are removed before assembly to avoid assembling the control sequence, but are **retained in fastp-processed reads** used for abundance quantification (see Step 18).

```yaml
Parameters:
  - very-sensitive: true
  - f2 flag: Extract both reads unmapped (flag 12)
  - F256: Exclude non-primary alignments
```

**Reference**: T7 phage genome (GCF_000844825.1)

### Phase 2: Assembly and Filtering

#### Step 3: Metagenomic Assembly
**Tool**: SPAdes v3.15.5 (metaSPAdes mode)
**Purpose**: Assemble metagenomic reads into contigs

```yaml
Parameters:
  - mode: --meta (metagenomic mode)
  - k-mer sizes: 21,33,55,77,101,121
  - memory: 256 GB
  - threads: 16
  - only-assembler: true (skip read correction)
```

**Output**: `scaffolds.fasta` with assembled contigs

#### Step 4: Sequence Filtering and Renaming
**Tool**: seqkit v2.3.1
**Purpose**: Filter contigs ≥5kb and rename with systematic IDs

```yaml
Parameters:
  - min_length: 5000 bp
  - naming scheme: {sample}_{number:010d}
  - nr_width: 10 (number width with leading zeros)
```

**Rationale**: 5kb minimum captures complete viral genomes while reducing false positives

### Phase 3: Viral Identification

#### Step 5: GeNomad Identification
**Tool**: GeNomad v1.5.2
**Purpose**: ML-based identification of viruses and plasmids

```yaml
Mode: end-to-end
Features:
  - Gene calling and annotation
  - Viral/plasmid classification
  - Taxonomic assignment
  - Protein annotation
Database: geNomad database v1.5
```

**Key outputs**:
- `*_virus_summary.tsv`: Viral predictions with scores
- `*_virus.fna`: Identified viral sequences
- `*_genes.tsv`: Annotated genes

#### Step 6: VirSorter2 Identification
**Tool**: VirSorter2 v2.2.4
**Purpose**: Neural network-based viral classification

```yaml
Parameters:
  - min_length: 5000 bp
  - min_score: 0.5
  - groups: dsDNAphage, NCLDV, lavidaviridae
  - keep-original-seq: true
```

**Algorithm**: Multi-classifier approach with viral hallmark genes

**Key outputs**:
- `final-viral-score.tsv`: Scores for all predictions
- `final-viral-combined.fa`: Identified sequences

#### Step 7: DeepVirFinder Identification
**Tool**: DeepVirFinder v1.0
**Purpose**: Deep learning-based viral sequence detection

```yaml
Parameters:
  - min_length: 3000 bp (uses k-mer based approach)
  - model: Pre-trained CNN model
  - score_cutoff: 0.9 (high confidence)
  - pvalue_cutoff: 0.05
```

**Algorithm**: Convolutional neural network trained on k-mer frequencies

**Key output**: `*_dvfpred.txt` with predictions and p-values

#### Step 8: Initial CheckV Analysis
**Tool**: CheckV v1.0.1
**Purpose**: Quality assessment of viral contigs before filtering

```yaml
Mode: end-to-end
Assessments:
  - Completeness estimation
  - Contamination detection
  - Provirus boundary detection
  - Host region identification
```

### Phase 4: Result Merging and Quality Control

#### Step 9: Merge Viral Predictions
**Script**: `merge_3_viral_identification.py`
**Purpose**: Union-based integration of three tool predictions

```yaml
Strategy:
  - UNION approach: Include sequences identified by ANY tool
  - Apply tool-specific quality thresholds:
    - DeepVirFinder: score ≥ 0.8 AND p-value < 0.05
    - VirSorter2: max_score ≥ 0.8 AND full sequences only (not partial)
    - GeNomad: Exclude proviruses (topology != "Provirus")
  - No weighting or confidence scoring
  - All predictions merged into single list
```

**Merging approach**:
1. Filter each tool's predictions using quality thresholds
2. Take UNION of all filtered sequences
3. Label each sequence with which tool(s) detected it (tag1=dvf, tag2=vs2, tag3=genomad)
4. All sequences pass to next step (no filtering by number of tools)

**Output**:
- `*_merged_results.csv`: All predictions with tool labels (columns: name, tag1, tag2, tag3)
- `*_merge3_list.txt`: List of all unique sequence IDs from any tool

#### Step 10: Extract Viral Contigs
**Tool**: seqkit grep
**Purpose**: Extract sequences identified as viral by ANY tool

**Input**: Merged viral ID list (union of all three tools)
**Output**: `*_filterd_vircontig.fasta` (all sequences from union)

#### Step 11: Second CheckV Analysis
**Tool**: CheckV v1.0.1
**Purpose**: Detailed quality assessment of merged viral contigs

**Metrics assessed**:
- **Completeness**: Percentage of complete genome
- **Contamination**: Non-viral sequence content
- **Quality tiers**: Complete, High-quality, Medium-quality, Low-quality
- **Provirus detection**: Integrated prophages
- **Host genes**: Potential host contamination

#### Step 12: Quality-based Filtering
**Script**: `checkv_sec_filter.py`
**Purpose**: Apply hallmark gene-based filtering using CheckV and original tool predictions

**Filtering strategy**: UNION of three filters (sequences passing ANY filter are kept)

```yaml
Filter 1 - VirSorter2 re-check:
  - hallmark > 0 (at least 1 hallmark gene)
  - full/partial = "full" (complete sequences only)

Filter 2 - GeNomad re-check:
  - topology != "Provirus" (exclude proviruses)
  - n_hallmarks > 0 (at least 1 hallmark gene)

Filter 3 - CheckV quality:
  - viral_genes > 0 (at least 1 viral gene detected)
```

**Approach**:
1. Re-check original VirSorter2 predictions for hallmark genes
2. Re-check original GeNomad predictions for hallmark genes
3. Check CheckV viral gene count
4. Take UNION of sequences passing any filter
5. Final list includes sequences with viral hallmark evidence

**Rationale**: Ensures retained sequences have viral hallmark genes or viral gene signatures, reducing false positives while maintaining broad sensitivity.

#### Step 13: Extract Final Viral Sequences
**Output**: High-quality viral contigs with hallmark gene support for downstream analysis

### Phase 5: Clustering and Representative Selection

#### Step 14: Per-Sample Clustering
**Method**: CheckV ANI-based clustering
**Purpose**: Dereplicate viral genomes within each sample

```yaml
Parameters:
  - min_ani: 95% (species-level threshold)
  - min_tcov: 85% (alignment coverage of target)
  - min_qcov: 0% (alignment coverage of query - no minimum)
Algorithm:
  1. All-vs-all BLAST (megablast)
  2. ANI calculation from BLAST results
  3. Single-linkage clustering
```

**Outputs**:
- `*_blast.tsv`: All-vs-all BLAST results
- `*_ani.tsv`: Pairwise ANI values
- `*_cluster.tsv`: Cluster assignments (representative, members)

#### Step 15: All-Samples Clustering
**Method**: Same as per-sample, applied to combined sequences
**Purpose**: Identify viral populations shared across samples

**Process**:
1. Concatenate all per-sample viral sequences
2. All-vs-all BLAST of combined dataset
3. ANI-based clustering at 95%/85%
4. Assign cluster representatives

**Advantage**: Enables cross-sample viral population tracking

#### Step 16: Extract Cluster Representatives
**Critical fix**: Use column 1 (representative ID), not column 2 (all members)

**Cluster file format**:
```
Column 1: Representative sequence ID (one per cluster)
Column 2: All member sequences (comma-separated)
```

**Command**: `cut -f 1 cluster.tsv | seqkit grep -f - input.fasta`

### Phase 6: Index Building with T7 Integration

#### Step 17a & 17b: Build Combined Bowtie2 Indices

**Innovation**: Combine T7 reference with viral cluster representatives

**Process**:
```bash
# Concatenate T7 reference with cluster representatives
cat T7_reference.fna cluster_representatives.fasta > combined_with_T7.fasta

# Build Bowtie2 index
bowtie2-build combined_with_T7.fasta cluster_index
```

**Rationale**:
- Simultaneously quantify T7 spike-in control and viral contigs
- T7 abundance serves as internal standard for cross-sample normalization
- Monitor T7 recovery to assess DNA extraction efficiency
- Single alignment step for both control and target sequences

**Outputs**:
- Per-sample indices: One per sample with sample-specific representatives
- All-samples index: Single index with cross-sample representatives

### Phase 7: Abundance Quantification

#### Step 18a & 18b: Read Alignment

**Tool**: Bowtie2 v2.5.0
**Input**: **fastp-processed reads** (retain T7 sequences)
**Reference**: Combined T7+viral indices

```yaml
Parameters:
  - mode: --sensitive
  - threads: 8
  - flags: -f 2 (properly paired reads only)
Post-processing:
  - samtools view: Convert SAM to BAM
  - samtools sort: Sort by coordinate
  - samtools index: Create BAI index
```

**Why fastp reads (not T7-removed reads)?**
- **Includes T7 sequences**: Essential for quantifying the spike-in control
- **Internal standard**: T7 abundance enables normalization across samples
- **Quality controlled**: Adapter-trimmed and quality-filtered
- **Experimental control**: T7 levels indicate DNA extraction/sequencing consistency

#### Step 19a & 19b: Abundance Calculation

**Tool**: CoverM v0.6.1
**Purpose**: Calculate multiple abundance metrics

**Metrics**:

1. **TPM (Transcripts Per Million)**:
   ```
   TPM = (mapped_reads / contig_length) × 1,000,000 / sum(mapped_reads / length)
   ```
   - Normalized for contig length and sequencing depth
   - Comparable across samples

2. **Mean coverage**:
   ```
   Mean = total_bases_aligned / contig_length
   ```
   - Average depth across contig
   - Indicates abundance and evenness

3. **Read count**:
   ```
   Count = number of reads mapped to contig
   ```
   - Raw abundance measure
   - Useful for low-abundance viruses

**CoverM parameters**:
```yaml
min_read_percent_identity: 95% (ANI threshold)
min_read_aligned_percent: 90% (read coverage)
contig_end_exclusion: 0 (use full contig)
no_zeros: true (exclude zero-abundance contigs)
```

#### Step 20: Taxonomic Annotation

**Tool**: GeNomad annotate
**Purpose**: Assign taxonomy and functional annotation to cluster representatives

**Annotations provided**:
- Taxonomic lineage (family, genus)
- Protein functions
- Marker genes
- Host predictions
- Lifestyle (lytic vs lysogenic)

**Outputs**:
- Per-sample annotations: `results/12_taxonomy/genomad/per_sample/{sample}/`
- All-samples annotations: `results/12_taxonomy/genomad/all_samples/`

## Key Features

### 1. Combined T7+Viral Index Strategy
**Innovation**: Integrates T7 phage reference with cluster representatives for comprehensive abundance profiling

**Benefits**:
- Simultaneous quantification of spike-in control and viral targets
- T7 as internal standard for cross-sample normalization
- Quality control through T7 recovery monitoring
- Assess DNA extraction and sequencing efficiency
- Enables accurate abundance comparison across samples

### 2. Multi-Tool Viral Identification
**Approach**: Union-based identification using three complementary tools with quality thresholds

**Strategy**:
- **Broad sensitivity**: UNION approach includes sequences detected by ANY tool
- **Quality filtering**: Each tool has stringent thresholds (score ≥ 0.8, p < 0.05)
- **Complementary detection**: Each tool has different strengths
  - GeNomad: Gene-based, good for novel viruses, taxonomic assignment
  - VirSorter2: Hallmark gene detection, group-specific models, provirus detection
  - DeepVirFinder: K-mer based, sequence composition, fast screening
- **Comprehensive capture**: Maximizes viral discovery by combining all tool outputs
- **Tool tagging**: Track which tool(s) detected each sequence for downstream analysis

### 3. Dual Clustering Strategy
**Two-level analysis**: Per-sample AND cross-sample clustering

**Per-sample clustering**:
- Identifies within-sample viral diversity
- Reduces computational burden
- Sample-specific viral populations

**Cross-sample clustering**:
- Tracks viral populations across samples
- Identifies core vs. variable viruses
- Temporal/spatial dynamics analysis

### 4. Corrected Representative Selection
**Fix**: Extract cluster representatives (column 1), not all members (column 2)

**Impact**:
- One representative per cluster (correct)
- Reduces index size dramatically
- Improves mapping accuracy
- Faster alignment

### 5. Quality-Controlled Input for Alignment
**Approach**: Use fastp-processed reads rather than T7-removed reads

**Rationale**:
- Retain T7 sequences for spike-in control quantification
- High-quality reads (adapter-trimmed, quality-filtered)
- Enables normalization using T7 as internal standard
- Consistent input for both control and viral mapping

### 6. SLURM-Optimized Workflow
**Features**:
- Resource specifications per rule
- Partition selection (compute vs. memory)
- Automatic job submission
- Parallel execution up to 50 jobs
- Runtime limits and memory allocation

## Requirements

### System Requirements

**Minimum**:
- CPU: 16 cores
- RAM: 64 GB
- Storage: 500 GB (for databases + results)
- OS: Linux (tested on Ubuntu 20.04, CentOS 7)

**Recommended** (for SLURM cluster):
- CPU: 32+ cores
- RAM: 256 GB (especially for assembly)
- Storage: 1 TB SSD
- Scheduler: SLURM workload manager

### Software Dependencies

**Required**:
- Snakemake ≥7.0
- Conda or Mamba (Mamba recommended for faster solving)
- Git

**Python version**: 3.9+

### Required Databases

All databases are available as a single package through Zenodo:

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.16737615.svg)](https://doi.org/10.5281/zenodo.16737615)

**Database package contents** (3.2 GB total):
1. **CheckV database v1.5** (1.2 GB)
   - Viral HMM profiles
   - Genome database
   - DTR (direct terminal repeat) database

2. **GeNomad database** (1.8 GB)
   - Marker protein database
   - Taxonomic profiles
   - Neural network models

3. **T7 reference genome** (40 KB)
   - GCF_000844825.1 (T7 bacteriophage)

4. **VirSorter2 database** (200 MB)
   - Downloaded separately during setup

### Tool Versions (Conda-managed)

```yaml
fastp: 0.23.4
bowtie2: 2.5.0
spades: 3.15.5
seqkit: 2.3.1
checkv: 1.0.1
genomad: 1.5.2
virsorter2: 2.2.4
deepvirfinder: 1.0
coverm: 0.6.1
blast: 2.13.0
samtools: 1.17
python: 3.9
```

## Installation

### Step 1: Clone Repository

```bash
git clone https://github.com/XuhanDeng/freshwater_viral_conservation_snakemake.git
cd freshwater_viral_conservation_snakemake
```

### Step 2: Install Snakemake

**Option A: Install Snakemake with Mamba (Recommended - fastest)**
```bash
# Install Mamba first (if not already installed)
conda install -c conda-forge mamba

# Create Snakemake environment with Mamba
mamba create -c conda-forge -c bioconda -n snakemake snakemake
```

**Option B: Install Snakemake with Conda**
```bash
conda create -c conda-forge -c bioconda -n snakemake snakemake
```

**Activate Snakemake environment**:
```bash
conda activate snakemake
```

**Note**: Individual tool environments (fastp, bowtie2, spades, etc.) will be **automatically created** by Snakemake when you run the pipeline with `--use-conda` flag. You don't need to manually create them!

### Step 3: Download and Extract Databases

```bash
# Create database directory
mkdir -p database && cd database

# Download from Zenodo (3.2 GB)
wget https://zenodo.org/records/16737615/files/database.zip

# Extract
unzip database.zip

# Clean up
rm database.zip

# Return to main directory
cd ..
```

**Expected directory structure**:
```
database/
├── checkv-db-v1.5/
│   ├── genome_db/
│   ├── hmm_db/
│   └── README.txt
├── genomad_db/
│   ├── genomad_db/
│   ├── genomad_db_proteins.dmnd
│   └── version.txt
└── T7_ref/
    └── GCF_000844825.1_ViralProj14460_genomic.fna
```

### Step 4: Set Up VirSorter2 Database (Optional - can be done later)

VirSorter2 requires a separate database download (~10 GB). You can either:

**Option A: Set up now (recommended)**
```bash
# Snakemake will auto-create vs2 environment when needed
# Manually activate it to download database
conda activate snakemake
snakemake --snakefile Snakefile --use-conda --until virsorter2_identification --cores 1

# OR manually set up if you prefer
conda create -c bioconda -n vs2 virsorter=2
conda activate vs2
virsorter setup -d database/virsorter2-db -j 4
conda deactivate
```

**Option B: Let Snakemake handle it automatically**
- Skip this step
- VirSorter2 will download its database on first run
- This happens automatically when the rule executes

### Step 5: Test Installation

```bash
# Make sure Snakemake environment is activated
conda activate snakemake

# Dry run to check pipeline setup
snakemake --snakefile Snakefile --dry-run --cores 1 --use-conda

# This will show the execution plan without running anything
```

## Quick Start

### 1. Prepare Input Data

**Directory structure**:
```
input/
├── sample1.R1.raw.fastq.gz
├── sample1.R2.raw.fastq.gz
├── sample2.R1.raw.fastq.gz
└── sample2.R2.raw.fastq.gz
```

**Naming convention**: `{sample}.R1.raw.fastq.gz` and `{sample}.R2.raw.fastq.gz`

### 2. Configure Samples

Edit `profiles/config.yaml`:

```yaml
samples:
  - sample1
  - sample2
  - sample3
```

### 3. Run Pipeline

**Local execution** (single machine):
```bash
# Test with dry-run
snakemake --snakefile Snakefile --dry-run --cores 8

# Run pipeline
snakemake --snakefile Snakefile --cores 8 --use-conda
```

**SLURM cluster execution**:
```bash
# Configure SLURM settings in profiles/config.yaml first
snakemake --snakefile Snakefile --profile profiles/ --jobs 50
```

**Partial workflow** (from clustering onwards):
```bash
# If you have existing viral sequences
snakemake --snakefile Snakefile_runed.smk --cores 8 --use-conda
```

## Configuration Guide

### Sample Configuration

Edit `profiles/config.yaml`:

```yaml
samples:
  - day0_2
  - day14_4_1
  - day3_20_2
  # Comment out samples to skip
  # - day0_1  # This sample will be skipped
```

### Resource Allocation

#### General Resources
```yaml
runtime: 2880  # Maximum runtime in minutes (48 hours)
slurm_account: "your-slurm-account"
regular_memory: 3900  # MB per CPU for regular tasks
regular_partition: "compute"
high_memory: 16000  # MB per CPU for memory-intensive tasks
high_partition: "memory"
```

#### Per-Tool Thread Configuration
```yaml
fastp:
  threads: 8
  qualified_quality_phred: 20
  length_required: 20

bowtie2_build:
  threads: 2

bowtie2_removal:
  threads: 4

spades:
  threads: 16
  memory: 256  # GB total
  k_values: "21,33,55,77,101,121"

genomad:
  threads: 16

virsorter2:
  threads: 16
  min_length: 5000
  min_score: 0.5
  groups: "dsDNAphage,NCLDV,lavidaviridae"

checkv:
  threads: 4

coverm:
  threads: 4
  runtime: 1440
  mem_mb_per_cpu: 3900
  partition: "compute"
```

### Tool Parameters

#### Clustering Parameters
```yaml
# In aniclust.py script call
min_ani: 95      # Minimum ANI for clustering (species-level)
min_tcov: 85     # Minimum target coverage (alignment breadth)
min_qcov: 0      # Minimum query coverage
```

**Adjust for different taxonomic levels**:
- Species-level: 95% ANI, 85% coverage
- Genus-level: 70-80% ANI
- Family-level: 50-60% ANI

#### Alignment Parameters
```yaml
# In CoverM abundance calculation
min_read_percent_identity: 95
min_read_aligned_percent: 90
contig_end_exclusion: 0
```

#### Quality Filtering
```yaml
seqkit:
  min_length: 5000  # Minimum contig length (bp)
  nr_width: 10      # Number width in sequence IDs
```

### Database Paths

```yaml
bowtie2_index: "database/T7_ref/GCF_000844825.1_ViralProj14460_genomic.fna"
checkv_db: "database/checkv-db-v1.5/"
genomad_db: "database/genomad_db"
virsorter2_database: "database/virsorter2-db"

# Script paths
merge_script: "script/merge_3_viral_identification.py"
checkv_filter_script: "script/checkv_sec_filter.py"
checkv_scripts: "script/checkv"
```

### SLURM Profile Configuration

Edit `profiles/config.yaml` for SLURM:

```yaml
executor: slurm
jobs: 50  # Maximum parallel jobs
use-conda: true
conda-prefix: conda_envs
rerun-incomplete: true
printshellcmds: true
keep-going: true

default-resources:
  - slurm_account="your-account"
  - slurm_partition="compute"
  - runtime=7200
  - mem_mb_per_cpu=3900
  - cpus_per_task=8
```

## Output Structure

```
results/
├── 1_fastp/                    # Quality control results
│   └── {sample}/
│       ├── {sample}_1P.fq.gz   # Paired R1
│       ├── {sample}_2P.fq.gz   # Paired R2
│       ├── {sample}_U1.fq.gz   # Unpaired R1
│       ├── {sample}_U2.fq.gz   # Unpaired R2
│       ├── {sample}.fastp.html # QC report
│       └── {sample}.fastp.json # QC metrics
│
├── 2_T7_removal/               # T7 spike-in control removal results
│   └── 5_removed_sequence/
│       └── {sample}/
│           ├── {sample}_host_removed_R1.fastq.gz
│           └── {sample}_host_removed_R2.fastq.gz
│
├── 3_spades_result/            # Assembly results
│   └── {sample}/
│       └── {sample}_no_correction/
│           ├── scaffolds.fasta
│           ├── contigs.fasta
│           └── assembly_graph.fastg
│
├── 4_rename_assembly/          # Filtered assemblies (≥5kb)
│   └── rename_5000/
│       └── {sample}_scaffolds_rename_5000.fasta
│
├── 5_genomad/                  # GeNomad results
│   └── {sample}/
│       └── {sample}_scaffolds_rename_5000_summary/
│           ├── {sample}_scaffolds_rename_5000_virus_summary.tsv
│           ├── {sample}_scaffolds_rename_5000_virus.fna
│           └── {sample}_scaffolds_rename_5000_virus_genes.tsv
│
├── 6_virsorter2/               # VirSorter2 results
│   └── {sample}/
│       ├── final-viral-score.tsv
│       ├── final-viral-combined.fa
│       └── final-viral-boundary.tsv
│
├── 7_deepvirfinder/            # DeepVirFinder results
│   └── {sample}/
│       └── {sample}_scaffolds_rename_5000.fasta_gt3000bp_dvfpred.txt
│
├── 8_python_merge_filter/      # Merged viral predictions
│   ├── csv/
│   │   └── {sample}_merged_results.csv
│   ├── list/
│   │   └── {sample}_merge3_list.txt
│   └── seqkit_filter_vircontig/
│       └── {sample}_filterd_vircontig.fasta
│
├── 9_checkv/                   # CheckV quality assessment
│   ├── checkv_sample_figure/   # Initial CheckV
│   │   └── {sample}/
│   │       ├── quality_summary.tsv
│   │       ├── completeness.tsv
│   │       └── contamination.tsv
│   └── chekv_result/           # Second CheckV
│       └── {sample}/
│           └── quality_summary.tsv
│
├── 10_cluster/                 # Viral clustering results
│   ├── 1_filtered_seq/         # Final filtered sequences
│   │   └── {sample}_filterd_final_vircontig.fasta
│   │
│   ├── 2_combined/             # Combined for cross-sample clustering
│   │   └── all_samples_combined.fasta
│   │
│   ├── 3_clustered/            # Clustering results
│   │   ├── per_sample/
│   │   │   └── {sample}/
│   │   │       ├── {sample}_blast.tsv       # BLAST results
│   │   │       ├── {sample}_ani.tsv         # ANI calculations
│   │   │       ├── {sample}_cluster.tsv     # Cluster assignments
│   │   │       └── {sample}_cluster_representatives.fasta
│   │   │
│   │   ├── all_samples_blast.tsv
│   │   ├── all_samples_ani.tsv
│   │   ├── all_samples_cluster.tsv
│   │   └── all_samples_cluster_representatives.fasta
│   │
│   ├── 4_cluster_fasta/        # Cluster representatives
│   │   ├── {sample}_cluster_representatives.fasta
│   │   └── all_samples_cluster_representatives.fasta
│   │
│   └── 5_bowtie_index/         # Combined indices (T7 + representatives)
│       ├── per_sample/
│       │   └── {sample}/
│       │       ├── combined_with_T7.fasta
│       │       ├── cluster_index.1.bt2
│       │       ├── cluster_index.2.bt2
│       │       └── cluster_index.done
│       │
│       └── all_samples/
│           ├── combined_with_T7.fasta
│           └── cluster_index.*.bt2
│
├── 11_bowtie2/                 # Abundance estimation
│   ├── per_sample/
│   │   └── {sample}/
│   │       ├── {sample}.f2.sorted.bam
│   │       ├── {sample}.f2.sorted.bam.bai
│   │       ├── {sample}_TPM.tsv
│   │       ├── {sample}_mean.tsv
│   │       └── {sample}_count.tsv
│   │
│   └── all_samples/
│       └── {sample}/
│           ├── {sample}.f2.sorted.bam
│           ├── {sample}.f2.sorted.bam.bai
│           ├── {sample}_TPM.tsv
│           ├── {sample}_mean.tsv
│           └── {sample}_count.tsv
│
└── 12_taxonomy/                # Taxonomic annotation
    └── genomad/
        ├── per_sample/
        │   └── {sample}/
        │       ├── {sample}_cluster_representatives_summary/
        │       │   ├── {sample}_cluster_representatives_virus_summary.tsv
        │       │   └── {sample}_cluster_representatives_virus_taxonomy.tsv
        │       └── {sample}_cluster_representatives_annotate/
        │           └── {sample}_cluster_representatives_genes.tsv
        │
        └── all_samples/
            └── all_samples_cluster_representatives_summary/
                ├── all_samples_cluster_representatives_virus_summary.tsv
                └── all_samples_cluster_representatives_virus_taxonomy.tsv
```

### Key Output Files Explained

#### Abundance Files
- **TPM.tsv**: Transcripts per million - normalized for length and depth
- **mean.tsv**: Mean coverage depth across each contig
- **count.tsv**: Raw read counts mapping to each contig

Format:
```
Contig ID                       Sample_name
day0_2_0000000001              1523.45
day0_2_0000000002              342.12
```

#### Clustering Files
- **cluster.tsv**: Two-column format
  - Column 1: Representative sequence ID
  - Column 2: All member sequences (comma-separated)

```
Representative              Members
day0_2_0000000001          day0_2_0000000001,day0_2_0000000025,day0_2_0000000089
day0_2_0000000002          day0_2_0000000002,day0_2_0000000456
```

#### CheckV Quality Summary
Key columns:
- `contig_id`: Sequence identifier
- `checkv_quality`: Complete/High-quality/Medium-quality/Low-quality/Not-determined
- `completeness`: Estimated percentage (0-100%)
- `contamination`: Percentage of non-viral content
- `viral_genes`: Number of viral genes detected
- `host_genes`: Number of host genes (contamination indicator)

## Advanced Usage

### Running Specific Steps Only

**Example 1: Run only clustering and abundance (skip upstream steps)**
```bash
# Use Snakefile_runed.smk
snakemake --snakefile Snakefile_runed.smk \
    --cores 8 \
    --use-conda \
    --forceall
```

**Example 2: Re-run only abundance calculation**
```bash
# Delete abundance outputs
rm -rf results/11_bowtie2/

# Re-run
snakemake --snakefile Snakefile \
    --cores 8 \
    --use-conda
```

**Example 3: Run for subset of samples**
```bash
# Edit config.yaml to comment out samples
# Then run
snakemake --snakefile Snakefile \
    --cores 8 \
    --use-conda \
    --rerun-incomplete
```

### Generating Workflow Diagrams

```bash
# Rule graph (high-level overview)
snakemake --snakefile Snakefile --rulegraph | dot -Tpng > rulegraph.png

# DAG (detailed with all jobs)
snakemake --snakefile Snakefile --dag | dot -Tpng > dag.png

# File graph (shows input/output relationships)
snakemake --snakefile Snakefile --filegraph | dot -Tpng > filegraph.png
```

### Resume After Failure

Snakemake automatically resumes from where it left off:

```bash
# Just rerun the same command
snakemake --snakefile Snakefile --cores 8 --use-conda
```

**Force rerun specific rules**:
```bash
snakemake --snakefile Snakefile \
    --cores 8 \
    --use-conda \
    --forcerun genomad_annotate_all_samples
```

### Cluster Execution with Different Parameters

**Test on single sample first**:
```yaml
# config.yaml
samples:
  - test_sample  # Only one sample
```

```bash
snakemake --snakefile Snakefile --cores 8 --use-conda
```

**Scale up to all samples**:
```yaml
samples:
  - sample1
  - sample2
  # ... all samples
```

```bash
snakemake --snakefile Snakefile --profile profiles/ --jobs 50
```

## Troubleshooting

### Common Issues and Solutions

#### Issue 1: MissingInputException for cluster representative files
```
MissingInputException: Missing input files for rule extract_representatives_all_samples:
    results/10_cluster/1_filtered_seq/sample1_filterd_final_vircontig.fasta
```

**Solution**:
1. Check that upstream rules completed successfully
2. Verify file existence: `ls results/10_cluster/1_filtered_seq/`
3. If files missing, re-run from viral identification:
```bash
snakemake --snakefile Snakefile --cores 8 --use-conda --forcerun merge_viral_results
```

#### Issue 2: Empty abundance output files
```
Results in *_TPM.tsv, *_mean.tsv, *_count.tsv are empty or have only header
```

**Causes**:
- No reads mapped to viral contigs
- Bowtie2 index built incorrectly
- CoverM parameters too stringent

**Solutions**:
1. Check BAM file has alignments:
```bash
samtools view -c results/11_bowtie2/all_samples/sample1/sample1.f2.sorted.bam
```
2. Verify combined index includes T7 and viral sequences:
```bash
grep "^>" results/10_cluster/5_bowtie_index/all_samples/combined_with_T7.fasta | wc -l
```
3. Relax CoverM parameters:
```yaml
min_read_percent_identity: 90  # Instead of 95
min_read_aligned_percent: 75   # Instead of 90
```

#### Issue 3: `cut -f 2` extracting wrong sequences
```
Error: Sequence IDs not found when extracting representatives
```

**Solution**: This should already be fixed in the current version
- Verify using `cut -f 1` not `cut -f 2`
- Check in rules: `extract_representatives_per_sample` and `extract_representatives_all_samples`
```bash
grep "cut -f" Snakefile | grep extract_representatives
# Should show: cut -f 1
```

#### Issue 4: SLURM jobs killed for exceeding memory
```
slurmstepd: error: Detected 1 oom-kill event(s) in step
```

**Solution**:
1. Increase memory allocation in config.yaml:
```yaml
spades:
  memory: 512  # Increase from 256 GB

# Or for specific rule
high_memory: 32000  # Increase from 16000 MB per CPU
```

2. Use memory partition:
```yaml
resources:
  slurm_partition: "memory"  # Instead of "compute"
```

#### Issue 5: Conda environment solving takes forever
```
Solving environment: / (hangs)
```

**Solution**:
1. Use Mamba instead of Conda:
```bash
conda install -c conda-forge mamba
```

2. Create environments with Mamba:
```bash
mamba env create -f envs/genomad.yaml
```

3. Clear conda cache:
```bash
conda clean --all
```

#### Issue 6: GeNomad database not found
```
Error: GeNomad database not found at specified path
```

**Solution**:
1. Verify database path in config.yaml:
```bash
ls -la database/genomad_db
```

2. Download if missing:
```bash
genomad download-database database/
```

#### Issue 7: VirSorter2 fails with database error
```
[Errno 2] No such file or directory: 'db/group/**/model'
```

**Solution**:
1. Set up VirSorter2 database:
```bash
conda activate vs2
virsorter setup -d database/virsorter2-db -j 4
conda deactivate
```

2. Update path in config.yaml:
```yaml
virsorter2:
  database: "database/virsorter2-db"
```

#### Issue 8: SPAdes assembly fails with "Not enough memory"
```
== Error == system call for: "[python path]/spades_init.py" finished abnormally
```

**Solution**:
1. Reduce k-mer sizes:
```yaml
spades:
  k_values: "21,33,55,77"  # Remove larger k-mers
```

2. Increase memory:
```yaml
spades:
  memory: 512
```

3. Use --meta flag (should already be enabled)

#### Issue 9: CheckV completes but produces no output
```
CheckV finished but quality_summary.tsv is empty
```

**Causes**:
- Input fasta file is empty
- Sequences don't pass minimum length
- Database issues

**Solutions**:
1. Check input has sequences:
```bash
grep -c "^>" results/8_python_merge_filter/seqkit_filter_vircontig/sample1_filterd_vircontig.fasta
```

2. Verify CheckV database:
```bash
ls database/checkv-db-v1.5/genome_db/
```

3. Check CheckV log:
```bash
cat log/6_checkv/sample1_checkv.log
```

### Performance Issues

#### Slow BLAST all-vs-all for large datasets
**Problem**: All-vs-all BLAST takes hours for 1000+ sequences

**Solutions**:
1. Increase threads for clustering:
```yaml
# In cluster rules
threads: 32  # Instead of 16
```

2. Pre-filter sequences before clustering:
```yaml
seqkit:
  min_length: 10000  # Increase from 5000 to reduce dataset
```

3. Use `blastn` instead of `megablast` for higher sensitivity but slower speed
- Edit in Snakefile clustering rules

#### Multiple samples slow down pipeline
**Problem**: Running 100+ samples sequentially

**Solution**: Use SLURM parallel execution
```bash
snakemake --snakefile Snakefile \
    --profile profiles/ \
    --jobs 100  # Allow up to 100 parallel jobs
```

## Performance Optimization

### Recommended Resource Allocation

**For local machine (16 cores, 64 GB RAM)**:
```yaml
fastp: 4 threads
bowtie2: 4 threads
spades: 8 threads, 48 GB memory
genomad: 8 threads
virsorter2: 8 threads
checkv: 4 threads
clustering: 8 threads
```

**For HPC cluster (unlimited resources)**:
```yaml
fastp: 8 threads
bowtie2: 8 threads
spades: 16 threads, 256 GB memory
genomad: 16 threads
virsorter2: 16 threads
checkv: 8 threads
clustering: 16 threads
```

### Speed Optimizations

1. **Use Mamba**: 10x faster conda environment solving
2. **SSD storage**: Store databases on SSD for faster I/O
3. **Parallel jobs**: Set `--jobs` to number of samples for maximum parallelism
4. **Pre-filter**: Remove short contigs early to reduce downstream processing
5. **Checkpoint resume**: Snakemake automatically resumes, no need to restart

### Expected Runtime

**Per sample** (assuming 10 GB input, 16 cores):
- fastp: 20 minutes
- T7 removal: 30 minutes
- SPAdes assembly: 4-12 hours (varies greatly)
- Viral identification (3 tools): 2-6 hours combined
- CheckV: 30 minutes
- Clustering: 1-4 hours (depends on sequence count)
- Abundance: 30 minutes

**Total per sample**: ~12-24 hours

**For 20 samples on SLURM** (parallel execution): ~24-48 hours total

## Citations

If you use this pipeline, please cite the following tools:

### Pipeline and Workflow
- **Snakemake**: Mölder, F., et al. (2021). Sustainable data analysis with Snakemake. F1000Research, 10:33.

### Quality Control and Preprocessing
- **fastp**: Chen, S., et al. (2018). fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics, 34(17), i884-i890.
- **Bowtie2**: Langmead, B., & Salzberg, S. L. (2012). Fast gapped-read alignment with Bowtie 2. Nature Methods, 9(4), 357-359.

### Assembly
- **SPAdes**: Nurk, S., et al. (2017). metaSPAdes: a new versatile metagenomic assembler. Genome Research, 27(5), 824-834.

### Viral Identification
- **CheckV**: Nayfach, S., et al. (2021). CheckV assesses the quality and completeness of metagenome-assembled viral genomes. Nature Biotechnology, 39(5), 578-585.
- **GeNomad**: Camargo, A.P., et al. (2023). Identification of mobile genetic elements with geNomad. Nature Biotechnology, 41, 1303-1312.
- **VirSorter2**: Guo, J., et al. (2021). VirSorter2: a multi-classifier, expert-guided approach to detect diverse DNA and RNA viruses. Microbiome, 9(1), 1-13.
- **DeepVirFinder**: Ren, J., et al. (2020). Identifying viruses from metagenomic data using deep learning. Quantitative Biology, 8(1), 64-77.

### Abundance Quantification
- **CoverM**: Woodcroft, B. (2022). CoverM: Read coverage calculator for metagenomics. GitHub repository: https://github.com/wwood/CoverM
- **samtools**: Danecek, P., et al. (2021). Twelve years of SAMtools and BCFtools. GigaScience, 10(2), giab008.

### Sequence Analysis
- **seqkit**: Shen, W., et al. (2016). SeqKit: A cross-platform and ultrafast toolkit for FASTA/Q file manipulation. PLoS ONE, 11(10), e0163962.
- **BLAST**: Camacho, C., et al. (2009). BLAST+: architecture and applications. BMC Bioinformatics, 10, 421.

## Version History

### v1.1 (2025-10-20) - Current
**Major improvements:**
- ✅ Added combined T7+cluster representative index building
- ✅ Corrected cluster representative extraction (column 1 instead of column 2)
- ✅ Updated alignment to use fastp-processed reads
- ✅ Added GeNomad taxonomic annotation step
- ✅ Implemented partial workflow (Snakefile_runed.smk)
- ✅ Comprehensive README with troubleshooting guide
- ✅ Sample configuration with active/commented samples

**Bug fixes:**
- Fixed extract_representatives rules extracting all members instead of representatives
- Fixed bowtie2_alignment using T7-removed reads instead of fastp reads
- Fixed rule all to match active uncommented rules

**Pipeline changes:**
- `build_cluster_index_per_sample`: Now combines T7+representatives
- `build_cluster_index_all_samples`: Now combines T7+representatives
- `extract_representatives_per_sample`: Corrected to column 1
- `extract_representatives_all_samples`: Corrected to column 1
- `bowtie2_alignment_per_sample`: Uses fastp output
- `bowtie2_alignment_all_samples`: Uses fastp output

### v1.0 (2025-10-13) - Initial Release
- Complete viral identification pipeline
- Dual clustering strategy (per-sample and all-samples)
- Integrated abundance quantification
- SLURM cluster support
- Multi-tool viral validation

## License

This pipeline is provided under the MIT License for research and educational purposes.

## Contact

**Issues and Questions**:
- Open an issue on GitHub: https://github.com/XuhanDeng/freshwater_viral_conservation_snakemake/issues
- Provide detailed error messages and log files for faster troubleshooting

**Repository**:
https://github.com/XuhanDeng/freshwater_viral_conservation_snakemake

## Acknowledgments

This pipeline was developed for freshwater viral conservation research. We thank:
- The developers of all integrated bioinformatics tools
- The Snakemake development team
- The conda/bioconda communities
- The viral metagenomics research community

**Funding**: [Add your funding information here]

**Contributors**:
- Xuhan Deng - Pipeline development and maintenance
- [Add other contributors]

---

**Last Updated**: October 20, 2025
**Pipeline Version**: 1.1
**Documentation Version**: 1.1
