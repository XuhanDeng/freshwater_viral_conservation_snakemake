# Freshwater Conservation - Viral Metagenome Pipeline

A Snakemake pipeline for viral metagenome analysis from freshwater samples.

## Overview

This pipeline performs comprehensive viral identification, clustering, and abundance quantification from metagenomic sequencing data. The workflow includes quality control, host removal, assembly, viral identification using multiple tools, ANI-based clustering at both per-sample and cross-sample levels, and detailed abundance estimation.

## Pipeline Workflow

![Pipeline Workflow](pipeline_rulegraph.png)

The workflow diagram shows the complete pipeline process:
1. Initial quality control (fastp) and host removal (T7)
2. Assembly with SPAdes
3. Parallel viral identification using three tools:
   - GeNomad
   - VirSorter2
   - DeepVirFinder
4. Merging and filtering viral predictions
5. Quality assessment with CheckV
6. Clustering and abundance estimation

## Pipeline Steps

### Steps 1-9: Quality Control and Viral Identification
1. Quality control with fastp
2. Host sequence removal (T7 phage) with Bowtie2
3. Metagenomic assembly with SPAdes
4. Sequence filtering and renaming with seqkit
5. Viral identification using:
   - GeNomad
   - VirSorter2
   - DeepVirFinder
   - CheckV

### Steps 10-14: Result Processing
10. Merge viral identification results from multiple tools
11. Extract viral contigs
12. Second CheckV analysis on filtered sequences
13. Quality filtering based on CheckV results
14. Extract final filtered viral sequences

### Steps 15-16: Clustering
15. Per-sample clustering using CheckV ANI-based method (95% ANI, 85% coverage)
16. All-samples combined clustering
17. Extract cluster representative sequences

### Steps 17-19: Abundance Quantification
18. Build Bowtie2 indices for cluster representatives
19. Align reads to clusters (per-sample and all-samples)
20. Calculate abundance metrics (TPM, mean coverage, read counts) using CoverM

## Requirements

### Software Dependencies
- Snakemake
- Conda/Mamba

### Required Databases
All required databases are provided in a single package through Zenodo [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.16737615.svg)](https://doi.org/10.5281/zenodo.16737615):
- CheckV database (v1.5)
- GeNomad database
- T7 reference genome
- Pre-built indices and supporting files

> **Note**: The complete database package (3.2GB) must be downloaded and placed in the `./database` directory before running the pipeline. See the installation section for details.

### Required Tools (installed via conda)
- fastp
- Bowtie2
- SPAdes
- seqkit
- CheckV
- GeNomad
- VirSorter2
- DeepVirFinder
- CoverM
- BLAST

## Installation

1. Clone the repository:
```bash
git clone https://github.com/yourusername/sluge_phage_Extraction_snakemake.git
cd sluge_phage_Extraction_snakemake
```

2. Create conda environments:
```bash
# Create all required conda environments
for env in envs/*.yaml; do
    conda env create -f $env
done
```

3. Download and set up the required databases:
```bash
# Create database directory
mkdir -p database

# Download the database package from Zenodo
wget https://doi.org/10.5281/zenodo.16737615/files/database.zip

# Extract the database files to the correct location
unzip database.zip -d database/
```

The database package includes:
- CheckV database (v1.5)
- GeNomad database
- T7 reference genome
- All necessary index files

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.16737615.svg)](https://doi.org/10.5281/zenodo.16737615)

> **Important**: The database files (3.2GB) must be placed in the `./database` directory for the pipeline to work correctly. The directory structure should match what's provided in the Zenodo archive.

## Usage

1. Configure your samples in `profiles/config.yaml`:
```yaml
samples:
  - sample1
  - sample2
```

2. Place your input files in the `input/` directory following the naming convention:
```
input/
  sample1.R1.raw.fastq.gz
  sample1.R2.raw.fastq.gz
  sample2.R1.raw.fastq.gz
  sample2.R2.raw.fastq.gz
```

3. Run the pipeline:
```bash
# Dry run to check workflow
snakemake --dryrun

# Run the pipeline using multiple cores
snakemake --cores 8

# Run on a SLURM cluster
snakemake --profile profiles/slurm
```

## Output Structure

```
results/
├── 1_fastp/                    # Quality control results
├── 2_T7_removal/              # Host removal results
├── 3_spades_result/           # Assembly results
├── 4_rename_assembly/         # Filtered assemblies
├── 5_genomad/                 # GeNomad results
├── 6_virsorter2/             # VirSorter2 results
├── 7_deepvirfinder/          # DeepVirFinder results
├── 8_python_merge_filter/    # Merged viral predictions
├── 9_checkv/                 # CheckV quality assessment
├── 10_cluster/               # Viral clustering results
│   ├── 1_filtered_seq/        # Filtered viral sequences
│   ├── 2_combined/            # Combined sequences for all-samples clustering
│   ├── 3_clustered/           # Clustering results
│   │   ├── per_sample/        # Per-sample clustering (ANI, BLAST, cluster assignments)
│   │   └── all_samples*       # All-samples clustering files
│   ├── 4_cluster_fasta/       # Cluster representative sequences
│   └── 5_bowtie_index/        # Bowtie2 indices for cluster representatives
└── 11_bowtie2/               # Abundance estimation
    ├── per_sample/            # Per-sample alignment and abundance (BAM, TPM, mean, count)
    └── all_samples/           # All-samples alignment and abundance (BAM, TPM, mean, count)
```

## Configuration

Edit `profiles/config.yaml` to configure:

### Sample Configuration
```yaml
samples:
  - sample1
  - sample2
```

### Resource Allocation
- SLURM account and partition settings
- Memory and CPU requirements per tool
- Runtime limits

### Tool Parameters
- Clustering thresholds (ANI: 95%, coverage: 85%)
- Alignment parameters (identity: 95%, coverage: 90%)
- Quality filtering criteria

### Database Paths
- CheckV database
- GeNomad database
- VirSorter2 database
- DeepVirFinder models

## Scripts

### `script/anicalc.py`
Calculate Average Nucleotide Identity (ANI) from BLAST results. Part of CheckV clustering workflow.

### `script/aniclust.py`
Cluster sequences based on ANI values using single-linkage clustering with customizable thresholds.

### `script/merge_3_viral_identification.py`
Merge and consolidate results from three viral identification tools (DeepVirFinder, GeNomad, VirSorter2).

### `script/checkv_sec_filter.py`
Filter viral sequences based on CheckV quality metrics and tool predictions.

## Key Features

- **Dual Clustering Strategy**: Performs both per-sample and all-samples clustering for comprehensive analysis
- **Multi-tool Validation**: Integrates results from three independent viral identification tools
- **Automated Abundance Quantification**: Calculates TPM, mean coverage, and read counts for all clusters
- **SLURM Integration**: Optimized for high-performance computing clusters
- **Modular Design**: Easy to customize and extend individual steps

## Citations

If you use this pipeline, please cite the following tools:

- **Snakemake**: Köster, J., & Rahmann, S. (2012). Snakemake—a scalable bioinformatics workflow engine. Bioinformatics, 28(19), 2520-2522.
- **CheckV**: Nayfach, S., et al. (2021). CheckV assesses the quality and completeness of metagenome-assembled viral genomes. Nature Biotechnology, 39(5), 578-585.
- **GeNomad**: Camargo, A.P., et al. (2023). Identification of mobile genetic elements with geNomad. Nature Biotechnology.
- **VirSorter2**: Guo, J., et al. (2021). VirSorter2: a multi-classifier, expert-guided approach to detect diverse DNA and RNA viruses. Microbiome, 9(1), 1-13.
- **DeepVirFinder**: Ren, J., et al. (2020). Identifying viruses from metagenomic data using deep learning. Quantitative Biology, 8(1), 64-77.
- **CoverM**: Woodcroft, B. (2022). CoverM: Read coverage calculator. https://github.com/wwood/CoverM
- **SPAdes**: Nurk, S., et al. (2017). metaSPAdes: a new versatile metagenomic assembler. Genome Research, 27(5), 824-834.

## License

This pipeline is provided as-is for research purposes.

## Contact

For questions or issues, please open an issue on GitHub.