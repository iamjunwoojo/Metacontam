# Metacontam v0.0.1
<img width="434" height="144" alt="Metacontam" src="https://github.com/user-attachments/assets/273b5685-b46d-4c9d-ab2e-929d73d7dfec" />

**Metacontam** is a contaminant detection tool for shotgun metagenomic sequencing data.  
Metacontam identifies contaminant species using a blacklist-guided Louvain algorithm and conANI comparison. It is best suited for low-biomass metagenomic samples.

---

## Overview

Metacontam detects microbial contamination through a 12-step pipeline:

```
Kraken2 → Bracken → Count Matrix → Prevalence Threshold
→ Network Analysis (NetCoMi + Louvain) → Genome Retrieval
→ Bowtie2 Alignment → MASH Distance → inStrain → Final Prediction
```

---

## Kraken2 database requirement

Metacontam requires a **self-built** Kraken2 database.  
Pre-built databases do not include raw `library.fna` files, which are needed for genome retrieval.  
In addition, `bracken-build` must be run on the database **before** use:

```bash
# Build standard Kraken2 DB (bacteria + archaea + viral + human + plasmid + UniVec_Core)
kraken2-build --standard --db /path/to/kraken2_db --threads 16

# Run bracken-build — set -l to your actual read length (e.g. 100, 150, 250)
bracken-build -d /path/to/kraken2_db -t 16 -k 35 -l <READ_LENGTH>
```

> **Note**: If download fails due to rsync issues, add `--use-ftp`:
> ```bash
> kraken2-build --standard --db /path/to/kraken2_db --threads 16 --use-ftp
> ```

---

## Installation

### Option 1 — Install with mamba (recommended)

> Fastest way to get started right now — all external tools and Python dependencies are installed automatically via a single environment file.

**Step 1. Install Miniforge (skip if mamba is already installed)**
```bash
wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh
# Restart your terminal after installation
```

**Step 2. Clone and install Metacontam**
```bash
git clone https://github.com/Lifemining-lab/Metacontam.git
cd Metacontam
mamba env create -f environment.yml
conda activate metacontam
pip install .
```

**Step 3. Verify installation**
```bash
metacontam --help
```

---

### Option 2 — Install with conda (coming soon)

> Conda package is currently under review. We will announce when it becomes available.

---

### Option 3 — Install from source (manual)

> Most flexible option — install each external tool yourself. Useful if you already have some tools installed or need fine-grained control, but requires more manual steps.

**Step 1. Install external tools**

Ensure the following tools are installed and available in your `$PATH`:

| Tool | Tested version |
|------|----------------|
| [Kraken2](https://github.com/DerrickWood/kraken2) | 2.1.3 |
| [Bracken](https://github.com/jenniferlu717/Bracken) | 2.8 |
| [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2) | 2.5 |
| [Samtools](http://www.htslib.org/) | 1.10 |
| [MASH](https://github.com/marbl/Mash) | 2.1 |
| [inStrain](https://github.com/MrOlm/inStrain) | 1.9.1 |
| [NetCoMi](https://github.com/stefpeschel/NetCoMi) (R package) | 1.1.0 |

**Step 2. Clone and install Metacontam**
```bash
git clone https://github.com/Lifemining-lab/Metacontam.git
cd Metacontam
mamba env create -f environment_min.yml
conda activate metacontam
pip install instrain
pip install .
```

**Step 3. Verify installation**
```bash
metacontam --help
```

---

## Quick Start

```bash
metacontam \
  --metadata    metadata.tsv \
  --output      ./output \
  --DB          /path/to/kraken2_db \
  --read-length <READ_LENGTH>   # set to your actual read length (e.g. 150)
```

### Metadata format (`metadata.tsv`)
Tab-separated, no header:
```
Samplename_A    sampletype_1    /path/to/SampleA_R1.fastq.gz    /path/to/SampleA_R2.fastq.gz
Samplename_B    sampletype_2    /path/to/SampleB_R1.fastq.gz    /path/to/SampleB_R2.fastq.gz
```

> **Note**: The sample type column is reserved for future batch-aware decontamination and is **not used by the current pipeline**. Any placeholder value (e.g., `sample`) is accepted.

---

## Arguments

### Required
| Argument | Description |
|----------|-------------|
| `--metadata` | Metadata TSV (sample name, sample type, R1 path, R2 path) |
| `--output` | Output directory |
| `--DB` | Kraken2 database directory |
| `--read-length` | Read length for Bracken re-estimation |

### Optional
| Argument | Default | Description |
|----------|---------|-------------|
| `--numcore` | 1 | Number of threads |
| `--min-reads` | 2 | Minimum reads to call a species |
| `--min-cor` | 0.45 | Pearson correlation threshold for network edges |
| `--kraken-report` |  | Pre-computed Kraken2 report directory (skips Kraken2) |
| `--bracken-report` |  | Pre-computed Bracken report directory (skips Bracken) |
| `--dist-matrix` |  | Pre-computed MASH distance matrix (skips MASH) |
| `--candidate-genome` |  | Pre-fetched candidate FASTA (skips genome retrieval) |
| `--bam-dir` |  | Pre-made BAM directory (skips alignment) |

---

## Output

| File / Directory | Description |
|-----------------|-------------|
| `Kraken_dir/` | Kraken2 classification reports |
| `Bracken_dir/` | Bracken species-level reports |
| `Abundance.txt` | Mean abundance per species across samples |
| `Prevalence.txt` | Prevalence per species per sample type |
| `kraken_filtered_matrix.txt` | Filtered count matrix used for network analysis |
| `Network_Output/` | Edge list, Louvain partition, community figure |
| `Genome_dir/Candidate.fasta` | Candidate contaminant genome sequences |
| `Bamfiles/` | Sorted BAM files |
| `mash_sketches/` | MASH sketch files |
| `dist_matrix.txt` | Pairwise MASH distance matrix |
| `pair_output.tsv` | Stratified sample pairs for inStrain |
| `ISfiles/` | inStrain profile outputs |
| `IScompare/` | inStrain compare outputs |
| `merged_IS_compare_Table.tsv` | Merged pairwise ANI comparison table |
| `Final_prediction.txt` | Final contaminant/non-contaminant classification |

---

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## Citation


If you use Metacontam, please cite our preprint:

> Jo J, Lee H, Baek JW, Lee S, Singh V, Shoaie S, Mardinoglu A, Choi J, Lee S. 
> **Metacontam: A Negative Control-Free Decontamination Method for Metagenomic Analysis.** 
> *bioRxiv* 2026.04.26.720876; doi: https://doi.org/10.64898/2026.04.26.720876

> **Note**: This work is currently a preprint and has not yet been peer-reviewed. 
> Citation will be updated once the peer-reviewed version is published.
