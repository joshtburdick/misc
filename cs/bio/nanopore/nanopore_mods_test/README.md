# nanopore_mods Docker Source

This directory contains the Docker configuration and helper scripts to build and run the direct RNA sequencing modification pipeline from [`vivian-cheung-lab/nanopore_mods`](https://github.com/vivian-cheung-lab/nanopore_mods).

## Overview

The container automatically:
1. Clones the `vivian-cheung-lab/nanopore_mods` repository into `/app/nanopore_mods`.
2. Sets up a Conda environment (`nanopore_mods_1`) equipped with Python 3.9, `snakemake`, `pysradb`, NCBI `sra-tools` (`fasterq-dump`, `prefetch`), `pysam`, `pybedtools`, `samtools`, and `tabix`.
3. Includes an updated data retrieval script (`get_nanopore_data_sra.py`) that fetches sequence metadata and sequence runs directly from SRA accession **SRP439844**.

---

## Building the Docker Image

Build the Docker image using:

```bash
docker build -t nanopore_mods:latest .
```

---

## Usage

### 1. Download Data from SRA Only
To download the SRA metadata and sequence FASTQ files into a local directory:

```bash
mkdir -p ./data
docker run --rm -v $(pwd)/data:/app/nanopore_mods/IVT_AANCR nanopore_mods:latest get-data
```

### 2. Run the Full Workflow (Data Download + Snakemake Analysis)
To fetch sequence data and run the modification pipeline across 8 CPU cores:

```bash
mkdir -p ./data
docker run --rm -e CORES=8 -v $(pwd)/data:/app/nanopore_mods/IVT_AANCR nanopore_mods:latest all
```

### 3. Run Interactive Shell inside Container
To inspect the container environment interactively:

```bash
docker run --rm -it -v $(pwd)/data:/app/nanopore_mods/IVT_AANCR nanopore_mods:latest bash
```

---

## Key Files

- [`Dockerfile`](file:///home/josh_t_burdick/git/misc/cs/bio/nanopore/nanopore_mods_test/Dockerfile): Container specification building the environment and checking out `vivian-cheung-lab/nanopore_mods`.
- [`get_nanopore_data_sra.py`](file:///home/josh_t_burdick/git/misc/cs/bio/nanopore/nanopore_mods_test/get_nanopore_data_sra.py): Script utilizing `pysradb` and `fasterq-dump` to download sequence data for study SRP439844 from SRA.
- [`entrypoint.sh`](file:///home/josh_t_burdick/git/misc/cs/bio/nanopore/nanopore_mods_test/entrypoint.sh): Workflow entrypoint script controlling data download and Snakemake execution.
