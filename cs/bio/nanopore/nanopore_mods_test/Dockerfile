# Dockerfile for nanopore_mods analysis pipeline with SRA sequence integration

FROM condaforge/mambaforge:latest

LABEL maintainer="Josh T. Burdick"
LABEL description="Container for nanopore RNA modification analysis pipeline with SRA data fetching"

# Avoid interactive prompts during apt install
ENV DEBIAN_FRONTEND=noninteractive
ENV PYTHONUNBUFFERED=1

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    git \
    curl \
    wget \
    ca-certificates \
    gzip \
    pigz \
    procps \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# Clone the reference repository
RUN git clone https://github.com/vivian-cheung-lab/nanopore_mods.git /app/nanopore_mods

WORKDIR /app/nanopore_mods

# Create mamba/conda environment with required bioinformatic packages
RUN mamba create -n nanopore_mods_1 -c conda-forge -c bioconda -y \
    python=3.9 \
    snakemake \
    numpy \
    scipy \
    matplotlib \
    seaborn \
    pandas \
    pybedtools \
    pysam \
    samtools \
    tabix \
    sra-tools \
    pysradb \
    ont-fast5-api \
    && mamba clean --all -f -y

# Make sure conda environment is on PATH
ENV PATH=/opt/conda/envs/nanopore_mods_1/bin:$PATH

# Copy updated SRA download script and entrypoint script
COPY get_nanopore_data_sra.py /app/nanopore_mods/get_nanopore_data_sra.py
COPY entrypoint.sh /app/entrypoint.sh

RUN chmod +x /app/nanopore_mods/get_nanopore_data_sra.py /app/entrypoint.sh

# Default mount point for output data persistence
VOLUME /app/nanopore_mods/IVT_AANCR

ENTRYPOINT ["/app/entrypoint.sh"]
CMD ["get-data"]
