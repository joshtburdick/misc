#!/bin/bash
set -e

# Activate conda environment
source /opt/conda/etc/profile.d/conda.sh
conda activate nanopore_mods_1

MODE="${1:-get-data}"
CORES="${CORES:-8}"

case "$MODE" in
  get-data|download)
    echo "=== Running SRA sequence data retrieval ==="
    cd /app/nanopore_mods
    python3 get_nanopore_data.py
    ;;
    
  run-snakemake|snakemake)
    echo "=== Running Snakemake pipeline (${CORES} cores) ==="
    cd /app/nanopore_mods/src/seq/nanopore/RNA/polish/f5c/IVT/mods/Snakemake
    snakemake --cores "${CORES}"
    ;;

  all)
    echo "=== Step 1: Downloading Data ==="
    cd /app/nanopore_mods
    python3 get_nanopore_data.py
    
    echo "=== Step 2: Executing Snakemake Workflow (${CORES} cores) ==="
    cd /app/nanopore_mods/src/seq/nanopore/RNA/polish/f5c/IVT/mods/Snakemake
    snakemake --cores "${CORES}"
    ;;

  *)
    exec "$@"
    ;;
esac
