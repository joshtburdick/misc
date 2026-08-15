#!/usr/bin/env python3
"""
get_nanopore_data_sra.py

Fetches sample metadata and FASTQ sequence data from SRA for accession SRP439844
(associated with the Vivian Cheung Lab nanopore_mods study).
"""

import os
import sys
import shutil
import subprocess
import pandas as pd

try:
    from pysradb import SRAweb
except ImportError:
    print("Error: pysradb is required. Install with `pip install pysradb` or `conda install -c bioconda pysradb`.")
    sys.exit(1)

OUTPUT_DIR = "IVT_AANCR"
SRA_ACCESSION = "SRP439844"
SNAKEMAKE_CSV = "src/seq/nanopore/RNA/polish/f5c/IVT/mods/Snakemake/IVT_samples.csv"

def format_file_name(s):
    """Clean up sample names for safe filesystem filenames."""
    s = str(s)
    s = s.replace(':', '_to_')
    s = s.replace('w/', 'with')
    s = s.replace('%', 'pct')
    s = s.replace('(', '').replace(')', '')
    s = s.replace(' ', '_')
    s = s.replace('/', '_')
    return s

def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    print(f"=== Fetching SRA Metadata for {SRA_ACCESSION} ===")

    # 1. Run clone extraction script if present
    if os.path.exists("./chr19_AANCR_clone_only.py"):
        print("Extracting chr19 clone reference...")
        subprocess.run(["python3", "./chr19_AANCR_clone_only.py"], check=False)

    # 2. Get SRA metadata using pysradb
    sra_client = SRAweb()
    try:
        sample_info = sra_client.sra_metadata(SRA_ACCESSION, detailed=True)
    except Exception as e:
        print(f"Warning: pysradb detailed fetch failed ({e}), falling back to standard fetch...")
        sample_info = sra_client.sra_metadata(SRA_ACCESSION)

    if sample_info is None or sample_info.empty:
        print(f"Error: Could not retrieve SRA metadata for {SRA_ACCESSION}.")
        sys.exit(1)

    if "run_total_bases" in sample_info.columns:
        sample_info["run_total_bases"] = pd.to_numeric(sample_info["run_total_bases"], errors='coerce').fillna(0).astype(int)

    metadata_path = os.path.join(OUTPUT_DIR, f"{SRA_ACCESSION}_sample_info.csv")
    sample_info.to_csv(metadata_path, index=False)
    print(f"Saved metadata to {metadata_path}")

    # 3. Merge with original sample mapping if available
    if os.path.exists(SNAKEMAKE_CSV):
        print("Merging with Snakemake sample configuration...")
        original_samples = pd.read_csv(SNAKEMAKE_CSV)
        if "short_name" in original_samples.columns:
            original_samples["short_name"] = original_samples["short_name"].replace({"Y": "pseudouridine"})
        if "mod_concentration" in original_samples.columns:
            original_samples['mod_concentration'] = original_samples['mod_concentration'].astype(str).str.replace(":", " to ")
            original_samples["mod_concentration"] = original_samples["mod_concentration"].replace({"100%": "1"})

    # 4. Download FASTQ files from SRA for each run accession
    run_accessions = sample_info["run_accession"].dropna().unique()
    print(f"Found {len(run_accessions)} SRA run accessions: {', '.join(run_accessions)}")

    fasterq_dump_bin = shutil.which("fasterq-dump") or shutil.which("fastq-dump")
    
    for run_acc in run_accessions:
        print(f"\n--- Processing run {run_acc} ---")
        expected_fastq = os.path.join(OUTPUT_DIR, f"{run_acc}.fastq")
        expected_fastq_gz = os.path.join(OUTPUT_DIR, f"{run_acc}.fastq.gz")

        if os.path.exists(expected_fastq) or os.path.exists(expected_fastq_gz):
            print(f"Sequence file for {run_acc} already exists. Skipping download.")
            continue

        if fasterq_dump_bin:
            print(f"Downloading sequences using {os.path.basename(fasterq_dump_bin)}...")
            cmd = [fasterq_dump_bin, "--outdir", OUTPUT_DIR, "--split-files", run_acc]
            res = subprocess.run(cmd)
            if res.returncode != 0:
                print(f"Warning: {fasterq_dump_bin} failed for {run_acc}. Trying prefetch...")
                subprocess.run(["prefetch", run_acc], check=False)
                subprocess.run([fasterq_dump_bin, "--outdir", OUTPUT_DIR, run_acc], check=False)
        else:
            print(f"fasterq-dump not found in PATH. Using pysradb download for {run_acc}...")
            sra_client.download(df=sample_info[sample_info["run_accession"] == run_acc], out_dir=OUTPUT_DIR)

    print("\n=== SRA Data Download Complete ===")

if __name__ == "__main__":
    main()
