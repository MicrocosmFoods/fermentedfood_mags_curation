#!/usr/bin/env python3
import argparse
import os
import shutil
from pathlib import Path
import csv

def create_batches(input_dir, batch_size, output_base, samplesheet_dir):
    """
    Create batches of FASTA files from input directory.
    Also creates a sample sheet CSV for each batch.
    
    Args:
        input_dir (str): Path to directory containing FASTA files
        batch_size (int): Number of files per batch
        output_base (str): Base directory for output batches
        samplesheet_dir (str): Base directory for output sample sheets
    """
    # Get all FASTA files recursively
    fasta_files = list(Path(input_dir).rglob("*.fa"))
    
    # Separate GCA/GCF and other files
    gca_gcf_files = [f for f in fasta_files if f.stem.startswith("GCA") or f.stem.startswith("GCF")]
    other_files = [f for f in fasta_files if not (f.stem.startswith("GCA") or f.stem.startswith("GCF"))]
    
    # Helper to create batches from a list of files
    def make_batches(files, batch_prefix, start_batch_num):
        batches = []
        for i in range(0, len(files), batch_size):
            batch_files = files[i:i+batch_size]
            batches.append((f"{batch_prefix}_{start_batch_num}", batch_files))
            start_batch_num += 1
        return batches, start_batch_num
    
    batch_num = 1
    all_batches = []
    
    # GCA/GCF batches
    gca_gcf_batches, batch_num = make_batches(gca_gcf_files, "batch", batch_num)
    all_batches.extend(gca_gcf_batches)
    
    # Other batches
    other_batches, batch_num = make_batches(other_files, "batch", batch_num)
    all_batches.extend(other_batches)
    
    # Create samplesheet directory
    os.makedirs(samplesheet_dir, exist_ok=True)
    
    # Create batches and samplesheets
    for batch_name, batch_files in all_batches:
        batch_dir = os.path.join(output_base, batch_name)
        os.makedirs(batch_dir, exist_ok=True)
        
        # Copy files to batch directory, preserving structure
        for file_path in batch_files:
            rel_path = file_path.relative_to(input_dir)
            dest_dir = os.path.join(batch_dir, rel_path.parent)
            os.makedirs(dest_dir, exist_ok=True)
            shutil.copy2(file_path, os.path.join(batch_dir, rel_path))
        
        print(f"Created {batch_name} with {len(batch_files)} files")
        
        # Create the batch sample sheet
        samplesheet_path = os.path.join(samplesheet_dir, f"{batch_name}.csv")
        with open(samplesheet_path, "w", newline="") as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(["Data type: assembly; Columns: 4; Version: 1"])
            writer.writerow(["staging_file_subdir_path", "type", "min_contig_length", "assembly_name"])
            writer.writerow(["FASTA file path", "Assembly type", "Minimum contig length", "Assembly object name"])
            for file_path in batch_files:
                rel_path = file_path.relative_to(input_dir)
                if file_path.stem.startswith("GCA") or file_path.stem.startswith("GCF"):
                    assembly_type = "Draft Isolate"
                else:
                    assembly_type = "Metagenome"
                writer.writerow([
                    str(rel_path),
                    assembly_type,
                    "500",
                    file_path.stem
                ])

def main():
    parser = argparse.ArgumentParser(description='Split FASTA files into batches')
    parser.add_argument('input_dir', help='Input directory containing FASTA files')
    parser.add_argument('--batch-size', type=int, default=500,
                      help='Number of files per batch (default: 500)')
    parser.add_argument('--output-dir', default='fasta_batches',
                      help='Base directory for output batches (default: fasta_batches)')
    parser.add_argument('--samplesheet-dir', default='batch_samplesheets',
                      help='Directory for batch sample sheets (default: batch_samplesheets)')
    
    args = parser.parse_args()
    
    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)
    os.makedirs(args.samplesheet_dir, exist_ok=True)
    
    create_batches(args.input_dir, args.batch_size, args.output_dir, args.samplesheet_dir)

if __name__ == '__main__':
    main() 