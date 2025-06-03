#!/usr/bin/env python3
import argparse
import os
import csv
import shutil

def main():
    parser = argparse.ArgumentParser(
        description="Subset representative genomes (species/strain) from a genome set using metadata."
    )
    parser.add_argument("metadata_tsv", help="Path to metadata TSV file")
    parser.add_argument("all_genomes_dir", help="Directory containing all .fa genome files")
    parser.add_argument("output_dir", help="Directory to copy representative genomes to")
    parser.add_argument("--rep-column", choices=["rep_95id", "rep_99id"],
                        help="Column in metadata to use for representatives (e.g., rep_95id or rep_99id)")
    parser.add_argument("--id-column", default="mag_id", help="Column in metadata with genome file IDs (default: mag_id)")
    parser.add_argument("--dry-run", action="store_true", help="Only print what would be copied, don't actually copy")
    parser.add_argument("--genome-list", help="Optional: Path to file with list of genome IDs or filenames to subset (one per line)")
    args = parser.parse_args()

    # Error if both or neither rep-column and genome-list are provided
    if (args.rep_column and args.genome_list) or (not args.rep_column and not args.genome_list):
        parser.error("You must provide exactly one of --rep-column or --genome-list.")

    # If genome-list is provided, use it
    if args.genome_list:
        genome_list_ids = set()
        with open(args.genome_list) as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.endswith('.fa'):
                    genome_id = line[:-3]
                else:
                    genome_id = line
                genome_list_ids.add(genome_id)
        final_ids = genome_list_ids
        print(f"User provided genome list with {len(final_ids)} entries.")

    # If rep-column is provided, build nonredundant set from that column
    else:
        with open(args.metadata_tsv, newline='') as tsvfile:
            reader = csv.DictReader(tsvfile, delimiter='\t')
            if args.rep_column not in reader.fieldnames:
                raise ValueError(f"Column '{args.rep_column}' not found in metadata.")
            rep_ids = set()
            for row in reader:
                rep_id = row[args.rep_column].strip()
                if rep_id:
                    rep_ids.add(rep_id)
        final_ids = rep_ids
        print(f"Found {len(final_ids)} unique representative genomes in column '{args.rep_column}'.")

    # Prepare output directory
    if not args.dry_run:
        os.makedirs(args.output_dir, exist_ok=True)

    # Copy genome files
    missing = []
    copied = 0
    for genome_id in final_ids:
        fa_name = f"{genome_id}.fa"
        src_path = os.path.join(args.all_genomes_dir, fa_name)
        dest_path = os.path.join(args.output_dir, fa_name)
        if os.path.exists(src_path):
            if args.dry_run:
                print(f"Would copy: {src_path} -> {dest_path}")
            else:
                shutil.copy2(src_path, dest_path)
                print(f"Copied: {src_path} -> {dest_path}")
            copied += 1
        else:
            print(f"WARNING: Genome file not found for: {fa_name}")
            missing.append(fa_name)

    print(f"\nDone. {copied} genomes {'would be' if args.dry_run else 'were'} copied to {args.output_dir}.")
    if missing:
        print(f"{len(missing)} genome files were missing:")
        for m in missing:
            print(f"  {m}")

if __name__ == "__main__":
    main() 