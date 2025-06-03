#!/usr/bin/env python3
import argparse
import os
import csv

def main():
    parser = argparse.ArgumentParser(description="Check correspondence between metadata and .fa files in a directory.")
    parser.add_argument("metadata_csv", help="Path to metadata TSV file (first column is file ID)")
    parser.add_argument("fasta_dir", help="Directory containing .fa files")
    parser.add_argument("--log", default="fasta_metadata_check.log", help="Path to log file (default: fasta_metadata_check.log)")
    parser.add_argument("--recursive", action="store_true", help="Recursively search for .fa files in fasta_dir")
    parser.add_argument("--move-extra", metavar="DIR", help="If set, move extra .fa files to this directory")
    args = parser.parse_args()

    # Read IDs from metadata
    with open(args.metadata_csv, newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter='\t')
        header = next(reader)
        metadata_ids = set()
        metadata_duplicates = set()
        all_metadata_ids = []
        for row in reader:
            if row and row[0].strip():
                id_ = row[0].strip()
                if id_ in metadata_ids:
                    metadata_duplicates.add(id_)
                metadata_ids.add(id_)
                all_metadata_ids.append(id_)

    # Find .fa files in directory
    fasta_files = set()
    file_paths = dict()
    fasta_duplicates = set()
    all_fasta_ids = []
    if args.recursive:
        for root, dirs, files in os.walk(args.fasta_dir):
            for fname in files:
                if fname.endswith(".fa"):
                    file_id = os.path.splitext(fname)[0]
                    if file_id in fasta_files:
                        fasta_duplicates.add(file_id)
                    fasta_files.add(file_id)
                    file_paths[file_id] = os.path.join(root, fname)
                    all_fasta_ids.append(file_id)
    else:
        for fname in os.listdir(args.fasta_dir):
            if fname.endswith(".fa"):
                file_id = os.path.splitext(fname)[0]
                if file_id in fasta_files:
                    fasta_duplicates.add(file_id)
                fasta_files.add(file_id)
                file_paths[file_id] = os.path.join(args.fasta_dir, fname)
                all_fasta_ids.append(file_id)

    # Check for missing and extra files
    missing = sorted(metadata_ids - fasta_files)
    extra = sorted(fasta_files - metadata_ids)

    # Move extra files if requested
    if args.move_extra and extra:
        os.makedirs(args.move_extra, exist_ok=True)
        for file_id in extra:
            src = file_paths[file_id]
            # Preserve subdirectory structure if recursive
            if args.recursive:
                rel_path = os.path.relpath(src, args.fasta_dir)
                dest = os.path.join(args.move_extra, rel_path)
                os.makedirs(os.path.dirname(dest), exist_ok=True)
            else:
                dest = os.path.join(args.move_extra, os.path.basename(src))
            print(f"Moving extra file: {src} -> {dest}")
            os.rename(src, dest)

    with open(args.log, "w") as log:
        if metadata_duplicates:
            log.write("Duplicate IDs in metadata file:\n")
            for dup in sorted(metadata_duplicates):
                log.write(f"  {dup}\n")
            log.write("\n")
        if fasta_duplicates:
            log.write("Duplicate .fa files (same stem) in directory:\n")
            for dup in sorted(fasta_duplicates):
                log.write(f"  {dup}.fa\n")
            log.write("\n")

        if missing:
            log.write("Files in metadata but missing in directory:\n")
            for m in missing:
                log.write(f"  {m}.fa\n")
        else:
            log.write("No missing files (all metadata entries have corresponding .fa files).\n")

        log.write("\n")

        if extra:
            log.write("Files in directory but not in metadata:\n")
            for e in extra:
                log.write(f"  {e}.fa\n")
        else:
            log.write("No extra files (all .fa files are listed in metadata).\n")

    print(f"Check complete. See log: {args.log}")
    if metadata_duplicates:
        print(f"Duplicate IDs in metadata: {len(metadata_duplicates)}")
    if fasta_duplicates:
        print(f"Duplicate .fa files in directory: {len(fasta_duplicates)}")
    if missing or extra:
        print("Discrepancies found:")
        if missing:
            print(f"  Missing files: {len(missing)}")
        if extra:
            print(f"  Extra files: {len(extra)}")
    elif not (metadata_duplicates or fasta_duplicates):
        print("No discrepancies found.")

if __name__ == "__main__":
    main() 