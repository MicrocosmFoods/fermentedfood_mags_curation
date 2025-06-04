#!/usr/bin/env python3
import argparse
import csv
import re
import ast
from typing import List, Tuple

def extract_accession_prefix(accession: str) -> Tuple[str, int]:
    """Extract the prefix and number from an accession."""
    # Match pattern like SRR1234567 or ERR1234567
    match = re.match(r'([A-Z]+)(\d+)', accession)
    if match:
        return match.group(1), int(match.group(2))
    return accession, 0

def find_sequential_groups(accessions: List[str]) -> List[List[str]]:
    """Split accessions into groups of sequential and non-sequential accessions."""
    if not accessions:
        return []
    groups = []
    current_group = [accessions[0]]
    for prev, curr in zip(accessions, accessions[1:]):
        prev_prefix, prev_num = extract_accession_prefix(prev)
        curr_prefix, curr_num = extract_accession_prefix(curr)
        if prev_prefix == curr_prefix and prev_num + 1 == curr_num:
            current_group.append(curr)
        else:
            groups.append(current_group)
            current_group = [curr]
    groups.append(current_group)
    return groups

def convert_groups_to_ranges(accessions: List[str]) -> str:
    """Convert groups of sequential accessions to ranges, keep others as is."""
    if not accessions:
        return ''
    groups = find_sequential_groups(accessions)
    result = []
    for group in groups:
        if len(group) > 1:
            prefix, start_num = extract_accession_prefix(group[0])
            _, end_num = extract_accession_prefix(group[-1])
            result.append(f"{prefix}{start_num}-{prefix}{end_num}")
        else:
            result.append(group[0])
    return ','.join(result)

def process_metadata(input_file: str, output_file: str):
    """Process the metadata file to convert sequential accessions to ranges and rename mag_id to name."""
    with open(input_file, 'r', newline='') as infile, \
         open(output_file, 'w', newline='') as outfile:
        
        reader = csv.DictReader(infile, delimiter='\t')
        # Rename mag_id to name in the fieldnames
        fieldnames = ["name" if fn == "mag_id" else fn for fn in reader.fieldnames]
        writer = csv.DictWriter(outfile, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        
        for row in reader:
            for col in ['sample_accession', 'run_accession']:
                if col in row and row[col]:
                    # Try to parse as a Python list, fallback to comma split
                    try:
                        accessions = ast.literal_eval(row[col])
                        if isinstance(accessions, str):
                            accessions = [accessions]
                    except Exception:
                        accessions = [acc.strip() for acc in row[col].split(',')]
                    # Sort accessions for proper range detection
                    def acc_key(acc):
                        prefix, num = extract_accession_prefix(acc)
                        return (prefix, num)
                    accessions = sorted(accessions, key=acc_key)
                    row[col] = convert_groups_to_ranges(accessions)
            # Rename mag_id to name in the row
            if 'mag_id' in row:
                row['name'] = row.pop('mag_id')
            writer.writerow(row)

def main():
    parser = argparse.ArgumentParser(
        description='Convert sequential accessions in metadata to ranges and rename mag_id to name.'
    )
    parser.add_argument('input_file', help='Input metadata TSV file')
    parser.add_argument('output_file', help='Output metadata TSV file')
    parser.add_argument('--backup', action='store_true',
                      help='Create a backup of the input file before processing')
    
    args = parser.parse_args()
    
    # Create backup if requested
    if args.backup:
        import shutil
        backup_file = f"{args.input_file}.bak"
        shutil.copy2(args.input_file, backup_file)
        print(f"Created backup: {backup_file}")
    
    # Process the file
    process_metadata(args.input_file, args.output_file)
    print(f"Processing complete. Output written to: {args.output_file}")

if __name__ == '__main__':
    main() 