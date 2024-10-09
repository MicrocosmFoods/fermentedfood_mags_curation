import argparse
from Bio import SeqIO
import os
import sys

def parse_stb(stb_file):
    contig_to_genome = {}
    try:
        with open(stb_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                try:
                    contig, genome = line.strip().split('\t')
                    contig_to_genome[contig] = genome
                except ValueError:
                    print(f"Warning: Invalid format in STB file at line {line_num}: {line.strip()}")
        print(f"Parsed {len(contig_to_genome)} entries from STB file")
    except Exception as e:
        print(f"Error parsing STB file: {e}")
        sys.exit(1)
    return contig_to_genome

def split_fasta(fasta_file, stb_file, output_dir):
    # Parse the STB file
    contig_to_genome = parse_stb(stb_file)

    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Dictionary to store sequences for each genome
    genome_sequences = {}

    # Parse the FASTA file and split sequences
    try:
        fasta_records = list(SeqIO.parse(fasta_file, "fasta"))
        print(f"Parsed {len(fasta_records)} sequences from FASTA file")
        
        for record in fasta_records:
            contig_id = record.id
            if contig_id in contig_to_genome:
                genome = contig_to_genome[contig_id]
                if genome not in genome_sequences:
                    genome_sequences[genome] = []
                genome_sequences[genome].append(record)
            else:
                print(f"Warning: Contig {contig_id} not found in STB file")
        
        print(f"Grouped sequences into {len(genome_sequences)} genomes")
    except Exception as e:
        print(f"Error parsing FASTA file: {e}")
        sys.exit(1)

    # Write sequences to individual genome files
    for genome, sequences in genome_sequences.items():
        output_file = os.path.join(output_dir, f"{genome}")
        try:
            SeqIO.write(sequences, output_file, "fasta")
            print(f"Wrote {len(sequences)} sequences to {output_file}")
        except Exception as e:
            print(f"Error writing to {output_file}: {e}")

def main():
    parser = argparse.ArgumentParser(description="Split a FASTA file into individual genome files based on an STB file.")
    parser.add_argument("fasta_file", help="Input FASTA file")
    parser.add_argument("stb_file", help="Input STB file")
    parser.add_argument("output_dir", help="Output directory for genome files")
    args = parser.parse_args()

    print(f"Input FASTA file: {args.fasta_file}")
    print(f"Input STB file: {args.stb_file}")
    print(f"Output directory: {args.output_dir}")

    split_fasta(args.fasta_file, args.stb_file, args.output_dir)

if __name__ == "__main__":
    main()