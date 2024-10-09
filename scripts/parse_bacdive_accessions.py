import csv
import argparse

def parse_bacdive_csv(input_file, output_file):
    current_id = None
    current_species = None
    current_gca = None

    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        # Write the header
        outfile.write("ID\tSpecies\tGCA\n")

        # Skip the first two lines
        next(infile)
        next(infile)

        csv_reader = csv.reader(infile)
        for row in csv_reader:
            if row[0]:  # New entry
                # Write the previous entry if it exists
                if current_id and current_species and current_gca:
                    outfile.write(f"{current_id}\t{current_species}\t{current_gca}\n")

                # Start a new entry
                current_id = row[0]
                current_species = row[4]
                current_gca = next((acc for acc in row if acc.startswith('GCA_')), None)
            else:  # Additional entry for the current ID
                if not current_gca:
                    current_gca = next((acc for acc in row if acc.startswith('GCA_')), None)

        # Write the last entry
        if current_id and current_species and current_gca:
            outfile.write(f"{current_id}\t{current_species}\t{current_gca}\n")

def main():
    parser = argparse.ArgumentParser(description='Parse BacDive CSV file into a simple TSV.')
    parser.add_argument('input_file', help='Path to the input CSV file')
    parser.add_argument('output_file', help='Path to the output TSV file')
    args = parser.parse_args()

    parse_bacdive_csv(args.input_file, args.output_file)
    print(f"Parsing complete. Output written to {args.output_file}")

if __name__ == '__main__':
    main()