import csv
import sys
import os

def main(contig_paths, species, output_dir):
    output_file = os.path.join(output_dir, "samplesheet.csv")
    
    # Write the header
    with open(output_file, mode='w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerow(["sample", "contigs", "species"])

        # Write each row
        for sample, contig in contig_paths:
            csv_writer.writerow([sample, contig, species])
    
    print(f"samplesheet.csv has been generated in {output_dir}.")

if __name__ == "__main__":
    # The first argument is the species
    species = sys.argv[1]
    output_dir = sys.argv[2]

    # Collect all sample and contig pairs
    contig_paths = []
    for line in sys.stdin:
        sample, contig = line.strip().split(',')
        contig_paths.append((sample, contig))
    
    main(contig_paths, species, output_dir)

