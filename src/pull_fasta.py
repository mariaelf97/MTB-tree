import pandas as pd
from Bio import Entrez, SeqIO
import argparse

def fetch_fasta(accessions, output_file):
    Entrez.email = "mahmadi@scripps.edu"
    with Entrez.efetch(db="nucleotide", id=accessions, rettype="fasta", retmode="text") as handle:
        # Parse the fetched sequences
        records = list(SeqIO.parse(handle, "fasta"))
        # Modify the accession IDs by replacing dots with underscores
        for record in records:
            record.description = ""  # Clear the description if needed
            # Write the modified records to the output file
            with open(output_file, "a") as out_handle:
                SeqIO.write(record, out_handle, "fasta")

        

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description='Pull GenBank (gbk) files from a list of accessions')
    parser.add_argument('accession_list', type=str, help='Input TSV file containing accession list')
    parser.add_argument('output_fasta', type=str, help='Output fasta file for the downloaded sequences')
    parser.add_argument('output_fasta_ref', type=str, help='Output reference fasta file for the downloaded sequences')
    
    # Parse the arguments
    args = parser.parse_args()

    # Read the TSV file and extract the accession numbers as a list
    metadata_df = pd.read_csv(args.accession_list, sep='\t')
    accessions = metadata_df['accession'].tolist()
    fetch_fasta(accessions, args.output_fasta)
    fetch_fasta("NC_000962.3",args.output_fasta_ref)
if __name__ == '__main__':
    main()