import pandas as pd
from Bio import Entrez, SeqIO
import argparse

def fetch_gbk(accessions, output_file, ref):
    Entrez.email = "mahmadi@scripps.edu"  # Set your email address
    if ref:
        # Fetch the records from GenBank in GenBank format
        with Entrez.efetch(db="nucleotide", id=accessions, rettype="gbwithparts", retmode="text") as handle:
            # Parse the fetched sequences in GenBank format
            records = list(SeqIO.parse(handle, "genbank"))
            # Write the records to the output file in GenBank format
            with open(output_file, "a") as out_handle:
                SeqIO.write(records, out_handle, "genbank")
    else:
        # Fetch the records from GenBank in GenBank format
        with Entrez.efetch(db="nucleotide", id=accessions, rettype="gb", retmode="text") as handle:
            # Parse the fetched sequences in GenBank format
            records = list(SeqIO.parse(handle, "genbank"))
            # Write the records to the output file in GenBank format
            with open(output_file, "a") as out_handle:
                SeqIO.write(records, out_handle, "genbank")

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description='Pull GenBank (gbk) files from a list of accessions')
    parser.add_argument('accession_list', type=str, help='Input TSV file containing accession list')
    parser.add_argument('output_gbk', type=str, help='Output GBK file for the downloaded sequences')
    parser.add_argument('output_ref', type=str, help='Output GBK file for the downloaded sequences')
    # Parse the arguments
    args = parser.parse_args()

    # Read the TSV file and extract the accession numbers as a list
    metadata_df = pd.read_csv(args.accession_list, sep='\t')
    accessions = metadata_df['accession'].tolist()
    
    # Fetch and save the GenBank records
    fetch_gbk(accessions, args.output_gbk,ref=False)
    fetch_gbk("NC_000962.3",args.output_ref,ref=True)
    
if __name__ == '__main__':
    main()
