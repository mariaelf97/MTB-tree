import argparse
from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
import pandas as pd
import re


def extract_genes_from_gbk(gbk_file):
    """
    Extract genes, their coordinates, strands, and locus tags from a GenBank file.
    Args:
    - gbk_file (str): Path to the GenBank file.

    Returns:
    - DataFrame: A DataFrame containing locus tags, gene names, start positions, end positions, and strands.
    """
    genes = []
    
    # Parse the GenBank file
    with open(gbk_file, "r") as handle:
        for record in SeqIO.parse(handle, "genbank"):
            for feature in record.features:
                # Only consider 'gene' features
                if feature.type == "gene":
                    gene_name = feature.qualifiers.get("gene", ["Unknown"])[0]  # Get gene name or 'Unknown'
                    locus_tag = feature.qualifiers.get("locus_tag", ["Unknown"])[0]  # Get locus tag or 'Unknown'
                    start = int(feature.location.start) + 1  # Convert to 1-based indexing
                    end = int(feature.location.end)
                    strand = "+" if feature.location.strand == 1 else "-"
                    
                    genes.append({
                        "locus_tag": locus_tag,
                        "gene": gene_name,
                        "start": start,
                        "end": end,
                        "strand": strand
                    })
    
    # Convert the list of genes to a DataFrame
    genes_df = pd.DataFrame(genes)
    return genes_df


def filter_genes(genes_df):
    """
    Filter genes matching the pattern PEN or PPEN, where N is any number.
    Args:
    - genes_df (DataFrame): DataFrame containing extracted genes.

    Returns:
    - DataFrame: Filtered DataFrame.
    """
    # Regular expression for gene names like PEN or PPEN
    pattern = r"^P+E[N\d]+$"
    
    # Filter the DataFrame
    filtered_genes = genes_df[genes_df['gene'].str.match(pattern, na=False)]
    return filtered_genes


def modify_positions_in_fasta(input_file, output_file, genes_df):
    """
    Mask regions in a FASTA file based on gene coordinates from a DataFrame.
    Args:
    - input_file (str): Path to the input FASTA file.
    - output_file (str): Path to the output FASTA file.
    - genes_df (DataFrame): DataFrame with start and end positions to mask.

    Returns:
    - None
    """
    # Read the input FASTA file
    records = list(SeqIO.parse(input_file, "fasta"))
    
    for record in records:
        sequence = str(record.seq)
        masked_sequence = sequence
        
        # Mask each region with 'N'
        for _, row in genes_df.iterrows():
            start = row['start']
            end = row['end']
            replacement = "N" * (end - start + 1)
            masked_sequence = masked_sequence[:start - 1] + replacement + masked_sequence[end:]
        
        record.seq = Seq(masked_sequence)
    
    # Write the modified sequences to the output file
    SeqIO.write(records, output_file, "fasta")
    print(f"Modified sequences saved to {output_file}.")


def main():
    """
    Main function to extract gene information from a GenBank file, filter it based on a pattern,
    and modify a FASTA file by masking regions corresponding to the filtered genes.
    """
    # Set up argument parser
    parser = argparse.ArgumentParser(description='Modify positions in a FASTA file by replacing a range with Ns.')
    parser.add_argument('input_fasta', type=str, help='Path to the input FASTA file')
    parser.add_argument('input_gbk', type=str, help='Path to the input reference GenBank file')
    parser.add_argument('output_file', type=str, help='Path to the output FASTA file')
    
    # Parse the arguments
    args = parser.parse_args()
    
    # Extract genes from the GenBank file
    genes_df = extract_genes_from_gbk(args.input_gbk)
    
    # Filter genes by pattern (e.g., PEN, PPEN, where N is any number)
    filtered_genes_df = filter_genes(genes_df)
    
    # Save the filtered genes to a CSV for verification
    filtered_genes_df.to_csv("filtered_genes.csv", index=False)
    print("Filtered gene information saved to 'filtered_genes.csv'.")
    
    # Modify the FASTA file based on filtered gene coordinates
    modify_positions_in_fasta(args.input_fasta, args.output_file, filtered_genes_df)


if __name__ == '__main__':
    main()
