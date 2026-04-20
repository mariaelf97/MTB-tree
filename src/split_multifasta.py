from Bio import SeqIO
import os
import argparse

def split_multifasta(input_fasta, output_dir):
    """
    Splits a multi-FASTA file into individual FASTA files.
    
    Args:
        input_fasta (str): Path to the input multi-FASTA file.
        output_dir (str): Directory to save the individual FASTA files.
    """
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Iterate through each record in the multi-FASTA file
    for record in SeqIO.parse(input_fasta, "fasta"):
        # Create a unique filename for each sequence
        filename = f"{record.id}.fasta"
        output_path = os.path.join(output_dir, filename)
        
        # Write the individual FASTA record to a file
        with open(output_path, "w") as output_handle:
            SeqIO.write(record, output_handle, "fasta")
        
        print(f"Saved: {output_path}")

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Split a multi-FASTA file into individual FASTA files.")
    parser.add_argument("input_fasta", type=str, help="Path to the input multi-FASTA file")
    parser.add_argument("output_dir", type=str, help="Directory to save the individual FASTA files")
    args = parser.parse_args()
    
    # Run the split function
    split_multifasta(args.input_fasta, args.output_dir)

if __name__ == "__main__":
    main()
