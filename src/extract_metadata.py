import argparse
from Bio import SeqIO
import pandas as pd
from datetime import datetime

# Function to extract geolocation and sample ID from a GenBank record
def extract_metadata(record):
    # Initialize variables
    sample_id = record.id
    geolocation = None
    strain = None
    host = None

    # Access the date in the LOCUS line
    date = record.annotations.get("date")

    # Iterate through features in the GenBank record
    for feature in record.features:
        if "geo_loc_name" in feature.qualifiers:
            geolocation = "; ".join(feature.qualifiers.get("geo_loc_name"))
        if "strain" in feature.qualifiers:
            strain = "; ".join(feature.qualifiers.get("strain"))
    return strain, date, sample_id, geolocation

# Function to parse a multi-sequence GenBank file and return a DataFrame
def parse_genbank_file(file_path):
    data = []
    # Parse the GenBank file
    for record in SeqIO.parse(file_path, "genbank"):
        strain, date, sample_id, geolocation = extract_metadata(record)
        region = geolocation.split(":")[0].strip() if geolocation else None
        
        # Convert the date format from '18-MAY-2021' to '2021-05-18'
        try:
            formatted_date = datetime.strptime(date, "%d-%b-%Y").strftime("%Y-%m-%d")
        except (ValueError, TypeError):
            formatted_date = None  # Handle cases where date format is not as expected

        data.append({"strain": sample_id, "accession": sample_id, "date": formatted_date,
                     "country": region})
    # Convert list of dictionaries to DataFrame
    df = pd.DataFrame(data)
    return df

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description='Parse a GenBank file and extract metadata to a TSV file.')
    parser.add_argument('genbank_file', type=str, help='Path to the input GenBank file')
    parser.add_argument('genotype_file', type=str, help='Path to the input genotype file')
    parser.add_argument('output_tsv', type=str, help='Path to the output TSV file')

    # Parse the arguments
    args = parser.parse_args()

    # Parse the GenBank file and get the DataFrame
    df = parse_genbank_file(args.genbank_file)
    genotype = pd.read_csv(args.genotype_file, sep="\t")
    joined_df = pd.merge(df, genotype, on=['accession'], how='left')
    # Save the DataFrame to a TSV file
    joined_df.to_csv(args.output_tsv, index=False, sep="\t")
if __name__ == '__main__':
    main()
