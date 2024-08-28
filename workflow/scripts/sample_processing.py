#!/usr/bin/env python

'''A function to get the sample data needed for the pipeline.'''

# Reading Metadata Columns:

# The script reads all columns from the CSV except for the sample_names column, which is used for matching.
# Storing Metadata:

# The script extracts metadata for each sample and appends it to the list of rows that will be written to the output file.
# Writing to Output CSV:

# The script writes the sample name, sample number, lane, read type, file path, and any additional metadata to the output CSV file.

import pandas as pd
import glob
import os
import re

def get_sample_data(csv_file):
    """
    A function to get the sample data needed for the pipeline.
    """
    df = pd.read_csv(csv_file)

    # Assuming the sample names are in a column named 'sample'
    sample_names = df['sample'].astype(str).tolist()

    # Directory where your fastq files are located
    fastq_dir = os.getcwd() + '/samples/'

    # Dictionary to hold the fastq file paths for each sample
    sample_fastq_files = {}

    missing_files = []

    for sample in sample_names:
        sample_fastq_files[sample] = []
        for lane in range(1, 2):  # Loop through lanes 1 to 4, for deployment change to (1, 2)
            r1_pattern = os.path.join(fastq_dir, "{}_S*_L00{}_R1_*.fastq.gz".format(sample, lane))
            r2_pattern = os.path.join(fastq_dir, "{}_S*_L00{}_R2_*.fastq.gz".format(sample, lane))

            r1_files = glob.glob(r1_pattern)
            r2_files = glob.glob(r2_pattern)

            if not r1_files:
                missing_files.append("Missing R1 file for sample {}, lane {}".format(sample, lane))
            if not r2_files:
                missing_files.append("Missing R2 file for sample {}, lane {}".format(sample, lane))

            for file in r1_files:
                match = re.search(r'_S(\d+)_L00(\d)_R1_', file)
                if match:
                    sample_number = match.group(1)
                    matched_lane = match.group(2)  # Use a different variable name for the matched lane
                    sample_fastq_files[sample].append((sample_number, 'L00{}'.format(matched_lane), 'R1', file))
            
            for file in r2_files:
                match = re.search(r'_S(\d+)_L00(\d)_R2_', file)
                if match:
                    sample_number = match.group(1)
                    matched_lane = match.group(2)  # Use a different variable name for the matched lane
                    sample_fastq_files[sample].append((sample_number, 'L00{}'.format(matched_lane), 'R2', file))

    # Check for missing files
    if missing_files:
        for missing_file in missing_files:
            print(missing_file)
        raise FileNotFoundError("Some files are missing. Please check the above messages.")
    # else:
        # print("All files are present.")

    # Merge the metadata with the fastq file information
    metadata_columns = df.columns.tolist()
    metadata_columns.remove('sample')

    output_rows = []

    for index, row in df.iterrows():
        sample = str(row['sample'])
        metadata = [row[col] for col in metadata_columns]
        
        for sample_number, lane, read, file in sample_fastq_files[sample]:
            merged_sample = "{}_S{}".format(sample, sample_number)
            output_rows.append([merged_sample, lane, read, file] + metadata)

    # Save the results to a CSV file
    output_file = 'sample_data.csv'
    with open(output_file, 'w') as f:
        header = ['sample', 'lane', 'read', 'file'] + metadata_columns
        f.write(','.join(header) + '\n')
        
        for row in output_rows:
            f.write(','.join(map(str, row)) + '\n')
    
    # print("Results have been saved to {}".format(output_file))
    
    # Return the CSV file as a dataframe
    return pd.read_csv(output_file)