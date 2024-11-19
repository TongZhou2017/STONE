#!/bin/bash

# Set the paths for the input and output folders
input_folder="./example_data/SAM_files"  # SAM files input folder
output_folder="./RNAframework_folder/example_output"  # Output folder

# Reference files
reference_folder="./example_data/reference_files"  # Reference files folder
reference_rRNA="yeast_rRNA_tRNA_mtRNA.fa"  # rRNA/tRNA/mtRNA reference file

# Temporary folder
tmp_folder="./RNAframework_folder/tmp"  # Temporary folder for processing

# Check if the input SAM folder exists
if [[ ! -d "$input_folder" ]]; then
    echo "Error: Input SAM folder $input_folder does not exist."
    exit 1
fi

# Check if the output folder exists, if not, create it
if [[ ! -d "$output_folder" ]]; then
    echo "Output folder $output_folder does not exist. Creating it."
    mkdir -p "$output_folder"
fi

# Check if the reference folder exists
if [[ ! -d "$reference_folder" ]]; then
    echo "Error: Reference folder $reference_folder does not exist."
    exit 1
fi

# Check if the reference rRNA file exists
if [[ ! -f "$reference_folder/$reference_rRNA" ]]; then
    echo "Error: Reference rRNA file $reference_folder/$reference_rRNA does not exist."
    exit 1
fi

# Check if the temporary folder exists, if not, create it
if [[ ! -d "$tmp_folder" ]]; then
    echo "Temporary folder $tmp_folder does not exist. Creating it."
    mkdir -p "$tmp_folder"
fi

# Process each SAM file in the input folder
for sam_file in "$input_folder"/*.sam; do
    # Skip if no SAM files are found
    if [[ ! -f "$sam_file" ]]; then
        echo "No SAM files found in $input_folder. Exiting."
        exit 1
    fi

    # Extract the file name (without path and extension)
    filename=$(basename -- "$sam_file")
    filename_no_extension="${filename%.*}"

    # Construct the output file path
    output_file="${output_folder}/${filename_no_extension}_output"

    # Log the current file being processed
    echo "Processing SAM file: $sam_file"
    echo "Output will be saved to: $output_file"

    # Execute the rf-count command
    echo "Running: rf-count -p 70 -a -o $output_file -ow -t $tmp_folder -f $reference_folder/$reference_rRNA -ni -m --orc $sam_file"
    
    rf-count -p 70 -a -o "$output_file" -ow -t "$tmp_folder" -f "$reference_folder/$reference_rRNA" -ni -m --orc "$sam_file"
    
    # Check if the rf-count command was successful
    if [[ $? -ne 0 ]]; then
        echo "Error: rf-count failed for $sam_file."
        continue  # Skip this file and proceed to the next one
    fi

    echo "Processing completed for $sam_file"
done

echo "All processing completed."
