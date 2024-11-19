#!/bin/bash

# Define paths for directories and files
sam_directory="./example_data/SAM_files"       # Input SAM files directory
output_directory="./RNAframework_folder/example_output/processed_files"  # Output directory
reference_directory="./example_data/reference_files"  # Reference files directory

# Check if the input SAM directory exists
if [[ ! -d "$sam_directory" ]]; then
    echo "Error: SAM directory $sam_directory does not exist."
    exit 1
fi

# Check if the output directory exists, if not, create it
if [[ ! -d "$output_directory" ]]; then
    echo "Output directory $output_directory does not exist. Creating it."
    mkdir -p "$output_directory"
fi

# Check if necessary reference files exist
if [[ ! -f "$reference_directory/example.len" ]]; then
    echo "Error: Reference file example.len does not exist in $reference_directory."
    exit 1
fi

if [[ ! -f "$reference_directory/example.fa" ]]; then
    echo "Error: Reference file example.fa does not exist in $reference_directory."
    exit 1
fi

# Process each SAM file in the SAM directory
for sam_file in "$sam_directory"/*.sam; do
    # Skip if no SAM files are found
    if [[ ! -f "$sam_file" ]]; then
        echo "No SAM files found in $sam_directory. Exiting."
        exit 1
    fi
    
    # Get the SAM file name without the path
    sam_filename=$(basename "$sam_file")
    
    # Construct corresponding output file paths
    tab_file="$output_directory/${sam_filename%.sam}.tab"
    output_file="$output_directory/${sam_filename%.sam}_countRT.c"
    shape_file="$output_directory/${sam_filename%.sam}_shape.gTab"
    csv_file="${output_file%.c}.csv"

    # Log file paths
    echo "Processing SAM file: $sam_file"
    echo "Output TAB file: $tab_file"
    echo "Output countRT file: $output_file"
    echo "Converted CSV file: $csv_file"

    # Run the icSHAPE-pipe sam2tab command
    echo "Running: icSHAPE-pipe sam2tab -in $sam_file -out $tab_file"
    icSHAPE-pipe sam2tab -in "$sam_file" -out "$tab_file"
    if [[ $? -ne 0 ]]; then
        echo "Error: icSHAPE-pipe sam2tab failed for $sam_file"
        continue
    fi
    
    # Run the icSHAPE-pipe countRT command
    echo "Running: icSHAPE-pipe countRT -in $tab_file -size $reference_directory/example.len -out $output_file"
    icSHAPE-pipe countRT -in "$tab_file" -size "$reference_directory/example.len" -out "$output_file" -omc 0
    if [[ $? -ne 0 ]]; then
        echo "Error: icSHAPE-pipe countRT failed for $tab_file"
        continue
    fi

    # Convert the countRT output file to CSV
    echo "Converting countRT output to CSV: sed 's/\t/,/g' $output_file > $csv_file"
    sed 's/\t/,/g' "$output_file" > "$csv_file"
    if [[ $? -ne 0 ]]; then
        echo "Error: Conversion to CSV failed for $output_file"
        continue
    fi

    echo "Processing completed for $sam_filename"
done

echo "All processing complete."
