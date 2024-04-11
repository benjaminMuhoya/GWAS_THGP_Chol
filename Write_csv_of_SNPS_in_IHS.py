#!/bin/bash

# Define paths
GENOTYPE_DATA="./GEMMA_files" # The directory where your PLINK files (.bed, .bim, .fam) are stored
OUTPUT_DIR="./output_snp_regions"
REGIONS_FILE="./top_IHS_SNPs.txt"
PLINK_PATH="./plink" # Ensure PLINK is accessible and executable

# Ensure the output directory exists
mkdir -p "$OUTPUT_DIR"

# Read the regions from the file and loop through each
while IFS=$'\r' read -r REGION; do
    # Remove carriage return, if present
    REGION=$(echo "$REGION" | tr -d '\r')
    
    # Clean up REGION variable
    IFS='_' read -ra ADDR <<< "$REGION"
    CHR="${ADDR[0]}"
    FROM_BP="${ADDR[1]}"
    TO_BP="${ADDR[2]}"

    # Subset data for the specified region with PLINK and recode into a text file with space delimiters
    ${PLINK_PATH} --bfile ${GENOTYPE_DATA} --chr $CHR --from-bp $FROM_BP --to-bp $TO_BP \
                  --recode A --out "${OUTPUT_DIR}/${REGION}"

    # Get the headers from the .raw file and format them as CSV
    HEADERS=$(head -n 1 "${OUTPUT_DIR}/${REGION}.raw" | tr ' ' ',')
    
    # Create a new CSV file and write the headers
    echo "$HEADERS" > "${OUTPUT_DIR}/${REGION}.csv"
    
    # Skip the header line in the .raw file, then append the rest to the CSV
    tail -n +2 "${OUTPUT_DIR}/${REGION}.raw" | tr ' ' ',' >> "${OUTPUT_DIR}/${REGION}.csv"
    
    # Optionally, remove the intermediate .raw file
    rm "${OUTPUT_DIR}/${REGION}.raw"

    echo "Extracted SNPs for region: ${REGION} in CSV format"
done < "$REGIONS_FILE"

echo "All regions processed. SNP data saved in CSV format in the directory: ${OUTPUT_DIR}"

