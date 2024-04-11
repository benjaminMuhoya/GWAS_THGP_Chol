#!/bin/bash

# Define the path to your PLINK tool
PLINK_PATH="./plink" # Adjust if PLINK is not in your PATH

# Define the path to your main genotype data
GENOTYPE_PREFIX="/Genomics/ayroleslab2/benja/wholegenome_AllSamples_updatednames"

# Define an array of your phenotypes
PHENOTYPES=("HDL" "LDL" "Chol")

# Directory containing the phenotype files
PHENOTYPE_DIR="/Genomics/ayroleslab2/benja/subset_genotype"

# Output directory for the subset genotype data
OUTPUT_DIR="/Genomics/ayroleslab2/benja/subset_genotype"

# Ensure the output directory exists
mkdir -p $OUTPUT_DIR

# Loop through each phenotype
for PHENO in "${PHENOTYPES[@]}"
do
    # Define the phenotype file name, adjusted to match your naming convention
    PHENOTYPE_FILE="${PHENOTYPE_DIR}/${PHENO}_for_GEMMA_phenotypeIID.txt"

    # Since the FID and IID are the same in the new phenotype files, directly extract them for PLINK's --keep option
    awk '{if(NR>1) print $1, $1}' $PHENOTYPE_FILE > "${OUTPUT_DIR}/${PHENO}_ids_to_keep.txt"

    # Use PLINK to subset the genotype data based on the provided identifiers
    $PLINK_PATH --bfile $GENOTYPE_PREFIX --keep "${OUTPUT_DIR}/${PHENO}_ids_to_keep.txt" --make-bed --out "${OUTPUT_DIR}/subset_${PHENO}"

    echo "Subset genotype data for ${PHENO} has been created."
done

echo "All subset genotype data have been created."

