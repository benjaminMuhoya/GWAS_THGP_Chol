#!/bin/bash

# Define paths
GENOTYPE_DATA="/Genomics/ayroleslab2/benja/March24_BED"
OUTPUT_DIR="/Genomics/ayroleslab2/benja/GWAS_outputs"
REGIONS_FILE="/Genomics/ayroleslab2/benja/top_IHS_SNPs.csv"
PHENOTYPE_FILES=("/Genomics/ayroleslab2/benja/HDL.phen" "/Genomics/ayroleslab2/benja/LDL.phen" "/Genomics/ayroleslab2/benja/Chol.phen")

# Create the output directory if it does not exist
mkdir -p $OUTPUT_DIR

# Calculate whole-genome GRM using all SNPs present in the genotype data
./gcta-1.94.1 --bfile $GENOTYPE_DATA --make-grm --out ${OUTPUT_DIR}/grm_all

# Read the regions from the file
mapfile -t region < <(cat $REGIONS_FILE)

# Loop through each phenotype file
for PHENOTYPE_FILE in "${PHENOTYPE_FILES[@]}"
do
    PHENO_NAME=$(basename "$PHENOTYPE_FILE" .phen)
    COVAR_FILE="/Genomics/ayroleslab2/benja/${PHENO_NAME}covar.covar"

    # Loop through each region
    for REGION in "${region[@]}"
    do
        # Clean up REGION variable, removing quotes and commas
        REGION=$(echo $REGION | tr -d '"' | tr -d ','| tr -d ' ')

        CHR=$(echo $REGION | cut -d'_' -f1)
        FROM_BP=$(echo $REGION | cut -d'_' -f2)
        TO_BP=$(echo $REGION | cut -d'_' -f3)
        SUBSET_PREFIX="subset_${CHR}_${FROM_BP}_${TO_BP}_${PHENO_NAME}"

        # Subset data for the region with PLINK
        ./plink --bfile $GENOTYPE_DATA --chr $CHR --from-bp $FROM_BP --to-bp $TO_BP --make-bed --out ${OUTPUT_DIR}/${SUBSET_PREFIX}

        # Note: Now, the GRM used is the whole-genome GRM calculated before.
        # Run mixed linear model analysis with GCTA, accounting for relatedness using the whole-genome GRM
        ./gcta-1.94.1 --grm ${OUTPUT_DIR}/grm_all --pheno $PHENOTYPE_FILE --mpheno 3 --covar $COVAR_FILE --reml --out ${OUTPUT_DIR}/analysis_result_${CHR}_${FROM_BP}_${TO_BP}_${PHENO_NAME}

        echo "Region-based analysis completed for: ${CHR}_${FROM_BP}_${TO_BP}, Phenotype: $PHENO_NAME"
    done
done

