
:## First generate the PCA
./plink --bfile wholegenome_AllSamples_updatednames --pca 10 'header' 'tabs' 'var-wts' --out whole_PCA

##Then generate the kinship matrix
./gemma -bfile wholegenome_AllSamples_updatednames -gk 1 -o kinship_matrix_Apr


##Then run the GWAS
./gemma -bfile wholegenome_AllSamples_updatednames -k kinship_matrix_Apr.cXX.txt -lmm 2 -c Covariate_file_Wholesome.txt -o Chol.lmm

##Make manhattan plot
