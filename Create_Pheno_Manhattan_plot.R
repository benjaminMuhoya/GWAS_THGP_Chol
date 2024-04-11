library(tidyverse)
library(reshape2)
library(qqman)
library(ggplot2)
library(rio)
library(dplyr)
## Import the needed files
Clean_meta = import("/Users/Bm0211/RegEx/Water_Energy/THGP_database_full_corrected_merged_2024-04-02.txt")
HDL_data = dplyr::select(Clean_meta, c("Unique.ID","HDL.cholesterol.mg.dL",
                                       "LDL.cholesterol.mg.dL","Total.cholesterol.mg.dL.",
                                       "Body.fat.percentage","BMI","Standing.height.cm.", "Age", "Sex"))
colnames(HDL_data)
colnames(HDL_data) <- c("Unique.ID","HDL","LDL","Chol","Body_fat","BMI","Height", "Age", "Sex")
##
calculate_LDL_HDL_Ratio <- function(ldl, hdl) {
  ratio <- rep(NA, length(ldl))
  # Loop through each element
  for(i in 1:length(ldl)) {
    # Check if both LDL and HDL are numeric and not NA
    if(is.numeric(ldl[i]) && !is.na(ldl[i]) && is.numeric(hdl[i]) && !is.na(hdl[i])) {
      # Perform the calculation
      ratio[i] <- ldl[i] / hdl[i]
    } 
    # If either LDL or HDL is not numeric or is NA, ratio[i] will remain NA
  }
  return(ratio)
}
# Apply the function to the LDL and HDL columns of your data frame
colnames(HDL_data)
HDL_data[,2:8] <- lapply(HDL_data[,2:8], as.numeric)
HDL_data$LDL_HDL_Ratio <- calculate_LDL_HDL_Ratio(HDL_data$LDL, HDL_data$HDL)
head(HDL_data)

##IMport Wholesome .fam file
April_GENO <- read.table("/Users/bm0211/RegEx/Water_Energy/wholegenome_AllSamples_updatednames.fam", header = TRUE)
head(April_GENO)
names(April_GENO)[1] <- "FID"
names(April_GENO)[2] <- "IID"
#HDL_data[c('DATE','barcodes','ID')] <- str_split_fixed(HDL_data$Unique.ID, '_', 3)
#Merge April_GENO and HDL_data
combined_data <- merge(April_GENO, HDL_data, by.x = "FID", by.y = "Unique.ID", all.x = TRUE)
colnames(combined_data)
dim(combined_data)
##Rename columns and drop the 6th empty column
colnames(combined_data)
combined_data <-  dplyr::select(combined_data, -c(5, 6))
##Arrange the columns to have Sex as the fifth column
combined_data <- combined_data[, c(1,2,3,4,12,5,6,7,8,9,10,11,13)]
names(combined_data)[5] <- "Sex"
colnames(combined_data)
table(combined_data$Sex) ##14 people(Female\Male ambigous)
combined_data$Sex <- ifelse(combined_data$Sex == "Male", 1, ifelse(combined_data$Sex == "Female", 2, -9))
colnames(combined_data)
names(combined_data)[3] <- "PID"
names(combined_data)[4] <- "MID"
names(combined_data)[5] <- "Sex"
##Turn columns 6-12 to numeric
colnames(combined_data)
combined_data[,6:13] <- lapply(combined_data[,6:13], as.numeric)
##Scale the using scale default R then Identify outliers in each 6-12 column using IQR then filter them out
# Scale each column from 6 to 12
combined_data[, 6:13] <- lapply(combined_data[, 6:13], scale)
# Function to remove outliers based on the IQR method
remove_outliers <- function(x) {
  Q1 <- quantile(x, 0.25, na.rm = TRUE)
  Q3 <- quantile(x, 0.75, na.rm = TRUE)
  IQR <- Q3 - Q1
  x[x >= (Q1 - 1.5 * IQR) & x <= (Q3 + 1.5 * IQR)]
}
# Apply the remove_outliers function to each column from 6 to 12
# Note: This will result in columns of different lengths due to removed outliers
# So, create a list where each element is a column after outlier removal
cleaned_columns <- lapply(combined_data[, 6:13], remove_outliers)
# replace outliers with NA (keeping the dataframe structure consistent):
replace_outliers_with_NA <- function(x) {
  Q1 <- quantile(x, 0.25, na.rm = TRUE)
  Q3 <- quantile(x, 0.75, na.rm = TRUE)
  IQR <- Q3 - Q1
  x[x < (Q1 - 1.5 * IQR) | x > (Q3 + 1.5 * IQR)] <- NA
  x
}
combined_data[, 6:13] <- lapply(combined_data[, 6:13], replace_outliers_with_NA)
##Sanity Check
dim(combined_data)
table(is.na(combined_data$HDL))

##Read in the whole_PCA.eigenvec file and merge with the combined_data
PCA_data <- read.table("/Users/bm0211/RegEx/Water_Energy/whole_PCA.eigenvec", header = TRUE)
colnames(PCA_data)
head(PCA_data)
PCA_data <- dplyr::select(PCA_data, c("IID","PC1"))
##merge with combined_data
combined_data <- merge(combined_data, PCA_data, by = "IID", all.x = TRUE)
##Create Covariate File with first column as Sex, second as Age and third as PC1
Covariate_file <- dplyr::select(combined_data, c("IID", "Sex", "Age", "PC1"))
head(Covariate_file)
##Ensure, covariate file is in the same order as the combined_data based on IID
Covariate_file <- Covariate_file[order(match(Covariate_file$IID, combined_data$IID)),]
dim(Covariate_file)
dim(combined_data)

##Write Covariate file into a .txt file without the first column and no headers
write.table(Covariate_file[, -1], file = "/Users/bm0211/RegEx/Water_Energy/Covariate_file_Wholesome.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
##Write the combineddata into .fam file
colnames(combined_data)
##Remove columns as needed and write the .fam file
#remove choose columns 1:5 then 7 only
select_data <- dplyr::select(combined_data, c(1:5, 13))
write.table(select_data, file = "/Users/bm0211/RegEx/Water_Energy/wholegenome_AllSamples_updatednames.fam", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")


######Manhattan Plot after GWAS
# Ensure p-values are numeric
# Load the GWAS results
gwas_results <- read.table("/Users/bm0211/RegEx/Water_Energy/Based_on_Video_Ratio_LDL_HDL.lmm.assoc.txt", header = TRUE, sep = "")
# Correct renaming was done previously
head(gwas_results)
names(gwas_results)[names(gwas_results) == "ps"] <- "pos"
names(gwas_results)[names(gwas_results) == "p_lrt"] <- "p"
gwas_results$p <- as.numeric(as.character(gwas_results$p))

# Check again for non-numeric p-values just in case
sum(is.na(gwas_results$p))

# Choose a y-axis upper limit AND BONFERRONI THRESHOLD
upper_ylim <- 10
bonferroni_threshold <- 3.84e-07 ##given The Bonferroni corrected p-value threshold, given 130,315 tests, is approximately 3.84 × 10−7 3.84×10−7
manhattan(gwas_results,
          chr = "chr",
          bp = "pos",
          snp = "rs",
          p = "p",
          main = "Manhattan Plot for Ratio_LDL_HDL GWAS",
          col = c("blue4", "orange3"),
          annotatePval = bonferroni_threshold,
          ylim = c(0, upper_ylim))

# QQPLOT
qq(gwas_results$p, main = "QQ Plot for Ratio_LDL_HDL GWAS P-values")


