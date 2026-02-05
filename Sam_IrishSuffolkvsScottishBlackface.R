# echo -"FID IID Cluster" > headers.txt
# cut -d" " -f 1,2 Dataset.fam | grep "IrishSuffolk" | awk '{ print $0, "Grp1" }' > grp1.txt
# cut -d" " -f 1,2 Dataset.fam | grep "ScottishBlackface" | awk '{ print $0, "Grp2" }' > grp2.txt
# cat headers.txt grp1.txt grp2.txt > clusterfst.txt

library(hierfstat)
library(tidyverse)

###############
## 1. Data Prep
###############

# Keep two breeds (IrishSuffolk and ScottishBlackface) from the given dataset 
# and upload the bfile (including log file) of your own specific data subset 
# with the two breeds in compressed zip format. (2.5 points) (save file as 
# example: Yourfirstname_ScottishTexelvsScottishBlackface)

extract.breeds <- function(bfileBasename, clusters, newBfileBasename) {
  system(paste("plink --bfile", bfileBasename,
               "--within clusterfst.txt",
               "--make-bed --out", newBfileBasename,
               sep = " "))
}

extract.breeds("Dataset", "clusterfst.txt", "Sam_IrishSuffolkvsScottishBlackface")

#########
## 2. FST
#########

# Perform a standardized FST analysis between the selected breeds and upload the
# resulting plot with repetitive grey and black colour (use a y-axis of 6 and 
# select only significant region with more than 4 SNPs). (5 points) (save file 
# as example: Yourfirtsname_ScottishTexelvsScottishBlackface)

gen.fst <- function(bfileBasename, clusterFile, fstOutFileBasename) {
  system(paste("plink --bfile", bfileBasename,
               "--chr-set 26 --autosome",
               "--fst --within", clusterFile,
               "--out", fstOutFileBasename,
               sep = " "))
}

gen.fst("Sam_IrishSuffolkvsScottishBlackface",
        "clusterfst.txt",
        "Sam_IrishSuffolkvsScottishBlackface_FST")


#------------------------------
# Load and visualize Fst values
#------------------------------
fst.Balothers.weir <- read.delim("Sam_IrishSuffolkvsScottishBlackface_FST.fst") %>%
  drop_na() 
fst.sd.oth <- sd(fst.Balothers.weir$FST, na.rm = TRUE)
fst.mean.oth <- mean(fst.Balothers.weir$FST, na.rm = TRUE)


#----------
# Calc Zfst
#----------
fst.Balothers.weir$Zfst <- (fst.Balothers.weir$FST - fst.mean.oth)/fst.sd.oth


################
## 3. Top 20 SNP
################

# Report and upload the top 20 SNPs in an excel file format (xlsx) including 
# chromosome, SNP name, BP position and ZFST values) (5 points). (save file as 
# example: Yourfirstname_ScottishTexelvsScottishBlackface)

report.snps.fst <- function(numTopSnps, outfile) {
  errorCondition("not implemented")
  
  # writexl::write_xlsx(topSnps, outfile)
}

report.snps.fst(20, "Sam_IrishSuffolkvsScottishBlackface_FST_SNP20.xlsx")


###########
## 4. XPEHH
###########

# Perform a XPEHH analysis between the breeds and upload the resulting plot with
# repetitive grey and black colour (use a y-axis of +8 to -8). (5points). (save 
# file as example: Yourfirstname_ScottishTexelvsScottishBlackface)

xpehh <- function(bfileBasename, outfile) {
  errorCondition("not implemented")
}
data.xpehh <- xpehh("Sam_IrishSuffolkvsScottishBlackface", "Sam_IrishSuffolkvsScottishBlackface_XPEHH.tiff")


##################
## 5. XPEHH Top 20
##################

# Report and upload the top 20 SNPs in an excel file (xlsx) including 
# chromosome, SNP name, BP position and XPEHH score) (5 points). (save file as 
# example: Yourfirstname_ScottishTexelvsScottishBlackface)

report.snps.xpehh <- function(data, outfile) {
  errorCondition("not implemented")
  
  # writexl::write_xlsx(topSnps, outfile)
}
data.xpehh.top20 <- report.snps.xpehh(data.xpehh, "Sam_IrishSuffolkvsScottishBlackface_XPEHH_SNP20.xlsx")


####################################
## 6. Indicate Breed Under Selection
####################################

# For the top 20 SNPs observed in XPEHH indicate the breed under selection in a 
# new column (column name: Breed) in the same excel sheet. (2.5 points) Note: 
# -Breeds to be compared will differ for most people, therefore ensure you only 
# use the breeds allocated to you for analysis.

addColumn.breedUnderSelection <- function(snpTop20, outfile) {
  errorCondition("not implemented")
  
  # writexl::write_xlsx(snpTop20_wBreedCol, outfile)
}

# Add breed column and overwrite xlsx
addColumn.breedUnderSelection(data.xpehh.top20, "Sam_IrishSuffolkvsScottishBlackface_XPEHH_SNP20.xlsx")