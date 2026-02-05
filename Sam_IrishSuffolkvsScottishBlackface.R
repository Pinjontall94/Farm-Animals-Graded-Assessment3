# echo -"FID IID Cluster" > headers.txt
# cut -d" " -f 1,2 Dataset.fam | grep "IrishSuffolk" | awk '{ print $0, "Grp1" }' > grp1.txt
# cut -d" " -f 1,2 Dataset.fam | grep "ScottishBlackface" | awk '{ print $0, "Grp2" }' > grp2.txt
# cat headers.txt grp1.txt grp2.txt > clusterfst.txt

library(hierfstat)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(rehh)
library(data.table)

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
n2 <- 1000000 #look at window of 1MB
# calc mean zfst within non-overlapping window of 1Mb
fst.calc.Balothers.n <- function(fstData, nonOverlappingWindow) {
  fstData %>% 
    group_by(fstData$CHR, mean = (fstData$POS) %/% nonOverlappingWindow) %>% 
    mutate(freq = n()) %>%
    mutate(ZFST = mean(Zfst,na.rm = TRUE))
}

fst.calc.Balothers.n2 <- function(balothersN, minSNPs) {
  balothersN %>% filter(freq > minSNPs)
}

fst.Balothers.n <- fst.calc.Balothers.n2(fst.Balothers.weir, 1000000)

fst.Balothers.n2 <- fst.calc.Balothers.n2(fst.Balothers.n, 4)
fst.meanSNP <- mean(fst.Balothers.n2$freq)
fst.sdSNP <- sd(fst.Balothers.n2$freq)
fst.range <- range(fst.Balothers.n2$freq)

# calculate total number of info windows
# for each chromosome (here, 26)
for(i in 1:26) {
  fst.Balothers.n3 <- fst.Balothers.n2 %>% filter(CHR==i)
  fst.Balothers.n3$mean <- as.factor(fst.Balothers.n3$mean)
  fst.Wnd = str(length(levels(fst.Balothers.n3$mean)))
}

# cat << EOF > cols.txt
# <copy-paste previous output into vt-100 emulator window>
# EOF

# order by fst and check fst value for top 1% SNP
fst.Balothers.n22 <- fst.Balothers.n2[order(-fst.Balothers.n2$ZFST), ]
fst.Balothers.n22[392,10]

fst.Balothers.n2$CHR<- as.integer(fst.Balothers.n2$CHR)#change CHROM to interger

fdata_cum <- fst.Balothers.n2 %>% 
  group_by(CHR) %>% 
  summarise(max_bp = max(POS)) %>% 
  mutate(bp_add = lag(cumsum(as.numeric(max_bp)), default = 0)) %>% 
  select(CHR, bp_add)

fst_data <- fst.Balothers.n2 %>% 
  inner_join(fdata_cum, by = "CHR") %>% 
  mutate(bp_cum = POS + bp_add)

faxisdf = fst_data %>% group_by(CHR) %>% summarize(center=( max(bp_cum) + min(bp_cum) ) / 2 )

fsig <- fst.Balothers.n22[ceiling(nrow(fst.Balothers.n22) * 0.01), "ZFST"]$ZFST


fstothers <- ggplot(fst_data, aes(x=bp_cum, y= ZFST)) +
  
  # Show all points
  geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=0.1) +
  scale_color_manual(values = rep(c("grey", "black"), 26 )) +
  geom_hline(yintercept = fsig , color = "grey40", linetype = "dashed") +
  # custom X axis:
  scale_x_continuous( label = c(1,2,3,4,5,6,7,8,9,"",11,"",13,"",15,"",17,"",19,"","",22,"","","",26), breaks= faxisdf$center ) +
  scale_y_continuous(expand = c(0, 0),limits = c(0,6)) +     # remove space between plot area and x axis
  labs(
    x =expression(paste("Chromosome")),
    y =expression(paste("ZF"["ST"]))
  ) + 
  # Custom the theme:
  theme_bw() +
  theme_classic()+
  theme( 
    legend.position="none",plot.margin = margin(0.5,3,0.5,1,"cm"),
    axis.line.x.bottom = element_line(linewidth = 0.6,arrow = NULL), axis.line.y.left = element_line(linewidth = 0.6,arrow = NULL),axis.text = element_text(size = 7),axis.title = element_text(size = 8))

tiff(filename = "Sam_IrishSuffolkvsScottishBlackface_FST.tiff",
     width=15,height=5 , units = "cm",
     compression = "lzw",
     bg = "white", res = 1200, family = "",
     type ="cairo")
par(mar=c(5.1, 4.1, 4.1, 2.1)) 
ggarrange(fstothers,
          labels ="",
          ncol = 1, nrow = 1, hjust= -0.1, vjust = 1.0,widths = 1, heights = 0.2,font.label = list(size = 12, color = "black", face = "bold", family = NULL),legend = NULL)
dev.off()

################
## 3. Top 20 SNP
################

# Report and upload the top 20 SNPs in an excel file format (xlsx) including 
# chromosome, SNP name, BP position and ZFST values) (5 points). (save file as 
# example: Yourfirstname_ScottishTexelvsScottishBlackface)

report.snps.fst <- function(numTopSnps, fstData, outfile) {
  topSnps <- fstData[numTopSnps,]
  
  writexl::write_xlsx(topSnps, outfile)
}

# report.snps.fst(20, fst.Balothers.n3, "Sam_IrishSuffolkvsScottishBlackface_FST_SNP20.xlsx")

write.csv(x = fst.Balothers.n3, file = "Sam_IrishSuffolkvsScottishBlackface_FST_SNP20.csv")

# not working, use write.csv and convert in libreoffice calc to xlsx
# writexl::write_xlsx(fst.top20, path = "Sam_IrishSuffolkvsScottishBlackface_FST_SNP20.xlsx")

###########
## 4. XPEHH
###########

# Perform a XPEHH analysis between the breeds and upload the resulting plot with
# repetitive grey and black colour (use a y-axis of +8 to -8). (5points). (save 
# file as example: Yourfirstname_ScottishTexelvsScottishBlackface)

xpehh <- function(bfileBasename, outfile) {
  data <- bfileBasename
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

bim <- read.table("Sam_IrishSuffolkvsScottishBlackface.bim", header = FALSE, stringsAsFactors = FALSE)
colnames(bim) <- c("CHR", "SNP", "CM", "POSITION", "A1", "A2")
bim$CHR <- as.integer(bim$CHR)

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