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
  # TODO: add the analysis below into this function body
}
data.xpehh <- xpehh("Sam_IrishSuffolkvsScottishBlackface", "Sam_IrishSuffolkvsScottishBlackface_XPEHH.tiff")

extract.breeds("Dataset", "clusterIrishSuffolk.txt", "Sam_IrishSuffolk")
extract.breeds("Dataset", "clusterScottishBlackface.txt", "Sam_ScottishBlackface")

#create bfile for each chromosome per population
CHR= c(1:26)
f="Sam_IrishSuffolk"
i="Sam_IrishSuffolk"
library(foreach)
foreach (r=CHR) %do% {
  cmd <- paste0("plink2 --bfile ", f, " --chr-set 26 --chr ", r, " --make-bed --out ", i, "_chr", r, "_data")
  system(cmd)
}

f="Sam_ScottishBlackface"
i="Sam_ScottishBlackface"
library(foreach)
foreach (r=CHR) %do% {
  cmd <- paste0("plink2 --bfile ", f, " --chr-set 26 --chr ", r, " --make-bed --out ", i, "_chr", r, "_data")
  system(cmd)
}


#create vcf for each chromosome per population
i="Sam_IrishSuffolk"
library(foreach)
foreach (r=CHR) %do% {
  cmd <- paste0("plink2 --bfile ", f, "_chr", r, "_data --chr-set 26 --chr ", r, " --recode vcf --out ", i, "_chr", r, "_data")
  system(cmd)
}

i="Sam_ScottishBlackface"
library(foreach)
foreach (r=CHR) %do% {
  cmd <- paste0("plink2 --bfile ", f, "_chr", r, "_data --chr-set 26 --chr ", r, " --recode vcf --out ", i, "_chr", r, "_data")
  system(cmd)
}


## Run Beagle to phase each VCF file in linux for each population and all chromosomes
for z in {1..26};do java -Xmx25000m -jar beagle.jar gt=Sam_IrishSuffolk_chr${z}_data.vcf  out=Sam_IrishSuffolk_chr${z}_phased nthreads=128;done
for z in {1..26};do java -Xmx25000m -jar beagle.jar gt=Sam_ScottishBlackface_chr${z}_data.vcf  out=Sam_ScottishBlackface_chr${z}_phased nthreads=128;done
read.table

##################################################################################Baldata
for(i in 1:26) {
  # haplotype file name for each chromosome
  vcf_file = paste("Sam_IrishSuffolk_chr", i, "_phased.vcf.gz", sep = "")
  # create internal representation
  hh <-data2haplohh(hap_file = vcf_file, 
                    chr.name = i, 
                    polarize_vcf = FALSE,
                    vcf_reader = "data.table")  # Skip ancestral allele polarization
  # perform scan on a single chromosome (calculate iHH values)
  scan <- scan_hh(hh)
  # concatenate chromosome-wise data frames to
  # a data frame for the whole genome
  # (more efficient ways certainly exist...)
  if (i == 1) {
    wgscan_IrishSuffolk <- scan
  } else {
    wgscan_IrishSuffolk <- rbind(wgscan_bal, scan)
  }
}

#########################################Bardata
for(i in 1:26) {
  # haplotype file name for each chromosome
  vcf_file = paste("Sam_ScottishBlackface_chr", i, "_phased.vcf.gz", sep = "")
  # create internal representation
  hh <-data2haplohh(hap_file = vcf_file, 
                    chr.name = i, 
                    polarize_vcf = FALSE,
                    vcf_reader = "data.table")  # Skip ancestral allele polarization
  # perform scan on a single chromosome (calculate iHH values)
  scan <- scan_hh(hh)
  # concatenate chromosome-wise data frames to
  # a data frame for the whole genome
  # (more efficient ways certainly exist...)
  if (i == 1) {
    wgscan_ScottishBlackface <- scan
  } else {
    wgscan_ScottishBlackface <- rbind(wgscan_bar, scan)
  }
}


#################################################################################################################
xpehh.IRISHSUFFOLK_SCOTTISHBLACKFACE <- ies2xpehh(scan_pop1 =  wgscan_IrishSuffolk,
                           scan_pop2 =  wgscan_ScottishBlackface,
                           popname1 = "IrishSuffolk",
                           popname2 = "ScottishBlackface",
                           standardize = TRUE,
                           p.adjust.method = "BH")


distribplot(xpehh.IRISHSUFFOLK_SCOTTISHBLACKFACE$XPEHH_IRISHSUFFOLK_SCOTTISHBLACKFACE, xlab = "XPEHH")

xpehh.IRISHSUFFOLK_SCOTTISHBLACKFACE$LOGPVALUE2 <- xpehh.IRISHSUFFOLK_SCOTTISHBLACKFACE$LOGPVALUE
Pval.xpehh2n <- (1-2*(abs(pnorm(xpehh.IRISHSUFFOLK_SCOTTISHBLACKFACE$XPEHH_IRISHSUFFOLK_SCOTTISHBLACKFACE)-0.5)))
#xpehh.IRISHSUFFOLK_SCOTTISHBLACKFACE$LOGPVALUE3<- -log(pnorm(xpehh.IRISHSUFFOLK_SCOTTISHBLACKFACE$XPEHH_IRISHSUFFOLK_SCOTTISHBLACKFACE))
xpehh.IRISHSUFFOLK_SCOTTISHBLACKFACE$LOGPVALUE<- -log(Pval.xpehh2n)
write.table(xpehh.IRISHSUFFOLK_SCOTTISHBLACKFACE[,], file = "TOTAL_IRISHSUFFOLK_SCOTTISHBLACKFACE_xpehh.csv",sep="\t", row.names=TRUE, col.names=TRUE,  quote = T)
OrderIRISHSUFFOLK_SCOTTISHBLACKFACE<-xpehh.IRISHSUFFOLK_SCOTTISHBLACKFACE[order(-xpehh.IRISHSUFFOLK_SCOTTISHBLACKFACE$LOGPVALUE), ]
OrderIRISHSUFFOLK_SCOTTISHBLACKFACE[392,3]
write.table(OrderIRISHSUFFOLK_SCOTTISHBLACKFACE[c(1:392),], file = "IRISHSUFFOLK_SCOTTISHBLACKFACENew.csv",sep="\t", row.names=TRUE, col.names=TRUE,  quote = T)


xpehh.IRISHSUFFOLK_SCOTTISHBLACKFACE$CHR<- as.integer(xpehh.IRISHSUFFOLK_SCOTTISHBLACKFACE$CHR)#change CHROM to interger

fdata_cumxp2n <- xpehh.IRISHSUFFOLK_SCOTTISHBLACKFACE %>% 
  group_by(CHR) %>% 
  summarise(max_bp = max(POSITION)) %>% 
  mutate(bp_add = lag(cumsum(as.numeric(max_bp)), default = 0)) %>% 
  select(CHR, bp_add)

fst_dataxp2n <- xpehh.IRISHSUFFOLK_SCOTTISHBLACKFACE %>% 
  inner_join(fdata_cumxp2n, by = "CHR") %>% 
  mutate(bp_cum = POSITION + bp_add)

faxisdfxp2n = fst_dataxp2n %>% group_by(CHR) %>% summarize(center=( max(bp_cum) + min(bp_cum) ) / 2 )

fsigxp2n <-c(-2.57,2.57)

#rep(c("#E69F00", "#56B4E9", "#009E73","#F0E442", "#0072B2", "#D55E00", "#CC79A7")


xpebb<-ggplot(fst_dataxp2n, aes(x=bp_cum, y= XPEHH_IRISHSUFFOLK_SCOTTISHBLACKFACE)) +
  
  # Show all points
  geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=0.1) +
  scale_color_manual(values = rep(c("grey", "black"), 26 )) +
  geom_hline(yintercept = fsigxp2n , color = "red", linetype = "solid",linewidth = 0.2) +
  # custom X axis:
  scale_x_continuous( label = c(1,2,3,4,5,6,7,"",9,"",11,"",13,"",15,"",17,"",19,"","",22,"",24,"",26), breaks= faxisdfxp2n$center ) +
  scale_y_continuous(expand = c(0, 0),limits = c(-8,8),breaks = scales::pretty_breaks(n = 8)) +     # remove space between plot area and x axis
  labs(
    subtitle = expression(paste("XP-EHH")),
    x =expression(paste("Chromosome")),
    y =expression(paste("XP-EHH"))
  ) + 
  # Custom the theme:
  theme_bw() +
  theme_classic()+
  theme(
    legend.position="none",plot.margin = margin(0,0.1,0.2,0.1,unit = "mm"),line = element_line(linewidth = 0.2,arrow = NULL),
    axis.line.x.bottom = element_line(linewidth = 0.2,arrow = NULL), axis.line.y.left = element_line(linewidth = 0.2,arrow = NULL),axis.text = element_text(size = 9),axis.title = element_text(size = 8),
    #plot.title = element_text(family="Arial",face="plain", color="Black", 
    #                            size=15, angle=0,hjust= 1, vjust = -5),
    plot.subtitle = element_text(family="Arial",face="plain", color="Black", 
                                 size=12, angle=0,hjust= 1, vjust = -5),
    axis.text.x = element_text(family= "Arial",face="plain", color="Black", 
                               size=5, angle=0,hjust= 0.5, vjust = 1),
    axis.text.y = element_text(family= "Arial",face="plain", color="Black", 
                               size=5, angle=0,hjust= 1, vjust = 0.5),
    axis.title.x = element_text(family= "Arial",face="plain", color="Black", 
                                size=9, angle=0,hjust= 0.5, vjust = 1),
    axis.title.y = element_text(family= "Arial",face="plain", color="Black", 
                                size=9, angle=90,hjust= 0.5, vjust = 0.5),
    axis.ticks=element_line(linewidth = 0.2,arrow = NULL))

tiff(filename = "Sam_IrishSuffolkvsScottishBlackfaceXPEHH.tiff",
     width=9,height=5, units = "cm",
     compression = "lzw",
     bg = "white", res = 1000, family = "",
     type ="cairo")
par(mar=c(5.1, 4.1, 4.1, 2.1)) 
ggarrange(xpebb,
          labels ="",
          ncol = 1, nrow = 1, hjust= -0.1, vjust = 1.0,widths = 1, heights = 0.2,font.label = list(size = 12, color = "black", face = "bold", family = NULL),legend = NULL)
dev.off()

tiff(filename = "Sam_IrishSuffolkvsScottishBlackface_AllDataset1and2_Final-afterSM.tiff",
     width=18,height=15 , units = "cm",
     compression = "lzw",
     bg = "white", res = 1000, family = "",
     type ="cairo")
par(mar=c(5.1, 4.1, 4.1, 2.1)) 
ggarrange(fstothers,xpebb,fstothershap,
          labels =c("A","B","","","","",""),
          ncol = 2, nrow = 3, hjust= -0.5, vjust = 1.0,widths = 1, heights = 0.2,font.label = list(size = 12, color = "black", face = "bold", family = NULL),legend = NULL)
dev.off()

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