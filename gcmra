###Original Resource: https://www.biostars.org/p/85124/

#Set some data directories
tempdir="C:/Users/bas859/OneDrive - The University of Texas-Rio Grande Valley/Dr.Upal/Bryan/MD/MachineLearn"


#install the core bioconductor packages, if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install()
# install additional bioconductor libraries, if not already installed
BiocManager::install(c("GEOquery", "affy"))
BiocManager::install(c("gcrma", "org.Hs.eg.db"))


# Install custom CDF packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("hgu133acdf")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("hgu133aprobe")


if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("hgu133a.db")


# Load all necessary libraries
library(GEOquery)
library(affy)
library(gcrma)
library(hgu133ahsentrezgcdf) #cdfname="HGU133A_HS_ENTREZG"
library(hgu133acdf) 
library(hgu133aprobe)
library(hgu133a.db)

# Set a data directory, download a GEO dataset, unpack and gunzip, and create a list of files for processing 
setwd(tempdir)
getGEOSuppFiles("GSE2034")
setwd(paste(tempdir,"GSE2034", sep=""))
untar("GSE2034_RAW.tar", exdir="data")
cels = list.files("data/", pattern = "CEL")
sapply(paste("data", cels, sep="/"), gunzip)
cels = list.files("data/", pattern = "CEL")

setwd("C:/Users/bas859/OneDrive - The University of Texas-Rio Grande Valley/Dr.Upal/Bryan/MD/MachineLearn/GSE2034/data")
raw.data=ReadAffy(verbose=TRUE, filenames=cels, cdfname="hgu133a")  ##print(raw.data@cdfName)
data.gcrma.norm=gcrma(raw.data)


gcrma=exprs(data.gcrma.norm)
#gcrma=gcrma[1:12030,]
gcrma=gcrma[which(!grepl("AFFX", rownames(gcrma))),]

probes=row.names(gcrma)
symbol = unlist(mget(probes, hgu133aSYMBOL))
ID = unlist(mget(probes, hgu133aENTREZID))
gcrma=cbind(probes,ID,symbol,gcrma)
dim(gcrma)

###Removing genes which had NA values in their gene ID or Symbol
gcrma.new=gcrma[which(ID!="NA"|symbol!="NA"),]
dim(gcrma.new)

setwd(tempdir)
write.table(gcrma.new, file = "trainset_gcrma.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

GSE2034_clindata=getGSEDataTables("GSE2034")[[2]][1:286,]
write.table(GSE2034_clindata, "trainset_clindetails.txt", quote = FALSE, sep = "\t", row.names = FALSE)
