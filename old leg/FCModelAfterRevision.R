##DataFile: file:///C:/Users/18134/OneDrive%20-%20The%20University%20of%20Texas-Rio%20Grande%20Valley/data2/GEO%20DataSet%20Browser.html

##ftp://ftp.ncbi.nlm.nih.gov/geo/

###https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE35864 ##HIV Data
##https://nntc.org/content/gene_array
###SAM Analysis: https://mdozmorov.github.io/BIOS567/assets/presentation_diffexpression/DiffExpr_SAM.html
##Gene Filtering http://math.usu.edu/jrstevens/stat5570/3.2.Filtering.pdf

###Frontal Cortex Analysis
##File 2

# Clear plots
if(!is.null(dev.list())) dev.off()
# Clear console
cat("\014") 
# Clean workspace
rm(list=ls())

setwd("C:/Users/18134/OneDrive - The University of Texas-Rio Grande Valley/Dr.Upal/Bryan/MD/HIVData/ModelAfterRevision")
load("C:/Users/18134/OneDrive - The University of Texas-Rio Grande Valley/Dr.Upal/Bryan/MD/HIVData/ModelAfterRevision/HIVDataTrainTestDataCreation.RData")

load("C:/Users/18134/OneDrive - The University of Texas-Rio Grande Valley/Dr.Upal/Bryan/MD/HIVData/ModelAfterRevision/FrontalCortexAnalysis.RData")

###Uni
setwd("C:/Users/bas859/OneDrive - The University of Texas-Rio Grande Valley/Dr.Upal/Bryan/MD/HIVData/ModelAfterRevision")
load("HIVDataTrainTestDataCreation.RData")
load("FrontalCortexAnalysis.RData")

#Get potential predictor variables for Basal Ganglia
newpredictor_traindata_FC=t(filt_TrainData_FC)  ##dim(newpredictor_traindata_FC)=19 by 794
newpredictor_testdata_FC=t(filt_TestData_FC)  ##dim(newpredictor_testdata_FC)=5 by 794

#TrainXDem[1:3,]
#age.ch1 race.ch1   tissue.ch1
#2      53    White White matter
#3      63 Hispanic White matter
#4      58    White White matter

#TrainXDem[,1]=age.ch1
#TrainXDem[,2]=race.ch1
#is.numeric(TrainXDem$age.ch1)
#TRUE

newpredictor_traindata1_FC=data.frame(newpredictor_traindata_FC,as.numeric(TrainXDem[,1]),TrainXDem[,2],stringsAsFactors=FALSE)
newpredictor_testdata1_FC=data.frame(newpredictor_testdata_FC,as.numeric(TestXDem[,1]),TestXDem[,2],stringsAsFactors=FALSE)
##dim(newpredictor_traindata1_FC)= 19 by 796
##dim(newpredictor_testdata1_FC)= 5 by 796

newpredictor_traindata1GeneID_FC=data.frame(newpredictor_traindata_FC,as.numeric(TrainXDem[,1]),TrainXDem[,2], stringsAsFactors=FALSE)
newpredictor_testdata1GeneID_FC=data.frame(newpredictor_testdata_FC,as.numeric(TestXDem[,1]),TestXDem[,2],stringsAsFactors=FALSE)

############################
predictor_dataCombined_FC=data.frame(predictor_dataFC,clincData24Individual[,3],target, stringsAsFactors=FALSE)
colnames(predictor_dataCombined_FC)=c(colnames(predictor_dataCombined_FC[,-c(20409:20410)]),"Age","Target")
##20409:20410 are the varoavbels "Age" and "Target". In order to rename 20409 by Age and 20410 by Target we do this.
## colnames(predictor_dataCombined_FC[1:6,20405:20410])
##"X20416" "X20417" "X20418" "X20419" "Age"    "Target"
##########################################

#Get potential predictor variables


predictor_names_FC=c(as.vector(filt_genes_FC[,3])) #gene symbol for FC
predictor_namesID_FC=c(as.vector(filt_genes_FC[,2])) #gene ID for FC
geneAnnot_FC=cbind(predictor_names_FC,predictor_namesID_FC)

colnames(newpredictor_traindata1GeneID_FC)=c(predictor_namesID_FC,"Age","Ethnicity")
colnames(newpredictor_testdata1GeneID_FC)=c(predictor_namesID_FC,"Age","Ethnicity")

colnames(newpredictor_traindata1_FC)=c(predictor_names_FC,"Age","Ethnicity")
colnames(newpredictor_testdata1_FC)=c(predictor_names_FC,"Age","Ethnicity")


#> newpredictor_traindata1GeneID_FC[1:5,792:796]
#                         9935      9956       998 Age Ethnicity
#GSM876863_26G758.CEL 5.685668  9.523312  8.930271  53     White
#GSM876864_27G758.CEL 6.009926  8.118530  9.687163  63  Hispanic
#GSM876865_28G758.CEL 7.075740  8.381182  7.587669  58     White
#GSM876866_29G758.CEL 5.230102  7.219253  9.546612  34     White
#GSM876867_30G758.CEL 7.196796 10.787860 10.280818  48     White

#> newpredictor_traindata1_FC[1:5,792:796]
#                         MAFB    HS3ST2     CDC42 Age Ethnicity
#GSM876863_26G758.CEL 5.685668  9.523312  8.930271  53     White
#GSM876864_27G758.CEL 6.009926  8.118530  9.687163  63  Hispanic
#GSM876865_28G758.CEL 7.075740  8.381182  7.587669  58     White
#GSM876866_29G758.CEL 5.230102  7.219253  9.546612  34     White
#GSM876867_30G758.CEL 7.196796 10.787860 10.280818  48     White


dim(newpredictor_traindata1_FC) #19 by 796
dim(newpredictor_testdata1_FC)  #5 by 796
str(newpredictor_testdata1_FC[,792:796])##Checking whether Age variable is taken as numeric and Ethnicity as a fcator

##Testing whether the column names of testing data set is same as for training data set
which(colnames(newpredictor_testdata1_FC)!=colnames(newpredictor_traindata1_FC))
str(colnames(newpredictor_testdata1_FC[796]))

###Creating a combined dataset for FC with both training and Testing data sets
newpredictor_trainTestcombined_FC=rbind(newpredictor_traindata1_FC,newpredictor_testdata1_FC)
newpredictor_trainTestcombined_FC_SAMAnalysis=cbind(as.data.frame.array(t(colnames(newpredictor_traindata1GeneID_FC)),t(newpredictor_trainTestcombined_FC)))


save(newpredictor_traindata1_FC, file="C:/Users/18134/OneDrive - The University of Texas-Rio Grande Valley/Dr.Upal/Bryan/MD/HIVData/ModelAfterRevision/DataFilesnewpredictor_traindata1_FC")
write.table(newpredictor_traindata1_FC, "C:/Users/18134/OneDrive - The University of Texas-Rio Grande Valley/Dr.Upal/Bryan/MD/HIVData/ModelAfterRevision/newpredictor_traindata1_FC.txt", quote = FALSE, sep = "\t", row.names = FALSE)

save(newpredictor_traindata1GeneID_FC, file="C:/Users/18134/OneDrive - The University of Texas-Rio Grande Valley/Dr.Upal/Bryan/MD/HIVData/ModelAfterRevision/newpredictor_traindata1GeneID_FC")
write.table(newpredictor_traindata1GeneID_FC, "C:/Users/18134/OneDrive - The University of Texas-Rio Grande Valley/Dr.Upal/Bryan/MD/HIVData/ModelAfterRevision/newpredictor_traindata1GeneID_FC.txt", quote = FALSE, sep = "\t", row.names = FALSE)

save(newpredictor_testdata1_FC, file="C:/Users/18134/OneDrive - The University of Texas-Rio Grande Valley/Dr.Upal/Bryan/MD/HIVData/ModelAfterRevision/newpredictor_testdata1_FC")
write.table(newpredictor_testdata1_FC, "C:/Users/18134/OneDrive - The University of Texas-Rio Grande Valley/Dr.Upal/Bryan/MD/HIVData/ModelAfterRevision/newpredictor_testdata1_FC.txt", quote = FALSE, sep = "\t", row.names = FALSE)

save(newpredictor_testdata1GeneID_FC, file="C:/Users/18134/OneDrive - The University of Texas-Rio Grande Valley/Dr.Upal/Bryan/MD/HIVData/ModelAfterRevision/newpredictor_testdata1GeneID_FC")
write.table(newpredictor_testdata1GeneID_FC, "C:/Users/18134/OneDrive - The University of Texas-Rio Grande Valley/Dr.Upal/Bryan/MD/HIVData/ModelAfterRevision/newpredictor_testdata1GeneID_FC.txt", quote = FALSE, sep = "\t", row.names = FALSE)

save(newpredictor_trainTestcombined_FC, file="C:/Users/18134/OneDrive - The University of Texas-Rio Grande Valley/Dr.Upal/Bryan/MD/HIVData/ModelAfterRevision/newpredictor_trainTestcombined_FC")
write.table(newpredictor_trainTestcombined_FC,"C:/Users/18134/OneDrive - The University of Texas-Rio Grande Valley/Dr.Upal/Bryan/MD/HIVData/ModelAfterRevision/newpredictor_trainTestcombined_FC.txt", quote = FALSE, sep = "\t", row.names = FALSE)


save.image("C:/Users/18134/OneDrive - The University of Texas-Rio Grande Valley/Dr.Upal/Bryan/MD/HIVData/ModelAfterRevision/HIVDataTrainTestDataCreation.RData")

###################################################################################################
setwd("C:/Users/18134/OneDrive - The University of Texas-Rio Grande Valley/Dr.Upal/Bryan/MD/HIVData/ModelAfterRevision")
load("C:/Users/18134/OneDrive - The University of Texas-Rio Grande Valley/Dr.Upal/Bryan/MD/HIVData/ModelAfterRevision/HIVDataTrainTestDataCreation.RData")

###Multidimensional scaling plot of distances between gene expression profiles
BiocManager::install(version = "3.9")
BiocManager::install("limma")
library(limma)

namessamples_FC<-seq(838,861,by=1) ##Naming by GEO Accession Number
rownames(filteredDataAll_FC)=namessamples_FC

plotMDS(t(filteredDataAll_FC),col=c(rep("black",6), rep("orange",6),rep("purple",7),rep("red",5)), labels= c(rep("Control",6), rep("HIV",6),rep("HAD",7),rep("HIVE",5)))

###Log 2fold Calculation

targetAll<-as.factor(c(trainY,testY))
newpredictor_trainTestcombinedAll_FC=cbind(newpredictor_trainTestcombined_FC,targetAll)
colnames(newpredictor_trainTestcombinedAll_FC)=c(colnames(newpredictor_traindata1_FC),"target")

geneID_FC=c(predictor_namesID_FC)
genenames_FC=colnames(newpredictor_traindata1_FC)[-c(795,796)] ##Filtering only the gene names withoot columns 795:Age, 796:Ethnicity

####*********************************Significance Analysis of MocroArrays#################
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("impute")
BiocManager::install("samr")
library(samr)  #
library(DT)
library(limma)
# source('http://bioconductor.org/biocLite.R') # Import biocLite() function into R environment biocLite('limma')
##1_control
##2_HAD
##3_HIV
##4_HIVE
SAMMultiClass_data_FC <- list(x =t(newpredictor_trainTestcombined_FC[,-c(795,796)]), y = t(targetAll), geneid = geneID_FC, genenames = genenames_FC, logged2 = TRUE)

samr.objMultiClass_FC <- samr(SAMMultiClass_data_FC,resp.type="Multiclass",assay.type="array", nperms=100,center.arrays=TRUE)
names(samr.objMultiClass_FC)
delta.tableMultiClass_FC<-samr.compute.delta.table(samr.objMultiClass_FC)
delta.tableMultiClass_FC
del<- 0.1422
siggenes.table_FC<- samr.compute.siggenes.table(samr.objMultiClass_FC,del,SAMMultiClass_data_FC,delta.tableMultiClass_FC,min.foldchange=2)
siggenes.table_FC
samr.plot(samr.objMultiClass_FC, del, min.foldchange=2)
pv=samr.pvalues.from.perms(samr.objMultiClass_FC$tt,samr.objMultiClass_FC$ttstar)
pv
write.table(siggenes.table_FC$genes.up,"DataFiles/FCMultiClassup.txt")
write.table(siggenes.table_FC$genes.lo,"DataFiles/FCMultiClasslow.txt")


#######Pairwise SAM Analysis######
##Control vs HAD
FC_y12 <- c(rep(1, 5), rep(2, 5),rep(1, 1), rep(2, 2))
##Control vs HIV
##Actual FC_y13 <- c(rep(1, 5), rep(3, 5),rep(1, 1), rep(3, 1))
FC_y13 <- c(rep(1, 5), rep(2, 5),rep(1, 1), rep(2, 1))
##Control vs HIVE
##Actual FC_y14 <- c(rep(1, 5), rep(4, 4),rep(1, 5), rep(4, 1))
FC_y14 <- c(rep(1, 5), rep(2, 4),rep(1, 1), rep(2, 1))


##HIV vs HAD
#Actual FC_y32 <- c(rep(3, 5), rep(2, 5),rep(3, 1), rep(2, 2))
FC_y32 <- c(rep(1, 5), rep(2, 5),rep(1, 1), rep(2, 2))
##HAD vs HIVE
##Actual FC_y24 <- c(rep(2, 5), rep(4, 4),rep(2, 2), rep(4, 1))
FC_y24 <- c(rep(1, 5), rep(2, 4),rep(1, 2), rep(2, 1))
##HIV vs HIVE
##Actual FC_y34 <- c(rep(3, 5), rep(4, 4),rep(3, 1), rep(4, 1))
FC_y34 <- c(rep(1, 5), rep(2, 4),rep(1, 1), rep(2, 1))

SAM_data_FC12 <- list(x =t(newpredictor_trainTestcombined_FC[c(1:5,11:15,20,22:23),-c(795,796)]), y = t(FC_y12), geneid = geneID_FC, genenames = genenames_FC, logged2 = TRUE)
SAM_data_FC13 <- list(x =t(newpredictor_trainTestcombined_FC[c(1:5,6:10,20,21),-c(795,796)]), y = t(FC_y13),geneid = geneID_FC, genenames = genenames_FC, logged2 = TRUE)
SAM_data_FC14 <- list(x =t(newpredictor_trainTestcombined_FC[c(1:5,16:19,20,24),-c(795,796)]), y = t(FC_y14),geneid = geneID_FC, genenames = genenames_FC, logged2 = TRUE)

SAM_data_FC32 <- list(x =t(newpredictor_trainTestcombined_FC[c(6:10,11:15,21,22:23),-c(795,796)]), y = t(FC_y32),geneid = geneID_FC, genenames = genenames_FC, logged2 = TRUE)
SAM_data_FC24 <- list(x =t(newpredictor_trainTestcombined_FC[c(11:15,16:19,22:23,24),-c(795,796)]), y = t(FC_y24),geneid = geneID_FC, genenames = genenames_FC, logged2 = TRUE)
SAM_data_FC34 <- list(x =t(newpredictor_trainTestcombined_FC[c(6:10,16:19,21,24),-c(795,796)]), y = t(FC_y34),geneid = geneID_FC, genenames = genenames_FC, logged2 = TRUE)


##Checking for CIRBP geneID
which(genenames_FC=="BTN3A3") ##40
geneID[187]  ##1153

##1. Control vs HAD
samr.obj_FC12 <- samr(SAM_data_FC12, resp.type = "Two class unpaired",assay.type="array", nperms=100,center.arrays=TRUE,)
names(samr.obj_FC12)
delta.table_FC12 <- samr.compute.delta.table(samr.obj_FC12, min.foldchange = 2)  # Compute thresholds for different deltas
delta.table_FC12
del<- 0.3124
siggenes.table_FC12<- samr.compute.siggenes.table(samr.obj_FC12,del,SAM_data_FC12,delta.table_FC12,min.foldchange=2)
siggenes.table_FC12
write.table(siggenes.table_FC12$genes.up,"DataFiles/FC12up.txt")
write.table(siggenes.table_FC12$genes.lo,"DataFiles/FC12low.txt")
samr.plot(samr.obj_FC12, del, min.foldchange=2)

##2. Control vs HIV
samr.obj_FC13 <- samr(SAM_data_FC13, resp.type = "Two class unpaired",assay.type="array", nperms=100,center.arrays=TRUE)
names(samr.obj_FC13)
delta.table_FC13 <- samr.compute.delta.table(samr.obj_FC13, min.foldchange = 2)  # Compute thresholds for different deltas
delta.table_FC13
del<- 0.5175
siggenes.table_FC13<- samr.compute.siggenes.table(samr.obj_FC13,del,SAM_data_FC13,delta.table_FC13,min.foldchange=2)
siggenes.table_FC13
write.table(siggenes.table_FC13$genes.up,"DataFiles/FC13up.txt")
write.table(siggenes.table_FC13$genes.lo,"DataFiles/FC13low.txt")
samr.plot(samr.obj_FC13, del, min.foldchange=2)

##3. Control vs HIVE
samr.obj_FC14 <- samr(SAM_data_FC14, resp.type = "Two class unpaired",assay.type="array", nperms=100,center.arrays=TRUE)
names(samr.obj_FC14)
delta.table_FC14 <- samr.compute.delta.table(samr.obj_FC14, min.foldchange = 2)  # Compute thresholds for different deltas
delta.table_FC14
del<- 1.0047
siggenes.table_FC14<- samr.compute.siggenes.table(samr.obj_FC14,del,SAM_data_FC14,delta.table_FC14,min.foldchange=2)
siggenes.table_FC14
write.table(siggenes.table_FC14$genes.up,"DataFiles/FC14up.txt")
write.table(siggenes.table_FC14$genes.lo,"DataFiles/FC14low.txt")
samr.plot(samr.obj_FC14, del, min.foldchange=2)

##4.HIV vs HAD
samr.obj_FC32 <- samr(SAM_data_FC32, resp.type = "Two class unpaired",assay.type="array", nperms=100,center.arrays=TRUE)
names(samr.obj_FC32)
delta.table_FC32 <- samr.compute.delta.table(samr.obj_FC32, min.foldchange = 2)  # Compute thresholds for different deltas
delta.table_FC32
del<-0.2418
siggenes.table_FC32<- samr.compute.siggenes.table(samr.obj_FC32,del,SAM_data_FC32,delta.table_FC32,min.foldchange=2)
siggenes.table_FC32
write.table(siggenes.table_FC32$genes.up,"DataFiles/FC32up.txt")
write.table(siggenes.table_FC32$genes.lo,"DataFiles/FC32low.txt")
samr.plot(samr.obj_FC32, del, min.foldchange=2)

##5.HAD vs HIVE
samr.obj_FC24 <- samr(SAM_data_FC24, resp.type = "Two class unpaired",assay.type="array", nperms=100,center.arrays=TRUE)
names(samr.obj_FC24)
delta.table_FC24 <- samr.compute.delta.table(samr.obj_FC24, min.foldchange = 2)  # Compute thresholds for different deltas
delta.table_FC24
del<- 0.8030
siggenes.table_FC24<- samr.compute.siggenes.table(samr.obj_FC24,del,SAM_data_FC24,delta.table_FC24,min.foldchange=2)
siggenes.table_FC24
write.table(siggenes.table_FC24$genes.up,"DataFiles/FC24up.txt")
write.table(siggenes.table_FC24$genes.lo,"DataFiles/FC24low.txt")
samr.plot(samr.obj_FC24, del, min.foldchange=2)

##6.HIV vs HIVE
samr.obj_FC34 <- samr(SAM_data_FC34, resp.type = "Two class unpaired",assay.type="array", nperms=100,center.arrays=TRUE)
names(samr.obj_FC34)
delta.table_FC34 <- samr.compute.delta.table(samr.obj_FC34, min.foldchange = 2)  # Compute thresholds for different deltas
delta.table_FC34
del<-0.4250
siggenes.table_FC34<- samr.compute.siggenes.table(samr.obj_FC34,del,SAM_data_FC34,delta.table_FC34,min.foldchange=2)
siggenes.table_FC34
write.table(siggenes.table_FC34$genes.up,"DataFiles/FC34up.txt")
write.table(siggenes.table_FC34$genes.lo,"DataFiles/FC34low.txt")
samr.plot(samr.obj_FC34, del, min.foldchange=2)

save.image("C:/Users/18134/OneDrive - The University of Texas-Rio Grande Valley/Dr.Upal/Bryan/MD/HIVData/ModelAfterRevision/HIVDataTrainTestDataCreation.RData")

#######################################################
load("C:/Users/18134/OneDrive - The University of Texas-Rio Grande Valley/Dr.Upal/Bryan/MD/HIVData/ModelAfterRevision/HIVDataTrainTestDataCreation.RData")
###################
##http://rstudio-pubs-static.s3.amazonaws.com/14943_86bc7b5d3b7547b0a3a3a359b10afcf4.html
splom(t(newpredictor_trainTestcombinedAll), panel = panel.smoothScatter, raster = TRUE)

##Correlation Analysis for White MAtter
install.packages("gplots")
library(gplots)
install.packages("RColorBrewer")
library(RColorBrewer)

hmcol <- colorRampPalette(brewer.pal(10,"PuOr"))(342)
newpredictor_trainTestcombined_FC[c(1:5,16:19,20,24),-c(795,796)]
corrs <- cor(t(newpredictor_trainTestcombined_FC[,-c(795,796)]), method = "spearman")
heatmap.2(corrs, margins = c(6, 11),col = hmcol, trace = "none", main="Correlation Analysis Frontal Cortex",keysize = 1,Rowv=FALSE,Colv=FALSE,dendrogram='none', 
          srtCol = 56)


###Principal Component Analysis for Frontal Cortex
prcomps_FC<- prcomp(t(newpredictor_trainTestcombined_FC[,-c(795,796)]))
summary(prcomps_FC)
plot(prcomps_FC$rotation, type = "n",main="Principal Components for Frontal Cortex")
text(prcomps_FC$rotation, rownames(prcomps_FC$rotation))

install.packages("rgl")
library(rgl)
plot3d(prcomps_FC$rotation[,1:3],col=c(rep(1,5),rep(3,5),rep(2,5),rep(4,4),c(1,3,2,2,4)))###loadings

install.packages("factoextra")
library(factoextra)
fviz_pca_ind(prcomps_FC,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)
colors<-rainbow(4)
#colors <- c('#999999', '#E69F00', '#56B4E9','#1004F0')
colors <- colors[as.numeric(targetAll)]


###col=c(rep(1,5),rep(3,5),rep(2,5),rep(4,4),c(1,3,2,2,4)) , Here the colors are assigned in the order the groups apppear in TargetAll data vector

plot3d(prcomps_FC$rotation[,1:3],col=c(rep(1,5),rep(3,5),rep(2,5),rep(4,4),c(1,3,2,2,4)),type="s",image=TRUE,xlab="PC1",ylab="PC2",zlab="PC3")###loadings
legend3d("topright", legend = levels(targetAll), pch = 16, col = levels(targetAll), cex=1, inset=c(0.01))



snapshot3d(filename = "C:/Users/18134/OneDrive - The University of Texas-Rio Grande Valley/Dr.Upal/Bryan/MD/HIVData/ModelAfterRevision/Graphs/GraphsPAnalysis3DFC.png", fmt = 'png')


###Verimax Rotation on White Matter
rawLoadings_FC     <- prcomps_FC$rotation[,1:3] %*% diag(prcomps$sdev, 3, 3)
rotatedLoadings_FC<- varimax(rawLoadings_FC)$loadings
rotatedLoadings_FC

colors<-rainbow(4)
#colors <- c('#999999', '#E69F00', '#56B4E9','#1004F0')
colors <- colors[as.numeric(targetAll)]

plot3d(rotatedLoadings_FC[,1:3],col=c(rep(1,5),rep(3,5),rep(2,5),rep(4,4),c(1,3,2,2,4)),type="s",image=TRUE,xlab="R1",ylab="R2",zlab="R3")###loadings
legend3d("topright", legend = levels(targetAll), pch = 20, col = levels(targetAll), cex=1, inset=c(0.01))

snapshot3d(filename ="C:/Users/18134/OneDrive - The University of Texas-Rio Grande Valley/Dr.Upal/Bryan/MD/HIVData/ModelAfterRevision/Graphs/RotatedPAnalysis3DFC.png", fmt = 'png')

########################################################################################
##Analysing Differentially Expressed Genes in Frontal Cortex with Moderated t-test
########################################################################################
design1 <- model.matrix(~0+targetAll)
colnames(design1)<-c("Control","HAD","HIV","HIVE")
array_FC<-arrayWeights(t(newpredictor_trainTestcombined_FC[,-c(795,796)]),design=design1)
fit1_FC <- lmFit(t(newpredictor_trainTestcombined_FC[,-c(795,796)]), design1,weights=array_FC)
dim(design1)
dim(fit1_FC)
contrast.matrix<-makeContrasts(HAD-Control,HIV-Control,HIVE-Control,HAD-HIV,HIVE-HAD,HIVE-HIV,levels=design1)

fit1_2_FC<-contrasts.fit(fit1_FC,contrast.matrix)

fit1_2_FC <- eBayes(fit1_2_FC)
FCT1<-topTable(fit1_2_FC,coef=1,adjust="BH",n="Inf") # A list of top genes differential expressed in HAD vs control
write.table(FCT1,"FCT1.txt")
FCT2<-topTable(fit1_2_FC,coef=2,adjust="BH",n="Inf") # A list of top genes differential expressed in HIV vs control
write.table(FCT2,"FCT2.txt")
FCT3<-topTable(fit1_2_FC,coef=3,adjust="BH",n="Inf") # A list of top genes differential expressed in HIVE vs control
write.table(FCT3,"FCT3.txt")
FCT4<-topTable(fit1_2_FC,coef=4,adjust="BH",n="Inf") # A list of top genes differential expressed in HAD vs HIV
write.table(FCT4,"FCT4.txt")
FCT5<-topTable(fit1_2_FC,coef=5,adjust="BH",n="inf") # A list of top genes differential expressed in HIVE vs HAD
write.table(FCT5,"FCT5.txt")
FCT6<-topTable(fit1_2_FC,coef=6,adjust="BH",n="inf") # A list of top genes differential expressed in HIVE vs HIV
write.table(FCT6,"FCT6.txt")

##Don't use the follwing lfc cut off as it is not recomended
#https://www.researchgate.net/post/Cut-off_values_for_gene_expression_fold_change_when_performing_RNA_seq10
#results_FC = decideTests(fit1_2_FC, adj.P.Val=0.05, lfc=1.1)

##Instead we have used this in model after revision
results_FC = decideTests(fit1_2_FC, adj.P.Val=0.01)
write.table(results_FC,"FC_eBayes.txt")
##Venn Diagram Creation for Up and Down Genes
vennDiagram(results_FC[,c(1:3)],main="DE Genes in Frontal Cortex",include=c("up", "down"),
            counts.col=c("red", "blue"),
            circle.col = c("red", "blue", "green3"))
#################################################################################

volcanoplot(fit1_2_FC)

#results = decideTests(fit2_2, adj.P.Val=0.05)
#vennDiagram(results)

plotMD(fit1_2_FC,coef=3,status =results1[,3],values=c(-1,1),hl.col=c("red","blue"))


#***********************************************************************************


##Checking for RBM3 geneID
which(genenames=="RBM3")# 1229
geneID[1229] ## 5935

##Checking for CIRBP geneID
which(genenames_FC=="CIRBP") ##187
geneID[187]  ##1153

##To identify the DE genes in all three sectors/all tw in FC
which(results3_1[,1]==1 & results3_1[,2]==1 & results3_1[,3]==1)
which(results3_1[,1]==1 & results3_1[,3]==1) ## HAD - Control and HIVE - Control


save.image("C:/Users/18134/OneDrive - The University of Texas-Rio Grande Valley/Dr.Upal/Bryan/MD/HIVData/ModelAfterRevision/FrontalCortexAnalysis.RData")
