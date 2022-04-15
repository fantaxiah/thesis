##DataFile: file:///C:/Users/18134/OneDrive%20-%20The%20University%20of%20Texas-Rio%20Grande%20Valley/data2/GEO%20DataSet%20Browser.html
###To get custom  CDf file: http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/CDF_download.asp
##ftp://ftp.ncbi.nlm.nih.gov/geo/
###https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE2034 ##breast Data
###Train Data GSE 2034
###Test Data GSE 2990
##https://nntc.org/content/gene_array
###Path Local Computer
setwd("/Users/myrinearevalo/Desktop/MachineLearn")

load ("/Users/myrinearevalo/Desktop/MachineLearn/R Analysis/Breast(2)/BreastAnalysis2.RData")

#### Creating the Random Forest Model###
install.packages("randomForest")
install.packages("ROCR")
install.packages("Hmisc")
install.packages("ggplot2")
install.packages("checkmate")
BiocManager::install("genefilter")

library(randomForest)
library(ROCR)
library(genefilter)
library(ggplot2)
library(checkmate)
library(Hmisc)

setwd("G:/MD/MachineLearn")

datafile="trainset_gcrma(2).txt" 
clindatafile="trainset_clindetails(2).txt"

outfile="trainset_RFoutput.txt"
varimp_pdffile="trainset_varImps.pdf"
MDS_pdffile="trainset_MDS.pdf"
ROC_pdffile="trainset_ROC.pdf"
case_pred_outfile="trainset_CasePredictions.txt"
vote_dist_pdffile="trainset_vote_dist.pdf"

#Read in data (expecting a tab-delimited file with header line and rownames)
data_import=read.table(datafile, header = TRUE, na.strings = "NA", sep="\t") ##dim(data_import)=12030 by  289
clin_data_import=read.table(clindatafile, header = TRUE, na.strings = "NA", sep="\t") ##dim(clin_data_import)=286  by 7

#NEED TO GET CLINICAL DATA IN SAME ORDER AS GCRMA DATA
clin_data_order=order(clin_data_import[,"GEO.asscession.number"])
clindata=clin_data_import[clin_data_order,]
data_order=order(colnames(data_import)[4:length(colnames(data_import))])+3 #Order data without first three columns, then add 3 to get correct index in original file
rawdata=data_import[,c(1:3,data_order)] #grab first three columns, and then remaining columns in order determined above
header=colnames(rawdata)

rawdata=rawdata[which(!is.na(rawdata[,3])),] #Remove rows with missing gene symbol; dim(rawdata)= 12029  by 289

##We will assign the variables to a new data structure. 
##Extract just the expression values from the filtered data and transpose the matrix. 
##The latter is necessary because RandomForests expects the predictor variables (genes) to be represented as columns instead of rows. 
##Finally, assign gene symbol as the predictor name.

predictor_data=t(rawdata[,4:length(header)]) ##dim(predictor_data)=286 by 12029

predictor_names=c(as.vector(rawdata[,3])) #gene symbol , length(predictor_names)=12029
predictor_probenames=c(as.vector(rawdata[,1])) #Probe ID
#colnames(predictor_data)=predictor_names
colnames(predictor_data)=predictor_probenames
predictor_data[1:10,1:10]

#Get target variable and specify as factor/categorical
target= clindata[,"relapse..1.True."]
target[target==0]="NoRelapse"
target[target==1]="Relapse"
target=as.factor(target)

####Creating the training and Testing datasets
install.packages("ggplot2")
install.packages("caret")
library(ggplot2)
library(caret)
library("RColorBrewer")
set.seed(3456)
trainIndex <- createDataPartition(target, p = .80, list = FALSE, times = 1)
head(trainIndex)

Train2 <- predictor_data[ trainIndex,]
Test2  <- predictor_data[-trainIndex,]
dim(Train2)  #dim(Train2)=230 by 12029
dim(Test2) # dim(Test2)=56 by 12029

library(randomForest) 
set.seed(1234) 

trainX=t(Train2)#dim(trainX)=12029 by 230
testX=t(Test2) #dim(testX)=12029 by 56
trainY=target[trainIndex]  
testY=target[-trainIndex]

###Next we filter out any variables (genes) that are not expressed or 
##do not have enough variance to be informative in classification. 
##We will first take the values and un-log2 them, then filter out any genes according to following criteria (recommended in multtest/MTP documentation): 
##(1) At least 20% of samples should have raw intensity greater than 100; 
##(2) The coefficient of variation (sd/mean) is between 0.7 and 10.
genefilter::Anova
genefilter::AUC

ffun=filterfun(pOverA(p = 0.2, A = 100), cv(a = 0.7, b = 10))
filt=genefilter(2^trainX,ffun)
filt_TrainData=trainX[filt,]  # dim(filt_TrainData)= 1246 by  230
filt_TestData=testX[filt,]  # dim(filt_TestData)= 1246 by  56
filt_genes=rawdata[filt,1:3]
prbname=rownames(filt_TrainData)
filteredDataAll<-predictor_data[,filt] #dim 286 by 1246

save(filt_TrainData,file="FilteredTrainDataBreast")
write.table(filt_TrainData,"FilteredTrainDataBreast.txt")
save(filt_TestData,file="FilteredTestDataBreast")
write.table(filt_TestData,"FilteredTestDataBreast.txt")

#Get potential predictor variables
newpredictor_traindata=t(filt_TrainData)
newpredictor_testdata=t(filt_TestData)

##Take the un-logged not filtered data to see the effect of the filter
unloggedTrainX=2^trainX
gene.mean<-apply(unloggedTrainX,1,mean)
gene.sd<-apply(unloggedTrainX,1,sd)
gene.cv<-gene.sd/gene.mean

par(mfrow=c(1,2))
plot.new()
par(mar=c(1,1,1,1))
plot(gene.mean,gene.sd,log='xy',pch=16,cex=0.1,main="Look at Mean of Genes")
abline(v=100,lwd=3,col='red')
dev.off()
hist(log2(gene.cv),main="Look at CV")
abline(v=log(.7),lwd=3,col='red')


###Just to see how the unloged filtered data 
filt_TrainDataunloged=2^filt_TrainData
genefiltered.mean1<-apply(filt_TrainDataunloged,1,mean)
genefiltered.sd1<-apply(filt_TrainDataunloged,1,sd)
genefiltered.cv1 <-genefiltered.sd1/genefiltered.mean1

dev.off()
plot(genefiltered.mean1,genefiltered.sd1,log='xy',pch=16,cex=0.1,main="After Filtering Gene expression>100")
hist(log2(genefiltered.cv1),main="After Filtering CV 0.7-10")


newpredictor_traindata=t(filt_TrainData)
newpredictor_testdata=t(filt_TestData)

predictor_probes=c(as.vector(filt_genes[,1]))
predictor_names=c(as.vector(filt_genes[,3])) #gene symbol
predictor_namesID=c(as.vector(filt_genes[,2])) #gene ID
geneAnnot=cbind(predictor_names,predictor_namesID,predictor_probes)


###Multidimensional scaling plot of distances between gene expression profiles
BiocManager::install(version = "3.9")
BiocManager::install("limma")
library(limma)

table(testY)
table(trainY)
plotMDS(t(filteredDataAll),col=c(rep("black",179), rep("red",107)), labels= c(rep("NoRelapsed",179), rep("Relapsed",107)))


###Creating a combined dataset with both training and Testing data sets
newpredictor_trainTestcombined=rbind(newpredictor_traindata,newpredictor_testdata)
newpredictor_trainTestcombined_SAMAnalysis=cbind(as.data.frame.array(predictor_namesID),newpredictor_trainTestcombined)
###Log 2fold Calculation

targetAll<-as.factor(c(trainY,testY))
newpredictor_trainTestcombinedAll=cbind(newpredictor_trainTestcombined,targetAll)
colnames(newpredictor_trainTestcombinedAll)=c(colnames(newpredictor_traindata),"target")

geneID=c(predictor_namesID)
genenames=c(predictor_names) ##Filtering only the gene names 

#################Significance Analysis of MicroArrays#################
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("impute")
install.packages("samr")
library(samr)  #
library(DT)
library(limma)
# source('http://bioconductor.org/biocLite.R') # Import biocLite() function into R environment biocLite('limma')

SAMTwoClass_data<- list(x =t(newpredictor_trainTestcombined), y = t(targetAll), geneid = geneID, genenames = genenames, logged2 = TRUE)
samr.obj<- samr(SAMTwoClass_data, resp.type = "Two class unpaired",assay.type="array", nperms=100,center.arrays=TRUE)
names(samr.obj)
delta.table <- samr.compute.delta.table(samr.obj, min.foldchange = 2)  # Compute thresholds for different deltas
delta.table
del<- 0.006
siggenes.table<- samr.compute.siggenes.table(samr.obj,del,SAMTwoClass_data,delta.table,min.foldchange=2)
siggenes.table
samr.plot(samr.obj, del, min.foldchange=2, )
pv=samr.pvalues.from.perms(samr.obj$tt,samr.obj$ttstar)
pv
write.table(siggenes.table$genes.up,file="BREASTup.txt")
write.table(siggenes.table$genes.lo,file="BREASTlow.txt")


##Correlation Analysis for White MAtter
install.packages("gplots")
library(gplots)
install.packages("RColorBrewer")
library(RColorBrewer)

hmcol <- colorRampPalette(brewer.pal(10,"PuOr"))(342)
corrs <- cor(t(newpredictor_trainTestcombined), method = "spearman")
new.plot()
dev.off()
heatmap.2(corrs, margins = c(6, 11),col = hmcol, trace = "none", main="Correlation Analysis",keysize = 1,Rowv=FALSE,Colv=FALSE,dendrogram='none', 
          srtCol = 56)
###Pricipal Component Analysis
prcomps<- prcomp(t(newpredictor_trainTestcombined))
summary(prcomps)
plot(prcomps$rotation, type = "n",main="Principal Components")
text(prcomps$rotation, rownames(prcomps$rotation))


install.packages("rgl")
library(rgl)
plot3d(prcomps$rotation[,1:3],col=levels(targetAll))###loadings


###Do the file creation .phe,.ann. and .exp file creations based on the GXNA tool 
###Save the files at C:\Users\18134\Desktop\gxna
##load the mm vector from the R file called "HIVDataGXNAInputFilesR.R" saved at C:\Users\18134\Desktop\gxna

### Load the initial clusters from GXNA into an array called "mm[[]]"
load("G:/MD/HIVData/gxnaTrainClusters.RData")

##Finally we run the RF algorithm. NOTE: we use an ODD number for ntree. This is because when the forest/ensembl is used on test data, ties are broken randomly. 
##Having an odd number of trees avoids this issue and makes the model fully deterministic. 
##Also note, we will use down-sampling to attempt to compensate for unequal class-sizes (less ##relapses than non-relapses).
#Setting up inputs for random forest function

#tmp holds values of relapse status, 179 No Relapse 107 Relapse
tmp = as.vector(table(trainY)) #dimensions 179 by 107
#number of classes (NoRelapse and Relapse) is 2
num_classes = length(tmp)
#sets Relapse class as minimum size of class (107)
min_size = tmp[order(tmp,decreasing=FALSE)[1]]
#replicares the values of x vector (min_size) into an object same type as x (num_classes)
sampsizes = rep(min_size,num_classes)

cutoff=c(0.75,0.25)
rf_output=randomForest(x=newpredictor_traindata, y=trainY, importance = TRUE, ntree = 10001, proximity=TRUE, sampsize=sampsizes, na.action = na.omit, mtry=6)

######
confusion=rf_output$confusion #takes values from confusion matrix in rf_output

#divides Relapse*Relapse by (Relapse*Relapse + Relapse*NoRelapse) and multiplies by 100 to get percentage
sensitivity=(confusion[2,2]/(confusion[2,2]+confusion[2,1]))*100
#divides NoRelapse*NoRelapse by (NoRelapse*NoRelapse + NoRelapse*Relapse) and multiplies by 100 to get percentage
specificity=(confusion[1,1]/(confusion[1,1]+confusion[1,2]))*100
#OOB error 34% --> 33.47826%
overall_error=rf_output$err.rate[length(rf_output$err.rate[,1]),1]*100
overall_accuracy=1-overall_error #1-34%=66% --> 1-33.47826%=66.52174%
class1_error=paste(rownames(confusion)[1]," error rate= ",confusion[1,3], sep="")
class2_error=paste(rownames(confusion)[2]," error rate= ",confusion[2,3], sep="")
overall_accuracy=100-overall_error

sens_out=paste("sensitivity=",sensitivity, sep="")
spec_out=paste("specificity=",specificity, sep="")
err_out=paste("overall error rate=",overall_error,sep="")
acc_out=paste("overall accuracy=",overall_accuracy,sep="")
misclass_1=paste(confusion[1,2], rownames(confusion)[1],"misclassified as", colnames(confusion)[2], sep=" ")
misclass_2=paste(confusion[2,1], rownames(confusion)[2],"misclassified as", colnames(confusion)[1], sep=" ")
confusion_out=confusion[1:2,1:2]
confusion_out=cbind(rownames(confusion_out), confusion_out)

write.table(rf_importances[,4],file=outfile, sep="\t", quote=FALSE, col.names=FALSE)
write("confusion table", file=outfile, append=TRUE)
write.table(confusion_out,file=outfile, sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE, append=TRUE)
write(c(sens_out,spec_out,acc_out,err_out,class1_error,class2_error,misclass_1,misclass_2), file=outfile, append=TRUE)

pdf(file=varimp_pdffile)
varImpPlot(rf_output, type=2, n.var=30, scale=FALSE, main="Variable Importance (Gini) for top 30 predictors")

pdf(file=MDS_pdffile)
target_labels=as.vector(target)
target_labels[target_labels=="NoRelapse"]="N"
target_labels[target_labels=="Relapse"]="R"
MDSplot(rf_output, target, k=2, xlab="", ylab="", pch=target_labels, palette=c("red", "blue"), main="MDS plot")

predictions=as.vector(rf_output$votes[,2])
pred=prediction(predictions,target)
#First calculate the AUC value
perf_AUC=performance(pred,"auc")
AUC=perf_AUC@y.values[[1]]
#Then, plot the actual ROC curve
perf_ROC=performance(pred,"tpr","fpr")
pdf(file=ROC_pdffile)
plot(perf_ROC, main="ROC plot")
text(0.5,0.5,paste("AUC = ",format(AUC, digits=5, scientific=FALSE)))

options(digits=2)
pdf(file=vote_dist_pdffile)
out <- histbackback(split(rf_output$votes[,"Relapse"], target), probability=FALSE, xlim=c(-50,50), main = 'Vote distributions for patients classified by RF', axes=TRUE, ylab="Fraction votes (Relapse)")
barplot(-out$left, col="red" , horiz=TRUE, space=0, add=TRUE, axes=FALSE)
barplot(out$right, col="blue", horiz=TRUE, space=0, add=TRUE, axes=FALSE)


case_predictions=cbind(clindata,target,rf_output$predicted,rf_output$votes)
write.table(case_predictions,file=case_pred_outfile, sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)


########################################################################################
##Analysing Differentially Expressed Genes in Frontal Cortex with Moderated t-test
########################################################################################
limma::plotMA

design1 <- model.matrix(~0+targetAll)
colnames(design1)<-c("NoRelapse", "Relapse")
array<-arrayWeights(t(newpredictor_trainTestcombined),design=design1)
fit1 <- lmFit(t(newpredictor_trainTestcombined),design1,weights=array)
dim(design1) #dim 286 by 2
dim(fit1) #dim 1246 by 2
contrast.matrix<-makeContrasts(Relapse-NoRelapse,levels=design1)

fit1_2<-contrasts.fit(fit1,contrast.matrix)

fit1_2<- eBayes(fit1_2)
T1<-topTable(fit1_2,coef=1,adjust="BH",n="Inf") # A list of top genes differential expressed in Relapse vs NoRelapse
write.table(T1,"T1.txt")

##Don't use the follwing lfc cut off as it is not recomended
#https://www.researchgate.net/post/Cut-off_values_for_gene_expression_fold_change_when_performing_RNA_seq10
#results_FC = decideTests(fit1_2_FC, adj.P.Val=0.05, lfc=1.1)

##Instead we have used this in model after revision
results = decideTests(fit1_2, adj.P.Val=0.01)
tail(results)
write.table(results,"eBayes.txt")
##Venn Diagram Creation for Up and Down Genes
vennDiagram(results,main="DE Genes",include=c("up", "down"),
            counts.col=c("red", "blue"),
            circle.col = c("red", "blue", "green3"))
#################################################################################

volcanoplot(fit1_2)

#results = decideTests(fit2_2, adj.P.Val=0.05)
#vennDiagram(results)

plotMD(fit1_2,coef=3,status =results1[,3],values=c(-1,1),hl.col=c("red","blue"))


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

save.image("/Users/myrinearevalo/Desktop/MachineLearn/BreastAnalysis2.RData")
