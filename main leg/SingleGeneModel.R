load("/Users/myrinearevalo/Desktop/MachineLearn/BreastAnalysis2.RData")
# Look at variable importance 
SingleGenes_imp.temp <- abs(rf_output$importance[,3]) 
SingleGenes_t <- order(SingleGenes_imp.temp,decreasing=TRUE) 

dev.off()
plot(c(1:ncol(newpredictor_traindata)),SingleGenes_imp.temp[SingleGenes_t],log='x',cex.main=1.5, xlab='Gene Rank',ylab='Variable Importance (Mean Decrease Accuracy)',cex.lab=1, pch=16,main='Single Genes RF Model') 

# Get subset of expression values for 25 most 'important' genes 
SingleGenes_gn.imp <- names(SingleGenes_imp.temp)[SingleGenes_t] 
SingleGenes_gn.25 <- SingleGenes_gn.imp[1:25] # vector of top 25 genes, in order 
# [1] "22795_at" "4320_at"  "27250_at" "79679_at" "9474_at"  "64699_at" "79682_at" "55110_at" "81611_at" "10635_at"
# [11] "4312_at"  "9124_at"  "7832_at"  "4604_at"  "8519_at"  "51659_at" "51280_at" "10057_at" "6790_at"  "7041_at" 
# [21] "6491_at"  "1289_at"  "28908_at" "3507_at"  "4591_at" 
SingleGenes_t<- is.element(colnames(newpredictor_traindata),SingleGenes_gn.25)

SingleGenes_sig.eset <- newpredictor_traindata[,SingleGenes_t]##Dimension 230 by 25

##Creating a dataset with geneId as column names
newpredictor_traindataGeneId<-newpredictor_traindata
colnames(newpredictor_traindataGeneId)<-geneAnnot[,1]

SingleGenes_sig.eset2<-newpredictor_traindataGeneId[,SingleGenes_t]
SingleGenes_sig.esetNew<-t(SingleGenes_sig.eset2)##Dimenstion 25 by 230
str(SingleGenes_sig.esetNew)
# num [1:25, 1:230] 6.72 6.31 9.22 8.64 8.89 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:25] "10057_at" "10635_at" "1289_at" "22795_at" ...
# ..$ : chr [1:230] "GSM36777.CEL" "GSM36779.CEL" "GSM36783.CEL" "GSM36784.CEL" ...
dim(SingleGenes_sig.esetNew)##Dimension 25 by 230

# matrix of expression values, not necessarily in order ## Make a heatmap, with group differences obvious on plot 
library(RColorBrewer) 
hmcol <- colorRampPalette(brewer.pal(11,"PuOr"))(256) 

colnames(SingleGenes_sig.esetNew) <- trainY # This will label the heatmap columns 

csc <- rep(hmcol[50],50) 
csc[colnames(SingleGenes_sig.esetNew)=='Relapse'] <- hmcol[200] # column side color will be purple for T and orange for B 
SingleGenes_sig.esetNew2<-SingleGenes_sig.esetNew[ ,order(colnames(SingleGenes_sig.esetNew))]


#reorder dataset according to column names
SingleGenes_sig.esetNew

install.packages("gplots")
library(gplots)
stats::lowess

dev.off()
heatmap.2(SingleGenes_sig.esetNew2,scale="row",col=hmcol,ColSideColors=csc[1:230],main="Top 25 Genes Using Single Genes RF",trace='none',Rowv=FALSE,Colv=FALSE,dendrogram='none')


dev.off()
heatmap.2(SingleGenes_sig.esetNew2[,c(1:25,206:230)],scale="row",col=hmcol,ColSideColors=csc[1:50],main="Top 25 Genes Using Single Genes RF",trace='none',Rowv=FALSE,Colv=FALSE,dendrogram='none')

save.image("/Users/myrinearevalo/Desktop/MachineLearn/SingleGeneHeatmapdata.RData")
