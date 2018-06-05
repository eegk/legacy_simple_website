##############################
### CURSO RNAseq
### Edgar & Gabi
### May - June 2018
##############################
#DATASET: GSE107218 #per.blood #illumina hiseq250
##############################
### INSTALL BIOCONDUCTOR PACKAGES
##############################
source("https://bioconductor.org/biocLite.R")
biocLite()
biocLite("ROTS")
biocLite("limma")
biocLite("pheatmap")
biocLite("gplots")
biocLite("edgeR")
biocLite("RColorBrewer")
##############################
### LOAD PACKAGES
##############################
library(pheatmap)
library(ROTS)
library(limma)
library(gplots) 
library(edgeR)
library(RColorBrewer)
##############################
### LOAD DATA
##############################
dados<-read.table("GSE107218_CBPB-hg19-counts.txt",sep="\t",header=TRUE)
### groups
group <- as.factor(c(rep("CB_CD34",3),rep("CB_BFUE",3),rep("CB_CFUE",3),rep("CB_PRO",3),
                     rep("CB_EBASO",3),rep("CB_LB",3),rep("CB_POLY",3),rep("CB_ORTHO",3),
                     rep("PB_CD34",3),rep("PB_BFU",3),rep("PB_CFU",3),rep("PB_PRO",3),
                     rep("PB_EBASO",3),rep("PB_LB",3),rep("PB_POLY",3),rep("PB_ORTHO",3)))
##############################
### CPM https://doi.org/10.1186/gb-2014-15-2-r29
##############################
x<-dados[,7:54]
rownames(x)<-dados$Geneid
cpm <- cpm(x)
lcpm <- cpm(x, log=TRUE)
### filter
keep.exprs <- rowSums(cpm>1)>=3
x <- x[keep.exprs,]
##############################
### plot to compare to unfiltered data
##############################
nsamples <- ncol(x)
col <- brewer.pal(nsamples, "Paired")
par(mfrow=c(1,2))
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.21), las=2, 
     main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", legend=group, text.col=col, bty="n")
###
lcpm <- cpm(x, log=TRUE)
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.21), las=2, 
     main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", legend=group, text.col=col, bty="n")
##############################
### normalize modified data
##############################
d.cpm.x <- DGEList(counts=x,group=group) #Create a DGEList object
d.cpm.x <- calcNormFactors(d.cpm.x, method = "TMM") #method TMM
#d.cpm.x$samples$norm.factors # to see the factors
##############################
### comparison to unormalized data
##############################
### makes a not normalized dataset
d.cpm.x2 <- d.cpm.x
d.cpm.x2$samples$norm.factors <- 1
d.cpm.x2$counts[,1] <- ceiling(d.cpm.x2$counts[,1]*0.05)
d.cpm.x2$counts[,2] <- d.cpm.x2$counts[,2]*5
##############################
### plot the unnormalized data and the normalized together
##############################
par(mfrow=c(1,2))
lcpm <- cpm(d.cpm.x2, log=TRUE)
boxplot(lcpm, las=2, col=col, main="")
title(main="A. Example: Unnormalised data",ylab="Log-cpm")
d.cpm.x2 <- calcNormFactors(d.cpm.x2,method = "TMM")
d.cpm.x2$samples
lcpm <- cpm(d.cpm.x2, log=TRUE)
boxplot(lcpm, las=2, col=col, main="")
title(main="B. Example: Normalised data",ylab="Log-cpm")
##############################
### MDS
##############################
lcpm <- cpm(d.cpm.x, log=TRUE)
plotMDS(lcpm, labels=group, col=as.numeric(group))
title(main="MDS - Sample groups")
##############################
### Voom
##############################
design = model.matrix( ~ 0 + group, data=d.cpm.x$samples)
colnames(design) <- levels(group)
d.cpm.x = estimateCommonDisp(d.cpm.x, verbose=TRUE)
d.cpm.x = estimateTagwiseDisp(d.cpm.x)
par(mfrow=c(1,2))
v <- voom(d.cpm.x, design, plot=TRUE)
##############################
### CONTRAST MATRIX & DE
##############################
contr.matrix <- makeContrasts(
  CB_CD34vsPB_CD34 = CB_CD34 - PB_CD34, #1
  CB_BFUEvsPB_BFU = CB_BFUE - PB_BFU, #2
  CB_CD34vsCB_ORTHO = CB_CD34 - CB_ORTHO, #3
  #CB_CFUEvsPB_CFU = CB_CFUE - PB_CFU, #3
  #CB_PROvsPB_PRO = CB_CFUE - PB_CFU, #4
  #CB_EBASOvsPB_EBASO = CB_EBASO - PB_EBASO, #5
  #CB_LBvsPB_LB = CB_LB - PB_LB, #6
  #CB_POLYvsPB_POLY = CB_POLY - PB_POLY, #7
  #CB_ORTHOvsPB_ORTHO = CB_ORTHO - PB_ORTHO, #8
  levels = colnames(design))
##############################
### STATISTICS
##############################
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
plotSA(efit, main="Final model: Mean Variance Trend")
summary(decideTests(efit))
tfit <- treat(vfit, lfc=1)
dt <- decideTests(tfit)
summary(dt)
##############################
### write tfit to file
##############################
write.fit(tfit, dt, file="results.csv", sep = ";")
##############################
### venn diagram
##############################
vennDiagram(dt[,1:3], circle.col=c("turquoise", "red","green"))
#de.common<-which(dt[,1]!=0 & dt[,2]!=1 & dt[,3]!=1 & dt[,2]!=-1 & dt[,3]!=-1)
#length(de.common)
#pickupGenes<-rownames(tfit)[de.common]
##############################
### topgenes
##############################
CB_CD34.vs.CB_ORTHO <- topTreat(tfit, coef=3, n=Inf, adjust.method = "fdr", p.value = 0.05, lfc = 1)
##############################
### MA plot
##############################
status<-rownames(v) %in% rownames(CB_CD34.vs.CB_ORTHO)
attr(status,"values") <- c("TRUE", "FALSE")
attr(status,"col") <- c("red", "black")
plotMA(tfit, coef = 3, xlab = "Average log-expression", ylab = "log-fold-change", main = "CB_CD34.vs.CB_ORTHO", status = status, cex = 0.5)
##############################
### #pheatmap
##############################
pheatmap(v$E[status,], color = colorRampPalette(c("navy", "white","firebrick4"))(255), 
         cluster_cols = F, cluster_rows=T, show_colnames = TRUE, 
         show_rownames = FALSE,clustering_distance_rows ="euclidean", scale="row")
##############################
### #gplots heatmap.2
##############################
heatmap.2(as.matrix(v$E[i,]),col=colorRampPalette(c("navy", "white","firebrick4"))(255), 
          trace="none",main="CB_CD34 vs CB_ORTHO Top Genes",
          scale="row", density.info="none",
          margin=c(14,10),colsep=c(1:48),rowsep=c(1:48),  
          sepwidth=c(0.01,0.01),sepcolor="white",
          #Rowv=FALSE,#Colv=FALSE,   
          dendrogram="row")
###################################################
### volcano plot
###################################################
all_CB_CD34.vs.CB_ORTHO <- topTreat(tfit, coef=3, n=Inf, adjust.method = "fdr")
volcanoData <- cbind(all_CB_CD34.vs.CB_ORTHO$logFC, -log10(all_CB_CD34.vs.CB_ORTHO$adj.P.Val))
colnames(volcanoData) <- c("logFC", "negLogPval")
DEGs <- all_CB_CD34.vs.CB_ORTHO$adj.P.Val < 0.05 & abs(all_CB_CD34.vs.CB_ORTHO$logFC) > 1
point.col <- ifelse(DEGs, "red", "black")
plot(volcanoData, pch = 16, col = point.col, cex = 0.5)
###################################################
### COMPARAR DOIS GRUPOS COM ROTS
###################################################
groups = as.numeric(group[c(1:3,4:6)])
### groups 2 = CB_CD34, 10 = CD_CBFUE
input = v$E[,c(1:3,25:28)]
###################################################
### statistics
###################################################
results.rots = ROTS(data = input, groups = groups , B = 100 , K = NULL , seed = 1234,progress = TRUE)
names(results.rots) 
###################################################
### summary
###################################################
summary(results.rots, fdr = 0.05)
###################################################
### volcano plot
###################################################
plot(results.rots, fdr = 0.05, type = "volcano")
###################################################
### heatmap
###################################################
plot(results.rots, fdr = 0.05, type = "heatmap")
##############################
### #save session as:
##############################
save.image("analysis.RNAseq.2018.Rdata")

