#setwd("/Users/stalenygard/ruvseq")
#install.packages("pROC")

print("ok")

library(EDASeq)
library(pROC)

files<-c("GSE47774_SEQC_ILM_AGR.txt",
"GSE47774_SEQC_ILM_BGI.txt",
"GSE47774_SEQC_ILM_CNL.txt",
"GSE47774_SEQC_ILM_COH.txt",
"GSE47774_SEQC_ILM_MAY.txt",
"GSE47774_SEQC_ILM_NVS.txt",
"GSE47774_SEQC_ROC_SQW.txt",
"GSE47774_SEQC_ROC_NYU.txt",
"GSE47774_SEQC_ROC_MGP.txt",
"GSE47774_SEQC_LIF_SQW.txt",
"GSE47774_SEQC_LIF_PSU.txt",
"GSE47774_SEQC_LIF_NWU.txt",
"GSE47774_SEQC_LIF_LIV.txt")

files<-files[-(7:9)]

D<-read.table(files[1],sep="\t",header=TRUE)
#for (i in 2:length(files)){
#    print(i)
#    data<-read.table(files[i],sep="\t",header=TRUE)
#    D<-cbind(D,data[,-1])
#}

#hkg<-read.table("hkg-short-2.txt",sep="\t",header=TRUE)
#hkg.navn<-hkg[,2]
navn<-as.character(D[,1])
data<-read.table("cms_095046.txt",sep="\t",header=TRUE)
mid<-as.character(data[which(data[,3]=="B"),2])
ix.neg<-match(gsub("-","_",mid),navn)
v<-sample(length(ix.neg))
ix.neg.train<-ix.neg[v[1:(length(v)/2)]]
ix.neg.test<-ix.neg[v[(length(v)/2+1):length(v)]]
mid<-as.character(data[which(data[,3]!="B"),2])
ix.pos<-match(gsub("-","_",mid),navn)
ix.ctrl<-c(ix.neg.test,ix.pos)
#install.packages("pROC")
library(pROC)
truth<-c(rep(0,length(ix.neg.test)),rep(1,length(ix.pos)))


print("ok")

ll<-length(files)
auc.ruvseq<-rep(0,ll)
auc.nonruvseq<-rep(0,ll)
auc.ruvseq.emp<-rep(0,ll)
#library(BiocParallel)
#register(MulticoreParam(7))
for (i in (1:ll)){
    D<-read.table(files[i],sep="\t",header=TRUE)
    dim(D)
    print("ok1")
    ix.1<-grep("_A_",names(D))
    ix.2<-grep("_B_",names(D))
    #source("https://bioconductor.org/biocLite.R")
    #biocLite("DESeq2")
    library(DESeq2)
    ix<-c(ix.1,ix.2)
    X<-D[,ix]
    condition=factor(c(rep("A",length(ix.1)),rep("B",length(ix.2))))
    rm(D)
    coldata<-data.frame(condition)
    row.names(coldata)<-colnames(X)
    dds <- DESeqDataSetFromMatrix(countData = as.matrix(X),
                              colData = as.matrix(coldata),
                              design= ~ condition)
    #dds <- dds[ rowSums(counts(dds)) > 20, ]
    dds<-DESeq(dds,parallel=FALSE)
    res <- results(dds,parallel=FALSE)
    pred1<-abs(res$log2FoldChange[ix.ctrl])
    roc.res.1<-roc(response=truth,predictor=pred1)
    auc.nonruvseq[i]<-auc(roc.res.1)
    rm(dds)
    V<-data.frame(condition)
    row.names(V)<-colnames(X)
    colnames(V)<-"condition"
    #rm(X)
    library(edgeR)
    filter <- apply(X,1, function(x) length(x[x>5])>=2)
    filtered <- X[filter,]

    set <- newSeqExpressionSet(as.matrix(filtered),phenoData =coldata)

    #establishing empirical
    design <- model.matrix(~condition, data=pData(set))
    y <- DGEList(counts=counts(set), group=coldata$condition)
    y <- calcNormFactors(y, method="upperquartile")
    y <- estimateGLMCommonDisp(y, design)
    fit <- glmFit(y, design)
    lrt <- glmLRT(fit, coef=2)
    top <- topTags(lrt, n=nrow(set))$table
    empirical <- as.numeric(rownames(set)[which(!(rownames(set) %in% rownames(top)[which(top$P<0.8)]))])
    L<-5
    auc.vec<-auc.vec.emp<-rep(0,L)
    library(RUVSeq)
    print(auc.nonruvseq)
    for (l in 1:L){
	seqRUVg <- RUVg(as.matrix(X),ix.neg.train, k=l)
	#V<-data.frame(condition)
	#row.names(V)<-colnames(X)
	col2<-cbind(V,seqRUVg$W)
	#col2<-as.matrix(col2)
	dds2 <- DESeqDataSetFromMatrix(countData = as.matrix(X),colData = col2,design=~condition)
	design(dds2)=as.formula(paste(c("~",paste(colnames(seqRUVg$W),rep("+" ,l)),"condition")))
	#dds2 <- dds2[ rowSums(counts(dds2)) > 20, ]
	dds2 <- DESeq(dds2,parallel=FALSE)
	res2 <- results(dds2,parallel=FALSE)
	pred2<-abs(res2$log2FoldChange[ix.ctrl])
	roc.res.2<-roc(response=truth,predictor=pred2)
        auc.vec[l]<-auc(roc.res.2)
	#print(auc.vec)
	seqRUVg <- RUVg(as.matrix(X),empirical,k=l)
        #V<-data.frame(condition)
        #row.names(V)<-colnames(X)
        col2<-cbind(V,seqRUVg$W)
        #col2<-as.matrix(col2)
        dds2 <- DESeqDataSetFromMatrix(countData = as.matrix(X),colData = col2,design=~condition)
        design(dds2)=as.formula(paste(c("~",paste(colnames(seqRUVg$W),rep("+" ,l)),"condition")))
        #dds2 <- dds2[ rowSums(counts(dds2)) > 20, ]
        dds2 <- DESeq(dds2,parallel=FALSE)
        res2 <- results(dds2,parallel=FALSE)
        pred2<-abs(res2$log2FoldChange[ix.ctrl])
        roc.res.2<-roc(response=truth,predictor=pred2)
        auc.vec.emp[l]<-auc(roc.res.2)
        print(auc.vec)
	print(auc.vec.emp)
    }
    #print(auc.vec)
    #print(auc.vec.emp)
    auc.ruvseq[i]<-max(auc.vec)
    auc.ruvseq.emp[i]<-max(auc.vec.emp)
    print(auc.nonruvseq)
    print(auc.ruvseq)
    print(auc.ruvseq.emp)
}
print(auc.nonruvseq)
print(auc.ruvseq)
print(auc.ruvseq.emp)