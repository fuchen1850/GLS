######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")


library(limma)       #引用包
expFile="diffGeneExp.txt"       #表达数据文件
methyFile="methyMatrix.txt"     #甲基化数据文件
setwd("D:\\1.毕业数据整理\\8.甲基化\\分析\\3.merge")     #设置工作目录

#读取表达数据文件
RNA=read.table(expFile, header=T, sep="\t", check.names=F)
RNA=as.matrix(RNA)
rownames(RNA)=RNA[,1]
exp=RNA[,2:ncol(RNA)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
RNA=data[rowMeans(data)>0,]

#读取甲基化数据文件
methy=read.table(methyFile, header=T, sep="\t", check.names=F, row.names=1)
methy=methy[rowMeans(methy)>0,]
methy=normalizeBetweenArrays(as.matrix(methy))

#提取正常样品和肿瘤样品表达数据
group=sapply(strsplit(colnames(RNA),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
mRNAnormal=RNA[,group==1]
RNA=RNA[,group==0]

#提取正常样品和肿瘤样品甲基化数据
group=sapply(strsplit(colnames(methy),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
Noraml=methy[,group==1]
methy=methy[,group==0]

#数据取交集
sameGene=intersect(rownames(methy),rownames(RNA))
colnames(methy)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*","\\1\\-\\2\\-\\3\\-\\4",colnames(methy))
colnames(Noraml)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*","\\1\\-\\2\\-\\3\\-\\4",colnames(Noraml))
colnames(RNA)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*","\\1\\-\\2\\-\\3\\-\\4",colnames(RNA))
colnames(mRNAnormal)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*","\\1\\-\\2\\-\\3\\-\\4",colnames(mRNAnormal))
sameSample=intersect(colnames(methy),colnames(RNA))

#输出整理的结果
mRNAnormal=rbind(id=colnames(mRNAnormal),mRNAnormal[sameGene,])
mRNAcancer=rbind(id=sameSample,RNA[sameGene,sameSample])
METnormal=rbind(id=colnames(Noraml),Noraml[sameGene,])
METcancer=rbind(id=sameSample,methy[sameGene,sameSample])
write.table(mRNAnormal,file="mRNA.normal.txt",sep="\t",quote=F,col.names=F)
write.table(mRNAcancer,file="mRNA.cancer.txt",sep="\t",quote=F,col.names=F)
write.table(METnormal,file="MET.normal.txt",sep="\t",quote=F,col.names=F)
write.table(METcancer,file="MET.cancer.txt",sep="\t",quote=F,col.names=F)


######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

