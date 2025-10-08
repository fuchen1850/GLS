######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("pheatmap")


#引用包
library(limma)
library(pheatmap)

logFCfilter=0.585    #logFC过滤条件(logFC=0.585,差异倍数1.5倍;logFC=1,差异2倍;logFC=2,差异4倍)
fdrFilter=0.05        #fdr过滤条件
expFile="TCGA.normalize.txt"    #表达数据文件
setwd("D:\\1.毕业数据整理\\8.甲基化\\分析\\2.差异")      #设置工作目录

#读取输入文件
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]

#正常和肿瘤数目
group=sapply(strsplit(colnames(data),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2", "1", group)
conNum=length(group[group==1])       #正常组样品数目
treatNum=length(group[group==0])     #肿瘤组样品数目
grade=c(rep(1,conNum), rep(2,treatNum))

#对基因进行循环，查找差异基因
outTab=data.frame()
for(i in row.names(data)){
	geneName=unlist(strsplit(i,"\\|",))[1]
	geneName=gsub("\\/", "_", geneName)
	rt=rbind(expression=data[i,], grade=grade)
	rt=as.matrix(t(rt))
	wilcoxTest=wilcox.test(expression ~ grade, data=rt)
	pvalue=wilcoxTest$p.value
	conGeneMeans=mean(data[i,1:conNum])
	treatGeneMeans=mean(data[i,(conNum+1):ncol(data)])
	logFC=treatGeneMeans-conGeneMeans
	conMed=median(data[i,1:conNum])
	treatMed=median(data[i,(conNum+1):ncol(data)])
	diffMed=treatMed-conMed
	if( ((logFC>0) & (diffMed>0)) | ((logFC<0) & (diffMed<0)) ){  
		  outTab=rbind(outTab,cbind(gene=i,conMean=conGeneMeans,treatMean=treatGeneMeans,logFC=logFC,pValue=pvalue))
	 }
}
pValue=outTab[,"pValue"]
fdr=p.adjust(as.numeric(as.vector(pValue)), method="fdr")
outTab=cbind(outTab, fdr=fdr)

#输出所有基因的差异情况
write.table(outTab,file="all.txt",sep="\t",row.names=F,quote=F)

#输出满足条件的差异基因的表格
outDiff=outTab[( abs(as.numeric(as.vector(outTab$logFC)))>logFCfilter & as.numeric(as.vector(outTab$fdr))<fdrFilter),]
write.table(outDiff,file="diff.txt",sep="\t",row.names=F,quote=F)

#输出差异基因的表达量
heatmap=rbind(ID=colnames(data[as.vector(outDiff[,1]),]),data[as.vector(outDiff[,1]),])
write.table(heatmap,file="diffGeneExp.txt",sep="\t",col.names=F,quote=F)

#绘制差异基因的热图
geneNum=50     #设置展示基因的数目
outDiff=outDiff[order(as.numeric(as.vector(outDiff$logFC))),]
diffGeneName=as.vector(outDiff[,1])
diffLength=length(diffGeneName)
hmGene=c()
if(diffLength>(2*geneNum)){
    hmGene=diffGeneName[c(1:geneNum,(diffLength-geneNum+1):diffLength)]
}else{
    hmGene=diffGeneName
}
hmExp=data[hmGene,]
Type=c(rep("Normal",conNum),rep("Tumor",treatNum))
names(Type)=colnames(data)
Type=as.data.frame(Type)
pdf(file="heatmap.pdf", width=10, height=7)
pheatmap(hmExp, 
         annotation=Type, 
         color = colorRampPalette(c(rep("blue",5), "white", rep("red",5)))(50),
         cluster_cols =F,
         show_colnames = F,
         scale="row",
         fontsize = 8,
         fontsize_row=5,
         fontsize_col=8)
dev.off()

#绘制差异基因的火山图
pdf(file="vol.pdf", width=5, height=5)
xMax=5
yMax=max(-log10(outTab$fdr))+1
plot(as.numeric(as.vector(outTab$logFC)), -log10(outTab$fdr), xlab="logFC",ylab="-log10(fdr)",
     main="Volcano", ylim=c(0,yMax),xlim=c(-xMax,xMax),yaxs="i",pch=20, cex=1.2)
diffSub=subset(outTab, fdr<fdrFilter & as.numeric(as.vector(logFC))>logFCfilter)
points(as.numeric(as.vector(diffSub$logFC)), -log10(diffSub$fdr), pch=20, col="red",cex=1.5)
diffSub=subset(outTab, fdr<fdrFilter & as.numeric(as.vector(logFC))<(-logFCfilter))
points(as.numeric(as.vector(diffSub$logFC)), -log10(diffSub$fdr), pch=20, col="green",cex=1.5)
abline(v=0,lty=2,lwd=3)
dev.off()


######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

