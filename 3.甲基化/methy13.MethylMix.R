######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

#if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("MethylMix")


library(MethylMix)    #引用包
pFilter=0.05          #pvalue过滤条件
logFCfilter=0.585         #logFC过滤条件(logFC=0.585,差异倍数1.5倍;logFC=1,差异2倍;logFC=2,差异4倍)
corFilter=-0.2        #相关系数过滤条件
GEfile="mRNA.cancer.txt"      #肿瘤样品的表达数据文件
METfile="MET.cancer.txt"      #肿瘤样品的甲基化数据文件
METNfile="MET.normal.txt"     #正常样品的甲基化数据文件
setwd("D:\\1.毕业数据整理\\8.甲基化\\分析\\4.甲基化")      #设置工作目录

#读取输入文件
GEcancer = read.table(GEfile, header=T, sep="\t", check.names=F, row.names=1)
METcancer = read.table(METfile, header=T, sep="\t", check.names=F, row.names=1)
METnormal = read.table(METNfile, header=T, sep="\t", check.names=F, row.names=1)
GEcancer=as.matrix(GEcancer)
METcancer=as.matrix(METcancer)
METnormal=as.matrix(METnormal)

#查找驱动基因
MethylMixResults=MethylMix(METcancer, GEcancer, METnormal)

#对基因进行过滤
outTab=data.frame()
for (gene in MethylMixResults$MethylationDrivers) {
    wilcoxTest=wilcox.test(METnormal[gene,], METcancer[gene,])
    wilcoxP=wilcoxTest$p.value
    normalGeneMeans=mean(METnormal[gene,])
    tumorGeneMeans=mean(METcancer[gene,])
    logFC=log2(tumorGeneMeans)-log2(normalGeneMeans)
    betaDiff=tumorGeneMeans-normalGeneMeans
    normalMed=median(METnormal[gene,])
    tumorMed=median(METcancer[gene,])
    diffMed=tumorMed-normalMed
    
    x=as.numeric(METcancer[gene,])
    y=as.numeric(GEcancer[gene,])
    corT=cor.test(x,y)
	cor=corT$estimate
	cor=round(cor,3)
	pvalue=corT$p.value
	if(cor<corFilter){
		if( ((logFC>0) & (diffMed>0)) | ((logFC<0) & (diffMed<0)) ){
		    if((wilcoxP<pFilter) & (abs(logFC)>logFCfilter)){
				if(pvalue<0.001){
					pval=signif(pvalue,4)
					pval=format(pval, scientific = TRUE)
				}else{pval=round(pvalue,3)}
				#输出甲基化差异的图形
				pdf(file=paste0("diff.",gene,".pdf"), width=7, height=5)
				plots=MethylMix_PlotModel(gene,MethylMixResults,METcancer,GEcancer,METnormal)
				print(plots$MixtureModelPlot)
				dev.off()
				#输出相关性图形
				pdf(file=paste0("cor.",gene,".pdf"), width=5, height=5)
				plot(x,y, type="p",pch=16,main=paste("Cor= ",cor," (p=",pval,")",sep=""),
								cex=0.9, cex.lab=1.05, cex.main=1.1, cex.axis=1,     #依次是点、坐标轴、标题、刻度字体大小
								xlab=paste0(gene," methylation"),ylab=paste0(gene," expression") )
				z=lm(y~x)
				lines(x,fitted(z),col=2)
				dev.off()
				#保存信息
				outTab=rbind(outTab,cbind(gene=gene,
				                    normalMean=normalGeneMeans,
				                    TumorMean=tumorGeneMeans,
				                    logFC=logFC,
				                    pValue=wilcoxP,
				                    cor=corT$estimate,
				                    corPavlue=corT$p.value))
	        }
	    }
	}
}
#输出驱动基因的结果
outTab=outTab[order(as.numeric(as.vector(outTab$pValue))),]
write.table(outTab,file="drivenGene.xls",sep="\t",row.names=F,quote=F)
write.table(outTab[,1],file="drivenGene.txt",sep="\t",row.names=F,quote=F,col.names=F)

#输出驱动基因的甲基化水平
drivenGeneMethy=METcancer[as.vector(outTab[,1]),]
drivenGeneMethy=rbind(ID=colnames(drivenGeneMethy),drivenGeneMethy)
write.table(drivenGeneMethy,file="drivenGeneMethy.txt",sep="\t",col.names=F,quote=F)

#输出驱动基因的表达水平
drivenGeneExp=GEcancer[as.vector(outTab[,1]),]
drivenGeneExp=rbind(ID=colnames(drivenGeneExp),drivenGeneExp)
write.table(drivenGeneExp,file="drivenGeneExp.txt",sep="\t",col.names=F,quote=F)


######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

