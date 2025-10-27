######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")
#BiocManager::install("GSEABase")
#BiocManager::install("GSVA")

#install.packages("reshape2")
#install.packages("ggplot2")


#引用包
library(limma)
library(GSEABase)
library(GSVA)
library(reshape2)
library(ggplot2)

expFile="symbol.txt"      #表达数据文件
riskFile="risk.TCGA.txt"    #风险文件
gmtFile="KEGG.gmt"          #基因集文件
setwd("C:\\Users\\周小珠\\Desktop\\结果较好\\19.风险评分&模型基因与GSVA相关分析")      #设置工作目录

#读取表达输入文件,并对输入文件整理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)

#GSVA分析
geneSets=getGmt(gmtFile, geneIdType=SymbolIdentifier())
gsvaResult=gsva(data, 
                geneSets, 
                min.sz=10, 
                max.sz=500, 
                verbose=TRUE,
                parallel.sz=1)
#gsvaOut=rbind(id=colnames(gsvaResult), gsvaResult)
#write.table(gsvaOut, file="gsvaOut.txt", sep="\t", quote=F, col.names=F)

#删除正常样品
data=t(gsvaResult)
group=sapply(strsplit(row.names(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
data=data[group==0,]
row.names(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", row.names(data))
data=avereps(data)

#读取风险文件
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
sameSample=intersect(row.names(data), row.names(risk))
data=data[sameSample,,drop=F]
risk=risk[sameSample,3:(ncol(risk)-1),drop=F]

#相关性分析
outTab=data.frame()
for(Geneset in colnames(data)){
	for(gene in colnames(risk)){
		x=as.numeric(data[,Geneset])
		y=as.numeric(risk[,gene])
		corT=cor.test(x,y,method="spearman")
		cor=corT$estimate
		pvalue=corT$p.value
		text=ifelse(pvalue<0.001,"***",ifelse(pvalue<0.01,"**",ifelse(pvalue<0.05,"*","")))
		outTab=rbind(outTab,cbind(Gene=gene, Geneset=Geneset, cor, text, pvalue))
	}
}

#绘制相关性热图
outTab$Gene=factor(outTab$Gene, levels=colnames(risk))
outTab$cor=as.numeric(outTab$cor)
pdf(file="GSVAcor.pdf", width=10, height=7)
ggplot(outTab, aes(Gene, Geneset)) + 
	geom_tile(aes(fill = cor), colour = "grey", size = 1)+
	scale_fill_gradient2(low = "green2", mid = "white", high = "red2") + 
	geom_text(aes(label=text),col ="black",size = 3) +
	theme_minimal() +    #去掉背景
	theme(axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y=element_blank(),
	      axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),   #x轴字体
	      axis.text.y = element_text(size = 10, face = "bold")) +       #y轴字体
	labs(fill =paste0("***  p<0.001","\n", "**  p<0.01","\n", " *  p<0.05","\n", "\n","Correlation")) +   #设置图例
	scale_x_discrete(position = "bottom")      #X轴名称显示位置
dev.off()


######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

