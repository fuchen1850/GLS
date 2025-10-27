######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("ggplot2")
#install.packages("ggpubr")


#引用包
library(limma)
library(ggplot2)
library(ggpubr)

pFilter=0.001                      #pvalue过滤条件
riskFile="risk.TCGA.txt"           #风险文件
drugFile="DrugPredictions.csv"     #药物敏感性文件
setwd("C:\\Users\\周小珠\\Desktop\\结果较好\\28.boxplot")     #设置工作目录

#读取风险输入文件
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)

#读取药物敏感性文件
senstivity=read.csv(drugFile, header=T, sep=",", check.names=F, row.names=1)
colnames(senstivity)=gsub("(.*)\\_(\\d+)", "\\1", colnames(senstivity))

#数据合并
sameSample=intersect(row.names(risk), row.names(senstivity))
risk=risk[sameSample, "Risk",drop=F]
senstivity=senstivity[sameSample,,drop=F]
rt=cbind(risk, senstivity)

#设置比较组
rt$Risk=factor(rt$Risk, levels=c("low", "high"))
type=levels(factor(rt[,"Risk"]))
comp=combn(type, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

#对药物进行循环, 绘制箱线图
for(drug in colnames(rt)[2:ncol(rt)]){
	rt1=rt[,c(drug, "Risk")]
	colnames(rt1)=c("Drug", "Risk")
	rt1=na.omit(rt1)
	rt1$Drug=log2(rt1$Drug+1)
	#差异分析
	test=wilcox.test(Drug ~ Risk, data=rt1)
	diffPvalue=test$p.value
	#对满足条件的药物绘制箱线图
	if(diffPvalue<pFilter){
		boxplot=ggboxplot(rt1, x="Risk", y="Drug", fill="Risk",
					      xlab="Risk",
					      ylab=paste0(drug, " senstivity"),
					      legend.title="Risk",
					      palette=c("green", "red")
					     )+ 
			stat_compare_means(comparisons=my_comparisons)
		#输出图形
		pdf(file=paste0("drugSenstivity.", drug, ".pdf"), width=5, height=4.5)
		print(boxplot)
		dev.off()
	}
}


######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

