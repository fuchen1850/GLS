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
#install.packages("ggExtra")


#引用包
library(limma)
library(ggplot2)
library(ggpubr)
library(ggExtra)
riskFile="risk.TCGA.txt"       #风险文件
RNAssFile="StemnessScores_DNAmeth_20170210.tsv"           #干细胞打分文件
setwd("C:\\Users\\周小珠\\Desktop\\结果较好\\17.风险评分与Stemness(干细胞）相关性分析")     #设置工作目录

#读取风险文件
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)

#读取RNA干细胞打分文件
RNAss=read.table(RNAssFile, header=T, sep="\t",check.names=F, row.names=1)
RNAss=t(RNAss[1,,drop=F])
rownames(RNAss)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(RNAss))
RNAss=avereps(RNAss)

#样品取交集
sameSample=intersect(row.names(risk), row.names(RNAss))
risk=risk[sameSample,"riskScore",drop=F]
RNAss=RNAss[sameSample,,drop=F]
data=cbind(RNAss, risk)

#绘制相关性散点图
xlab="riskScore"
ylab="DNAss"
outFile="DNAss.cor.pdf"
x=as.numeric(data[,xlab])
x[x>quantile(x,0.99)]=quantile(x,0.99)
y=as.numeric(data[,ylab])
df1=as.data.frame(cbind(x,y))
p1=ggplot(df1, aes(x, y)) + 
		  xlab("Risk score") + ylab(ylab)+ ylim(0,0.7)+
		  geom_point() + geom_smooth(method="lm",formula = y ~ x) + theme_bw()+
		  stat_cor(method = 'spearman', aes(x =x, y =y))
p2=ggMarginal(p1, type="density", xparams=list(fill = "orange"), yparams=list(fill = "blue"))
#相关性图形
pdf(file=outFile, width=5.2, height=5)
print(p2)
dev.off()


######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

