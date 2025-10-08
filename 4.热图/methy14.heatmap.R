######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

#install.packages("pheatmap")


library(pheatmap)      #引用包
setwd("D:\\1.毕业数据整理\\8.甲基化\\分析\\5.热图")      #设置工作目录
gene=read.table("drivenGene.txt", header=F, sep="\t", check.names=F)     #读取驱动基因文件

#绘制基因表达的热图
normal="mRNA.normal.txt"      #正常样品的表达文件
tumor="mRNA.cancer.txt"       #肿瘤样品的表达文件
rt1=read.table(normal, header=T, sep="\t", check.names=F, row.names=1)
rt2=read.table(tumor, header=T, sep="\t", check.names=F, row.names=1)
rt=cbind(rt1,rt2)
rt=rt[as.vector(gene[,1]),]
Type=c( rep("normal",ncol(rt1)),rep("tumor",ncol(rt2)) )
names(Type)=colnames(rt)
Type=as.data.frame(Type)
pdf(file="mRNA.heatmap.pdf", width=7, height=6)
pheatmap(rt,
         annotation=Type,
         cluster_cols=F,
         show_colnames=F,
         scale="row",
         color=colorRampPalette(c(rep("blue",2), "white", rep("red",2)))(50),
		 fontsize=6,
		 fontsize_row=4,
		 fontsize_col=6)
dev.off()

#绘制甲基化数据的热图
normal="MET.normal.txt"       #正常样品的甲基化文件
tumor="MET.cancer.txt"        #肿瘤样品的甲基化文件
rt1=read.table(normal, header=T, sep="\t", check.names=F, row.names=1)
rt2=read.table(tumor, header=T, sep="\t", check.names=F, row.names=1)
rt=cbind(rt1,rt2)
rt=rt[as.vector(gene[,1]),]
Type=c( rep("normal",ncol(rt1)),rep("tumor",ncol(rt2)) )
names(Type)=colnames(rt)
Type=as.data.frame(Type)
pdf(file="MET.heatmap.pdf", width=7, height=6)
pheatmap(rt,
         annotation=Type,
         cluster_cols=F,
         show_colnames=F,
         scale="row",
         color=colorRampPalette(c(rep("blue",2), "white", rep("red",2)))(50),
         fontsize=6,
		 fontsize_row=4,
		 fontsize_col=6)
dev.off()


######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

