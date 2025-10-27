######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("survival")
#install.packages("survminer")


#引用包
library(limma)
library(survival)
library(survminer)

riskFile="risk.TCGA.txt"        #风险文件
methyFile="MET.cancer.txt"      #甲基化数据文件
setwd("C:\\Users\\周小珠\\Desktop\\结果较好\\24.methySur")     #设置工作目录

#读取风险文件
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
gene=colnames(risk)[3:(ncol(risk)-2)]

#读取甲基化数据文件
methy=read.table(methyFile, header=T, sep="\t", check.names=F, row.names=1)
methy=t(methy[gene,,drop=F])
row.names(methy)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", row.names(methy))
methy=avereps(methy)

#数据合并
sameSample=intersect(row.names(risk), row.names(methy))
risk=risk[sameSample,c("futime","fustat"),drop=F]
methy=methy[sameSample,,drop=F]
rt=cbind(risk, methy)

#对模型基因进行循环
outTab=data.frame()
for(gene in colnames(rt)[3:ncol(rt)]){
	if(sd(rt[,gene])<0.001){next}
	data=rt[,c("futime", "fustat", gene)]
	colnames(data)=c("futime", "fustat", "gene")
	
	#获取最优cutoff
	res.cut=surv_cutpoint(data, time = "futime", event = "fustat", variables =c("gene"))
	res.cat=surv_categorize(res.cut)
	fit=survfit(Surv(futime, fustat) ~gene, data = res.cat)
	
	#比较高低表达组的生存差异
	diff=survdiff(Surv(futime, fustat) ~gene,data =res.cat)
	pValue=1-pchisq(diff$chisq, df=1)
	outVector=cbind(gene, res.cut$cutpoint[1], pValue)
	outTab=rbind(outTab,outVector)
	if(pValue<0.001){
		pValue="p<0.001"
	}else{
		pValue=paste0("p=",sprintf("%.03f",pValue))
	}
		
	#绘制生存曲线
	surPlot=ggsurvplot(fit,
					data=res.cat,
					pval=pValue,
					pval.size=6,
					legend.title=paste0(gene, " methylation"),
					legend.labs=c("high","low"),
					xlab="Time(years)",
					break.time.by=1,
					palette=c("red", "blue"),
					conf.int=F,
					risk.table=F,
					risk.table.title="",
					risk.table.height=.25)
	
	#输出图形
	pdf(file=paste0("Survival.",gene,".pdf"), width=5, height=4.5, onefile=FALSE)
	print(surPlot)
	dev.off()
}

#输出模型基因名称和pvalue的表格
write.table(outTab,file="survival.result.txt",sep="\t",row.names=F,quote=F)


######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

