######Video source: http://ke.biowolf.cn
######生信自学网: http://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：2749657388@qq.com
######答疑微信: 18520221056

#install.packages('survival')
#install.packages("survminer")

library(survival)
library(survminer)

setwd("C:\\Users\\周小珠\\Desktop\\4.风险评分")                 #设置工作目录
rt=read.table("lassoSigExp.txt",header=T,sep="\t",check.names=F,row.names=1)  #读取train输入文件
rt[,"futime"]=rt[,"futime"]/365

#使用train组构建COX模型
multiCox=coxph(Surv(futime, fustat) ~ ., data = rt)
multiCox=step(multiCox,direction = "both")
multiCoxSum=summary(multiCox)

#输出模型相关信息
outTab=data.frame()
outTab=cbind(
             coef=multiCoxSum$coefficients[,"coef"],
             HR=multiCoxSum$conf.int[,"exp(coef)"],
             HR.95L=multiCoxSum$conf.int[,"lower .95"],
             HR.95H=multiCoxSum$conf.int[,"upper .95"],
             pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
outTab=cbind(id=row.names(outTab),outTab)
write.table(outTab,file="multiCox.xls",sep="\t",row.names=F,quote=F)

#绘制森林图
pdf(file="forest.pdf",
       width = 8,             #图片的宽度
       height = 5,            #图片的高度
       )
ggforest(multiCox,
         main = "Hazard ratio",
         cpositions = c(0.02,0.22, 0.4), 
         fontsize = 0.7, 
         refLabel = "reference", 
         noDigits = 2)
dev.off()

#输出train组风险文件
riskScore=predict(multiCox,type="risk",newdata=rt)           #利用train得到模型预测train样品风险
coxGene=rownames(multiCoxSum$coefficients)
coxGene=gsub("`","",coxGene)
outCol=c("futime","fustat",coxGene)
medianTrainRisk=median(riskScore)
risk=as.vector(ifelse(riskScore>medianTrainRisk,"high","low"))
write.table(cbind(id=rownames(cbind(rt[,outCol],riskScore,risk)),cbind(rt[,outCol],riskScore,risk)),
    file="riskTrain.txt",
    sep="\t",
    quote=F,
    row.names=F)

#输出test组风险文件
rtTest=read.table("test.txt",header=T,sep="\t",check.names=F,row.names=1)          #读取test输入文件
rtTest[,"futime"]=rtTest[,"futime"]/365
riskScoreTest=predict(multiCox,type="risk",newdata=rtTest)      #利用train得到模型预测test样品风险
riskTest=as.vector(ifelse(riskScoreTest>medianTrainRisk,"high","low"))
write.table(cbind(id=rownames(cbind(rtTest[,outCol],riskScoreTest,riskTest)),cbind(rtTest[,outCol],riskScore=riskScoreTest,risk=riskTest)),
    file="riskTest.txt",
    sep="\t",
    quote=F,
    row.names=F)


######Video source: http://ke.biowolf.cn
######生信自学网: http://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：2749657388@qq.com
######答疑微信: 18520221056