######Video source: http://ke.biowolf.cn
######生信自学网: http://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：2749657388@qq.com
######答疑微信: 18520221056

#install.packages("survival")

setwd("C:\\Users\\周小珠\\Desktop\\结果较好\\5.生存")   #工作目录（需修改）

library(survival)

#绘制train组生存曲线
rt=read.table("riskTrain.txt",header=T,sep="\t")
diff=survdiff(Surv(futime, fustat) ~risk,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,4)
pValue=format(pValue, scientific = TRUE)
fit <- survfit(Surv(futime, fustat) ~ risk, data = rt)
summary(fit)    #查看五年生存率
pdf(file="survivalTrain.pdf",width=5.5,height=5)
plot(fit, 
     lwd=2,
     col=c("red","blue"),
     xlab="Time (year)",
     ylab="Survival rate",
     main=paste("Survival train (p=", pValue ,")",sep=""),
     mark.time=T)
legend("topright", 
       c("high risk", "low risk"),
       lwd=2,
       col=c("red","blue"))
dev.off()

#绘制test组生存曲线
rt=read.table("riskTest.txt",header=T,sep="\t")
diff=survdiff(Surv(futime, fustat) ~risk,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,4)
pValue=format(pValue, scientific = TRUE)
fit <- survfit(Surv(futime, fustat) ~ risk, data = rt)
summary(fit)    #查看五年生存率
pdf(file="survivalTest.pdf",width=5.5,height=5)
plot(fit, 
     lwd=2,
     col=c("red","blue"),
     xlab="Time (year)",
     ylab="Survival rate",
     main=paste("Survival test (p=", pValue ,")",sep=""),
     mark.time=T)
legend("topright", 
       c("high risk", "low risk"),
       lwd=2,
       col=c("red","blue"))
dev.off()

######Video source: http://ke.biowolf.cn
######生信自学网: http://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：2749657388@qq.com
######答疑微信: 18520221056