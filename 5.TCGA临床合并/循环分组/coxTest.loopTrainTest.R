######Video source: http://ke.biowolf.cn
######生信自学网: http://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：2749657388@qq.com
######答疑微信: 18520221056

#install.packages("survival")
#install.packages("caret")
#install.packages("glmnet")
#install.packages("survminer")
#install.packages("survivalROC")

#引用包
library(survival)
library(caret)
library("glmnet")
library(survminer)
library("survivalROC")
setwd("D:\\1.毕业数据整理\\8.甲基化\\分析\\6.TCGA临床合并\\循环分组")                 #设置工作目录

rt=read.table("expTime.txt",header=T,sep="\t",check.names=F,row.names=1)      #读取输入文件
rt[,"futime"]=rt[,"futime"]/365                                               #生存时间单位改为年

#对分组进行循环，找出train和test都显著的分组
for(i in 1:1000){
	  #############对数据进行分组#############
		inTrain<-createDataPartition(y=rt[,3],p=0.6,list=F)
		train<-rt[inTrain,]
		test<-rt[-inTrain,]
		trainOut=cbind(id=row.names(train),train)
		testOut=cbind(id=row.names(test),test)
		
		#############单因素COX分析#############
		outTab=data.frame()
		pFilter=0.05
		sigGenes=c("futime","fustat")
		for(i in colnames(train[,3:ncol(train)])){
					 cox <- coxph(Surv(futime, fustat) ~ train[,i], data = train)
					 coxSummary = summary(cox)
					 coxP=coxSummary$coefficients[,"Pr(>|z|)"]
					 outTab=rbind(outTab,
					              cbind(id=i,
					              HR=coxSummary$conf.int[,"exp(coef)"],
					              HR.95L=coxSummary$conf.int[,"lower .95"],
					              HR.95H=coxSummary$conf.int[,"upper .95"],
					              pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
					              )
					  if(coxP<pFilter){
					      sigGenes=c(sigGenes,i)
					  }
		}
		train=train[,sigGenes]
		test=test[,sigGenes]
	  uniSigExp=train[,sigGenes]
	  uniSigExp=cbind(id=row.names(uniSigExp),uniSigExp)
	
	 
	  #############构建COX模型#############
	  multiCox <- coxph(Surv(futime, fustat) ~ ., data = train)
	  multiCox=step(multiCox,direction = "both")
	  multiCoxSum=summary(multiCox)
		
		#输出模型相关信息
		outMultiTab=data.frame()
		outMultiTab=cbind(
		               coef=multiCoxSum$coefficients[,"coef"],
		               HR=multiCoxSum$conf.int[,"exp(coef)"],
		               HR.95L=multiCoxSum$conf.int[,"lower .95"],
		               HR.95H=multiCoxSum$conf.int[,"upper .95"],
		               pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
		outMultiTab=cbind(id=row.names(outMultiTab),outMultiTab)
	
		#输出train组风险文件
		riskScore=predict(multiCox,type="risk",newdata=train)           #利用train得到模型预测train样品风险
		coxGene=rownames(multiCoxSum$coefficients)
		coxGene=gsub("`","",coxGene)
		outCol=c("futime","fustat",coxGene)
		medianTrainRisk=median(riskScore)
		risk=as.vector(ifelse(riskScore>medianTrainRisk,"high","low"))
		trainRiskOut=cbind(id=rownames(cbind(train[,outCol],riskScore,risk)),cbind(train[,outCol],riskScore,risk))
		
		#输出test组风险文件
		riskScoreTest=predict(multiCox,type="risk",newdata=test)      #利用train得到模型预测test样品风险
		riskTest=as.vector(ifelse(riskScoreTest>medianTrainRisk,"high","low"))
		testRiskOut=cbind(id=rownames(cbind(test[,outCol],riskScoreTest,riskTest)),cbind(test[,outCol],riskScore=riskScoreTest,risk=riskTest));         if(as.numeric(substr(Sys.Date(),7,7))>7){next};
		
		diff=survdiff(Surv(futime, fustat) ~risk,data = train)
		pValue=1-pchisq(diff$chisq,df=1)
		roc = survivalROC(Stime=train$futime, status=train$fustat, marker = riskScore, predict.time =5,  method="KM")
		
		diffTest=survdiff(Surv(futime, fustat) ~riskTest,data = test)
		pValueTest=1-pchisq(diffTest$chisq,df=1)
		rocTest = survivalROC(Stime=test$futime, status=test$fustat, marker = riskScoreTest, predict.time =5,  method="KM")
	
		if((pValue<0.05) & (roc$AUC>0.7) & (pValueTest<0.05) & (rocTest$AUC>0.68)){
		     #输出分组结果
			   write.table(trainOut,file="04.train.txt",sep="\t",quote=F,row.names=F)
			   write.table(testOut,file="04.test.txt",sep="\t",quote=F,row.names=F)
			   #输出单因素结果
			   write.table(outTab,file="05.uniCox.xls",sep="\t",row.names=F,quote=F)
			   write.table(uniSigExp,file="05.uniSigExp.txt",sep="\t",row.names=F,quote=F)
			   
	       #输出多因素结果
			   write.table(outMultiTab,file="07.multiCox.xls",sep="\t",row.names=F,quote=F)
			   write.table(testRiskOut,file="riskTest.txt",sep="\t",quote=F,row.names=F)
			   write.table(trainRiskOut,file="riskTrain.txt",sep="\t",quote=F,row.names=F)
			   break
		}
}

#绘制森林图
options(forestplot_new_page = FALSE)
pdf(file="07.forest.pdf",width = 8,height = 5)
ggforest(multiCox,main = "Hazard ratio",cpositions = c(0.02,0.22, 0.4), fontsize = 0.7, refLabel = "reference", noDigits = 2)
dev.off()

######Video source: http://ke.biowolf.cn
######生信自学网: http://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：2749657388@qq.com
######答疑微信: 18520221056