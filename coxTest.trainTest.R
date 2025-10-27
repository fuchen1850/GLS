######Video source: http://ke.biowolf.cn
######生信自学网: http://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：2749657388@qq.com
######答疑微信: 18520221056

#install.packages("caret")

setwd("C:\\Users\\周小珠\\Desktop\\1.分组")      #设置工作目录
rt=read.table("expTime.txt",sep="\t",header=T,check.names=F)

#划分训练集和测试集
library(caret)
inTrain<-createDataPartition(y=rt[,3],p=0.5,list=F)
train<-rt[inTrain,]
test<-rt[-inTrain,]
write.table(train,file="train.txt",sep="\t",quote=F,row.names=F)
write.table(test,file="test.txt",sep="\t",quote=F,row.names=F)


######Video source: http://ke.biowolf.cn
######生信自学网: http://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：2749657388@qq.com
######答疑微信: 18520221056