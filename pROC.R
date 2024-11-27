#307686155@qq.com
#18666983305
#更多课程请关注生信碱移
#af

#引用包
library(pROC)
library(ggplot2)

#正常对照名称
C="Normal"
#roc图片颜色
afcol="#0073C2FF"
hub="keycluster.txt" #核心基因文件名称
expfile="1.rawexp_GSE54236.txt"  

exp=read.table(expfile,sep = "\t",header = F,check.names = F)
exp=t(exp)
colnames(exp)=exp[1,]
exp=exp[2:nrow(exp),]
colnames(exp)[1]="ID"
cli=read.table("sample.txt",sep = "\t",header = F,check.names = F)         #载入临床性状数据
cli[,2]=ifelse(cli[,2]==C,0,1)
colnames(cli)=c("ID","P")
af=merge(cli,exp,by = "ID")
af[,3:ncol(af)]=lapply(af[,3:ncol(af)],as.numeric)

genes=read.table(hub,sep = " ",header = F,check.names = F)[,1]
genes=intersect(colnames(af),genes)

par(mfrow=c(4,6))  #ROC图呈4行6列摆放
for (i in c(genes)) {
  for (j in c("P")) {              
afroc=roc(af[,j],af[,i])
print(paste0(auc(afroc),"_",i,"_",j))
#auc.polygon.col参数为颜色修改
#pdf(file=paste0(i,"_",j,".pdf"), width=6, height=6)
plot(afroc, print.auc=TRUE, auc.polygon=TRUE,axes=FALSE, grid.col=c("green", "red"), max.auc.polygon=TRUE,auc.polygon.col=afcol, print.thres=F,main=i)
#dev.off()

#ggroc(list(s100b=afroc))

}
}


#307686155@qq.com
#18666983305
#更多课程请关注生信碱移
#af
