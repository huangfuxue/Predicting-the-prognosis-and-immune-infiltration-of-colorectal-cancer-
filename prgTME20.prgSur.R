

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages('survival')
#install.packages("survminer")


#引用包
library(limma)
library(survival)
library(survminer)

expFile="nrgGeneExp.txt"     #表达数据文件
cliFile="time.txt"           #生存数据文件
setwd("C:\\Users\\ASUS\\Desktop\\结直肠癌坏死性凋亡TME\\20.nrgSur")     #设置工作目录

#读取表达文件，并对输入文件整理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]
data=t(data)
rownames(data)=gsub("(.*?)\\_(.*?)", "\\2", rownames(data))

#读取生存数据
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
cli$futime=cli$futime/365

#数据合并
sameSample=intersect(row.names(data), row.names(cli))
data=data[sameSample,]
cli=cli[sameSample,]
rt=cbind(cli, data)

#对基因进行循环，找出预后相关的基因
outTab=data.frame()
km=c()
for(i in colnames(rt[,3:ncol(rt)])){
	#cox分析
	cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
	coxSummary = summary(cox)
	coxP=coxSummary$coefficients[,"Pr(>|z|)"]
	outTab=rbind(outTab,
			         cbind(id=i,
			         HR=coxSummary$conf.int[,"exp(coef)"],
			         HR.95L=coxSummary$conf.int[,"lower .95"],
			         HR.95H=coxSummary$conf.int[,"upper .95"],
			         pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
			        )
	#km分析
	data=rt[,c("futime", "fustat", i)]
	colnames(data)=c("futime", "fustat", "gene")
	#获取最优cutoff
	res.cut=surv_cutpoint(data, time = "futime", event = "fustat", variables =c("gene"))
	res.cat=surv_categorize(res.cut)
	fit=survfit(Surv(futime, fustat) ~gene, data = res.cat)
	#print(paste0(i, " ", res.cut$cutpoint[1]))
	#比较高低表达生存差异
	diff=survdiff(Surv(futime, fustat) ~gene,data =res.cat)
	pValue=1-pchisq(diff$chisq, df=1)
	km=c(km, pValue)
	#对pvalue<0.05的基因绘制生存曲线
	if(pValue<0.05){
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
						   legend.title=i,
						   legend.labs=c("high","low"),
						   xlab="Time(years)",
						   ylab="Overall survival",
						   palette=c("red", "blue"),
						   break.time.by=1,
						   conf.int=F,
						   risk.table=F,
						   risk.table.title="",
						   risk.table.height=.25)
		pdf(file=paste0("sur.", i, ".pdf"),onefile = FALSE,
			width = 5,         #图片的宽度
			height =4.5)         #图片的高度
		print(surPlot)
		dev.off()
	}
}

#输出单因素的结果
outTab=cbind(outTab, km)
write.table(outTab,file="uniCox.txt",sep="\t",row.names=F,quote=F)




