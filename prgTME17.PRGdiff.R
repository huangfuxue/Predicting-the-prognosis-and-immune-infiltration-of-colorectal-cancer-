######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("reshape2")
#install.packages("ggpubr")


#引用包
library(limma)
library(reshape2)
library(ggpubr)

expFile="TCGA.TPM.txt"      #表达数据矩阵
geneFile="gene.txt"         #基因列表
setwd("C:\\Users\\ASUS\\Desktop\\140prgTME\\17.PRGdiff")       #设置工作目录

#读取输入文件
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)

#提取细胞焦亡基因的表达量
gene=read.table(geneFile, header=T, sep="\t", check.names=F)
sameGene=intersect(as.vector(gene[,1]), row.names(data))
data=data[sameGene,]

#正常和肿瘤数目
group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
conNum=length(group[group==1])       #正常组样品数目
treatNum=length(group[group==0])     #肿瘤组样品数目
sampleType=c(rep(1,conNum), rep(2,treatNum))

#提取显著差异的基因
exp=log2(data+1)
exp=as.data.frame(t(exp))
exp=cbind(exp, Type=sampleType)
exp$Type=ifelse(exp$Type==1, "Normal", "Tumor")
sigGene=c()
for(i in colnames(exp)[1:(ncol(exp)-1)]){
	if(sd(exp[,i])<0.001){next}
	wilcoxTest=wilcox.test(exp[,i] ~ exp[,"Type"])
	pvalue=wilcoxTest$p.value
	if(wilcoxTest$p.value<0.05){
		sigGene=c(sigGene, i)
	}
}
sigGene=c(sigGene, "Type")
exp=exp[,sigGene]

#把数据转换成ggplot2输入文件
data=melt(exp, id.vars=c("Type"))
colnames(data)=c("Type", "Gene", "Expression")

#绘制箱线图
p=ggboxplot(data, x="Gene", y="Expression", color = "Type", 
	     ylab="Gene expression",
	     xlab="",
	     legend.title="Type",
	     palette = c("blue", "red"),
	     width=1)
p=p+rotate_x_text(60)
p1=p+stat_compare_means(aes(group=Type),
	      method="wilcox.test",
	      symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
	      label = "p.signif")

#输出箱线图
pdf(file="boxplot.pdf", width=9, height=6)
print(p1)
dev.off()


######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

