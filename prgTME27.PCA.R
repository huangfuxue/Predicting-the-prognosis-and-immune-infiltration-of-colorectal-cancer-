######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("limma")

install.packages("ggplot2")


#引用包
library(limma)
library(ggplot2)

expFile="prgGeneExp.txt"         #表达数据文件
clusterFile="PRGcluster.txt"     #分型结果文件
setwd("G:\\吉大\\放疗科数据\\A生信分析相关\\14dead\\27.PCA")      #设置工作目录

#读取输入文件,并对输入文件进行整理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]
data=t(data)

#PCA分析
data.pca=prcomp(data, scale. = TRUE)
pcaPredict=predict(data.pca)
write.table(pcaPredict, file="newTab.xls", quote=F, sep="\t")

#读取分型文件
cluster=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)
PRGcluster=as.vector(cluster[,1])

#设置颜色
bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
prgCluCol=bioCol[1:length(levels(factor(PRGcluster)))]

#可视化
PCA=data.frame(PC1=pcaPredict[,1], PC2=pcaPredict[,2], PRGcluster=PRGcluster)
PCA.mean=aggregate(PCA[,1:2], list(PRGcluster=PCA$PRGcluster), mean)
pdf(file="PCA.pdf", width=6.5, height=5)
ggplot(data = PCA, aes(PC1, PC2)) + geom_point(aes(color = PRGcluster)) +
	scale_colour_manual(name="PRGcluster", values =prgCluCol)+
    theme_bw()+
    theme(plot.margin=unit(rep(1.5,4),'lines'))+
    annotate("text",x=PCA.mean$PC1, y=PCA.mean$PC2, label=PCA.mean$PRGcluster, cex=7)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()


######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

