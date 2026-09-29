######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")
#BiocManager::install("GSEABase")
#BiocManager::install("GSVA")

#install.packages("ggpubr")


#引用包
library(reshape2)
library(ggpubr)
library(limma)
library(GSEABase)
library(GSVA)
library(BiocParallel)

expFile="merge.txt"                #表达数据文件
clusterFile="PRGcluster.txt"       #分型结果文件
gmtFile="immune.gmt"               #免疫基因集文件
setwd("G:\\吉大\\放疗科数据\\A生信分析相关\\14dead\\26.ssGSEA")     #设置工作目录

#读取表达输入文件,并对输入文件整理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)

#读取基因集文件
geneSets <- getGmt(gmtFile)

#ssGSEA分析
param <- gsvaParam(
  expr = data,
  geneSets = geneSets,
  minSize = 1,
  maxSize = 500
)

ssgseaScore <- gsva(
  param,
  verbose = TRUE,
  BPPARAM = SerialParam()  # 单线程运行，避免报错
)
#对ssGSEA打分进行矫正
normalize=function(x){
  return((x-min(x))/(max(x)-min(x)))}
ssgseaScore=normalize(ssgseaScore)
#输出ssGSEA打分结果
ssgseaOut=rbind(id=colnames(ssgseaScore), ssgseaScore)
write.table(ssgseaOut,file="ssGSEA.result.txt",sep="\t",quote=F,col.names=F)

#读取分型文件
cluster=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)

#数据合并
ssgseaScore=t(ssgseaScore)
sameSample=intersect(row.names(ssgseaScore), row.names(cluster))
ssgseaScore=ssgseaScore[sameSample,,drop=F]
cluster=cluster[sameSample,,drop=F]
scoreCluster=cbind(ssgseaScore, cluster)

#把数据转换成ggplot2输入文件
data=melt(scoreCluster, id.vars=c("PRGcluster"))
colnames(data)=c("PRGcluster", "Immune", "Fraction")

#绘制箱线图
bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(levels(factor(data[,"PRGcluster"])))]
p=ggboxplot(data, x="Immune", y="Fraction", color="PRGcluster", 
     ylab="Immune infiltration",
     xlab="",
     legend.title="PRGcluster",
     palette=bioCol)
p=p+rotate_x_text(50)

#保存图形
pdf(file="boxplot.pdf", width=8, height=6.5)
p+stat_compare_means(aes(group=PRGcluster),symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "")),label = "p.signif")
dev.off()


######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

