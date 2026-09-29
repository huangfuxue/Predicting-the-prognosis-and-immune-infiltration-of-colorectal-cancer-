######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

#install.packages("pheatmap")


library(pheatmap)         #引用包
expFile="prgGeneExp.txt"         #表达输入文件
clusterFile="PRGcluster.txt"     #分型的结果文件
cliFile="clinical.txt"           #临床数据文件
setwd("H:\\生信分析相关\\A生信分析相关\\14dead\\24.heatmap")     #设置工作目录

#读取输入文件
exp=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
exp=t(exp)
cluster=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)

#合并表达和分型数据
sameSample=intersect(row.names(exp), row.names(cluster))
exp=exp[sameSample, , drop=F]
cluster=cluster[sameSample, , drop=F]
expCluster=cbind(exp, cluster)
Project=gsub("(.*?)\\_.*", "\\1", rownames(expCluster))
rownames(expCluster)=gsub("(.*?)\\_(.*?)", "\\2", rownames(expCluster))
expCluster=cbind(expCluster, Project)

#合并临床数据
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
cli[,"Age"]=ifelse(cli[,"Age"]=="unknow", "unknow", ifelse(cli[,"Age"]>65,">65","<=65"))
sameSample=intersect(row.names(expCluster), row.names(cli))
expCluster=expCluster[sameSample,,drop=F]
cli=cli[sameSample,,drop=F]
data=cbind(expCluster, cli)

#提取热图数据
data=data[order(data$PRGcluster),]
Type=data[,((ncol(exp)+1):ncol(data))]
data=t(data[,1:ncol(exp)])

#聚类颜色
bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
ann_colors=list()
prgCluCol=bioCol[1:length(levels(factor(Type$PRGcluster)))]
names(prgCluCol)=levels(factor(Type$PRGcluster))
ann_colors[["PRGcluster"]]=prgCluCol

#热图可视化
pdf("heatmap.pdf", width=6, height=5.5)
pheatmap(data,
         annotation=Type,
         annotation_colors = ann_colors,
         color = colorRampPalette(c(rep("blue",5), "white", rep("red",5)))(100),
         cluster_cols =F,
         cluster_rows =F,
         scale="row",
         show_colnames=F,
         fontsize=6,
         fontsize_row=4,
         fontsize_col=6)
dev.off()


######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

