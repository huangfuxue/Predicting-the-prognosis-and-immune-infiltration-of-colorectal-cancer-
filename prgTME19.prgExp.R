

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")


library(limma)           #引用包
expFile="merge.txt"      #表达数据文件
geneFile="gene.txt"      #基因列表文件
setwd("C:\\Users\\ASUS\\Desktop\\结直肠癌坏死性凋亡TME\\19.nrgExp")    #设置工作目录

#读取输入文件，并对数据进行处理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]

#获取坏死性凋亡基因的表达量
gene=read.table(geneFile, header=T, sep="\t", check.names=F)
sameGene=intersect(as.vector(gene[,1]), rownames(data))
geneExp=data[sameGene,]

#输出结果
out=rbind(ID=colnames(geneExp),geneExp)
write.table(out,file="nrgGeneExp.txt",sep="\t",quote=F,col.names=F)




