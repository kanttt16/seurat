
library(Seurat)
library(SeuratData)
library(SingleR)
library(harmony)
library(tidyverse)
library(tibble)
library(ggplot2)
##########################################################################################
#多文件读取
setwd("D:/R Windows/RStudio program/program/scRNA/data/GSE157645")

samples <- list.files("GSE157645/")
samples


seurat_list <- list()
for (sample in samples){
  data.path <- paste0("GSE157645/",sample)
  
  #读取10x
  seurat_data <- Read10X(data.dir = data.path)
  seurat_obj <- CreateSeuratObject(counts = seurat_data,project = sample,
                                   min.features=200,min.cell=3)
  seurat_list <- append(seurat_list,seurat_obj)
}


seurat_combined <- merge(seurat_list[[1]],
                         y=seurat_list[-1],
                         add.cell.ids = samples)
#layers融合
Seurat_data1 <- JoinLayers(seurat_combined)

#######################################################################################
#10x类型读入示例(单)
data1 = Read10X("D:/R Windows/RStudio program/program/scRNA/data/GSE291735_RAW/")
#转化为seurat类型
Seurat_data1 = CreateSeuratObject(counts=data1,project="GSE291735",
                                  min.features=200,
                                  min.cell=3)
#CreateSeuratObject是一个初步筛选
#首先counts=data1是把data1转入到Seurat_data1@assays$RNA$counts
#project就是样本名传入到meta.data，min.features,min.cell是质控
#另外在做多个测序比对时，meta.data的行名可能重复，所以要标注
Seurat_data1 <- merge(Seurat_data1,add.cell.ids="GSE291735")

##########################################################################################
#seurat函数详解，总的来说是一个r4类
#assay是表达矩阵
#meta_data有细胞类型，count等
#seurat实际上就是一个大数据库，每个条目相当于一个超链接
#这个超链接用@、$依次取下，当然也有[[]]
#而data@assays$RNA$counts是储存原始矩阵的地方
#meta.data储存细胞注释，相当于临床信息
#标准化方式，scaledata，sctdata

dim(Seurat_data1)#基因数*细胞数

#查看一下
table(Seurat_data1@meta.data$orig.ident)

table(str_split(colnames(Seurat_data1),"-",simplify=T)[,2])

#以上这两步是为了防止，在实验时会用不同细胞做多次重复
#于是会发生在meta.data的基因片段（列名）出现后缀1，2，3，4等等
#先前的merge函数，会定向匹配，即data.n赋名于ids，从而达到GSE291735_AAACCCAAGGAAGTGA-1的区分效果
#table(Seurat_data1@meta.data$orig.ident)是查看了meta.data的第一列，统计诸如GSE291735 有多少
#table(str_split(colnames(Seurat_data1),"-",simplify=T)[,2])，是把列名提取最后一个数再统计

colnames(Seurat_data1)[Seurat_data1$orig.ident == "GSM4771998"] <- sub("-1$", "-2", colnames(Seurat_data1)[Seurat_data1$orig.ident == "GSM4771998"])

Seurat_data1 <- AddMetaData(object = Seurat_data1,
                            metadata = str_split(colnames(Seurat_data1),"-",simplify=T)[,2],
                            col.name = "orig.ident")

#以上在质控之前的工作流：
#首先，我们预先知道这个数据做了很多实验，它的meta.data的列名会呈现如AAACCCAAGGAAGTGA-1，AAACCCAAGGAAGTGA-2，等，
#那么我们先用read10x读入，然后对他进行creatseuratobject，
#然后用joinlayers融合，这一步会把之前说的像-1，-2放到同一列名里面，
#但是-1，-2不足以区分他们的“身份”，所以Seurat_data1 <- merge(Seurat_data1,Seurat_data2，add.cell.ids=c（"GSE291735"，"GSE291736"）)这样来增强识别，
#然后呢，在meta.data$orig.ident中也没有储存身份，所以要Seurat_data1 <- AddMetaData(object = Seurat_data1,
#                                                                metadata = str_split(colnames(Seurat_data1),"-",simplify=T)[,2],
#                                                               col.name = "orig.ident")，这样列名后缀和orig.ident身份就一致了

#质控
#首先寻找线粒体基因，由于细胞衰老会导致线粒体暴增，所以以数量角度看，
#线粒体的过多和过少都不可取
Seurat_data1[["percent.mt"]] <- PercentageFeatureSet(Seurat_data1,
                                                     pattern = "^MT-")
#寻找“MT-”开头基因（即线粒体）

#然后找红细胞
HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB.genes <- CaseMatch(HB.genes,rownames(Seurat_data1))
view(HB.genes)#匹配存在的红细胞

Seurat_data1[["percent.HB"]] <- PercentageFeatureSet(Seurat_data1,features = HB.genes)

#相关性
FeatureScatter(Seurat_data1,"nCount_RNA","percent.mt",group.by = "orig.ident")#数越小，代表数据不受线粒体污染越小，质控越好
FeatureScatter(Seurat_data1,"nCount_RNA","nFeature_RNA",group.by = "orig.ident")#越大越好

theme.set1 <- theme(axis.title.x = element_blank())
plot.features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.HB")
group = "orig.ident"
#在质控（数据切割）前先用小提琴图探路
plots = list()
for(i in c(1:length(plot.features))){
  plots[[i]] = VlnPlot(Seurat_data1,group.by = group,
                       pt.size = 0,
                       features = plot.features[i])+theme.set1
  
}
library(patchwork)
violin <- wrap_plots(plots=plots,norw=2)
violin


#接下来用quantile查看前10%后10%的跨度，若每1%跨度过大要适当考虑删去
quantile(Seurat_data1$nFeature_RNA,seq(0.01,0.1,0.01))
quantile(Seurat_data1$nFeature_RNA,seq(0.9,1,0.01))

quantile(Seurat_data1$nCount_RNA,seq(0.01,0.1,0.01))
quantile(Seurat_data1$nCount_RNA,seq(0.9,1,0.01))

quantile(Seurat_data1$percent.mt,seq(0.9,1,0.01))
quantile(Seurat_data1$percent.HB,seq(0.9,1,0.01))

#然后根据目前结果来设置质控标准
minGene=200
maxGene=2500
minUMI=500
pctMT=8

#筛选
Seurat_data1 <- subset(Seurat_data1,
                       nFeature_RNA>minGene&
                         nFeature_RNA<maxGene&
                         nCount_RNA>minUMI&
                         percent.mt<pctMT)
dim(Seurat_data1)

for(i in c(1:length(plot.features))){
  plots[[i]] = VlnPlot(Seurat_data1,group.by = group,
                       pt.size = 0,
                       features = plot.features[i])+theme.set1
  
}
library(patchwork)
violin <- wrap_plots(plots=plots,norw=2)
violin


#接下来是双细胞质控
#由于双细胞会导致某个片段出现高counts所以原则上删去
#实现的方法是k邻近找到它
library(DoubletFinder)

Seurat_data1 <- Seurat_data1%>%
  NormalizeData()%>%#标准化
  FindVariableFeatures()%>%#寻找高变基因
  ScaleData()#z-scale
##双细胞没写完，不过先不写了因为暂时用不上


#去除细胞周期（即生长周期）的影响
#正常来说，我们希望能够分开不同种类和亚群的细胞
#而不希望对一个细胞来说，在它的不同生长周期被划分成了不同的亚群
#这会导致结果可信度的受到影响

#细胞周期评分
Seurat_data1 <- NormalizeData(Seurat_data1)

#生长周期有G1，S，G2,M期
#获取G2M相关基因
g2m_genes <- cc.genes$g2m.genes
g2m_genes <- CaseMatch(search = g2m_genes,
                       match = rownames(Seurat_data1))


s_genes <- cc.genes$s.genes
s_genes <- CaseMatch(search = s_genes,
                     match = rownames(Seurat_data1))

#对每期进行评分
Seurat_data1 <- CellCycleScoring(Seurat_data1,
                                 g2m.features = g2m_genes,
                                 s.features = s_genes)

#看看每期多少
table(Seurat_data1$Phase)

#PCA分析
Seurat_data1 <- RunPCA(Seurat_data1)#主成分
Seurat_data1 <- RunTSNE(Seurat_data1)#前十个维度

setwd("D:/R Windows/RStudio program/program/scRNA/GSE157645_plot")
DimPlot(Seurat_data1,group.by = "Phase",reduction ="tsne")

#标准化，高变，放缩
#每个细胞最初包含相同数量的rna分子假设
#结果储存在layer$data
Seurat_data1 <- NormalizeData(Seurat_data1,normalization.method = "LogNormalize")#标准化

#在counts_matrix[1:100, 1:100]可以知道，大多数细胞对多数的基因表达是0
#那么就要找大量的高变基因
Seurat_data1 <- FindVariableFeatures(Seurat_data1,normalization.method = "vst")#寻找高变基因

#给每个基因进行权重的中心配比，基于同等权重
Seurat_data1 <- ScaleData(Seurat_data1,vars.to.regress = c("S.Score","G2M.Score"))#z-scale


#PCA,聚类，去除批次效应（harmony算法，不同样本均匀分布）
install.packages("Rogue")
install.packages("clustree")
install.packages("harmony")

library(Rogue)
library(clustree)
library(harmony)

#先把高变基因画出来
Seurat_data1 <- FindVariableFeatures(Seurat_data1)
pdf(file = "1.pdf",width = 7,height = 6)
VariableFeaturePlot(object = Seurat_data1)
dev.off()

#再画一个标注
top10 <- head(VariableFeatures(Seurat_data1),10)
pdf(file = "2.pdf",width = 7,height = 6)
LabelPoints(plot= VariableFeaturePlot(object = Seurat_data1),
            points = top10,repel = T)
dev.off()

#之前已经跑过RUNPCA，不用再跑一遍
#可以直接画PCA
pdf(file = "3.pdf",width = 7,height = 6)
DimPlot(object = Seurat_data1,reduction = "pca")
dev.off()

#PCA基因相关
pdf(file = "4.pdf",width = 7,height = 6)
VizDimLoadings(object = Seurat_data1,dims = 1:4,reduction = "pca",nfeatures = 20)
dev.off()

#热图
pdf(file = "5.pdf",width = 10,height = 9)
DimHeatmap(object = Seurat_data1,dims = 1:4,cells = 500,reduction = "pca",nfeatures = 30,ncol = 2)
dev.off()

#累计贡献
pdf(file = "6.pdf",width = 7,height = 6)
ElbowPlot(Seurat_data1,ndims = 50)
dev.off()

#确定每个pc百分比
pct <- Seurat_data1[["pca"]]@stdev / sum(Seurat_data1[["pca"]]@stdev) * 100

#做累加和
cumu <- cumsum(pct)
cumu

pcs=40

#harmony,消除来自不同人，不同样本和技术差异导致的偏差
Seurat_data1 <- RunHarmony(Seurat_data1,group.by.vars="orig.ident",
                           assay.use="RNA",max.iter=20)


#选取分辨率
seq=seq(0.1,2,by=0.1)
Seurat_data1 <- FindNeighbors(Seurat_data1,dims = pcs)

for (res in seq){
  Seurat_data1 = FindClusters(Seurat_data1,resolution = res)
}

#画图
p1 = clustree(Seurat_data1,prefix = "RNA_snn_res.")+coord_flip()
p=p1+plot_layout(widths = c(3,1))
ggsave("RNA_snn_res.png",p,width=30,height = 14)

#降维聚类，这里harmony和pca都可
Seurat_data1 <- FindNeighbors(Seurat_data1,
                              reduction = "harmony",
                              dims = pcs)%>%FindClusters(resolution = 1)
Seurat_data1 <- RunUMAP(Seurat_data1,reduction = "harmony",dims =1: pcs)%>%RunTSNE(dims=1:pcs,reduction = "harmony")
                                                                                 

#画图
pdf(file = "7.pdf",width = 7,height = 6)
DimPlot(Seurat_data1,reduction = "umap",label = T)
dev.off()


pdf(file = "8.pdf",width = 7,height = 6)
DimPlot(Seurat_data1,reduction = "umap",label = F,group.by = "orig.ident")
dev.off()

pdf(file = "9.pdf",width = 7,height = 6)
DimPlot(Seurat_data1,reduction = "umap",label = F,group.by = "Is_Double")
dev.off()

pdf(file = "10.pdf",width = 7,height = 6)
DimPlot(Seurat_data1,reduction = "tsne",label = T)
dev.off()


pdf(file = "11.pdf",width = 7,height = 6)
DimPlot(Seurat_data1,reduction = "tsne",label = F,group.by = "orig.ident")
dev.off()

pdf(file = "12.pdf",width = 7,height = 6)
DimPlot(Seurat_data1,reduction = "tsne",label = F,group.by = "Is_Double")
dev.off()









































