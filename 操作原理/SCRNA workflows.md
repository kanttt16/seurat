Main packages: seurat (supporting cran)
URL :[seurat offical](https://satijalab.org/seurat/articles/pbmc3k_tutorial)


# 导入数据（后续会添加txt，h5等格式读取）
They use the pbmc3k(10x) from GEO,The `[Read10X()](https://satijalab.org/seurat/reference/read10x)` function reads in the output of the [cellranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger) pipeline from 10X, returning a unique molecular identified (UMI) count matrix. The values in this matrix represent the number of molecules for each feature (i.e. gene; row) that are detected in each cell (column). Note that more recent versions of cellranger now also output using the [h5 file format](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/h5_matrices), which can be read in using the Read10X_h5() function in Seurat.
```r
#easy read
data1 = Read10X("D:/R Windows/RStudio program/program/scRNA/data/GSE291735_RAW/")
```

# 创建一个s4类的seurat格式
and then ,create a seuart
```r
library(Seurat)
library(SeuratData)
library(SingleR)
library(harmony)
library(tidyverse)
library(tibble)
library(ggplot2)

Seurat_data1 = CreateSeuratObject(counts=data1,project="GSE291735",
                                  min.features=200,
                                  min.cell=3)

```
这里着重要考虑导入的seurat结构，
![seurat 结构](C:/Users/kantt/OneDrive/Pictures/Screenshots/seurat结构.png)
assays的目录下有 layers $ count & data
他们都是稀疏矩阵，是细胞于对应基因的表达情况，
区别在于count是源数据，而data是标准化数据（这有什么区别我目前没搞清楚）

[为什么采用稀疏矩阵，原因在于储存结构，如下](https://zhuanlan.zhihu.com/p/714294089)
```r
# 密集矩阵
dense.size <- object.size(as.matrix(pbmc.data))
dense.size
# 709591472 bytes

# 稀疏矩阵
sparse.size <- object.size(pbmc.data)
sparse.size
# 29905192 bytes

# 两者相差倍数
dense.size/sparse.size
# 23.7 bytes
```
可以看到密集与稀疏矩阵的储存倍数差距相当大

# 筛选凋亡细胞（线粒体）与红细胞
接下来是scRNA的中游workflows
标准：
Seurat allows you to easily explore QC metrics and filter cells based on any user-defined criteria. A few QC metrics [commonly used](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4758103/) by the community include

- The number of unique genes detected in each cell.
    - Low-quality cells or empty droplets will often have very few genes
    - Cell doublets or multiplets may exhibit an aberrantly high gene count
- Similarly, the total number of molecules detected within a cell (correlates strongly with unique genes)
- The percentage of reads that map to the mitochondrial genome
    - Low-quality / dying cells often exhibit extensive mitochondrial contamination
    - We calculate mitochondrial QC metrics with the `[PercentageFeatureSet()](https://satijalab.org/seurat/reference/percentagefeatureset)` function, which calculates the percentage of counts originating from a set of features
    - We use the set of all genes starting with `MT-` as a set of mitochondrial genes
                                                                 --from [url](https://satijalab.org/seurat/articles/pbmc3k_tutorial)

这指向了两个问题：
1.在上游测序时定义了“低质量细胞”和空泡细胞，前者指的是衰老或者濒死细胞，这类细胞的明显特征是线粒体的表达明显增高，因此我们可以根据线粒体的百分位数来选择去除线粒体高表达；相反的是，后者细胞往往是无表达基因或者很少，而且表达意愿不高
2.双细胞或多细胞问题，细胞双胞体或多胞体可能表现出异常高的基因计数，在==原则上==去除

```r
Seurat_data1[["percent.mt"]] <- PercentageFeatureSet(Seurat_data1,
                                                     pattern = "^MT-")
#寻找“MT-”开头基因（即线粒体）

#然后找红细胞
HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB.genes <- CaseMatch(HB.genes,rownames(Seurat_data1))
view(HB.genes)#匹配存在的红细胞

Seurat_data1[["percent.HB"]] <- PercentageFeatureSet(Seurat_data1,features = HB.genes)

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
```

# 初步质控
```r
#接下来用quantile查看前10%后10%的跨度，若每1%跨度过大要适当考虑删去
quantile(Seurat_data1$nFeature_RNA,seq(0.01,0.1,0.01))
quantile(Seurat_data1$nFeature_RNA,seq(0.9,1,0.01))

quantile(Seurat_data1$nCount_RNA,seq(0.01,0.1,0.01))
quantile(Seurat_data1$nCount_RNA,seq(0.9,1,0.01))

quantile(Seurat_data1$percent.mt,seq(0.9,1,0.01))
quantile(Seurat_data1$percent.HB,seq(0.9,1,0.01))

#QC
#然后根据目前结果来设置质控标准
minGene=200
maxGene=9000
minUMI=500
pctMT=8

#筛选
Seurat_data1 <- subset(Seurat_data1,
                       nFeature_RNA>minGene&
                         nFeature_RNA<maxGene&
                         nCount_RNA>minUMI&
                         percent.mt<pctMT)

```

这一步是质控全过程，简单说的说就是先拿小提琴图探路，找n_feature,n_count,mt,HB的表达
`nFeature_RNA`是对所有细胞检测到的基因数量的小提琴图，横坐标为样品名，纵坐标为每个细胞中包含的基因数量，图中的每个点代表一个细胞；
`nCount_RNA`为样品所有细胞内检测到的分子总数的小提琴图，横坐标为样品名，纵坐标为每个细胞内检测到的分子总数，图中的每个点代表一个细胞
`percent.mt`为样品所有细胞的线粒体基因比例的小提琴图，横坐标为样品名，纵坐标为每个细胞的线粒体基因比例，图中的每个点代表一个细胞；HB同理
由此可以知道在哪些数值上属于高表达

```r
#然后要查看相关性
FeatureScatter(Seurat_data1,"nCount_RNA","percent.mt",group.by = "orig.ident")#数越小，代表数据不受线粒体污染越小，质控越好
FeatureScatter(Seurat_data1,"nCount_RNA","nFeature_RNA",group.by = "orig.ident")#越大越好
#FeatureScatter函数会直接画好mt~n_count,n_count~n_feature
```
随着测序深度的增加，单细胞所检测到的基因数量也在增加，两者是有一定关联性的。
对含有线粒体基因的细胞进行筛选是因为线粒体基因比例过高的细胞通常为低质量的细胞，这些数据杂糅进来会干扰细胞分群。
`nCount_RNA vs percent.mt`，是线粒体基因比例与细胞中检测到的分子总数的关系的散点图，通常两者没有什么相关关系，我们通常对这两个指标的标准是：==线粒体基因比例应相对较低以排除低质量/濒死的细胞==，细胞中检测到的分子总数也不宜过高以排除双胞和多胞细胞的影响。
`nCount_RNA vs nFeature_RNA`，是基因数量与细胞中检测到的分子总数的关系的散点图，即测序深度与基因数量的关系。高质量的测序数据中两者==基本处于正相关的关系==，**但要排除由于双胞和多胞造成的分子数量过大的部分数据，即右上方离群点**

然后根据以上结果自适应过滤就行
官网推荐过滤指标：
In the example below, we visualize QC metrics, and use these to filter cells.
- We filter cells that have unique feature counts over 2,500 or less than 200
- We filter cells that have >5% mitochondrial counts


# 去除生长周期的影响
```r
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

```

# 数据进一步标准化
```r
#标准化，高变，放缩
#每个细胞最初包含相同数量的rna分子假设
#结果储存在layer$data
Seurat_data1 <- NormalizeData(Seurat_data1,normalization.method = "LogNormalize")#标准化

#在counts_matrix[1:100, 1:100]可以知道，大多数细胞对多数的基因表达是0
#那么就要找大量的高变基因
Seurat_data1 <- FindVariableFeatures(Seurat_data1,normalization.method = "vst")#寻找高变基因

#给每个基因进行权重的中心配比，基于同等权重
Seurat_data1 <- ScaleData(Seurat_data1,vars.to.regress = c("S.Score","G2M.Score"))#z-scale


```


# 查看高变基因，并且用pca寻找最佳dim
```r
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
#pcs记为dim，是滚石图所得结果
pcs=35
```
其中ElbowPlot作用见 ：[[PCA原理以及在scRNA的应用]]


# Harmony
```r
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

```


# 降维聚类可视化
这里注意，先前我们跑过了PCA和harmony，这里的图片呈现只能选择一种类型

```r
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


```
umap和tsne的意义见：[[PCA原理以及在scRNA的应用]]
