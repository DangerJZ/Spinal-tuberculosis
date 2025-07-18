library(Seurat)
library(data.table)
library(dplyr)
library(tidyverse)
library(cowplot)
library(clustree)
library(harmony)
library(patchwork)
library(rhdf5)
library(GENIE3)
library(AUCell)
library(devtools)

setwd("D:/1/Macrophage")

Idents(sce.all)

Monocytes <- subset(sce.all, idents = "Monocytes")

Idents(Monocytes) <- 'orig.ident'

all.genes=rownames(Monocytes)

g2m_genes <- cc.genes$g2m.genes
g2m_genes <- CaseMatch(search = g2m_genes, match = rownames(Monocytes))

s_genes <- cc.genes$s.genes
s_genes <- CaseMatch(search = s_genes, match = rownames(Monocytes))

options(future.globals.maxSize = 60000 * 1024^3)
Monocytes <- NormalizeData(Monocytes) %>% FindVariableFeatures() %>% ScaleData(features = all.genes) %>% RunPCA()

Monocytes <- CellCycleScoring(Monocytes, g2m.features = g2m_genes, s.features = s_genes)

colnames(Monocytes@meta.data)
table(Monocytes$Phase)
DimPlot(Monocytes, group.by = "Phase")

Monocytes <- NormalizeData(Monocytes, normalization.method = "LogNormalize", scale.factor = 10000)
Monocytes <- FindVariableFeatures(Monocytes, selection.method = "vst", nfeatures = 2000)
Monocytes <- ScaleData(Monocytes, vars.to.regress = c("S.Score", "G2M.Score"))
DimPlot(Monocytes, group.by = "Phase")

top10 <- head(VariableFeatures(Monocytes), 10)
pdf(file = "y1.pdf", width = 7, height = 6)
VariableFeaturePlot(object = Monocytes)
dev.off()

pdf(file = "y2.pdf", width = 7, height = 6)
LabelPoints(plot= VariableFeaturePlot(object = Monocytes), points = top10, repel = TRUE)
dev.off()

Monocytes <- RunPCA(Monocytes, verbose = F)
pdf(file = "y3.pdf", width = 7, height = 6)
DimPlot(object = Monocytes, reduction = "pca")
dev.off()

pdf(file = "y4.pdf", width = 10, height = 9)
VizDimLoadings(object = Monocytes, dims = 1:4, reduction = "pca", nfeatures = 20)
dev.off()

pdf(file = "y5.pdf", width = 10, height = 9)
DimHeatmap(object = Monocytes, dims = 1:4, cells = 500, balanced = TRUE, nfeatures = 30, ncol = 2)
dev.off()

pdf(file = "y6.pdf", width = 10, height = 9)
ElbowPlot(Monocytes, ndims = 50)
dev.off()

pct <- Monocytes [["pca"]]@stdev / sum(Monocytes [["pca"]]@stdev) * 100
pct
cumu <- cumsum(pct)
cumu
co1 = which(cumu > 90 & pct < 5)[1]
co2 = sort(which((pct[1:length(pct)-1]-pct[2:length(pct)]) > 0.1), decreasing = T)[1]+1
pcs = min(co1, co2)
pcs
pcs = 1:15

Monocytes <- RunHarmony(Monocytes, group.by.vars="orig.ident", assay.use="SCT", max.iter.harmony = 20)
table(Monocytes@meta.data$orig.ident)

library(clustree)
seq = seq(0.1,1,by=0.1)
Monocytes <- FindNeighbors(Monocytes, dims = pcs)
for (res in seq){
  Monocytes = FindClusters(Monocytes, resolution = res)
}

library(patchwork)
p8 = clustree(Monocytes, prefix = "RNA_snn_res.")+coord_flip()
p = p8+plot_layout(widths = c(3,1))
p
ggsave("y7_T_RNA_snn_res-.pdf", p, width = 30, height = 14)

for (i in c(10,11,12,13,14,15,16,17,18,19,20)) {
  Monocytes <- FindNeighbors(Monocytes,reduction = "harmony", dims = 1:i) %>% FindClusters(resolution = 0.1)
  Monocytes <- RunUMAP(Monocytes, reduction = "harmony", dims = 1:i)
  plot_i <- print(DimPlot(Monocytes,reduction = "umap", label = TRUE, pt.size = 1, repel = TRUE, label.size = 5))
  plotname <- paste("plot_", i, sep = "")
  assign(plotname, plot_i)
  print(plot_i)
}
p = plot_10+plot_11+plot_12+plot_13+plot_14+plot_15+plot_16+plot_17+plot_18+plot_19+plot_20+plot_layout(ncol = 3)
ggsave(p,filename="dim.10-20.pdf",width = 15,height = 16)

for (i in seq) {
  Monocytes <- FindNeighbors(Monocytes,reduction = "harmony", dims = pcs) %>% FindClusters(resolution = i)
  Monocytes <- RunUMAP(Monocytes, reduction = "harmony", dims = pcs)
  plot_i <- print(DimPlot(Monocytes,reduction = "umap", label = TRUE, pt.size = 1, repel = TRUE, label.size = 5))
  plotname <- paste("plot_", i, sep = "")
  assign(plotname, plot_i)
  print(plot_i)
}
p = plot_0.1+plot_0.2+plot_0.3+plot_0.4+plot_0.5+plot_0.6+plot_0.7+plot_0.8+plot_0.9+plot_1+plot_layout(ncol = 4)
ggsave(p,filename="dim.15,res.0.1-1.pdf",width = 20,height = 12)


Monocytes <- FindNeighbors(Monocytes, reduction = "harmony", dims = pcs) %>% FindClusters(resolution = 0.1)
Monocytes <- RunUMAP(Monocytes, reduction = "harmony", dims = pcs) %>% RunTSNE(dims = pcs, reduction = "harmony")
DimPlot(Monocytes,reduction = "umap",cols = Biocols)

pdf(file = "y8-.pdf", width = 7, height = 6)
DimPlot(Monocytes, reduction = "umap", label = T,cols = Biocols)
dev.off()

pdf(file = "y9.pdf", width = 7, height = 6)
DimPlot(Monocytes, reduction = "umap", label = F,group.by = "orig.ident",cols = Biocols)
dev.off()

pdf(file = "y10.pdf", width = 7, height = 6)
DimPlot(Monocytes, reduction = "tsne", label = T,cols = Biocols)
dev.off()

pdf(file = "y11.pdf", width = 7, height = 6)
DimPlot(Monocytes, reduction = "tsne", label = F,group.by = "orig.ident",cols = Biocols)
dev.off()

Monocytes.markers1 <- FindAllMarkers(Monocytes,only.pos = TRUE,logfc.threshold = 0.25)
write.csv(Monocytes.markers1,file="Monocytes.1.RNA.csv")

markers <- c("ITM2C","HLA-DQB2","HLA-DPB1","BIRC3","CLEC9A","CD1C","LAMP3","LILRA4","CST3","LYZ","FCGR3A","CD83","CCL3","VCAN","C1QA","C1QB","CD68","CD163","C1QC","CCL3L1","EREG","DNAJB1","WDR74","CD14","HLA-DRA","CTSK","NFATC1","FOS","MMP9","ACP5","ATP6V0D2","OCSTAMP","DCSTAMP","SIGLEC15","TNFRSF11A")

p11 <- FeaturePlot(Monocytes,features = markers,ncol = 3)
p11
ggsave("y14.pdf",p11,width = 12,height = 12)

p12 <- DotPlot(Monocytes,features = markers) + RotatedAxis()
p12
ggsave("y15.pdf",p12,width = 10,height = 4)

library(ComplexHeatmap)
library(pheatmap)
library(circlize)
group1_means <- AggregateExpression(
  object = Monocytes,
  features = markers,
  group.by = "RNA_snn_res.0.1",
  slot = "data"
)
group1_means <- as.data.frame(group1_means$RNA)
df <- data.frame(colnames(group1_means))
colnames(df) <- 'class'
top_anno = HeatmapAnnotation(df = df,
                             border = T,
                             show_annotation_name = F,
                             gp = gpar(col = 'black'),
                             col = list(class = c('g0' = "#AB3282",
                                                  'g1' = "#53A85F",
                                                  'g2' = "#F1BB72",
                                                  'g3' = "#D6E7A3",
                                                  'g4' = "#57C3F3",
                                                  'g5' = "#E95C59",
                                                  'g6' = "#E59CC4",
                                                  'g7' = "#9FA3A8")))#颜色设置
marker_exp <- t(scale(t(group1_means),scale = T,center = F))
p=Heatmap(marker_exp,
          cluster_rows = T,
          cluster_columns = T,
          show_column_names = F,
          show_row_names = T,
          column_title = NULL,
          heatmap_legend_param = list(
            title=' '),
          col = colorRampPalette(c("#ADD8E6", "#F0E68C", "#FF4500"))(50),
          border = 'black',
          rect_gp = gpar(col = "black", lwd = 1),
          row_names_gp = gpar(fontsize = 10),
          column_names_gp = gpar(fontsize = 10),
          top_annotation = top_anno)
p
pdf("y12.pdf", width = 5, height = 4)
print(p)
dev.off()

p13 <- VlnPlot(Monocytes,features = markers,stack = T,flip = T) + NoLegend()
p13
ggsave("y16.pdf",p13,width = 8,height = 12)

Monocytes$celltype <- recode(Monocytes@meta.data$seurat_clusters,
                             "0" = "Macro_C1QC",
                             "1" = "Mono_VCAN",
                             "2" = "Macro_CCL3L1",
                             "3" = "Macro_DNAJB1",
                             "4" = "Mono_CD1C",
                             "5" = "DC_CLEC9A",
                             "6" = "Macro_ITM2C",
                             "7" = "DC_LAMP3")
table(Monocytes@meta.data$celltype)

p=DimPlot(Monocytes, split.by = "group1", label = TRUE,cols = Biocols)
p
ggsave("A13.pdf",p,width = 10,height = 4)

p=DimPlot(Monocytes, split.by = "group1", label = F,cols = Biocols)
p
ggsave("A13'.pdf",p,width = 10,height = 4)

p=DimPlot(Monocytes, split.by = "group2", label = T,cols = Biocols)
p
ggsave("A14.pdf",p,width = 8,height = 4)

p=DimPlot(Monocytes, split.by = "group2", label = F,cols = Biocols)
p
ggsave("A14'.pdf",p,width = 8,height = 4)

pdf(file = "A7.pdf", width = 6, height = 4)
DimPlot(Monocytes, reduction = "umap", label = T,cols = Biocols)
dev.off()

pdf(file = "A7'.pdf", width = 6, height = 4)
DimPlot(Monocytes, reduction = "umap", label = F,cols = Biocols)
dev.off()

pdf(file = "A8.pdf", width = 5, height = 4)
DimPlot(Monocytes, reduction = "umap", label = F,group.by = "orig.ident",cols = Biocols)
dev.off()

pdf(file = "A9.pdf", width = 5, height = 4)
DimPlot(Monocytes, reduction = "umap", label = F,group.by = "group1",cols = Biocols)
dev.off()

pdf(file = "A9'.pdf", width = 5.5, height = 4)
DimPlot(Monocytes, reduction = "umap", label = F,group.by = "group2",cols = Biocols)
dev.off()

pdf(file = "A10.pdf", width = 6, height = 4)
DimPlot(Monocytes, reduction = "tsne", label = T,cols = Biocols)
dev.off()

pdf(file = "A10'.pdf", width = 6, height = 4)
DimPlot(Monocytes, reduction = "tsne", label = F,cols = Biocols)
dev.off()

pdf(file = "A11.pdf", width = 5, height = 4)
DimPlot(Monocytes, reduction = "tsne", label = F,group.by = "orig.ident",cols = Biocols)
dev.off()

pdf(file = "A11.pdf", width = 5, height = 4)
DimPlot(Monocytes, reduction = "tsne", label = F,group.by = "group1",cols = Biocols)
dev.off()

pdf(file = "A11'.pdf", width = 5.5, height = 4)
DimPlot(Monocytes, reduction = "tsne", label = F,group.by = "group2",cols = Biocols)
dev.off()

Monocytes.markers3 <- FindAllMarkers(Monocytes, only.pos = TRUE, logfc.threshold =1)
write.csv(Monocytes.markers3, file = "Monocytes.celltype.RNA.csv")

top10 <- Monocytes.markers3 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

markers = as.data.frame(top10[,"gene"])
Monocytes <- ScaleData(Monocytes, features = as.character(unique(markers$gene)))
p = DoHeatmap(Monocytes,
              features = as.character(unique(markers$gene)),
              group.by = "celltype")
p
ggsave("y19.pdf", p, width = 10,height = 10)

allCells = names(Idents(Monocytes))
allType = levels(Idents(Monocytes))
choose_Cells = unlist(lapply(allType, function(x){
  cgCells = allCells[Idents(Monocytes) == x]
  n = min(table(Monocytes@meta.data$celltype))
  if (length(cgCells) >= n) {
    cg = sample(cgCells, n)
  } else {
    cg = sample(cgCells, length(cgCells))
  }
  cg
}))

cg_sce = Monocytes[, allCells %in% choose_Cells]
table(Idents(cg_sce))

p = DoHeatmap(cg_sce,
              features = as.character(unique(markers$gene)),
              group.by = "celltype")
p
ggsave("y20.pdf", p, width = 10,height = 10)

cell.prop <- as.data.frame(prop.table(table(Monocytes@meta.data$celltype,Monocytes@meta.data$orig.ident)))
colnames(cell.prop)<-c("cluster", "group", "proportion")

p= ggplot(cell.prop, aes(group,proportion,fill=cluster))+
  geom_bar(stat = "identity",position = "fill")+
  ggtitle("")+
  theme_bw()+
  theme(axis.ticks.length = unit(0.5, "cm"))+
  guides(fill=guide_legend(title = NULL))
p
ggsave("y21.pdf", p, width = 7,height = 5)

cell.prop <- as.data.frame(prop.table(table(Monocytes@meta.data$celltype,Monocytes@meta.data$group1)))
colnames(cell.prop)<-c("cluster", "group", "proportion")

p= ggplot(cell.prop, aes(group,proportion,fill=cluster))+
  geom_bar(stat = "identity",position = "fill")+
  ggtitle("")+
  theme_bw()+
  theme(axis.ticks.length = unit(0.5, "cm"))+
  guides(fill=guide_legend(title = NULL))
p
ggsave("y22.pdf", p, width = 4,height = 4)

cell.prop <- as.data.frame(prop.table(table(Monocytes@meta.data$celltype,Monocytes@meta.data$group2)))
colnames(cell.prop)<-c("cluster", "group", "proportion")

p= ggplot(cell.prop, aes(group,proportion,fill=cluster))+
  geom_bar(stat = "identity",position = "fill")+
  ggtitle("")+
  theme_bw()+
  theme(axis.ticks.length = unit(0.5, "cm"))+
  guides(fill=guide_legend(title = NULL))
p
ggsave("y23.pdf", p, width = 4,height = 4)

DefaultAssay(Monocytes) <- "RNA"
options(future.globals.maxSize = 20 * 1024^3)
Monocytes <- NormalizeData(Monocytes)
dim(Monocytes)

expr = as.matrix(GetAssayData(object = Monocytes@assays$RNA, layer="data"))
dim(expr)
expr <- expr[rowSums(expr)>0,]
dim(expr)
library(readxl)
geneSets2 = read_excel("genesets.xls",col_names=F)
write.table(geneSets2,file = "my_genesets.gmt",sep="\t",row.names = F,col.names = F,quote = F)
library(GSEABase)
geneSets=getGmt("my_genesets.gmt",geneIdType=SymbolIdentifier())

library(rhdf5)
library(GSVA)

ssgsea_params <- ssgseaParam(exprData = expr, geneSets = geneSets, 
                             minSize = 1, maxSize = Inf, alpha = 0.25, normalize = TRUE)
gsvaResult <- gsva(ssgsea_params)

normalize=function(x){return((x-min(x))/(max(x)-min(x)))}
gsvaResult=normalize(gsvaResult)
gsvaOut=rbind(id=colnames(gsvaResult),gsvaResult)
write.table(gsvaOut, file="ssgseaOut2.txt",sep="\t",col.names = F,quote = F)
Monocytes <- AddMetaData(Monocytes, metadata = t(gsvaResult))

p <- FeaturePlot(Monocytes, split.by = "group1", features = rownames(gsvaResult), reduction = "umap", ncol=2)
p
ggsave("z49_Monocytes_OC_group1.pdf",p,width = 10,height = 3)

p <- FeaturePlot(Monocytes, split.by = "group2", features = rownames(gsvaResult), reduction = "umap", ncol=2)
p
ggsave("z49_Monocytes_OC_group2.pdf",p,width = 6,height = 3)

data <- data.frame(group1 = Monocytes$group1,celltype = Monocytes$celltype,Migration_scores = Monocytes$`OC scores`)
data$group1 <- factor(data$group1, levels = c("Health", "Mild", "Severe"))
p <- ggviolin(data, x = "group1", y = "Migration_scores", fill = "group1",
              palette = c("#00AFBB", "#E7B800", "#FC4E07"),  # 设置颜色
              ylab = "Average Expression", xlab = "Group",
              add = "boxplot", add.params = list(fill = "white"))  # 添加箱线图
p
p <- p + stat_compare_means(method = "wilcox.test", comparisons = my_comparisons, label = "p.signif")
p
ggsave("z49_OC-group1-violin.pdf",plot = p,width = 3.5,height = 4.5)

p <- ggviolin(data, x = "group1", y = "Migration_scores", fill = "group1",
              palette = c("#00AFBB", "#E7B800", "#FC4E07"),  # 设置颜色
              ylab = "Average Expression", xlab = "Group",
              add = "boxplot", add.params = list(fill = "white")) +  # 添加箱线图
  facet_wrap(~ celltype, scales = "free",ncol=4)  # 按celltype分面，并允许每个分面有不同的y轴范围
p
# 假设你想在每个细胞类型中比较group1的分组差异
p <- p + stat_compare_means(method = "wilcox.test", label = "p.signif", comparisons = my_comparisons)
p
ggsave("z49_OC-group1-violin-celltype.pdf",plot = p,width = 10,height = 8)

data2 <- data.frame(group1 = Monocytes$group2,celltype = Monocytes$celltype,Migration_scores = Monocytes$`OC scores`)
data2$group1 <- factor(data2$group1, levels = c("Health", "Tuberculosis"))
p <- ggviolin(data2, x = "group1", y = "Migration_scores", fill = "group1",
              palette = c("#00AFBB", "#E76980"),  # 设置颜色
              ylab = "Average Expression", xlab = "Group",
              add = "boxplot", add.params = list(fill = "white"))  # 添加箱线图
p <- p + stat_compare_means(method = "wilcox.test", label = "p.signif", comparisons = comparisons)
p
ggsave("z49_OC-group2.pdf",plot = p,width = 3.2,height = 4)

p <- ggviolin(data2, x = "group1", y = "Migration_scores", fill = "group1",
              palette = c("#00AFBB", "#E76980"),  # 设置颜色
              ylab = "Average Expression", xlab = "Group",
              add = "boxplot", add.params = list(fill = "white")) +  # 添加箱线图
  facet_wrap(~ celltype, scales = "free",ncol=4)  # 按celltype分面，并允许每个分面有不同的y轴范围
p <- p + stat_compare_means(method = "wilcox.test", label = "p.signif", comparisons = comparisons)
p
ggsave("z49_OC-group2-celltype.pdf",plot = p,width = 10,height = 8)

geneSets2 = read_excel("genesets5.xls",col_names=F)
write.table(geneSets2,file = "my_genesets5.gmt",sep="\t",row.names = F,col.names = F,quote = F)
library(GSEABase)
geneSets=getGmt("my_genesets5.gmt",geneIdType=SymbolIdentifier())

library(rhdf5)
library(GSVA)

ssgsea_params <- ssgseaParam(exprData = expr, geneSets = geneSets, 
                             minSize = 1, maxSize = Inf, alpha = 0.25, normalize = TRUE)
gsvaResult <- gsva(ssgsea_params)

normalize=function(x){return((x-min(x))/(max(x)-min(x)))}
gsvaResult=normalize(gsvaResult)
gsvaOut=rbind(id=colnames(gsvaResult),gsvaResult)
write.table(gsvaOut, file="ssgseaOut5.txt",sep="\t",col.names = F,quote = F)
Monocytes <- AddMetaData(Monocytes, metadata = t(gsvaResult))

p <- FeaturePlot(Monocytes, split.by = "group2", features = rownames(gsvaResult), reduction = "umap", ncol=2)
p
ggsave("z73_Monocytes_Cytokine_group2.pdf",p,width = 6,height = 3)

p <- FeaturePlot(Monocytes, split.by = "group1", features = rownames(gsvaResult), reduction = "umap", ncol=2)
p
ggsave("z73_Monocytes8_Cytokine_group1.pdf",p,width = 10,height = 3)

data <- data.frame(group1 = Monocytes$group1,celltype = Monocytes$celltype,Cytokine_scores = Monocytes$`Cytokine scores`)
data$group1 <- factor(data$group1, levels = c("Health", "Mild", "Severe"))
p <- ggviolin(data, x = "group1", y = "Cytokine_scores", fill = "group1",
              palette = c("#00AFBB", "#E7B800", "#FC4E07"),  # 设置颜色
              ylab = "Average Expression", xlab = "Group",
              add = "boxplot", add.params = list(fill = "white"))  # 添加箱线图
p
p <- p + stat_compare_means(method = "wilcox.test", comparisons = my_comparisons, label = "p.signif")
p
ggsave("z73_Cytokine-group1-violin.pdf",plot = p,width = 3.5,height = 4.5)

p <- ggviolin(data, x = "group1", y = "Cytokine_scores", fill = "group1",
              palette = c("#00AFBB", "#E7B800", "#FC4E07"),  # 设置颜色
              ylab = "Average Expression", xlab = "Group",
              add = "boxplot", add.params = list(fill = "white")) +  # 添加箱线图
  facet_wrap(~ celltype, scales = "free",ncol=4)  # 按celltype分面，并允许每个分面有不同的y轴范围
p
# 假设你想在每个细胞类型中比较group1的分组差异
p <- p + stat_compare_means(method = "wilcox.test", label = "p.signif", comparisons = my_comparisons)
p
ggsave("z73_Cytokine-group1-violin-celltype.pdf",plot = p,width = 10,height = 8)

data2 <- data.frame(group1 = Monocytes$group2,celltype = Monocytes$celltype,Cytokine_scores = Monocytes$`Cytokine scores`)
data2$group1 <- factor(data2$group1, levels = c("Health", "Tuberculosis"))
p <- ggviolin(data2, x = "group1", y = "Cytokine_scores", fill = "group1",
              palette = c("#00AFBB", "#E76980"),  # 设置颜色
              ylab = "Average Expression", xlab = "Group",
              add = "boxplot", add.params = list(fill = "white"))  # 添加箱线图
p <- p + stat_compare_means(method = "wilcox.test", label = "p.signif", comparisons = comparisons)
p
ggsave("z73_Cytokine-group2.pdf",plot = p,width = 3.2,height = 4)

p <- ggviolin(data2, x = "group1", y = "Cytokine_scores", fill = "group1",
              palette = c("#00AFBB", "#E76980"),  # 设置颜色
              ylab = "Average Expression", xlab = "Group",
              add = "boxplot", add.params = list(fill = "white")) +  # 添加箱线图
  facet_wrap(~ celltype, scales = "free",ncol=4)  # 按celltype分面，并允许每个分面有不同的y轴范围
p <- p + stat_compare_means(method = "wilcox.test", label = "p.signif", comparisons = comparisons)
p
ggsave("z73_Cytokine-group2-celltype.pdf",plot = p,width = 10,height = 8)

library(ComplexHeatmap)
library(pheatmap)
library(circlize)

OC <- c("CTSK","NFATC1","FOS","MMP9","ACP5","ATP6V0D2","OCSTAMP","DCSTAMP","SIGLEC15","TNFRSF11A")
group1_means <- AggregateExpression(
  object = Monocytes,
  features = OC,  # 指定耗竭基因列表
  group.by = "group1",          # 按`group1`分组
  slot = "data"                 # 使用归一化后的数据
)
group1_means <- as.data.frame(group1_means$RNA)
df <- data.frame(colnames(group1_means))
colnames(df) <- 'class'
top_anno = HeatmapAnnotation(df = df,#细胞名/cluster
                             border = T,
                             show_annotation_name = F,
                             gp = gpar(col = 'black'),
                             col = list(class = c( 'Health'="#00AFBB",
                                                   'Mild'="#E7B800",
                                                   'Severe'="#FC4E07")))#颜色设置
marker_exp <- t(scale(t(group1_means),scale = T,center = F))
p=Heatmap(marker_exp,
          cluster_rows = F,
          cluster_columns = F,
          show_column_names = F,
          show_row_names = T,
          column_title = NULL,
          heatmap_legend_param = list(
            title=' '),
          col = colorRampPalette(c("white", "#F9D46C","#c85d4d"))(50),
          border = 'black',
          rect_gp = gpar(col = "black", lwd = 1),
          row_names_gp = gpar(fontsize = 10),
          column_names_gp = gpar(fontsize = 10),
          top_annotation = top_anno)
p
pdf("z49_OC-heatmap-group1.pdf", width = 3, height = 2.5)
print(p)
dev.off()

Idents(Monocytes) <- 'seurat_clusters'
Macrophages <- subset(Monocytes, idents = c("0", "2", "3", "6"))

Idents(Macrophages) <- 'orig.ident'

all.genes=rownames(Macrophages)

g2m_genes <- cc.genes$g2m.genes
g2m_genes <- CaseMatch(search = g2m_genes, match = rownames(Macrophages))

s_genes <- cc.genes$s.genes
s_genes <- CaseMatch(search = s_genes, match = rownames(Macrophages))

options(future.globals.maxSize = 60000 * 1024^3)
Macrophages <- NormalizeData(Macrophages) %>% FindVariableFeatures() %>% ScaleData(features = all.genes) %>% RunPCA()

Macrophages <- CellCycleScoring(Macrophages, g2m.features = g2m_genes, s.features = s_genes)

colnames(Macrophages@meta.data)
table(Macrophages$Phase)
DimPlot(Macrophages, group.by = "Phase")

Macrophages <- NormalizeData(Macrophages, normalization.method = "LogNormalize", scale.factor = 10000)
Macrophages <- FindVariableFeatures(Macrophages, selection.method = "vst", nfeatures = 2000)
Macrophages <- ScaleData(Macrophages, vars.to.regress = c("S.Score", "G2M.Score"))
DimPlot(Macrophages, group.by = "Phase")

top10 <- head(VariableFeatures(Macrophages), 10)
pdf(file = "y1.pdf", width = 7, height = 6)
VariableFeaturePlot(object = Macrophages)
dev.off()

pdf(file = "y2.pdf", width = 7, height = 6)
LabelPoints(plot= VariableFeaturePlot(object = Macrophages), points = top10, repel = TRUE)
dev.off()

Macrophages <- RunPCA(Macrophages, verbose = F)
pdf(file = "y3.pdf", width = 7, height = 6)
DimPlot(object = Macrophages, reduction = "pca")
dev.off()

pdf(file = "y4.pdf", width = 10, height = 9)
VizDimLoadings(object = Macrophages, dims = 1:4, reduction = "pca", nfeatures = 20)
dev.off()

pdf(file = "y5.pdf", width = 10, height = 9)
DimHeatmap(object = Macrophages, dims = 1:4, cells = 500, balanced = TRUE, nfeatures = 30, ncol = 2)
dev.off()

pdf(file = "y6.pdf", width = 10, height = 9)
ElbowPlot(Macrophages, ndims = 50)
dev.off()

pct <- Macrophages [["pca"]]@stdev / sum(Macrophages [["pca"]]@stdev) * 100
pct
cumu <- cumsum(pct)
cumu
co1 = which(cumu > 90 & pct < 5)[1]
co2 = sort(which((pct[1:length(pct)-1]-pct[2:length(pct)]) > 0.1), decreasing = T)[1]+1
pcs = min(co1, co2)
pcs
pcs = 1:14

Macrophages <- RunHarmony(Macrophages, group.by.vars="orig.ident", assay.use="SCT", max.iter.harmony = 20)
table(Macrophages@meta.data$orig.ident)

library(clustree)
seq = seq(0.1,1,by=0.1)
Macrophages <- FindNeighbors(Macrophages, dims = pcs)
for (res in seq){
  Macrophages = FindClusters(Macrophages, resolution = res)
}

library(patchwork)
p8 = clustree(Macrophages, prefix = "RNA_snn_res.")+coord_flip()
p = p8+plot_layout(widths = c(3,1))
p
ggsave("y7_T_RNA_snn_res-.pdf", p, width = 30, height = 14)

for (i in c(10,11,12,13,14,15,16,17,18,19,20)) {
  Macrophages <- FindNeighbors(Macrophages,reduction = "harmony", dims = 1:i) %>% FindClusters(resolution = 0.1)
  Macrophages <- RunUMAP(Macrophages, reduction = "harmony", dims = 1:i)
  plot_i <- print(DimPlot(Macrophages,reduction = "umap", label = TRUE, pt.size = 1, repel = TRUE, label.size = 5))
  plotname <- paste("plot_", i, sep = "")
  assign(plotname, plot_i)
  print(plot_i)
}
p = plot_10+plot_11+plot_12+plot_13+plot_14+plot_15+plot_16+plot_17+plot_18+plot_19+plot_20+plot_layout(ncol = 3)
ggsave(p,filename="dim.10-20.pdf",width = 15,height = 16)

for (i in seq) {
  Macrophages <- FindNeighbors(Macrophages,reduction = "harmony", dims = pcs) %>% FindClusters(resolution = i)
  Macrophages <- RunUMAP(Macrophages, reduction = "harmony", dims = pcs)
  plot_i <- print(DimPlot(Macrophages,reduction = "umap", label = TRUE, pt.size = 1, repel = TRUE, label.size = 5))
  plotname <- paste("plot_", i, sep = "")
  assign(plotname, plot_i)
  print(plot_i)
}
p = plot_0.1+plot_0.2+plot_0.3+plot_0.4+plot_0.5+plot_0.6+plot_0.7+plot_0.8+plot_0.9+plot_1+plot_layout(ncol = 4)
ggsave(p,filename="dim.16,res.0.1-1.pdf",width = 20,height = 12)


Macrophages <- FindNeighbors(Macrophages, reduction = "harmony", dims = pcs) %>% FindClusters(resolution = 0.1)
Macrophages <- RunUMAP(Macrophages, reduction = "harmony", dims = pcs) %>% RunTSNE(dims = pcs, reduction = "harmony")
DimPlot(Macrophages,reduction = "umap",cols = Biocols)

pdf(file = "y8.pdf", width = 7, height = 6)
DimPlot(Macrophages, reduction = "umap", label = T,cols = Biocols)
dev.off()

pdf(file = "y9.pdf", width = 7, height = 6)
DimPlot(Macrophages, reduction = "umap", label = F,group.by = "orig.ident",cols = Biocols)
dev.off()

pdf(file = "y10.pdf", width = 7, height = 6)
DimPlot(Macrophages, reduction = "tsne", label = T,cols = Biocols)
dev.off()

pdf(file = "y11.pdf", width = 7, height = 6)
DimPlot(Macrophages, reduction = "tsne", label = F,group.by = "orig.ident",cols = Biocols)
dev.off()

Macrophages.markers1 <- FindAllMarkers(Macrophages,only.pos = TRUE,logfc.threshold = 0.25)
write.csv(Macrophages.markers1,file="Macrophages.1.RNA.csv")

markers.OC <- c("CST3","CCL3","CCL5","IL1B","TNF","VCAN","CD274","TGFB1","CD68","C1QC","CCL3L1","DNAJB1","CD74","FCGR3A","CD14","LYZ","CTSK","NFATC1","FOS","MMP9","ACP5","ATP6V0D2","OCSTAMP","DCSTAMP","SIGLEC15","TNFRSF11A")

p11 <- FeaturePlot(Macrophages,features = markers.OC,ncol = 3)
p11
ggsave("y14-OC.pdf",p11,width = 12,height = 12)

p12 <- DotPlot(Macrophages,features = markers.OC) + RotatedAxis()
p12
ggsave("y15-OC.pdf",p12,width = 8,height = 4)

library(ComplexHeatmap)
library(pheatmap)
library(circlize)
group1_means <- AggregateExpression(
  object = Macrophages,
  features = markers.OC,
  group.by = "RNA_snn_res.0.1",
  slot = "data"
)
group1_means <- as.data.frame(group1_means$RNA)
df <- data.frame(colnames(group1_means))
colnames(df) <- 'class'
top_anno = HeatmapAnnotation(df = df,
                             border = T,
                             show_annotation_name = F,
                             gp = gpar(col = 'black'),
                             col = list(class = c('g0' = "#AB3282",
                                                  'g1' = "#53A85F",
                                                  'g2' = "#F1BB72",
                                                  'g3' = "#D6E7A3",
                                                  'g4' = "#57C3F3",
                                                  'g5' = "#E95C59",
                                                  'g6' = "#E59CC4")))#颜色设置
marker_exp <- t(scale(t(group1_means),scale = T,center = F))
p=Heatmap(marker_exp,
          cluster_rows = T,
          cluster_columns = T,
          show_column_names = F,
          show_row_names = T,
          column_title = NULL,
          heatmap_legend_param = list(
            title=' '),
          col = colorRampPalette(c("#ADD8E6", "#F0E68C", "#FF4500"))(50),
          border = 'black',
          rect_gp = gpar(col = "black", lwd = 1),
          row_names_gp = gpar(fontsize = 10),
          column_names_gp = gpar(fontsize = 10),
          top_annotation = top_anno)
p
pdf("y12-OC.pdf", width = 5, height = 4)
print(p)
dev.off()

p13 <- VlnPlot(Macrophages,features = markers.OC,stack = T,flip = T) + NoLegend()
p13
ggsave("y16-OC.pdf",p13,width = 8,height = 12)

Macrophages$celltype <- recode(Macrophages@meta.data$seurat_clusters,
                               "0" = "Macro_C1QC",
                               "1" = "pOC_CCL3",
                               "2" = "pOC",
                               "3" = "Macro_TNF",
                               "4" = "Macro_VCAN",
                               "5" = "Macro_DNAJB1",
                               "6" = "mOC")
table(Macrophages@meta.data$seurat_clusters)

Idents(Macrophages) <- 'celltype'
table(Macrophages@meta.data$celltype)

p=DimPlot(Macrophages, split.by = "group1", label = TRUE,cols = Biocols)
p
ggsave("A13.pdf",p,width = 10,height = 4)

p=DimPlot(Macrophages, split.by = "group1", label = F,cols = Biocols)
p
ggsave("A13'.pdf",p,width = 10,height = 4)

p=DimPlot(Macrophages, split.by = "group2", label = T,cols = Biocols)
p
ggsave("A14.pdf",p,width = 8,height = 4)

p=DimPlot(Macrophages, split.by = "group2", label = F,cols = Biocols)
p
ggsave("A14'.pdf",p,width = 8,height = 4)

pdf(file = "A7.pdf", width = 6, height = 4)
DimPlot(Macrophages, reduction = "umap", label = T,cols = Biocols)
dev.off()

pdf(file = "A7'.pdf", width = 6, height = 4)
DimPlot(Macrophages, reduction = "umap", label = F,cols = Biocols)
dev.off()

pdf(file = "A8.pdf", width = 5, height = 4)
DimPlot(Macrophages, reduction = "umap", label = F,group.by = "orig.ident",cols = Biocols)
dev.off()

pdf(file = "A9.pdf", width = 5, height = 4)
DimPlot(Macrophages, reduction = "umap", label = F,group.by = "group1",cols = Biocols)
dev.off()

pdf(file = "A9'.pdf", width = 5.5, height = 4)
DimPlot(Macrophages, reduction = "umap", label = F,group.by = "group2",cols = Biocols)
dev.off()

pdf(file = "A10.pdf", width = 6, height = 4)
DimPlot(Macrophages, reduction = "tsne", label = T,cols = Biocols)
dev.off()

pdf(file = "A10'.pdf", width = 6, height = 4)
DimPlot(Macrophages, reduction = "tsne", label = F,cols = Biocols)
dev.off()

pdf(file = "A11.pdf", width = 5, height = 4)
DimPlot(Macrophages, reduction = "tsne", label = F,group.by = "orig.ident",cols = Biocols)
dev.off()

pdf(file = "A11.pdf", width = 5, height = 4)
DimPlot(Macrophages, reduction = "tsne", label = F,group.by = "group1",cols = Biocols)
dev.off()

pdf(file = "A11'.pdf", width = 5.5, height = 4)
DimPlot(Macrophages, reduction = "tsne", label = F,group.by = "group2",cols = Biocols)
dev.off()

dim(Macrophages)
DimPlot(Macrophages,group.by = "celltype",label = T)
Macrophages$celltype <- Macrophages$celltype
colnames(Macrophages@meta.data)

data <- GetAssayData(Macrophages,assay = "RNA",slot = "counts")

pd <- new('AnnotatedDataFrame',data = Macrophages@meta.data[,c(1,11,22,24)])

fData <- data.frame(gene_short_name = row.names(data),row.names = row.names(data))
fd <- new("AnnotatedDataFrame", data = fData)
library(monocle)
dim(Macrophages)
mycds <- newCellDataSet(as.matrix(data),phenoData = pd,featureData = fd,lowerDetectionLimit = 0.5,expressionFamily = negbinomial.size())

mycds <- estimateSizeFactors(mycds)
mycds <- estimateDispersions(mycds,cores=6)

disp_table <- dispersionTable(mycds)
order.genes <- subset(disp_table, mean_expression >= 0.005 & dispersion_empirical >= 1 * dispersion_fit) %>% pull(gene_id) %>%as.character()


mycds <- setOrderingFilter(mycds,order.genes)

plot_ordering_genes(mycds)
p <- plot_ordering_genes(mycds)
ggsave("z0_OrderGenes.pdf",p,width = 8,height = 6)

mycds <- reduceDimension(mycds,max_components = 2,reduction_method = 'DDRTree',residualModelFormulaStr = "~orig.ident")
mycds <- orderCells(mycds, root_state = 5)   

p1 <- plot_cell_trajectory(mycds,color_by = "State")+
  scale_color_manual(values = Biocols)
p1
ggsave("z1_Trajectory_State.pdf",plot = p1,width = 5,height = 4)

p2 <- plot_cell_trajectory(mycds,color_by = "Pseudotime")
p2
ggsave("z2_Trajectory_Pseudotime.pdf",plot = p2,width = 5,height = 4)

p3 <- plot_cell_trajectory(mycds,color_by = "celltype")+
  scale_color_manual(values = Biocols)
p3
ggsave("z3_Trajectory_celltype.pdf",plot = p3,width = 5,height = 4)

p4 <- plot_cell_trajectory(mycds, color_by = "celltype") +
  facet_wrap(~group1, nrow = 1) +
  scale_color_manual(values = Biocols)
p4
ggsave("z4_Trajectory_celltype_group1.pdf",plot = p4,width = 10,height = 4)

p4 <- plot_cell_trajectory(mycds, color_by = "group1") +
  scale_color_manual(values = Biocols)
p4
ggsave("z5_Trajectory_group1.pdf",plot = p4,width = 5,height = 4)

p4 <- plot_cell_trajectory(mycds, color_by = "group1") +
  facet_wrap(~group1, nrow = 1) +
  scale_color_manual(values = Biocols)
p4
ggsave("z6_Trajectory_group1_group1.pdf",plot = p4,width = 10,height = 4)

p5 <- plot_complex_cell_trajectory(mycds,x = 1,y = 2,
                                   color_by = "celltype")+
  scale_color_manual(values = Biocols)
p5
ggsave("z7_Trajectory_dendrogram.pdf",plot = p5,width = 6,height = 4)

p6 <- ggplot(pData(mycds),aes(Pseudotime,colour = celltype,fill = celltype))+
  geom_density(bw = 0.5, size = 1,alpha = 0.5)+theme_classic()
p6
ggsave("z8_Trajectory_density.pdf",plot = p6,width = 6,height = 4)

genes = c("CD14","CTSK","FOS","MMP9","ACP5")
p1 = plot_genes_in_pseudotime(mycds[genes],color_by = "State")
p2 = plot_genes_in_pseudotime(mycds[genes],color_by = "celltype")
p3 = plot_genes_in_pseudotime(mycds[genes],color_by = "Pseudotime")
ggsave("z9_Trajectory_pseudotime.pdf",plot = p1|p2|p3,width = 10,height = 5)

p1 = plot_genes_jitter(mycds[genes],grouping="State",color_by = "State")
p2 = plot_genes_violin(mycds[genes],grouping="State",color_by = "State")
p3 = plot_genes_in_pseudotime(mycds[genes],color_by = "State")
ggsave("z10_Trajectory_jitter.pdf",plot = p1|p2|p3,width = 10,height = 5)

p1 = plot_genes_jitter(mycds[genes],grouping="celltype",color_by = "celltype")
p2 = plot_genes_violin(mycds[genes],grouping="celltype",color_by = "celltype")
p3 = plot_genes_in_pseudotime(mycds[genes],color_by = "celltype")
ggsave("z11_Trajectory_jitter.pdf",plot = p1|p2|p3,width = 10,height = 5)

p=plot_cells(cds,genes="CD14",
             cell_size=1,
             show_trajectory_graph=FALSE,
             label_cell_groups=FALSE,
             label_leaves=FALSE)
p
ggsave("z16_trajectory-umap-CD14.pdf",p,width = 5,height = 4)

p=plot_cells(cds,genes="ACP5",
             cell_size=1,
             show_trajectory_graph=FALSE,
             label_cell_groups=FALSE,
             label_leaves=FALSE)
p
ggsave("z16_trajectory-umap-ACP5.pdf",p,width = 5,height = 4)

p=plot_cells(cds,genes="CTSK",
             cell_size=1,
             show_trajectory_graph=FALSE,
             label_cell_groups=FALSE,
             label_leaves=FALSE)
p
ggsave("z16_trajectory-umap-CTSK.pdf",p,width = 5,height = 4)

p=plot_cells(cds,genes="MMP9",
             cell_size=1,
             show_trajectory_graph=FALSE,
             label_cell_groups=FALSE,
             label_leaves=FALSE)
p
ggsave("z16_trajectory-umap-MMP9.pdf",p,width = 5,height = 4)

table(sce.all$celltype)

DefaultAssay(Macrophages) <- "RNA"
options(future.globals.maxSize = 20 * 1024^3)
Macrophages <- NormalizeData(Macrophages)
dim(Macrophages)

expr = as.matrix(GetAssayData(object = Macrophages@assays$RNA, layer="data"))
dim(expr)
expr <- expr[rowSums(expr)>0,]
dim(expr)
library(readxl)
geneSets2 = read_excel("genesets2.xls",col_names=F)
write.table(geneSets2,file = "my_genesets2.gmt",sep="\t",row.names = F,col.names = F,quote = F)
library(GSEABase)
geneSets=getGmt("my_genesets2.gmt",geneIdType=SymbolIdentifier())

library(rhdf5)
library(GSVA)

ssgsea_params <- ssgseaParam(exprData = expr, geneSets = geneSets, 
                             minSize = 1, maxSize = Inf, alpha = 0.25, normalize = TRUE)
gsvaResult <- gsva(ssgsea_params)

normalize=function(x){return((x-min(x))/(max(x)-min(x)))}
gsvaResult=normalize(gsvaResult)
gsvaOut=rbind(id=colnames(gsvaResult),gsvaResult)
write.table(gsvaOut, file="ssgseaOut2.txt",sep="\t",col.names = F,quote = F)
Macrophages <- AddMetaData(Macrophages, metadata = t(gsvaResult))

p <- FeaturePlot(Macrophages, split.by = "group1", features = rownames(gsvaResult), reduction = "umap", ncol=2)
p
ggsave("z49_Macrophages_OC_group1.pdf",p,width = 10,height = 3)

p <- FeaturePlot(Macrophages, split.by = "group2", features = rownames(gsvaResult), reduction = "umap", ncol=2)
p
ggsave("z49_Macrophages_OC_group2.pdf",p,width = 6,height = 3)

data <- data.frame(group1 = Macrophages$group1,celltype = Macrophages$celltype,Migration_scores = Macrophages$`OC scores`)
data$group1 <- factor(data$group1, levels = c("Health", "Mild", "Severe"))
p <- ggviolin(data, x = "group1", y = "Migration_scores", fill = "group1",
              palette = c("#00AFBB", "#E7B800", "#FC4E07"),  # 设置颜色
              ylab = "Average Expression", xlab = "Group",
              add = "boxplot", add.params = list(fill = "white"))  # 添加箱线图
p
p <- p + stat_compare_means(method = "wilcox.test", comparisons = my_comparisons, label = "p.signif")
p
ggsave("z49_OC-group1-violin.pdf",plot = p,width = 3.5,height = 4.5)

p <- ggviolin(data, x = "group1", y = "Migration_scores", fill = "group1",
              palette = c("#00AFBB", "#E7B800", "#FC4E07"),  # 设置颜色
              ylab = "Average Expression", xlab = "Group",
              add = "boxplot", add.params = list(fill = "white")) +  # 添加箱线图
  facet_wrap(~ celltype, scales = "free",ncol=4)  # 按celltype分面，并允许每个分面有不同的y轴范围
p
# 假设你想在每个细胞类型中比较group1的分组差异
p <- p + stat_compare_means(method = "wilcox.test", label = "p.signif", comparisons = my_comparisons)
p
ggsave("z49_OC-group1-violin-celltype.pdf",plot = p,width = 10,height = 8)

data2 <- data.frame(group1 = Macrophages$group2,celltype = Macrophages$celltype,Migration_scores = Macrophages$`OC scores`)
data2$group1 <- factor(data2$group1, levels = c("Health", "Tuberculosis"))
p <- ggviolin(data2, x = "group1", y = "Migration_scores", fill = "group1",
              palette = c("#00AFBB", "#E76980"),  # 设置颜色
              ylab = "Average Expression", xlab = "Group",
              add = "boxplot", add.params = list(fill = "white"))  # 添加箱线图
p <- p + stat_compare_means(method = "wilcox.test", label = "p.signif", comparisons = comparisons)
p
ggsave("z49_OC-group2.pdf",plot = p,width = 3.2,height = 4)

p <- ggviolin(data2, x = "group1", y = "Migration_scores", fill = "group1",
              palette = c("#00AFBB", "#E76980"),  # 设置颜色
              ylab = "Average Expression", xlab = "Group",
              add = "boxplot", add.params = list(fill = "white")) +  # 添加箱线图
  facet_wrap(~ celltype, scales = "free",ncol=4)  # 按celltype分面，并允许每个分面有不同的y轴范围
p <- p + stat_compare_means(method = "wilcox.test", label = "p.signif", comparisons = comparisons)
p
ggsave("z49_OC-group2-celltype.pdf",plot = p,width = 10,height = 8)

data3 <- data.frame(group1 = Macrophages$celltype,celltype = Macrophages$celltype,Migration_scores = Macrophages$`OC scores`)
data3$group1 <- factor(data3$group1, levels = c("Macro_C1QC","Macro_DNAJB1", "Macro_VCAN","Macro_TNF","pOC","pOC_CCL3","mOC"))
p <- ggviolin(data3, x = "group1", y = "Migration_scores", fill = "group1",
              palette = c("#AB3282", "#53A85F","#F1BB72","#D6E7A3","#57C3F3","#E95C59","#E59CC4"),  # 设置颜色
              ylab = "Average Expression", xlab = "Group",
              add = "boxplot", add.params = list(fill = "white"))  # 添加箱线图
p
ggsave("z49_OC-celltype.pdf",plot = p,width = 6,height = 4)

geneSets2 = read_excel("genesets3.xls",col_names=F)
write.table(geneSets2,file = "my_genesets3.gmt",sep="\t",row.names = F,col.names = F,quote = F)
geneSets=getGmt("my_genesets3.gmt",geneIdType=SymbolIdentifier())


ssgsea_params <- ssgseaParam(exprData = expr, geneSets = geneSets, 
                             minSize = 1, maxSize = Inf, alpha = 0.25, normalize = TRUE)
gsvaResult <- gsva(ssgsea_params)

normalize=function(x){return((x-min(x))/(max(x)-min(x)))}
gsvaResult=normalize(gsvaResult)
gsvaOut=rbind(id=colnames(gsvaResult),gsvaResult)
write.table(gsvaOut, file="ssgseaOut3.txt",sep="\t",col.names = F,quote = F)
Macrophages <- AddMetaData(Macrophages, metadata = t(gsvaResult))

p <- FeaturePlot(Macrophages, split.by = "group1", features = rownames(gsvaResult), reduction = "umap", ncol=2)
p
ggsave("z57_Macrophages_Apoptosis_group1.pdf",p,width = 10,height = 3)

p <- FeaturePlot(Macrophages, split.by = "group2", features = rownames(gsvaResult), reduction = "umap", ncol=2)
p
ggsave("z57_Macrophages_Apoptosis_group2.pdf",p,width = 6,height = 3)

data <- data.frame(group1 = Macrophages$group1,celltype = Macrophages$celltype,Apoptosis_scores = Macrophages$`Apoptosis scores`)
data$group1 <- factor(data$group1, levels = c("Health", "Mild", "Severe"))
p <- ggviolin(data, x = "group1", y = "Apoptosis_scores", fill = "group1",
              palette = c("#00AFBB", "#E7B800", "#FC4E07"),  # 设置颜色
              ylab = "Average Expression", xlab = "Group",
              add = "boxplot", add.params = list(fill = "white"))  # 添加箱线图
p
p <- p + stat_compare_means(method = "wilcox.test", comparisons = my_comparisons, label = "p.signif")
p
ggsave("z57_Apoptosis-group1-violin.pdf",plot = p,width = 3.5,height = 4.5)

p <- ggviolin(data, x = "group1", y = "Apoptosis_scores", fill = "group1",
              palette = c("#00AFBB", "#E7B800", "#FC4E07"),  # 设置颜色
              ylab = "Average Expression", xlab = "Group",
              add = "boxplot", add.params = list(fill = "white")) +  # 添加箱线图
  facet_wrap(~ celltype, scales = "free",ncol=4)  # 按celltype分面，并允许每个分面有不同的y轴范围
p
# 假设你想在每个细胞类型中比较group1的分组差异
p <- p + stat_compare_means(method = "wilcox.test", label = "p.signif", comparisons = my_comparisons)
p
ggsave("z57_Apoptosis-group1-violin-celltype.pdf",plot = p,width = 10,height = 8)

data2 <- data.frame(group1 = Macrophages$group2,celltype = Macrophages$celltype,Apoptosis_scores = Macrophages$`Apoptosis scores`)
data2$group1 <- factor(data2$group1, levels = c("Health", "Tuberculosis"))
p <- ggviolin(data2, x = "group1", y = "Apoptosis_scores", fill = "group1",
              palette = c("#00AFBB", "#E76980"),  # 设置颜色
              ylab = "Average Expression", xlab = "Group",
              add = "boxplot", add.params = list(fill = "white"))  # 添加箱线图
p <- p + stat_compare_means(method = "wilcox.test", label = "p.signif", comparisons = comparisons)
p
ggsave("z57_Apoptosis-group2.pdf",plot = p,width = 3.2,height = 4)

p <- ggviolin(data2, x = "group1", y = "Apoptosis_scores", fill = "group1",
              palette = c("#00AFBB", "#E76980"),  # 设置颜色
              ylab = "Average Expression", xlab = "Group",
              add = "boxplot", add.params = list(fill = "white")) +  # 添加箱线图
  facet_wrap(~ celltype, scales = "free",ncol=4)  # 按celltype分面，并允许每个分面有不同的y轴范围
p <- p + stat_compare_means(method = "wilcox.test", label = "p.signif", comparisons = comparisons)
p
ggsave("z57_Apoptosis-group2-celltype.pdf",plot = p,width = 10,height = 8)

data3 <- data.frame(group1 = Macrophages$celltype,celltype = Macrophages$celltype,Migration_scores = Macrophages$`Apoptosis scores`)
data3$group1 <- factor(data3$group1, levels = c("Macro_C1QC","Macro_DNAJB1", "Macro_VCAN","Macro_TNF","pOC","pOC_CCL3","mOC"))
p <- ggviolin(data3, x = "group1", y = "Migration_scores", fill = "group1",
              palette = c("#AB3282", "#53A85F","#F1BB72","#D6E7A3","#57C3F3","#E95C59","#E59CC4"),  # 设置颜色
              ylab = "Average Expression", xlab = "Group",
              add = "boxplot", add.params = list(fill = "white"))  # 添加箱线图
p
ggsave("z57_Apoptosis-celltype.pdf",plot = p,width = 6,height = 4)

geneSets2 = read_excel("genesets5.xls",col_names=F)
write.table(geneSets2,file = "my_genesets5.gmt",sep="\t",row.names = F,col.names = F,quote = F)
geneSets=getGmt("my_genesets5.gmt",geneIdType=SymbolIdentifier())

ssgsea_params <- ssgseaParam(exprData = expr, geneSets = geneSets, 
                             minSize = 1, maxSize = Inf, alpha = 0.25, normalize = TRUE)
gsvaResult <- gsva(ssgsea_params)

#gsvaResult=gsva(expr,geneSets,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)

normalize=function(x){return((x-min(x))/(max(x)-min(x)))}
gsvaResult=normalize(gsvaResult)
gsvaOut=rbind(id=colnames(gsvaResult),gsvaResult)
write.table(gsvaOut, file="ssgseaOut5.txt",sep="\t",col.names = F,quote = F)
Macrophages <- AddMetaData(Macrophages, metadata = t(gsvaResult))

p <- FeaturePlot(Macrophages, split.by = "group2", features = rownames(gsvaResult), reduction = "umap", ncol=2)
p
ggsave("z73_Macrophages_Cytokine_group2.pdf",p,width = 6,height = 3)

p <- FeaturePlot(Macrophages, split.by = "group1", features = rownames(gsvaResult), reduction = "umap", ncol=2)
p
ggsave("z73_Macrophages_Cytokine_group1.pdf",p,width = 10,height = 3)

data <- data.frame(group1 = Macrophages$group1,celltype = Macrophages$celltype,Cytokine_scores = Macrophages$`Cytokine scores`)
data$group1 <- factor(data$group1, levels = c("Health", "Mild", "Severe"))
p <- ggviolin(data, x = "group1", y = "Cytokine_scores", fill = "group1",
              palette = c("#00AFBB", "#E7B800", "#FC4E07"),  # 设置颜色
              ylab = "Average Expression", xlab = "Group",
              add = "boxplot", add.params = list(fill = "white"))  # 添加箱线图
p
p <- p + stat_compare_means(method = "wilcox.test", comparisons = my_comparisons, label = "p.signif")
p
ggsave("z73_Cytokine-group1-violin.pdf",plot = p,width = 3.5,height = 4.5)

p <- ggviolin(data, x = "group1", y = "Cytokine_scores", fill = "group1",
              palette = c("#00AFBB", "#E7B800", "#FC4E07"),  # 设置颜色
              ylab = "Average Expression", xlab = "Group",
              add = "boxplot", add.params = list(fill = "white")) +  # 添加箱线图
  facet_wrap(~ celltype, scales = "free",ncol=4)  # 按celltype分面，并允许每个分面有不同的y轴范围
p
# 假设你想在每个细胞类型中比较group1的分组差异
p <- p + stat_compare_means(method = "wilcox.test", label = "p.signif", comparisons = my_comparisons)
p
ggsave("z73_Cytokine-group1-violin-celltype.pdf",plot = p,width = 10,height = 8)

data2 <- data.frame(group1 = Macrophages$group2,celltype = Macrophages$celltype,Cytokine_scores = Macrophages$`Cytokine scores`)
data2$group1 <- factor(data2$group1, levels = c("Health", "Tuberculosis"))
p <- ggviolin(data2, x = "group1", y = "Cytokine_scores", fill = "group1",
              palette = c("#00AFBB", "#E76980"),  # 设置颜色
              ylab = "Average Expression", xlab = "Group",
              add = "boxplot", add.params = list(fill = "white"))  # 添加箱线图
p <- p + stat_compare_means(method = "wilcox.test", label = "p.signif", comparisons = comparisons)
p
ggsave("z73_Cytokine-group2.pdf",plot = p,width = 3.2,height = 4)

p <- ggviolin(data2, x = "group1", y = "Cytokine_scores", fill = "group1",
              palette = c("#00AFBB", "#E76980"),  # 设置颜色
              ylab = "Average Expression", xlab = "Group",
              add = "boxplot", add.params = list(fill = "white")) +  # 添加箱线图
  facet_wrap(~ celltype, scales = "free",ncol=4)  # 按celltype分面，并允许每个分面有不同的y轴范围
p <- p + stat_compare_means(method = "wilcox.test", label = "p.signif", comparisons = comparisons)
p
ggsave("z73_Cytokine-group2-celltype.pdf",plot = p,width = 10,height = 8)

data3 <- data.frame(group1 = Macrophages$celltype,celltype = Macrophages$celltype,Migration_scores = Macrophages$`Cytokine scores`)
data3$group1 <- factor(data3$group1, levels = c("Macro_C1QC","Macro_DNAJB1", "Macro_VCAN","Macro_TNF","pOC","pOC_CCL3","mOC"))
p <- ggviolin(data3, x = "group1", y = "Migration_scores", fill = "group1",
              palette = c("#AB3282", "#53A85F","#F1BB72","#D6E7A3","#57C3F3","#E95C59","#E59CC4"),  # 设置颜色
              ylab = "Average Expression", xlab = "Group",
              add = "boxplot", add.params = list(fill = "white"))  # 添加箱线图
p
ggsave("z73_Cytokine-celltype.pdf",plot = p,width = 6,height = 4)