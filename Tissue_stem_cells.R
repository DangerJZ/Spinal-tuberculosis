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

setwd("D:/1/Tissue_stem_cells")

Idents(sce.all)

Tissue_stem_cells <- subset(sce.all, idents = "Tissue_stem_cells")

Idents(Tissue_stem_cells) <- 'orig.ident'

all.genes=rownames(Tissue_stem_cells)

g2m_genes <- cc.genes$g2m.genes
g2m_genes <- CaseMatch(search = g2m_genes, match = rownames(Tissue_stem_cells))

s_genes <- cc.genes$s.genes
s_genes <- CaseMatch(search = s_genes, match = rownames(Tissue_stem_cells))

options(future.globals.maxSize = 60000 * 1024^3)
Tissue_stem_cells <- NormalizeData(Tissue_stem_cells) %>% FindVariableFeatures() %>% ScaleData(features = all.genes) %>% RunPCA()

Tissue_stem_cells <- CellCycleScoring(Tissue_stem_cells, g2m.features = g2m_genes, s.features = s_genes)

colnames(Tissue_stem_cells@meta.data)
table(Tissue_stem_cells$Phase)
DimPlot(Tissue_stem_cells, group.by = "Phase")

Tissue_stem_cells <- NormalizeData(Tissue_stem_cells, normalization.method = "LogNormalize", scale.factor = 10000)
Tissue_stem_cells <- FindVariableFeatures(Tissue_stem_cells, selection.method = "vst", nfeatures = 2000)
Tissue_stem_cells <- ScaleData(Tissue_stem_cells, vars.to.regress = c("S.Score", "G2M.Score"))
DimPlot(Tissue_stem_cells, group.by = "Phase")

top10 <- head(VariableFeatures(Tissue_stem_cells), 10)
pdf(file = "y1.pdf", width = 7, height = 6)
VariableFeaturePlot(object = Tissue_stem_cells)
dev.off()

pdf(file = "y2.pdf", width = 7, height = 6)
LabelPoints(plot= VariableFeaturePlot(object = Tissue_stem_cells), points = top10, repel = TRUE)
dev.off()

Tissue_stem_cells <- RunPCA(Tissue_stem_cells, verbose = F)
pdf(file = "y3.pdf", width = 7, height = 6)
DimPlot(object = Tissue_stem_cells, reduction = "pca")
dev.off()

pdf(file = "y4.pdf", width = 10, height = 9)
VizDimLoadings(object = Tissue_stem_cells, dims = 1:4, reduction = "pca", nfeatures = 20)
dev.off()

pdf(file = "y5.pdf", width = 10, height = 9)
DimHeatmap(object = Tissue_stem_cells, dims = 1:4, cells = 500, balanced = TRUE, nfeatures = 30, ncol = 2)
dev.off()

pdf(file = "y6.pdf", width = 10, height = 9)
ElbowPlot(Tissue_stem_cells, ndims = 50)
dev.off()

pct <- Tissue_stem_cells [["pca"]]@stdev / sum(Tissue_stem_cells [["pca"]]@stdev) * 100
pct
cumu <- cumsum(pct)
cumu
co1 = which(cumu > 90 & pct < 5)[1]
co2 = sort(which((pct[1:length(pct)-1]-pct[2:length(pct)]) > 0.1), decreasing = T)[1]+1
pcs = min(co1, co2)
pcs
pcs = 1:18

Tissue_stem_cells <- RunHarmony(Tissue_stem_cells, group.by.vars="orig.ident", assay.use="SCT", max.iter.harmony = 20)
table(Tissue_stem_cells@meta.data$orig.ident)

library(clustree)
seq = seq(0.1,1,by=0.1)
Tissue_stem_cells <- FindNeighbors(Tissue_stem_cells, dims = pcs)
for (res in seq){
  Tissue_stem_cells = FindClusters(Tissue_stem_cells, resolution = res)
}

library(patchwork)
p8 = clustree(Tissue_stem_cells, prefix = "RNA_snn_res.")+coord_flip()
p = p8+plot_layout(widths = c(3,1))
p
ggsave("y7_T_RNA_snn_res-.pdf", p, width = 30, height = 14)

for (i in c(10,11,12,13,14,15,16,17,18,19,20)) {
  Tissue_stem_cells <- FindNeighbors(Tissue_stem_cells,reduction = "harmony", dims = 1:i) %>% FindClusters(resolution = 0.1)
  Tissue_stem_cells <- RunUMAP(Tissue_stem_cells, reduction = "harmony", dims = 1:i)
  plot_i <- print(DimPlot(Tissue_stem_cells,reduction = "umap", label = TRUE, pt.size = 1, repel = TRUE, label.size = 5))
  plotname <- paste("plot_", i, sep = "")
  assign(plotname, plot_i)
  print(plot_i)
}
p = plot_10+plot_11+plot_12+plot_13+plot_14+plot_15+plot_16+plot_17+plot_18+plot_19+plot_20+plot_layout(ncol = 3)
ggsave(p,filename="dim.10-20.pdf",width = 15,height = 16)

for (i in seq) {
  Tissue_stem_cells <- FindNeighbors(Tissue_stem_cells,reduction = "harmony", dims = pcs) %>% FindClusters(resolution = i)
  Tissue_stem_cells <- RunUMAP(Tissue_stem_cells, reduction = "harmony", dims = pcs)
  plot_i <- print(DimPlot(Tissue_stem_cells,reduction = "umap", label = TRUE, pt.size = 1, repel = TRUE, label.size = 5))
  plotname <- paste("plot_", i, sep = "")
  assign(plotname, plot_i)
  print(plot_i)
}
p = plot_0.1+plot_0.2+plot_0.3+plot_0.4+plot_0.5+plot_0.6+plot_0.7+plot_0.8+plot_0.9+plot_1+plot_layout(ncol = 4)
ggsave(p,filename="dim.16,res.0.1-1.pdf",width = 20,height = 12)


Tissue_stem_cells <- FindNeighbors(Tissue_stem_cells, reduction = "harmony", dims = pcs) %>% FindClusters(resolution = 0.1)
Tissue_stem_cells <- RunUMAP(Tissue_stem_cells, reduction = "harmony", dims = pcs) %>% RunTSNE(dims = pcs, reduction = "harmony")
DimPlot(Tissue_stem_cells,reduction = "umap",cols = Biocols)

pdf(file = "y8.pdf", width = 7, height = 6)
DimPlot(Tissue_stem_cells, reduction = "umap", label = T,cols = Biocols)
dev.off()

pdf(file = "y9.pdf", width = 7, height = 6)
DimPlot(Tissue_stem_cells, reduction = "umap", label = F,group.by = "orig.ident",cols = Biocols)
dev.off()

pdf(file = "y10.pdf", width = 7, height = 6)
DimPlot(Tissue_stem_cells, reduction = "tsne", label = T,cols = Biocols)
dev.off()

pdf(file = "y11.pdf", width = 7, height = 6)
DimPlot(Tissue_stem_cells, reduction = "tsne", label = F,group.by = "orig.ident",cols = Biocols)
dev.off()

Tissue_stem_cells.markers1 <- FindAllMarkers(Tissue_stem_cells,only.pos = TRUE,logfc.threshold = 0.25)
write.csv(Tissue_stem_cells.markers1,file="Tissue_stem_cells.1.RNA.csv")

markers <- c("ENG","NT5E","THY1","PPARG","FBN1","LUM","COL1A2","COL1A1","COL3A1","RUNX2","SP7","BGLAP","IBSP","SPP1","PDPN","HAS1","DPT","CD164")

p11 <- FeaturePlot(Tissue_stem_cells,features = markers,ncol = 3)
p11
ggsave("y14.pdf",p11,width = 12,height = 12)

p12 <- DotPlot(Tissue_stem_cells,features = markers) + RotatedAxis()
p12
ggsave("y15.pdf",p12,width = 8,height = 4)

library(ComplexHeatmap)
library(pheatmap)
library(circlize)
group1_means <- AggregateExpression(
  object = Tissue_stem_cells,
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

p13 <- VlnPlot(Tissue_stem_cells,features = markers,stack = T,flip = T) + NoLegend()
p13
ggsave("y16.pdf",p13,width = 8,height = 12)

Tissue_stem_cells$celltype <- recode(Tissue_stem_cells@meta.data$seurat_clusters,
                                     "0" = "MSC_fibro",
                                     "1" = "MSC_adipo",
                                     "2" = "MSC_osteo",
                                     "3" = "MSC_osteo",
                                     "4" = "Osteoblasts",
                                     "5" = "MSC_osteo",
                                     "6" = "MSC_osteo",
                                     "7" = "MSC_fibro")
table(Tissue_stem_cells@meta.data$seurat_clusters)

Idents(Tissue_stem_cells) <- 'celltype'
table(Tissue_stem_cells@meta.data$celltype)

p=DimPlot(Tissue_stem_cells, split.by = "group1", label = TRUE,cols = Biocols)
p
ggsave("A13.pdf",p,width = 10,height = 4)

p=DimPlot(Tissue_stem_cells, split.by = "group1", label = F,cols = Biocols)
p
ggsave("A13'.pdf",p,width = 10,height = 4)

p=DimPlot(Tissue_stem_cells, split.by = "group2", label = T,cols = Biocols)
p
ggsave("A14.pdf",p,width = 8,height = 4)

p=DimPlot(Tissue_stem_cells, split.by = "group2", label = F,cols = Biocols)
p
ggsave("A14'.pdf",p,width = 8,height = 4)

pdf(file = "A7.pdf", width = 6, height = 4)
DimPlot(Tissue_stem_cells, reduction = "umap", label = T,cols = Biocols)
dev.off()

pdf(file = "A7'.pdf", width = 6, height = 4)
DimPlot(Tissue_stem_cells, reduction = "umap", label = F,cols = Biocols)
dev.off()

pdf(file = "A8.pdf", width = 5, height = 4)
DimPlot(Tissue_stem_cells, reduction = "umap", label = F,group.by = "orig.ident",cols = Biocols)
dev.off()

pdf(file = "A9.pdf", width = 5, height = 4)
DimPlot(Tissue_stem_cells, reduction = "umap", label = F,group.by = "group1",cols = Biocols)
dev.off()

pdf(file = "A9'.pdf", width = 5.5, height = 4)
DimPlot(Tissue_stem_cells, reduction = "umap", label = F,group.by = "group2",cols = Biocols)
dev.off()

pdf(file = "A10.pdf", width = 6, height = 4)
DimPlot(Tissue_stem_cells, reduction = "tsne", label = T,cols = Biocols)
dev.off()

pdf(file = "A10'.pdf", width = 6, height = 4)
DimPlot(Tissue_stem_cells, reduction = "tsne", label = F,cols = Biocols)
dev.off()

pdf(file = "A11.pdf", width = 5, height = 4)
DimPlot(Tissue_stem_cells, reduction = "tsne", label = F,group.by = "orig.ident",cols = Biocols)
dev.off()

pdf(file = "A12.pdf", width = 5, height = 4)
DimPlot(Tissue_stem_cells, reduction = "tsne", label = F,group.by = "group1",cols = Biocols)
dev.off()

pdf(file = "A12'.pdf", width = 5.5, height = 4)
DimPlot(Tissue_stem_cells, reduction = "tsne", label = F,group.by = "group2",cols = Biocols)
dev.off()

DefaultAssay(Tissue_stem_cells) <- "RNA"

bulk <- AggregateExpression(Tissue_stem_cells, return.seurat = T, slot = "counts", assays = "RNA",                             
                            group.by = c("celltype","orig.ident", "group2")) 

colnames(bulk)

library(DESeq2)
Idents(bulk) <- "group2"
de_markers <- FindMarkers(bulk, ident.1 ="Tuberculosis", ident.2 = "Health",
                          slot = "counts", test.use = "DESeq2", verbose = F)

de_markers$gene <- rownames(de_markers)


head(de_markers,n=100)

k1 = de_markers$avg_log2FC< -1 & de_markers$p_val <0.01
k2 = de_markers$avg_log2FC> 1 & de_markers$p_val <0.01
de_markers$change <- ifelse(k1,"down",ifelse(k2,"up","not"))

library(ggplot2)
library(RColorBrewer)
library(ggrepel)
library(RUnit)
library(ggforce)
library(tidyverse)
library(ggpubr)
library(ggprism)
library(paletteer)

up_genes <- de_markers[de_markers$p_val < 0.01 & de_markers$avg_log2FC > 1, ]
up_genes <- head(up_genes[order(up_genes$p_val), ], 10)

down_genes <- de_markers[de_markers$p_val < 0.01 & de_markers$avg_log2FC < -1, ]
down_genes <- head(down_genes[order(down_genes$p_val), ], 10)

top_de_markers <- rbind(up_genes, down_genes)

p=ggplot(de_markers, aes(avg_log2FC, -log10(p_val),color = change)) +   
  geom_point(size = 2, alpha = 0.5) +   
  geom_vline(xintercept = c(1,-1),linetype = 4)+  
  geom_hline(yintercept = -log10(0.01),linetype = 4)+  
  geom_text_repel(
    data = top_de_markers,
    aes(label = gene),
    size = 3,
    color = "black",
    segment.color = "black", show.legend = FALSE)+
  scale_color_manual(values = c("#2874C5", "grey", "#f87669"))+  
  theme_bw() + 
  ylab("-log10(unadjusted p-value)") 
p
ggsave("z16_plotdt_(TB VS H).pdf",p,width = 6,height = 4)

DefaultAssay(Tissue_stem_cells) <- "RNA"

bulk <- AggregateExpression(Tissue_stem_cells, return.seurat = T, slot = "counts", assays = "RNA",                             
                            group.by = c("celltype","orig.ident", "group1")) 

colnames(bulk)

library(DESeq2)
Idents(bulk) <- "group1"
de_markers <- FindMarkers(bulk, ident.1 ="Severe", ident.2 = "Mild",
                          slot = "counts", test.use = "DESeq2", verbose = F)

de_markers$gene <- rownames(de_markers)


head(de_markers,n=100)

k1 = de_markers$avg_log2FC< -1 & de_markers$p_val <0.01
k2 = de_markers$avg_log2FC> 1 & de_markers$p_val <0.01
de_markers$change <- ifelse(k1,"down",ifelse(k2,"up","not"))

up_genes <- de_markers[de_markers$p_val < 0.01 & de_markers$avg_log2FC > 1, ]
up_genes <- head(up_genes[order(up_genes$p_val), ], 10)

down_genes <- de_markers[de_markers$p_val < 0.01 & de_markers$avg_log2FC < -1, ]
down_genes <- head(down_genes[order(down_genes$p_val), ], 10)

top_de_markers <- rbind(up_genes, down_genes)

p=ggplot(de_markers, aes(avg_log2FC, -log10(p_val),color = change)) +   
  geom_point(size = 2, alpha = 0.5) +   
  geom_vline(xintercept = c(1,-1),linetype = 4)+  
  geom_hline(yintercept = -log10(0.01),linetype = 4)+  
  geom_text_repel(
    data = top_de_markers,
    aes(label = gene),
    size = 3,
    color = "black",
    segment.color = "black", show.legend = FALSE)+
  scale_color_manual(values = c("#2874C5", "grey", "#f87669"))+  
  theme_bw() + 
  ylab("-log10(unadjusted p-value)") 
p
ggsave("z16_plotdt_(S VS M) 207down51up.pdf",p,width = 6,height = 4)

Idents(Tissue_stem_cells)="celltype"
cell.markers <- FindAllMarkers(object = Tissue_stem_cells,
                               only.pos = FALSE,
                               test.use = "wilcox",
                               slot = "data",
                               min.pct = 0.25,
                               logfc.threshold = 0.25)
library(scRNAtoolVis)
p=jjVolcano(diffData = cell.markers,     
            log2FC.cutoff = 0.25,
            size=3.5,
            tile.col = c("#EE3B3B", "#76EEC6", "#F5F5DC", "#97FFFF", "#528B8B", "#9400D3", "#EE1289",  
                         "#00FF00", "#191970", "#FFFF00", "#4A708B", "#00FF7F", "#8B8B00",  
                         "#FF1493", "#FFA500", "#8B4513"))
p
ggsave("z20_jjVolcano.pdf", p, width = 6,height = 6)

DefaultAssay(Tissue_stem_cells) <- "RNA"
options(future.globals.maxSize = 20 * 1024^3)
Tissue_stem_cells <- NormalizeData(Tissue_stem_cells)
dim(Tissue_stem_cells)

expr = as.matrix(GetAssayData(object = Tissue_stem_cells@assays$RNA, layer="data"))
dim(expr)
expr <- expr[rowSums(expr)>0,]
dim(expr)
library(readxl)
geneSets2 = read_excel("genesets3.xls",col_names=F)
write.table(geneSets2,file = "my_genesets3.gmt",sep="\t",row.names = F,col.names = F,quote = F)
library(GSEABase)
geneSets=getGmt("my_genesets3.gmt",geneIdType=SymbolIdentifier())

library(rhdf5)
library(GSVA)

ssgsea_params <- ssgseaParam(exprData = expr, geneSets = geneSets, 
                             minSize = 1, maxSize = Inf, alpha = 0.25, normalize = TRUE)
gsvaResult <- gsva(ssgsea_params)

normalize=function(x){return((x-min(x))/(max(x)-min(x)))}
gsvaResult=normalize(gsvaResult)
gsvaOut=rbind(id=colnames(gsvaResult),gsvaResult)
write.table(gsvaOut, file="ssgseaOut3.txt",sep="\t",col.names = F,quote = F)
Tissue_stem_cells <- AddMetaData(Tissue_stem_cells, metadata = t(gsvaResult))

p <- FeaturePlot(Tissue_stem_cells, split.by = "group1", features = rownames(gsvaResult), reduction = "umap", ncol=2)
p
ggsave("z57_Tissue_stem_cells_Apoptosis_group1.pdf",p,width = 10,height = 3)

p <- FeaturePlot(Tissue_stem_cells, split.by = "group2", features = rownames(gsvaResult), reduction = "umap", ncol=2)
p
ggsave("z57_Tissue_stem_cells_Apoptosis_group2.pdf",p,width = 6,height = 3)

data <- data.frame(group1 = Tissue_stem_cells$group1,celltype = Tissue_stem_cells$celltype,Apoptosis_scores = Tissue_stem_cells$`Apoptosis scores`)
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
ggsave("z57_Apoptosis-group1-violin-celltype.pdf",plot = p,width = 8,height = 3.5)

data2 <- data.frame(group1 = Tissue_stem_cells$group2,celltype = Tissue_stem_cells$celltype,Apoptosis_scores = Tissue_stem_cells$`Apoptosis scores`)
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
ggsave("z57_Apoptosis-group2-celltype.pdf",plot = p,width = 8,height = 3.5)

data3 <- data.frame(group1 = Tissue_stem_cells$celltype,celltype = Tissue_stem_cells$celltype,Migration_scores = Tissue_stem_cells$`Apoptosis scores`)
p <- ggviolin(data3, x = "group1", y = "Migration_scores", fill = "group1",
              palette = c("#AB3282", "#53A85F","#F1BB72","#D6E7A3","#57C3F3","#E95C59","#E59CC4"),  # 设置颜色
              ylab = "Average Expression", xlab = "Group",
              add = "boxplot", add.params = list(fill = "white"))  # 添加箱线图
p
ggsave("z57_Apoptosis-celltype.pdf",plot = p,width = 4,height = 3)

geneSets2 = read_excel("genesets5.xls",col_names=F)
write.table(geneSets2,file = "my_genesets5.gmt",sep="\t",row.names = F,col.names = F,quote = F)
geneSets=getGmt("my_genesets5.gmt",geneIdType=SymbolIdentifier())

ssgsea_params <- ssgseaParam(exprData = expr, geneSets = geneSets, 
                             minSize = 1, maxSize = Inf, alpha = 0.25, normalize = TRUE)
gsvaResult <- gsva(ssgsea_params)

normalize=function(x){return((x-min(x))/(max(x)-min(x)))}
gsvaResult=normalize(gsvaResult)
gsvaOut=rbind(id=colnames(gsvaResult),gsvaResult)
write.table(gsvaOut, file="ssgseaOut5.txt",sep="\t",col.names = F,quote = F)
Tissue_stem_cells <- AddMetaData(Tissue_stem_cells, metadata = t(gsvaResult))

p <- FeaturePlot(Tissue_stem_cells, split.by = "group2", features = rownames(gsvaResult), reduction = "umap", ncol=2)
p
ggsave("z73_Tissue_stem_cells_Cytokine_group2.pdf",p,width = 6,height = 3)

p <- FeaturePlot(Tissue_stem_cells, split.by = "group1", features = rownames(gsvaResult), reduction = "umap", ncol=2)
p
ggsave("z73_Tissue_stem_cells_Cytokine_group1.pdf",p,width = 10,height = 3)

data <- data.frame(group1 = Tissue_stem_cells$group1,celltype = Tissue_stem_cells$celltype,Cytokine_scores = Tissue_stem_cells$`Cytokine scores`)
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
ggsave("z73_Cytokine-group1-violin-celltype.pdf",plot = p,width = 8,height = 3.5)

data2 <- data.frame(group1 = Tissue_stem_cells$group2,celltype = Tissue_stem_cells$celltype,Cytokine_scores = Tissue_stem_cells$`Cytokine scores`)
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
ggsave("z73_Cytokine-group2-celltype.pdf",plot = p,width = 8,height = 3.5)

data3 <- data.frame(group1 = Tissue_stem_cells$celltype,celltype = Tissue_stem_cells$celltype,Migration_scores = Tissue_stem_cells$`Cytokine scores`)
p <- ggviolin(data3, x = "group1", y = "Migration_scores", fill = "group1",
              palette = c("#AB3282", "#53A85F","#F1BB72","#D6E7A3","#57C3F3","#E95C59","#E59CC4"),  # 设置颜色
              ylab = "Average Expression", xlab = "Group",
              add = "boxplot", add.params = list(fill = "white"))  # 添加箱线图
p
ggsave("z73_Cytokine-celltype.pdf",plot = p,width = 6,height = 4)