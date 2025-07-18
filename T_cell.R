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

setwd("D:/1/T_cells/")

Idents(sce.all)
Idents(sce.all) <- 'celltype'

T_cells <- subset(sce.all, idents = "T_cells")

all.genes=rownames(T_cells)

g2m_genes <- cc.genes$g2m.genes
g2m_genes <- CaseMatch(search = g2m_genes, match = rownames(T_cells))

s_genes <- cc.genes$s.genes
s_genes <- CaseMatch(search = s_genes, match = rownames(T_cells))

options(future.globals.maxSize = 60000 * 1024^3)
T_cells <- NormalizeData(T_cells) %>% FindVariableFeatures() %>% ScaleData(features = all.genes) %>% RunPCA()

T_cells <- CellCycleScoring(T_cells, g2m.features = g2m_genes, s.features = s_genes)

colnames(T_cells@meta.data)
table(T_cells$Phase)
DimPlot(T_cells, group.by = "Phase")

T_cells <- NormalizeData(T_cells, normalization.method = "LogNormalize", scale.factor = 10000)
T_cells <- FindVariableFeatures(T_cells, selection.method = "vst", nfeatures = 2000)
T_cells <- ScaleData(T_cells, vars.to.regress = c("S.Score", "G2M.Score"))
DimPlot(T_cells, group.by = "Phase")

top10 <- head(VariableFeatures(T_cells), 10)
pdf(file = "y1.pdf", width = 7, height = 6)
VariableFeaturePlot(object = T_cells)
dev.off()

pdf(file = "y2.pdf", width = 7, height = 6)
LabelPoints(plot= VariableFeaturePlot(object = T_cells), points = top10, repel = TRUE)
dev.off()

T_cells <- RunPCA(T_cells, verbose = F)
pdf(file = "y3.pdf", width = 7, height = 6)
DimPlot(object = T_cells, reduction = "pca")
dev.off()

pdf(file = "y4.pdf", width = 10, height = 9)
VizDimLoadings(object = T_cells, dims = 1:4, reduction = "pca", nfeatures = 20)
dev.off()

pdf(file = "y5.pdf", width = 10, height = 9)
DimHeatmap(object = T_cells, dims = 1:4, cells = 500, balanced = TRUE, nfeatures = 30, ncol = 2)
dev.off()

pdf(file = "y6.pdf", width = 10, height = 9)
ElbowPlot(T_cells, ndims = 50)
dev.off()

pct <- sce.all [["pca"]]@stdev / sum(sce.all [["pca"]]@stdev) * 100
pct
cumu <- cumsum(pct)
cumu
co1 = which(cumu > 90 & pct < 5)[1]
co2 = sort(which((pct[1:length(pct)-1]-pct[2:length(pct)]) > 0.1), decreasing = T)[1]+1
pcs = min(co1, co2)
pcs
pcs = 1:15

T_cells <- RunHarmony(T_cells, group.by.vars="orig.ident", assay.use="SCT", max.iter.harmony = 20)
table(T_cells@meta.data$orig.ident)

library(clustree)
seq = seq(0.1,1,by=0.1)
T_cells <- FindNeighbors(T_cells, dims = pcs)
for (res in seq){
  T_cells = FindClusters(T_cells, resolution = res)
}

library(patchwork)
p8 = clustree(T_cells, prefix = "RNA_snn_res.")+coord_flip()
p = p8+plot_layout(widths = c(3,1))
p
ggsave("y7_T_RNA_snn_res.pdf", p, width = 30, height = 14)

for (i in c(10,11,12,13,14,15,16,17,18,19,20)) {
  T_cells <- FindNeighbors(T_cells,reduction = "harmony", dims = 1:i) %>% FindClusters(resolution = 0.1)
  T_cells <- RunUMAP(T_cells, reduction = "harmony", dims = 1:i)
  plot_i <- print(DimPlot(T_cells,reduction = "umap", label = TRUE, pt.size = 1, repel = TRUE, label.size = 5))
  plotname <- paste("plot_", i, sep = "")
  assign(plotname, plot_i)
  print(plot_i)
}
p = plot_10+plot_11+plot_12+plot_13+plot_14+plot_15+plot_16+plot_17+plot_18+plot_19+plot_20+plot_layout(ncol = 3)
ggsave(p,filename="dim.10-20.pdf",width = 15,height = 16)

for (i in seq) {
  T_cells <- FindNeighbors(T_cells,reduction = "harmony", dims = pcs) %>% FindClusters(resolution = i)
  T_cells <- RunUMAP(T_cells, reduction = "harmony", dims = pcs)
  plot_i <- print(DimPlot(T_cells,reduction = "umap", label = TRUE, pt.size = 1, repel = TRUE, label.size = 5))
  plotname <- paste("plot_", i, sep = "")
  assign(plotname, plot_i)
  print(plot_i)
}
p = plot_0.1+plot_0.2+plot_0.3+plot_0.4+plot_0.5+plot_0.6+plot_0.7+plot_0.8+plot_0.9+plot_1+plot_layout(ncol = 4)
ggsave(p,filename="dim.13,res.0.1-1.pdf",width = 20,height = 12)


T_cells <- FindNeighbors(T_cells, reduction = "harmony", dims = pcs) %>% FindClusters(resolution = 0.1)
T_cells <- RunUMAP(T_cells, reduction = "harmony", dims = pcs) %>% RunTSNE(dims = pcs, reduction = "harmony")
DimPlot(T_cells,reduction = "umap",cols = Biocols)

pdf(file = "y8.pdf", width = 7, height = 6)
DimPlot(T_cells, reduction = "umap", label = T,cols = Biocols)
dev.off()

pdf(file = "y9.pdf", width = 7, height = 6)
DimPlot(T_cells, reduction = "umap", label = F,group.by = "orig.ident",cols = Biocols)
dev.off()

pdf(file = "y10.pdf", width = 7, height = 6)
DimPlot(T_cells, reduction = "tsne", label = T,cols = Biocols)
dev.off()

pdf(file = "y11.pdf", width = 7, height = 6)
DimPlot(T_cells, reduction = "tsne", label = F,group.by = "orig.ident",cols = Biocols)
dev.off()

T_cells.markers1 <- FindAllMarkers(T_cells,only.pos = TRUE,logfc.threshold = 0.25)
write.csv(T_cells.markers1,file="T_markers.1.RNA.csv")

markers <- c("CD40LG","CD4","CD3D","CD3E","CD8A","CD8B","SLC4A10","TRAV1-2","TRDV2","TRGV9","KLRF1","NKG7","TYROBP","FGFBP2","FCGR3A","CX3CR1")

p11 <- FeaturePlot(T_cells,features = markers,ncol = 3)
p11
ggsave("y14.pdf",p11,width = 12,height = 12)

p12 <- DotPlot(T_cells,features = markers) + RotatedAxis()
p12
ggsave("y15.pdf",p12,width = 8,height = 4)

library(ComplexHeatmap)
library(pheatmap)
library(circlize)
group1_means <- AggregateExpression(
  object = T_cells,
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
                                                  'g7' = "#9FA3A8",
                                                  'g8' = "#CCE0F5")))
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
pdf("y15'.pdf", width = 5, height = 4)
print(p)
dev.off()

p13 <- VlnPlot(T_cells,features = markers,stack = T,flip = T) + NoLegend()
p13
ggsave("y16.pdf",p13,width = 8,height = 12)

T_cells$celltype <- recode(T_cells@meta.data$seurat_clusters,
                           "0" = "CD8_T",
                           "1" = "CD4_T",
                           "2" = "CD8_T",
                           "3" = "CD4_T",
                           "4" = "CD4_T",
                           "5" = "CD8_T",
                           "6" = "NK_cells",
                           "7" = "Neutrophils",
                           "8" = "B_cells")
table(T_cells@meta.data$celltype)
