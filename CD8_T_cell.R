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

setwd("D:/1/CD8_T_cells/")

Idents(T_cells)
Idents(T_cells) <- "celltype"

CD8_T_cells <- subset(T_cells, idents = "CD8_T")

all.genes=rownames(CD8_T_cells)

g2m_genes <- cc.genes$g2m.genes
g2m_genes <- CaseMatch(search = g2m_genes, match = rownames(CD8_T_cells))

s_genes <- cc.genes$s.genes
s_genes <- CaseMatch(search = s_genes, match = rownames(CD8_T_cells))

CD8_T_cells <- NormalizeData(CD8_T_cells) %>% FindVariableFeatures() %>% ScaleData(features = all.genes) %>% RunPCA()

CD8_T_cells <- CellCycleScoring(CD8_T_cells, g2m.features = g2m_genes, s.features = s_genes)

colnames(CD8_T_cells@meta.data)
table(CD8_T_cells$Phase)
DimPlot(CD8_T_cells, group.by = "Phase")

options(future.globals.maxSize = 60000 * 1024^3)
CD8_T_cells <- NormalizeData(CD8_T_cells, normalization.method = "LogNormalize", scale.factor = 10000)
CD8_T_cells <- FindVariableFeatures(CD8_T_cells, selection.method = "vst", nfeatures = 2000)
CD8_T_cells <- ScaleData(CD8_T_cells, vars.to.regress = c("S.Score", "G2M.Score"))
DimPlot(CD8_T_cells, group.by = "Phase")

top10 <- head(VariableFeatures(CD8_T_cells), 10)
pdf(file = "y1.pdf", width = 7, height = 6)
VariableFeaturePlot(object = CD8_T_cells)
dev.off()

pdf(file = "y2.pdf", width = 7, height = 6)
LabelPoints(plot= VariableFeaturePlot(object = CD8_T_cells), points = top10, repel = TRUE)
dev.off()

CD8_T_cells <- RunPCA(CD8_T_cells, verbose = F)
pdf(file = "y3.pdf", width = 7, height = 6)
DimPlot(object = CD8_T_cells, reduction = "pca")
dev.off()

pdf(file = "y4.pdf", width = 10, height = 9)
VizDimLoadings(object = CD8_T_cells, dims = 1:4, reduction = "pca", nfeatures = 20)
dev.off()

pdf(file = "y5.pdf", width = 10, height = 9)
DimHeatmap(object = CD8_T_cells, dims = 1:4, cells = 500, balanced = TRUE, nfeatures = 30, ncol = 2)
dev.off()

pdf(file = "y6.pdf", width = 10, height = 9)
ElbowPlot(CD8_T_cells, ndims = 50)
dev.off()

pct <- CD8_T_cells [["pca"]]@stdev / sum(CD8_T_cells [["pca"]]@stdev) * 100
pct
cumu <- cumsum(pct)
cumu
pcs = 1:15

seq = seq(0.1,1,by=0.1)
CD8_T_cells <- FindNeighbors(CD8_T_cells, dims = pcs)
for (res in seq){
  CD8_T_cells = FindClusters(CD8_T_cells, resolution = res)
}

library(patchwork)
p8 = clustree(CD8_T_cells, prefix = "RNA_snn_res.")+coord_flip()
p = p8+plot_layout(widths = c(3,1))
p
ggsave("y7_CD8_T_RNA_snn_res.pdf", p, width = 30, height = 14)

CD8_T_cells <- RunHarmony(CD8_T_cells, group.by.vars="orig.ident", max.iter.harmony = 20)
table(CD8_T_cells@meta.data$orig.ident)

for (i in c(10,11,12,13,14,15,16,17,18,19,20)) {
  CD8_T_cells <- FindNeighbors(CD8_T_cells,reduction = "harmony", dims = 1:i) %>% FindClusters(resolution = 0.1)
  CD8_T_cells <- RunUMAP(CD8_T_cells, reduction = "harmony", dims = 1:i)
  plot_i <- print(DimPlot(CD8_T_cells,reduction = "umap", label = TRUE, pt.size = 1, repel = TRUE, label.size = 5))
  plotname <- paste("plot_", i, sep = "")
  assign(plotname, plot_i)
  print(plot_i)
}
p = plot_10+plot_11+plot_12+plot_13+plot_14+plot_15+plot_16+plot_17+plot_18+plot_19+plot_20+plot_layout(ncol = 3)
ggsave(p,filename="dim.10-20.pdf",width = 15,height = 16)

for (i in c(25,30,35,40)) {
  CD8_T_cells <- FindNeighbors(CD8_T_cells,reduction = "harmony", dims = 1:i) %>% FindClusters(resolution = 0.1)
  CD8_T_cells <- RunUMAP(CD8_T_cells, reduction = "harmony", dims = 1:i)
  plot_i <- print(DimPlot(CD8_T_cells,reduction = "umap", label = TRUE, pt.size = 1, repel = TRUE, label.size = 5))
  plotname <- paste("plot_", i, sep = "")
  assign(plotname, plot_i)
  print(plot_i)
}
p = plot_25+plot_30+plot_35+plot_40+plot_layout(ncol = 2)
ggsave(p,filename="dim.25-40.pdf",width = 10,height = 10)

for (i in seq) {
  CD8_T_cells <- FindNeighbors(CD8_T_cells,reduction = "harmony", dims = pcs) %>% FindClusters(resolution = i)
  CD8_T_cells <- RunUMAP(CD8_T_cells, reduction = "harmony", dims = pcs)
  plot_i <- print(DimPlot(CD8_T_cells,reduction = "umap", label = TRUE, pt.size = 1, repel = TRUE, label.size = 5))
  plotname <- paste("plot_", i, sep = "")
  assign(plotname, plot_i)
  print(plot_i)
}
p = plot_0.1+plot_0.2+plot_0.3+plot_0.4+plot_0.5+plot_0.6+plot_0.7+plot_0.8+plot_0.9+plot_1+plot_layout(ncol = 4)
ggsave(p,filename="dim.15,res.0.1-1.pdf",width = 20,height = 12)


CD8_T_cells <- FindNeighbors(CD8_T_cells, reduction = "harmony", dims = pcs) %>% FindClusters(resolution = 0.3)
CD8_T_cells <- RunUMAP(CD8_T_cells, reduction = "harmony", dims = pcs) %>% RunTSNE(dims = pcs, reduction = "harmony")
DimPlot(CD8_T_cells,reduction = "umap",cols = Biocols)

pdf(file = "y8.pdf", width = 5, height = 4)
DimPlot(CD8_T_cells, reduction = "umap", label = T,cols = Biocols)
dev.off()

pdf(file = "y9.pdf", width = 5, height = 4)
DimPlot(CD8_T_cells, reduction = "umap", label = F,group.by = "orig.ident",cols = Biocols)
dev.off()

pdf(file = "y10.pdf", width = 5, height = 4)
DimPlot(CD8_T_cells, reduction = "tsne", label = T,cols = Biocols)
dev.off()

pdf(file = "y11.pdf", width = 5, height = 4)
DimPlot(CD8_T_cells, reduction = "tsne", label = F,group.by = "orig.ident",cols = Biocols)
dev.off()

CD8_T_cells.markers1 <- FindAllMarkers(CD8_T_cells,only.pos = TRUE,logfc.threshold = 0.25)
write.csv(CD8_T_cells.markers1,file="CD8_T_markers.1.RNA.csv")

markers <- c("CD3E","CD8A","CCR7","SELL","GPR183","S100A4","GZMA","GZMB","GZMH","GZMM","GZMK","GNLY","TYMS","MKI67","PDCD1","TIGIT","LAG3","CTLA4","HAVCR2","BTLA","CD80")

Idents(CD8_T_cells)

p11 <- FeaturePlot(CD8_T_cells,features = markers,ncol = 3)
p11
ggsave("y14.pdf",p11,width = 10,height = 20)

p12 <- DotPlot(CD8_T_cells,features = markers) + RotatedAxis()
p12
ggsave("y15.pdf",p12,width = 8,height = 4)

CD8_T_cells <- subset(CD8_T_cells, idents = c("0", "2", "5", "6"))

setwd("D:/1/CD8_T/0256/")

all.genes=rownames(CD8_T_cells)

g2m_genes <- cc.genes$g2m.genes
g2m_genes <- CaseMatch(search = g2m_genes, match = rownames(CD8_T_cells))

s_genes <- cc.genes$s.genes
s_genes <- CaseMatch(search = s_genes, match = rownames(CD8_T_cells))

CD8_T_cells <- NormalizeData(CD8_T_cells) %>% FindVariableFeatures() %>% ScaleData(features = all.genes) %>% RunPCA()

CD8_T_cells <- CellCycleScoring(CD8_T_cells, g2m.features = g2m_genes, s.features = s_genes)

colnames(CD8_T_cells@meta.data)
table(CD8_T_cells$Phase)
DimPlot(CD8_T_cells, group.by = "Phase")

options(future.globals.maxSize = 60000 * 1024^3)
CD8_T_cells <- NormalizeData(CD8_T_cells, normalization.method = "LogNormalize", scale.factor = 10000)
CD8_T_cells <- FindVariableFeatures(CD8_T_cells, selection.method = "vst", nfeatures = 2000)
CD8_T_cells <- ScaleData(CD8_T_cells, vars.to.regress = c("S.Score", "G2M.Score"))
DimPlot(CD8_T_cells, group.by = "Phase")

top10 <- head(VariableFeatures(CD8_T_cells), 10)
pdf(file = "y1.pdf", width = 7, height = 6)
VariableFeaturePlot(object = CD8_T_cells)
dev.off()

pdf(file = "y2.pdf", width = 7, height = 6)
LabelPoints(plot= VariableFeaturePlot(object = CD8_T_cells), points = top10, repel = TRUE)
dev.off()

CD8_T_cells <- RunPCA(CD8_T_cells, verbose = F)
pdf(file = "y3.pdf", width = 7, height = 6)
DimPlot(object = CD8_T_cells, reduction = "pca")
dev.off()

pdf(file = "y4.pdf", width = 10, height = 9)
VizDimLoadings(object = CD8_T_cells, dims = 1:4, reduction = "pca", nfeatures = 20)
dev.off()

pdf(file = "y5.pdf", width = 10, height = 9)
DimHeatmap(object = CD8_T_cells, dims = 1:4, cells = 500, balanced = TRUE, nfeatures = 30, ncol = 2)
dev.off()

pdf(file = "y6.pdf", width = 10, height = 9)
ElbowPlot(CD8_T_cells, ndims = 50)
dev.off()

pct <- CD8_T_cells [["pca"]]@stdev / sum(CD8_T_cells [["pca"]]@stdev) * 100
pct
cumu <- cumsum(pct)
cumu
pcs = 1:16

seq = seq(0.1,1,by=0.1)
CD8_T_cells <- FindNeighbors(CD8_T_cells, dims = pcs)
for (res in seq){
  CD8_T_cells = FindClusters(CD8_T_cells, resolution = res)
}

library(patchwork)
p8 = clustree(CD8_T_cells, prefix = "RNA_snn_res.")+coord_flip()
p = p8+plot_layout(widths = c(3,1))
p
ggsave("y7_CD8_T_cells_RNA_snn_res.pdf", p, width = 30, height = 14)

CD8_T_cells <- RunHarmony(CD8_T_cells, group.by.vars="orig.ident", max.iter.harmony = 20)
table(CD8_T_cells@meta.data$orig.ident)

for (i in c(10,11,12,13,14,15,16,17,18,19,20)) {
  CD8_T_cells <- FindNeighbors(CD8_T_cells,reduction = "harmony", dims = 1:i) %>% FindClusters(resolution = 0.1)
  CD8_T_cells <- RunUMAP(CD8_T_cells, reduction = "harmony", dims = 1:i)
  plot_i <- print(DimPlot(CD8_T_cells,reduction = "umap", label = TRUE, pt.size = 1, repel = TRUE, label.size = 5))
  plotname <- paste("plot_", i, sep = "")
  assign(plotname, plot_i)
  print(plot_i)
}
p = plot_10+plot_11+plot_12+plot_13+plot_14+plot_15+plot_16+plot_17+plot_18+plot_19+plot_20+plot_layout(ncol = 3)
ggsave(p,filename="dim.10-20.pdf",width = 15,height = 16)

for (i in c(25,30,35,40)) {
  CD8_T_cells <- FindNeighbors(CD8_T_cells,reduction = "harmony", dims = 1:i) %>% FindClusters(resolution = 0.1)
  CD8_T_cells <- RunUMAP(CD8_T_cells, reduction = "harmony", dims = 1:i)
  plot_i <- print(DimPlot(CD8_T_cells,reduction = "umap", label = TRUE, pt.size = 1, repel = TRUE, label.size = 5))
  plotname <- paste("plot_", i, sep = "")
  assign(plotname, plot_i)
  print(plot_i)
}
p = plot_25+plot_30+plot_35+plot_40+plot_layout(ncol = 2)
ggsave(p,filename="dim.25-40.pdf",width = 10,height = 10)

for (i in seq) {
  CD8_T_cells <- FindNeighbors(CD8_T_cells,reduction = "harmony", dims = pcs) %>% FindClusters(resolution = i)
  CD8_T_cells <- RunUMAP(CD8_T_cells, reduction = "harmony", dims = pcs)
  plot_i <- print(DimPlot(CD8_T_cells,reduction = "umap", label = TRUE, pt.size = 1, repel = TRUE, label.size = 5))
  plotname <- paste("plot_", i, sep = "")
  assign(plotname, plot_i)
  print(plot_i)
}
p = plot_0.1+plot_0.2+plot_0.3+plot_0.4+plot_0.5+plot_0.6+plot_0.7+plot_0.8+plot_0.9+plot_1+plot_layout(ncol = 4)
ggsave(p,filename="dim.16,res.0.1-1.pdf",width = 20,height = 12)


CD8_T_cells <- FindNeighbors(CD8_T_cells, reduction = "harmony", dims = pcs) %>% FindClusters(resolution = 0.3)
CD8_T_cells <- RunUMAP(CD8_T_cells, reduction = "harmony", dims = pcs) %>% RunTSNE(dims = pcs, reduction = "harmony")
DimPlot(CD8_T_cells,reduction = "umap",cols = Biocols)

pdf(file = "y8.pdf", width = 5, height = 4)
DimPlot(CD8_T_cells, reduction = "umap", label = T,cols = Biocols)
dev.off()

pdf(file = "y9.pdf", width = 5, height = 4)
DimPlot(CD8_T_cells, reduction = "umap", label = F,group.by = "orig.ident",cols = Biocols)
dev.off()

pdf(file = "y10.pdf", width = 5, height = 4)
DimPlot(CD8_T_cells, reduction = "tsne", label = T,cols = Biocols)
dev.off()

pdf(file = "y11.pdf", width = 5, height = 4)
DimPlot(CD8_T_cells, reduction = "tsne", label = F,group.by = "orig.ident",cols = Biocols)
dev.off()

CD8_T_cells.markers1 <- FindAllMarkers(CD8_T_cells,only.pos = TRUE,logfc.threshold = 0.25)
write.csv(CD8_T_cells.markers1,file="CD8_T_cells_markers.1.RNA.csv")

markers <- c("CD3E","CD8A","CCR7","TCF7","SATB1","GPR183","S100A4","GZMA","GZMB","GZMH","GZMM","GZMK","GNLY","PDCD1","TIGIT","LAG3","CTLA4","HAVCR2","BTLA","CD80")

Idents(CD8_T_cells)

p11 <- FeaturePlot(CD8_T_cells,features = markers,ncol = 3)
p11
ggsave("y14.pdf",p11,width = 10,height = 20)

p12 <- DotPlot(CD8_T_cells,features = markers) + RotatedAxis()
p12
ggsave("y15.pdf",p12,width = 8,height = 4)

library(ComplexHeatmap)
library(pheatmap)
library(circlize)
group1_means <- AggregateExpression(
  object = CD8_T_cells,
  features = markers,
  group.by = "RNA_snn_res.0.3",
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
                                                  'g4' = "#57C3F3")))#颜色设置
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
pdf("y15'.pdf", width = 4, height = 6)
print(p)
dev.off()

p13 <- VlnPlot(CD8_T_cells,features = markers,stack = T,flip = T) + NoLegend()
p13
ggsave("y16.pdf",p13,width = 5,height = 10)

CD8_T_cells$celltype <- recode(CD8_T_cells@meta.data$seurat_clusters,
                               "0" = "CD8_Tem",
                               "1" = "CD8_Te",
                               "2" = "CD8_Naive",
                               "3" = "CD8_Tex",
                               "4" = "CD8_Tem")
table(CD8_T_cells@meta.data$celltype)

p = DimPlot(CD8_T_cells, reduction = "umap", label = T, group.by = "celltype",cols = Biocols)
p
ggsave("y17.pdf", p, width = 6,height = 4)
p = DimPlot(CD8_T_cells, reduction = "tsne", label = T, group.by = "celltype",cols = Biocols)
p
ggsave("y18.pdf", p, width = 6,height = 4)

p = DimPlot(CD8_T_cells, reduction = "umap", label = F, group.by = "celltype",cols = Biocols)
p
ggsave("y17'.pdf", p, width = 6,height = 4)
p = DimPlot(CD8_T_cells, reduction = "tsne", label = F, group.by = "celltype",cols = Biocols)
p
ggsave("y18'.pdf", p, width = 6,height = 4)

table(Idents(CD8_T_cells))
Idents(CD8_T_cells)="celltype"
table(Idents(CD8_T_cells))

CD8_T_cells.markers3 <- FindAllMarkers(CD8_T_cells, only.pos = TRUE, logfc.threshold =1)
write.csv(CD8_T_cells.markers3, file = "CD8_T_markers.celltype.RNA.csv")

top10 <- CD8_T_cells.markers3 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

markers = as.data.frame(top10[,"gene"])
CD8_T_cells <- ScaleData(CD8_T_cells, features = as.character(unique(markers$gene)))
p = DoHeatmap(CD8_T_cells,
              features = as.character(unique(markers$gene)),
              group.by = "celltype")
p
ggsave("y19.pdf", p, width = 10,height = 10)

allCells = names(Idents(CD8_T_cells))
allType = levels(Idents(CD8_T_cells))
choose_Cells = unlist(lapply(allType, function(x){
  cgCells = allCells[Idents(CD8_T_cells) == x]
  n = min(table(CD8_T_cells@meta.data$celltype))
  if (length(cgCells) >= n) {
    cg = sample(cgCells, n)
  } else {
    cg = sample(cgCells, length(cgCells))
  }
  cg
}))

cg_sce = CD8_T_cells[, allCells %in% choose_Cells]
table(Idents(cg_sce))

p = DoHeatmap(cg_sce,
              features = as.character(unique(markers$gene)),
              group.by = "celltype")
p
ggsave("y20.pdf", p, width = 10,height = 10)

cell.prop <- as.data.frame(prop.table(table(CD8_T_cells@meta.data$celltype,CD8_T_cells@meta.data$orig.ident)))
colnames(cell.prop)<-c("cluster", "group", "proportion")

p= ggplot(cell.prop, aes(group,proportion,fill=cluster))+
  geom_bar(stat = "identity",position = "fill")+
  ggtitle("")+
  theme_bw()+
  theme(axis.ticks.length = unit(0.5, "cm"))+
  guides(fill=guide_legend(title = NULL))
p
ggsave("y21.pdf", p, width = 7,height = 5)

cell.prop <- as.data.frame(prop.table(table(CD8_T_cells@meta.data$celltype,CD8_T_cells@meta.data$group1)))
colnames(cell.prop)<-c("cluster", "group", "proportion")

p= ggplot(cell.prop, aes(group,proportion,fill=cluster))+
  geom_bar(stat = "identity",position = "fill")+
  ggtitle("")+
  theme_bw()+
  theme(axis.ticks.length = unit(0.5, "cm"))+
  guides(fill=guide_legend(title = NULL))
p
ggsave("y22.pdf", p, width = 4,height = 4)

cell.prop <- as.data.frame(prop.table(table(CD8_T_cells@meta.data$celltype,CD8_T_cells@meta.data$group2)))
colnames(cell.prop)<-c("cluster", "group", "proportion")

p= ggplot(cell.prop, aes(group,proportion,fill=cluster))+
  geom_bar(stat = "identity",position = "fill")+
  ggtitle("")+
  theme_bw()+
  theme(axis.ticks.length = unit(0.5, "cm"))+
  guides(fill=guide_legend(title = NULL))
p
ggsave("y23.pdf", p, width = 4,height = 4)

dim(CD8_T_cells)
DimPlot(CD8_T_cells,group.by = "celltype",label = T)
CD8_T_cells$celltype <- CD8_T_cells$celltype
colnames(CD8_T_cells@meta.data)

data <- GetAssayData(CD8_T_cells,assay = "RNA",slot = "counts")

pd <- new('AnnotatedDataFrame',data = CD8_T_cells@meta.data[,c(1,11,22,24)])

fData <- data.frame(gene_short_name = row.names(data),row.names = row.names(data))
fd <- new("AnnotatedDataFrame", data = fData)
library(monocle)
dim(CD8_T_cells)
mycds <- newCellDataSet(as.matrix(data),phenoData = pd,featureData = fd,lowerDetectionLimit = 0.5,expressionFamily = negbinomial.size())

mycds <- estimateSizeFactors(mycds)
mycds <- estimateDispersions(mycds,cores=6)

disp_table <- dispersionTable(mycds)
order.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit) %>% pull(gene_id) %>%as.character()

order.genes <- subset(CD8_T_cells.markers3, p_val_adj < 0.05)$gene

mycds <- setOrderingFilter(mycds,order.genes)

plot_ordering_genes(mycds)
p <- plot_ordering_genes(mycds)
ggsave("z0_OrderGenes.pdf",p,width = 8,height = 6)

mycds <- reduceDimension(mycds,max_components = 2,reduction_method = 'DDRTree',residualModelFormulaStr = "~orig.ident")

mycds <- orderCells(mycds)

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

genes = c(order.genes[1:4])
p1 = plot_genes_in_pseudotime(mycds[genes],color_by = "State")
p2 = plot_genes_in_pseudotime(mycds[genes],color_by = "celltype")
p3 = plot_genes_in_pseudotime(mycds[genes],color_by = "Pseudotime")
ggsave("z9_Trajectory_pseudotime.pdf",plot = p1|p2|p3,width = 7,height = 4)

p1 = plot_genes_jitter(mycds[genes],grouping="State",color_by = "State")
p2 = plot_genes_violin(mycds[genes],grouping="State",color_by = "State")
p3 = plot_genes_in_pseudotime(mycds[genes],color_by = "State")
ggsave("z10_Trajectory_jitter.pdf",plot = p1|p2|p3,width = 7,height = 4)

pData(mycds)$BTLA = log2(exprs(mycds)["BTLA",]+1)
p1 = plot_cell_trajectory(mycds,color_by = "BTLA") + 
  scale_color_continuous(type = "viridis")

pData(mycds)$TIGIT = log2(exprs(mycds)["TIGIT",]+1)
p2 = plot_cell_trajectory(mycds,color_by = "TIGIT") + 
  scale_color_continuous(type = "viridis")

pData(mycds)$PDCD1 = log2(exprs(mycds)["PDCD1",]+1)
p3 = plot_cell_trajectory(mycds,color_by = "PDCD1") + 
  scale_color_continuous(type = "viridis")

pData(mycds)$CTLA4 = log2(exprs(mycds)["CTLA4",]+1)
p4 = plot_cell_trajectory(mycds,color_by = "CTLA4") + 
  scale_color_continuous(type = "viridis")

pData(mycds)$CD80 = log2(exprs(mycds)["CD80",]+1)
p5 = plot_cell_trajectory(mycds,color_by = "CD80") + 
  scale_color_continuous(type = "viridis")

pData(mycds)$LAG3 = log2(exprs(mycds)["LAG3",]+1)
p6 = plot_cell_trajectory(mycds,color_by = "LAG3") + 
  scale_color_continuous(type = "viridis")

pData(mycds)$HAVCR2 = log2(exprs(mycds)["HAVCR2",]+1)
p7 = plot_cell_trajectory(mycds,color_by = "HAVCR2") + 
  scale_color_continuous(type = "viridis")

ggsave("z11_Trajectory_Expression.pdf",plot = p1|p2|p3|p4|p5|p6|p7,width = 18,height = 4)

pData(mycds)$BTLA = log2(exprs(mycds)["BTLA",]+1)
p1 = plot_cell_trajectory(mycds,color_by = "BTLA") + 
  scale_color_continuous(type = "viridis")+
  facet_wrap(~group1, nrow = 1)

pData(mycds)$TIGIT = log2(exprs(mycds)["TIGIT",]+1)
p2 = plot_cell_trajectory(mycds,color_by = "TIGIT") + 
  scale_color_continuous(type = "viridis")+
  facet_wrap(~group1, nrow = 1)

pData(mycds)$PDCD1 = log2(exprs(mycds)["PDCD1",]+1)
p3 = plot_cell_trajectory(mycds,color_by = "PDCD1") + 
  scale_color_continuous(type = "viridis")+
  facet_wrap(~group1, nrow = 1)

pData(mycds)$CTLA4 = log2(exprs(mycds)["CTLA4",]+1)
p4 = plot_cell_trajectory(mycds,color_by = "CTLA4") + 
  scale_color_continuous(type = "viridis")+
  facet_wrap(~group1, nrow = 1)

pData(mycds)$CD80 = log2(exprs(mycds)["CD80",]+1)
p5 = plot_cell_trajectory(mycds,color_by = "CD80") + 
  scale_color_continuous(type = "viridis")+
  facet_wrap(~group1, nrow = 1)

pData(mycds)$LAG3 = log2(exprs(mycds)["LAG3",]+1)
p6 = plot_cell_trajectory(mycds,color_by = "LAG3") + 
  scale_color_continuous(type = "viridis")+
  facet_wrap(~group1, nrow = 1)

pData(mycds)$HAVCR2 = log2(exprs(mycds)["HAVCR2",]+1)
p7 = plot_cell_trajectory(mycds,color_by = "HAVCR2") + 
  scale_color_continuous(type = "viridis")+
  facet_wrap(~group1, nrow = 1)
plot_layout <- p1 + p2 + p3 + p4 + p5 + p6 + p7 + 
  plot_layout(ncol = 1, nrow = 7)
ggsave("z12_Trajectory_Expression_group1.pdf",plot = plot_layout,width = 6,height = 24)

pdata <- Biobase::pData(mycds)
CD8_T_cells <- AddMetaData(CD8_T_cells,metadata = pdata[,c("Pseudotime","State")])

Time_diff <- differentialGeneTest(mycds,cores = 10,
                                  fullModelFormulaStr = "~sm.ns(Pseudotime)")
write.csv(Time_diff,"z10_Time_diff_all.csv",row.names = F)

Time_genes <-Time_diff[order(Time_diff$qval),"gene_short_name"][1:100]
p = plot_pseudotime_heatmap(mycds[Time_genes,],num_clusters = 3,
                            show_rownames = T,return_heatmap = T)
ggsave("z13_Time_heatmap.pdf",p,width = 8,height = 10)

hp.genes <- p$tree_row$labels[p$tree_row$order]
Time_diff_sig <- Time_diff[hp.genes,c("gene_short_name","pval","qval")]
write.csv(Time_diff_sig,"z14_Time_diff_sig.csv",row.names = F)

library(paletteer)
library(Seurat)
library(monocle3)
library(dplyr)
library(BiocParallel)
library(ggplot2)
library(ggsci)

DimPlot(CD8_T_cells,pt.size=0.8,group.by="celltype",label=T)
levels(Idents(CD8_T_cells))
Idents(CD8_T_cells)<-CD8_T_cells$celltype
expression_matrix<-GetAssayData(CD8_T_cells,assay='RNA',layer='counts')
cell_metadata<-CD8_T_cells@meta.data
gene_annotation<-data.frame(gene_short_name=rownames(expression_matrix))
rownames(gene_annotation)<-rownames(expression_matrix)

cds<-new_cell_data_set(expression_matrix,
                       cell_metadata=cell_metadata,
                       gene_metadata=gene_annotation)
cds<-preprocess_cds(cds,num_dim=50,norm_method=c("none"))
cds<-align_cds(cds,alignment_group="orig.ident")
cds<-reduce_dimension(cds,cores=5)
cds<-cluster_cells(cds,resolution=0.0000001)
cds<-learn_graph(cds)
myselect<-function(cds,select.classify,my_select){
  cell_ids<-which(colData(cds)[,select.classify]==my_select)
  closest_vertex<-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex<-as.matrix(closest_vertex[colnames(cds),])
  root_pr_nodes<-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]
  root_pr_nodes}
cds<-order_cells(cds,root_pr_nodes=myselect(cds,select.classify='celltype',my_select="CD8_Naive"))
cds.embed<-cds@int_colData$reducedDims$UMAP
int.embed<-Embeddings(CD8_T_cells,reduction="umap")
int.embed<-int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP<-int.embed
p=plot_cells(cds,color_cells_by="pseudotime",
             show_trajectory_graph=T)
p
ggsave("z14_trajectory-umap.pdf",p,width = 5,height = 4)
p=plot_cells(cds,
             color_cells_by="celltype",
             label_cell_groups=FALSE,
             label_leaves=FALSE,
             label_branch_points=FALSE,
             graph_label_size=1)+scale_color_manual(values=pal_jco("default",alpha=0.6)(9))
p
ggsave("z15_trajectory-umap-celltype.pdf",p,width = 5.5,height = 4)

Track_genes<-graph_test(cds,neighbor_graph="principal_graph",cores=6)

Track_genes_sig<-Track_genes%>%top_n(n=10,morans_I)%>%pull(gene_short_name)%>%as.character()
plot_genes_in_pseudotime(cds[Track_genes_sig,],color_cells_by="celltype",min_expr=0.5,ncol=2,cell_size=1.5)+scale_color_manual(values=pal_jco("default",alpha=0.6)(9))

genes<-c('PDCD1','HAVCR2','TIGIT','LAG3','CTLA4')
genes_cds<-cds[rowData(cds)$gene_short_name%in%genes,]
plot_genes_in_pseudotime(genes_cds,color_cells_by="celltype",min_expr=0.5,cell_size=1.5)+scale_color_manual(values=pal_jco("default",alpha=0.6)(9))

p=plot_cells(cds,genes="PDCD1",
             cell_size=1,
             show_trajectory_graph=FALSE,
             label_cell_groups=FALSE,
             label_leaves=FALSE)
p
ggsave("z16_trajectory-umap-PDCD1.pdf",p,width = 5,height = 4)

p=plot_cells(cds,genes="HAVCR2",
             cell_size=1,
             show_trajectory_graph=FALSE,
             label_cell_groups=FALSE,
             label_leaves=FALSE)
p
ggsave("z16_trajectory-umap-HAVCR2.pdf",p,width = 5,height = 4)

p=plot_cells(cds,genes="TIGIT",
             cell_size=1,
             show_trajectory_graph=FALSE,
             label_cell_groups=FALSE,
             label_leaves=FALSE)
p
ggsave("z16_trajectory-umap-TIGIT.pdf",p,width = 5,height = 4)

p=plot_cells(cds,genes="LAG3",
             cell_size=1,
             show_trajectory_graph=FALSE,
             label_cell_groups=FALSE,
             label_leaves=FALSE)
p
ggsave("z16_trajectory-umap-LAG3.pdf",p,width = 5,height = 4)

p=plot_cells(cds,genes="CTLA4",
             cell_size=1,
             show_trajectory_graph=FALSE,
             label_cell_groups=FALSE,
             label_leaves=FALSE)
p
ggsave("z16_trajectory-umap-CTLA4.pdf",p,width = 5,height = 4)

CD8_T_cells_DEGs <- FindMarkers(CD8_T_cells,
                                ident.1 ="Tuberculosis",
                                ident.2 ="Health",
                                min.pct =0.15,
                                logfc.threshold =0.25,
                                test.use ="wilcox")
filtered_DEGs <- CD8_T_cells_DEGs[CD8_T_cells_DEGs$p_val_adj <0.05 & abs(CD8_T_cells_DEGs$avg_log2FC) >1, ]
filtered_DEGs <- filtered_DEGs[order(filtered_DEGs$avg_log2FC, decreasing =TRUE), ]
up_genes <- sum(filtered_DEGs$avg_log2FC > 1)
down_genes <- sum(filtered_DEGs$avg_log2FC < -1)
cat("上调基因数量:", up_genes, "\n")
cat("下调基因数量:", down_genes, "\n")
cat("总差异基因数量:", nrow(filtered_DEGs), "\n")
write.csv(filtered_DEGs,"CD8_T_cells_Tuberculosis_vs_Health_filtered_DEGs.csv")

library(ggprism)
library(tidyverse)
library(org.Hs.eg.db)
library(clusterProfiler)

pal <- c('#7bc4e2','#acd372','#fbb05b','#ed6ca4')
pal <- c('#eaa052','#b74147','#90ad5b','#23929c')
pal <- c('#c3e1e6','#f3dfb7','#dcc6dc','#96c38e')

ego_readable <- setReadable(ego_all, OrgDb = org.Hs.eg.db, keyType ="ENTREZID")
GO <- as.data.frame(ego_readable)

ekegg_readable <- setReadable(ekegg_all, OrgDb = org.Hs.eg.db, keyType ="ENTREZID")
KEGG <- as.data.frame(ekegg_readable)

use_pathway <- group_by(GO, ONTOLOGY) %>%
  top_n(5, wt = -p.adjust) %>%
  group_by(p.adjust) %>%
  top_n(1, wt = Count) %>%
  rbind(
    top_n(KEGG,5, -p.adjust) %>%
      group_by(p.adjust) %>%
      top_n(1, wt = Count) %>%
      mutate(ONTOLOGY ='KEGG')
  ) %>%
  ungroup() %>%
  mutate(ONTOLOGY = factor(ONTOLOGY,
                           levels = rev(c('BP','CC','MF','KEGG')))) %>%
  dplyr::arrange(ONTOLOGY, p.adjust) %>%
  mutate(Description = factor(Description, levels = Description)) %>%
  tibble::rowid_to_column('index')
width <-0.5
xaxis_max <- max(-log10(use_pathway$p.adjust)) + 1
rect.data <- group_by(use_pathway, ONTOLOGY) %>%
  summarize(n = n()) %>%
  ungroup() %>%
  mutate(
    xmin = -3* width,
    xmax = -2* width,
    ymax = cumsum(n),
    ymin = lag(ymax, default =0) +0.6,
    ymax = ymax +0.4
  )

plot_enrichment <-function() {
  p <- use_pathway %>%
    ggplot(aes(-log10(p.adjust), y = index, fill = ONTOLOGY)) +
    geom_col(
      aes(y = Description), width =0.6, alpha =0.8
    ) +
    geom_text(
      aes(x =0.05, label = Description),
      hjust =0, size =5
    ) +
    geom_text(
      aes(x =0.1, label = geneID, colour = ONTOLOGY),
      hjust =0, vjust =2.6, size =3.5, fontface ='italic',
      show.legend =FALSE
    ) +
    geom_point(
      aes(x = -width, size = Count),
      shape =21
    ) +
    geom_text(
      aes(x = -width, label = Count)
    ) +
    scale_size_continuous(name ='Count', range = c(5,12)) +
    geom_rect(
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
          fill = ONTOLOGY),
      data = rect.data,
      radius = unit(2,'mm'),
      inherit.aes =FALSE
    ) +
    geom_text(
      aes(x = (xmin + xmax) /2, y = (ymin + ymax) /2, label = ONTOLOGY),
      data = rect.data,
      inherit.aes =FALSE
    ) +
    geom_segment(
      aes(x =0, y =0, xend = xaxis_max, yend =0),
      linewidth =1.5,
      inherit.aes =FALSE
    ) +
    labs(y =NULL) +
    scale_fill_manual(name ='Category', values = pal) +
    scale_colour_manual(values = pal) +
    scale_x_continuous(
      breaks = seq(0, xaxis_max,2),
      expand = expansion(c(0,0))
    ) +
    theme_prism() +
    theme(
      axis.text.y = element_blank(),
      axis.line = element_blank(),
      axis.ticks.y = element_blank(),
      legend.title = element_text()
    )
  return(p)
}

enrichment_plot <- plot_enrichment()

p=enrichment_plot
p
ggsave("B24 plot_enrichment T vs H.pdf",p,width = 20,height = 8)

Idents(CD8_T_cells) <- CD8_T_cells$group1
CD8_T_cells_DEGs <- FindMarkers(CD8_T_cells,
                                ident.1 ="Severe",
                                ident.2 ="Mild",
                                min.pct =0.15,
                                logfc.threshold =0.25,
                                test.use ="wilcox")
filtered_DEGs <- CD8_T_cells_DEGs[CD8_T_cells_DEGs$p_val <0.05 & abs(CD8_T_cells_DEGs$avg_log2FC) >1, ]
filtered_DEGs <- filtered_DEGs[order(filtered_DEGs$avg_log2FC, decreasing =TRUE), ]
up_genes <- sum(filtered_DEGs$avg_log2FC > 1)
down_genes <- sum(filtered_DEGs$avg_log2FC < -1)
cat("上调基因数量:", up_genes, "\n")
cat("下调基因数量:", down_genes, "\n")
cat("总差异基因数量:", nrow(filtered_DEGs), "\n")
write.csv(filtered_DEGs,"CD8_T_cells_Severe_vs_Mild_filtered_DEGs.csv")

ego_readable <- setReadable(ego_all, OrgDb = org.Hs.eg.db, keyType ="ENTREZID")
GO <- as.data.frame(ego_readable)

ekegg_readable <- setReadable(ekegg_all, OrgDb = org.Hs.eg.db, keyType ="ENTREZID")
KEGG <- as.data.frame(ekegg_readable)

use_pathway <- group_by(GO, ONTOLOGY) %>%
  top_n(5, wt = -p.adjust) %>%
  group_by(p.adjust) %>%
  top_n(1, wt = Count) %>%
  rbind(
    top_n(KEGG,5, -p.adjust) %>%
      group_by(p.adjust) %>%
      top_n(1, wt = Count) %>%
      mutate(ONTOLOGY ='KEGG')
  ) %>%
  ungroup() %>%
  mutate(ONTOLOGY = factor(ONTOLOGY,
                           levels = rev(c('BP','CC','MF','KEGG')))) %>%
  dplyr::arrange(ONTOLOGY, p.adjust) %>%
  mutate(Description = factor(Description, levels = Description)) %>%
  tibble::rowid_to_column('index')
width <-0.5
xaxis_max <- max(-log10(use_pathway$p.adjust)) + 1
rect.data <- group_by(use_pathway, ONTOLOGY) %>%
  summarize(n = n()) %>%
  ungroup() %>%
  mutate(
    xmin = -3* width,
    xmax = -2* width,
    ymax = cumsum(n),
    ymin = lag(ymax, default =0) +0.6,
    ymax = ymax +0.4
  )

plot_enrichment <-function() {
  p <- use_pathway %>%
    ggplot(aes(-log10(p.adjust), y = index, fill = ONTOLOGY)) +
    geom_col(
      aes(y = Description), width =0.6, alpha =0.8
    ) +
    geom_text(
      aes(x =0.05, label = Description),
      hjust =0, size =5
    ) +
    geom_text(
      aes(x =0.1, label = geneID, colour = ONTOLOGY),
      hjust =0, vjust =2.6, size =3.5, fontface ='italic',
      show.legend =FALSE
    ) +
    geom_point(
      aes(x = -width, size = Count),
      shape =21
    ) +
    geom_text(
      aes(x = -width, label = Count)
    ) +
    scale_size_continuous(name ='Count', range = c(5,12)) +
    geom_rect(
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
          fill = ONTOLOGY),
      data = rect.data,
      radius = unit(2,'mm'),
      inherit.aes =FALSE
    ) +
    geom_text(
      aes(x = (xmin + xmax) /2, y = (ymin + ymax) /2, label = ONTOLOGY),
      data = rect.data,
      inherit.aes =FALSE
    ) +
    geom_segment(
      aes(x =0, y =0, xend = xaxis_max, yend =0),
      linewidth =1.5,
      inherit.aes =FALSE
    ) +
    labs(y =NULL) +
    scale_fill_manual(name ='Category', values = pal) +
    scale_colour_manual(values = pal) +
    scale_x_continuous(
      breaks = seq(0, xaxis_max,2),
      expand = expansion(c(0,0))
    ) +
    theme_prism() +
    theme(
      axis.text.y = element_blank(),
      axis.line = element_blank(),
      axis.ticks.y = element_blank(),
      legend.title = element_text()
    )
  return(p)
}

enrichment_plot <- plot_enrichment()

p=enrichment_plot
p
ggsave("B24 plot_enrichment S vs M.pdf",p,width = 12,height = 7)

cluster3 <- unique(markers$cluster)
gene_lists <- lapply(cluster3, function(cluster3) {
  input_gene <- as.character(subset(markers, cluster == cluster3)$gene)
  input_gene <- input_gene[!is.na(input_gene)]
  entrezIDs <- mget(input_gene, org.Hs.egSYMBOL2EG, ifnotfound = NA)
  entrezIDs <- as.character(entrezIDs)
  gene <- entrezIDs[entrezIDs != "NA"]
  gene <- gsub("c\\(\"(\\d+)\".*", "\\1", gene)
  return(gene)
})

names(gene_lists) <- cluster3

safe_compareCluster <- safely(compareCluster)

go_comparison <- compareCluster(gene_lists, fun = "enrichGO", 
                                OrgDb = org.Hs.eg.db, 
                                ont = "ALL", 
                                pvalueCutoff = 0.05)
summary(go_comparison)
p=dotplot(go_comparison, showCategory = 10, title = "Cluster Comparison")
p
ggsave("B18_GO_dotplot.pdf",p,width = 8,height = 12)

p=cnetplot(go_comparison, showCategory = 10, title = "GO Enrichment Comparison")
p
ggsave("B19_GO_cnet.pdf",p,width = 20,height = 20)

kegg_comparison <- compareCluster(gene_lists, fun = "enrichKEGG", 
                                  organism = "hsa", 
                                  pvalueCutoff = 0.05)
summary(kegg_comparison)
p=dotplot(kegg_comparison, showCategory = 10, title = "KEGG Enrichment Comparison")
p
ggsave("B20_KEGG_dotplot.pdf",p,width = 6,height = 8)

p=cnetplot(kegg_comparison, showCategory = 10, title = "KEGG Enrichment Comparison")
p
ggsave("B21_GO_cnet.pdf",p,width = 20,height = 20)

DefaultAssay(CD8_T_cells) <- "RNA"

bulk <- AggregateExpression(CD8_T_cells, return.seurat = T, slot = "counts", assays = "RNA",                             
                            group.by = c("celltype","group1", "group2")) 

colnames(bulk)

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
library(openxlsx)

up_genes <- de_markers[de_markers$p_val < 0.01 & de_markers$avg_log2FC > 1, ]
write.xlsx(up_genes, file = "B22_up_genes(TB vs Health).xlsx")
up_genes <- head(up_genes[order(up_genes$p_val), ], 10)

down_genes <- de_markers[de_markers$p_val < 0.05 & de_markers$avg_log2FC < -1, ]
write.xlsx(down_genes, file = "B23_down_genes(TB vs Health).xlsx")
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
ggsave("B24_plotdt_(TB VS H).pdf",p,width = 6,height = 4)

Idents(bulk) <- "group1"
de_markers2 <- FindMarkers(bulk, ident.1 ="Severe", ident.2 = "Mild",
                           slot = "counts", test.use = "DESeq2", verbose = F)

de_markers2$gene <- rownames(de_markers2)

head(de_markers2,n=100)

k1 = de_markers2$avg_log2FC< -1 & de_markers2$p_val <0.01
k2 = de_markers2$avg_log2FC> 1 & de_markers2$p_val <0.01
de_markers2$change <- ifelse(k1,"down",ifelse(k2,"up","not"))

# 提取上调基因并选择前10个
up_genes2 <- de_markers2[de_markers2$p_val < 0.01 & de_markers2$avg_log2FC > 1, ]
write.xlsx(up_genes2, file = "B22_up_genes(Severe vs Mild).xlsx")
up_genes2 <- head(up_genes2[order(up_genes2$p_val), ], 10)

down_genes2 <- de_markers2[de_markers2$p_val < 0.01 & de_markers2$avg_log2FC < -1, ]
write.xlsx(down_genes2, file = "B23_down_genes(Severe vs Mild).xlsx")
down_genes2 <- head(down_genes2[order(down_genes2$p_val), ], 10)

top_de_markers2 <- rbind(up_genes2, down_genes2)

p=ggplot(de_markers2, aes(avg_log2FC, -log10(p_val),color = change)) +   
  geom_point(size = 2, alpha = 0.5) +   
  geom_vline(xintercept = c(1,-1),linetype = 4)+  
  geom_hline(yintercept = -log10(0.01),linetype = 4)+  
  geom_text_repel(
    data = top_de_markers2,
    aes(label = gene),
    size = 3,
    color = "black",
    segment.color = "black", show.legend = FALSE)+
  scale_color_manual(values = c("#2874C5", "grey", "#f87669"))+  
  theme_bw() + 
  ylab("-log10(unadjusted p-value)") 
p
ggsave("B24_plotdt_(Severe vs Mild).pdf",p,width = 6,height = 4)

library(scRNAtoolVis)
p=jjVolcano(diffData = cell.markers,     
            log2FC.cutoff = 0.25,
            size=3.5,
            legend.position = "right",
            tile.col = Biocols)
p
ggsave("z21_jjVolcano.pdf", p, width = 12,height = 10)

Idents(CD8_T_cells) <- "celltype"
table(Idents(CD8_T_cells))
dim(CD8_T_cells)

p=DimPlot(CD8_T_cells, group.by = "group1", cols = Biocols)
p
ggsave("y24_UMAP-group1.pdf",p,width = 5,height = 4)

p=DimPlot(CD8_T_cells, group.by = "group2", cols = Biocols)
p
ggsave("y25_UMAP-group2.pdf",p,width = 5,height = 4)

p=DimPlot(CD8_T_cells, split.by = "group1", cols = Biocols)
p
ggsave("y26_UMAP-group1.pdf",p,width = 10,height = 4)

p=DimPlot(CD8_T_cells, split.by = "group2", cols = Biocols)
p
ggsave("y27_UMAP-group2.pdf",p,width = 8,height = 4)

table(Idents(CD8_T_cells))
Idents(CD8_T_cells) <- "celltype"
library(CellChat)
data.input <- GetAssayData(CD8_T_cells, assay = "RNA", slot = "data")
meta <- CD8_T_cells@meta.data
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "celltype")
CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)
cellchat@DB
cellchat@DB <- CellChatDB
cellchat <- subsetData(cellchat)
future::plan("multisession", workers=6)
options(future.globals.maxSize = 2 * 1024^3) 
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat,raw.use = F)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
df.net <- subsetCommunication(cellchat)
write.csv(df.net, "I1_Gene.csv")
df.netP <- subsetCommunication(cellchat, slot.name="netP")
write.csv(df.netP, "I2_Pathway.csv")
cellchat <- aggregateNet(cellchat)

groupSize <- as.numeric(table(cellchat@idents))
pdf("I3_NetVisual_overview_all.pdf", width = 6, height = 4)
par(xpd=TRUE)
netVisual_circle(cellchat@net$count,vertex.weight = groupSize, weight.scale = T,
                 label.edge = F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight,vertex.weight = groupSize, weight.scale = T,
                 label.edge = F, title.name = "Interaction weights/strength")
dev.off()

pdf("I4_NetVisual_overview_split.pdf", width = 6, height = 4)
mat <- cellchat@net$weight
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i,] <- mat[i,]
  par(xpd=TRUE)
  netVisual_circle(mat2,vertex.weight = groupSize,weight.scale = T,
                   edge.weight.max = max(mat), title.name = rownames(mat)[i])
  
}
dev.off()

mypathways <- cellchat@netP$pathways
mypathways <- mypathways[1:10]
mypathways

pdf("I5_NetVisual_pathways_circle.pdf", width = 6, height = 4)
for (pathways.show in mypathways) {
  par(xpd=TRUE)
  netVisual_aggregate(cellchat,signaling = pathways.show,layout = "circle")
  
}
dev.off()

pdf("I6_NetVisual_pathways_chord.pdf", width = 6, height = 4)
for (pathways.show in mypathways) {
  netVisual_aggregate(cellchat,signaling = pathways.show,layout = "chord")
}
dev.off()

pdf("I7_NetVisual_pathways_heatmap.pdf", width = 4, height = 3)
for (pathways.show in mypathways) {
  par(xpd=TRUE)
  p <- netVisual_heatmap(cellchat,signaling = pathways.show,color.heatmap = "Reds")
  plot(p)
}
dev.off()

pdf("I8_Pathways.pdf",width = 6, height = 4)
for (pathways.show in mypathways) {
  netAnalysis_contribution(cellchat,signaling = pathways.show)
  pairLR <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)$interaction_name
  for (LR.show in pairLR) {
    netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")
  }
  for (LR.show in pairLR) {
    netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "chord") 
  }
  dev.off()
}

levels(cellchat@idents)
p <- netVisual_bubble(cellchat, sources.use = 1:length(levels(cellchat@idents)),
                      targets.use = 1:length(levels(cellchat@idents)), remove.isolate = FALSE)
ggsave("I9_CCI_all.pdf",p,width = 6,height = 16,limitsize = F)

p <- netVisual_bubble(cellchat, signaling = c("CCL"), remove.isolate = FALSE)
ggsave("I11_CCI_subcell_subpathway.pdf",p,width = 4,height = 4,limitsize = F)

p <- netVisual_bubble(cellchat, signaling = c("MIF"), remove.isolate = FALSE)
ggsave("I11_CCI_subcell_subpathway-MIF.pdf",p,width = 4,height = 4,limitsize = F)

p <- netVisual_bubble(cellchat, signaling = c("TGFb"), remove.isolate = FALSE)
ggsave("I11_CCI_subcell_subpathway-TGFb.pdf",p,width = 6,height = 4,limitsize = F)

p <- netVisual_bubble(cellchat, signaling = c("CLEC"), remove.isolate = FALSE)
ggsave("I11_CCI_subcell_subpathway-CLEC.pdf",p,width = 3.5,height = 3,limitsize = F)

p <- netVisual_bubble(cellchat, signaling = c("PARs"), remove.isolate = FALSE)
ggsave("I11_CCI_subcell_subpathway-PARs.pdf",p,width = 6,height = 4,limitsize = F)

pairLR.use <- c("CCL3_CCR1","CCL4_CCR5","CCL5_CCR3","TNF_TNFRSF1A","TNF_TNFRSF1B")
pairLR.use <- data.frame(interaction_name=pairLR.use)
p <- netVisual_bubble(cellchat, sources.use = 1:3,targets.use = 4:5, pairLR.use = pairLR.use, remove.isolate = FALSE)
ggsave("I12_CCI_subcell_subLR.pdf",p,width = 4,height = 4,limitsize = F)

p<-plotGeneExpression(cellchat,signaling = "CCL")
ggsave("I13_GeneExpressionn_violin_sig.pdf",p,width = 4,height = 5,limitsize = F)

p<-plotGeneExpression(cellchat,signaling = "MIF")
ggsave("I13_GeneExpressionn_violin_sig-MIF.pdf",p,width = 4,height = 5,limitsize = F)

p<-plotGeneExpression(cellchat,signaling = "TGFb")
ggsave("I13_GeneExpressionn_violin_sig-TGFb.pdf",p,width = 4,height = 5,limitsize = F)

p<-plotGeneExpression(cellchat,signaling = "CLEC")
ggsave("I13_GeneExpressionn_violin_sig-CLEC.pdf",p,width = 4,height = 4,limitsize = F)

p<-plotGeneExpression(cellchat,signaling = "PARs")
ggsave("I13_GeneExpressionn_violin_sig-PARs.pdf",p,width = 4,height = 4.5,limitsize = F)

p<-plotGeneExpression(cellchat,signaling = "CCL", enriched.only = FALSE)
ggsave("I14_GeneExpressionn_violin_all.pdf",p,width = 4,height = 9,limitsize = F)

cellchat2 <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
pdf("I15_SignalingRole.pdf",width = 5,height = 6)
for (pathways.show in mypathways) {
  netAnalysis_signalingRole_network(cellchat2,signaling = pathways.show,
                                    width = 8,height = 2.5,font.size = 10)
}                                          
dev.off()

DefaultAssay(CD8_T_cells) <- "RNA"
options(future.globals.maxSize = 20 * 1024^3)
CD8_T_cells <- NormalizeData(CD8_T_cells)
dim(CD8_T_cells)

expr = as.matrix(GetAssayData(object = CD8_T_cells@assays$RNA, layer="data"))
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
write.table(gsvaOut, file="ssgseaOut.txt",sep="\t",col.names = F,quote = F)
CD8_T_cells <- AddMetaData(CD8_T_cells, metadata = t(gsvaResult))

p <- FeaturePlot(CD8_T_cells, split.by = "group1", features = rownames(gsvaResult), reduction = "umap", ncol=2)
p
ggsave("z32_CD8_exhaust_group1.pdf",p,width = 10,height = 3)

p <- FeaturePlot(CD8_T_cells, split.by = "group2", features = rownames(gsvaResult), reduction = "umap", ncol=2)
p
ggsave("z32_CD8_exhaust_group2.pdf",p,width = 6,height = 3)

my_comparisons <- list(c("Health", "Mild"), c("Health", "Severe"), c("Mild", "Severe"))
comparisons <- list(c("Health", "Tuberculosis"))

data <- data.frame(group1 = CD8_T_cells$group1,celltype = CD8_T_cells$celltype,Genes_of_cell_exhaustion_response = CD8_T_cells$`Exhaustion scores`)
data$group1 <- factor(data$group1, levels = c("Health", "Mild", "Severe"))
p <- ggviolin(data, x = "group1", y = "Genes_of_cell_exhaustion_response", fill = "group1",
              palette = c("#00AFBB", "#E7B800", "#FC4E07"),  # 设置颜色
              ylab = "Average Expression", xlab = "Group",
              add = "boxplot", add.params = list(fill = "white"))  # 添加箱线图
p
p <- p + stat_compare_means(method = "wilcox.test", comparisons = my_comparisons, label = "p.signif")
p
ggsave("z33_exhaustion-group1-violin.pdf",plot = p,width = 3.5,height = 4.5)

p <- ggviolin(data, x = "group1", y = "Genes_of_cell_exhaustion_response", fill = "group1",
              palette = c("#00AFBB", "#E7B800", "#FC4E07"),  # 设置颜色
              ylab = "Average Expression", xlab = "Group",
              add = "boxplot", add.params = list(fill = "white")) +  # 添加箱线图
  facet_wrap(~ celltype, scales = "free",ncol=4)  # 按celltype分面，并允许每个分面有不同的y轴范围
p

p <- p + stat_compare_means(method = "wilcox.test", label = "p.signif", comparisons = my_comparisons)
p
ggsave("z33_exhaustion-group1-violin-celltype.pdf",plot = p,width = 10,height = 5)

data2 <- data.frame(group1 = CD8_T_cells$group2,celltype = CD8_T_cells$celltype,Genes_of_cell_exhaustion_response = CD8_T_cells$`Exhaustion scores`)
data2$group1 <- factor(data2$group1, levels = c("Health", "Tuberculosis"))
p <- ggviolin(data2, x = "group1", y = "Genes_of_cell_exhaustion_response", fill = "group1",
              palette = c("#00AFBB", "#E76980"),  # 设置颜色
              ylab = "Average Expression", xlab = "Group",
              add = "boxplot", add.params = list(fill = "white"))  # 添加箱线图
p <- p + stat_compare_means(method = "wilcox.test", comparisons = comparisons, label = "p.signif")
p
ggsave("z33_exhaustion-group2.pdf",plot = p,width = 3.2,height = 4)

p <- ggviolin(data2, x = "group1", y = "Genes_of_cell_exhaustion_response", fill = "group1",
              palette = c("#00AFBB", "#E76980"),  # 设置颜色
              ylab = "Average Expression", xlab = "Group",
              add = "boxplot", add.params = list(fill = "white")) +  # 添加箱线图
  facet_wrap(~ celltype, scales = "free",ncol=4)  # 按celltype分面，并允许每个分面有不同的y轴范围
p <- p + stat_compare_means(method = "wilcox.test", label = "p.signif", comparisons = comparisons)
p
ggsave("z33_exhaustion-group2-celltype.pdf",plot = p,width = 10,height = 5)

data3 <- data.frame(group1 = CD8_T_cells$celltype,celltype = CD8_T_cells$celltype,Migration_scores = CD8_T_cells$`Exhaustion scores`)
data3$group1 <- factor(data3$group1, levels = c("CD8_Naive","CD8_Te", "CD8_Tem","CD8_Tem_exh"))
p <- ggviolin(data3, x = "group1", y = "Migration_scores", fill = "group1",
              palette = c("#AB3282", "#53A85F","#F1BB72","#D6E7A3"),  # 设置颜色
              ylab = "Average Expression", xlab = "Group",
              add = "boxplot", add.params = list(fill = "white"))  # 添加箱线图
p
ggsave("z33_exhaustion-celltype.pdf",plot = p,width = 5,height = 4)

geneSets2 = read_excel("genesets1.xls",col_names=F)
write.table(geneSets2,file = "my_genesets1.gmt",sep="\t",row.names = F,col.names = F,quote = F)
library(GSEABase)
geneSets=getGmt("my_genesets1.gmt",geneIdType=SymbolIdentifier())

library(rhdf5)
library(GSVA)

ssgsea_params <- ssgseaParam(exprData = expr, geneSets = geneSets, 
                             minSize = 1, maxSize = Inf, alpha = 0.25, normalize = TRUE)
gsvaResult <- gsva(ssgsea_params)

normalize=function(x){return((x-min(x))/(max(x)-min(x)))}
gsvaResult=normalize(gsvaResult)
gsvaOut=rbind(id=colnames(gsvaResult),gsvaResult)
write.table(gsvaOut, file="ssgseaOut1.txt",sep="\t",col.names = F,quote = F)
CD8_T_cells <- AddMetaData(CD8_T_cells, metadata = t(gsvaResult))

p <- FeaturePlot(CD8_T_cells, split.by = "group1", features = rownames(gsvaResult), reduction = "umap", ncol=2)
p
ggsave("z41_CD8_Cytotoxicity_group1.pdf",p,width = 10,height = 3)

p <- FeaturePlot(CD8_T_cells, split.by = "group2", features = rownames(gsvaResult), reduction = "umap", ncol=2)
p
ggsave("z41_CD8_Cytotoxicity_group2.pdf",p,width = 6,height = 3)

data <- data.frame(group1 = CD8_T_cells$group1,celltype = CD8_T_cells$celltype,Cytotoxicity_scores = CD8_T_cells$`Cytotoxicity scores`)
data$group1 <- factor(data$group1, levels = c("Health", "Mild", "Severe"))
p <- ggviolin(data, x = "group1", y = "Cytotoxicity_scores", fill = "group1",
              palette = c("#00AFBB", "#E7B800", "#FC4E07"),  # 设置颜色
              ylab = "Average Expression", xlab = "Group",
              add = "boxplot", add.params = list(fill = "white"))  # 添加箱线图
p
p <- p + stat_compare_means(method = "wilcox.test", comparisons = my_comparisons, label = "p.signif")
p
ggsave("z43_Cytotoxicity-violin.pdf",plot = p,width = 3.5,height = 4.5)

p <- ggviolin(data, x = "group1", y = "Cytotoxicity_scores", fill = "group1",
              palette = c("#00AFBB", "#E7B800", "#FC4E07"),  # 设置颜色
              ylab = "Average Expression", xlab = "Group",
              add = "boxplot", add.params = list(fill = "white")) +  # 添加箱线图
  facet_wrap(~ celltype, scales = "free",ncol=4)  # 按celltype分面，并允许每个分面有不同的y轴范围
p

p <- p + stat_compare_means(method = "wilcox.test", label = "p.signif", comparisons = my_comparisons)
p
ggsave("z43_Cytotoxicity-violin-celltype.pdf",plot = p,width = 10,height = 5)

data2 <- data.frame(group1 = CD8_T_cells$group2,celltype = CD8_T_cells$celltype,Cytotoxicity_scores = CD8_T_cells$`Cytotoxicity scores`)
data2$group1 <- factor(data2$group1, levels = c("Health", "Tuberculosis"))
p <- ggviolin(data2, x = "group1", y = "Cytotoxicity_scores", fill = "group1",
              palette = c("#00AFBB", "#E76980"),  # 设置颜色
              ylab = "Average Expression", xlab = "Group",
              add = "boxplot", add.params = list(fill = "white"))  # 添加箱线图
p <- p + stat_compare_means(method = "wilcox.test", comparisons = comparisons, label = "p.signif")
p
ggsave("z43_Cytotoxicity-group2.pdf",plot = p,width = 3.2,height = 4)

p <- ggviolin(data2, x = "group1", y = "Cytotoxicity_scores", fill = "group1",
              palette = c("#00AFBB", "#E76980"),  # 设置颜色
              ylab = "Average Expression", xlab = "Group",
              add = "boxplot", add.params = list(fill = "white")) +  # 添加箱线图
  facet_wrap(~ celltype, scales = "free",ncol=4)  # 按celltype分面，并允许每个分面有不同的y轴范围
p <- p + stat_compare_means(method = "wilcox.test", label = "p.signif", comparisons = comparisons)
p
ggsave("z43_Cytotoxicity-group2-celltype.pdf",plot = p,width = 10,height = 5)

data3 <- data.frame(group1 = CD8_T_cells$celltype,celltype = CD8_T_cells$celltype,Migration_scores = CD8_T_cells$`Cytotoxicity scores`)
data3$group1 <- factor(data3$group1, levels = c("CD8_Naive","CD8_Te", "CD8_Tem","CD8_Tem_exh"))
p <- ggviolin(data3, x = "group1", y = "Migration_scores", fill = "group1",
              palette = c("#AB3282", "#53A85F","#F1BB72","#D6E7A3"),  # 设置颜色
              ylab = "Average Expression", xlab = "Group",
              add = "boxplot", add.params = list(fill = "white"))  # 添加箱线图
p
ggsave("z33_exhaustion-celltype.pdf",plot = p,width = 5,height = 4)

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
CD8_T_cells <- AddMetaData(CD8_T_cells, metadata = t(gsvaResult))

p <- FeaturePlot(CD8_T_cells, split.by = "group1", features = rownames(gsvaResult), reduction = "umap", ncol=2)
p
ggsave("z49_CD8_Migration_group1.pdf",p,width = 10,height = 3)

p <- FeaturePlot(CD8_T_cells, split.by = "group2", features = rownames(gsvaResult), reduction = "umap", ncol=2)
p
ggsave("z49_CD8_Migration_group2.pdf",p,width = 6,height = 3)

data <- data.frame(group1 = CD8_T_cells$group1,celltype = CD8_T_cells$celltype,Migration_scores = CD8_T_cells$`Migration scores`)
data$group1 <- factor(data$group1, levels = c("Health", "Mild", "Severe"))
p <- ggviolin(data, x = "group1", y = "Migration_scores", fill = "group1",
              palette = c("#00AFBB", "#E7B800", "#FC4E07"),  # 设置颜色
              ylab = "Average Expression", xlab = "Group",
              add = "boxplot", add.params = list(fill = "white"))  # 添加箱线图
p
p <- p + stat_compare_means(method = "wilcox.test", comparisons = my_comparisons, label = "p.signif")
p
ggsave("z49_Migration-group1-violin.pdf",plot = p,width = 3.5,height = 4.5)

p <- ggviolin(data, x = "group1", y = "Migration_scores", fill = "group1",
              palette = c("#00AFBB", "#E7B800", "#FC4E07"),  # 设置颜色
              ylab = "Average Expression", xlab = "Group",
              add = "boxplot", add.params = list(fill = "white")) +  # 添加箱线图
  facet_wrap(~ celltype, scales = "free",ncol=4)  # 按celltype分面，并允许每个分面有不同的y轴范围
p

p <- p + stat_compare_means(method = "wilcox.test", label = "p.signif", comparisons = my_comparisons)
p
ggsave("z49_Migration-group1-violin-celltype.pdf",plot = p,width = 10,height = 5)

data2 <- data.frame(group1 = CD8_T_cells$group2,celltype = CD8_T_cells$celltype,Migration_scores = CD8_T_cells$`Migration scores`)
data2$group1 <- factor(data2$group1, levels = c("Health", "Tuberculosis"))
p <- ggviolin(data2, x = "group1", y = "Migration_scores", fill = "group1",
              palette = c("#00AFBB", "#E76980"),  # 设置颜色
              ylab = "Average Expression", xlab = "Group",
              add = "boxplot", add.params = list(fill = "white"))  # 添加箱线图
p <- p + stat_compare_means(method = "wilcox.test", label = "p.signif", comparisons = comparisons)
p
ggsave("z49_Migration-group2.pdf",plot = p,width = 3.2,height = 4)

p <- ggviolin(data2, x = "group1", y = "Migration_scores", fill = "group1",
              palette = c("#00AFBB", "#E76980"),  # 设置颜色
              ylab = "Average Expression", xlab = "Group",
              add = "boxplot", add.params = list(fill = "white")) +  # 添加箱线图
  facet_wrap(~ celltype, scales = "free",ncol=4)  # 按celltype分面，并允许每个分面有不同的y轴范围
p <- p + stat_compare_means(method = "wilcox.test", label = "p.signif", comparisons = comparisons)
p
ggsave("z49_Migration-group2-celltype.pdf",plot = p,width = 10,height = 5)

data3 <- data.frame(group1 = CD8_T_cells$celltype,celltype = CD8_T_cells$celltype,Migration_scores = CD8_T_cells$`Migration scores`)
data3$group1 <- factor(data3$group1, levels = c("CD8_Naive","CD8_Te", "CD8_Tem","CD8_Tem_exh"))
p <- ggviolin(data3, x = "group1", y = "Migration_scores", fill = "group1",
              palette = c("#AB3282", "#53A85F","#F1BB72","#D6E7A3"),  # 设置颜色
              ylab = "Average Expression", xlab = "Group",
              add = "boxplot", add.params = list(fill = "white"))  # 添加箱线图
p
ggsave("z49_Migration-celltype.pdf",plot = p,width = 5,height = 4)

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
CD8_T_cells <- AddMetaData(CD8_T_cells, metadata = t(gsvaResult))

p <- FeaturePlot(CD8_T_cells, split.by = "group1", features = rownames(gsvaResult), reduction = "umap", ncol=2)
p
ggsave("z57_CD8_Apoptosis_group1.pdf",p,width = 10,height = 3)

p <- FeaturePlot(CD8_T_cells, split.by = "group2", features = rownames(gsvaResult), reduction = "umap", ncol=2)
p
ggsave("z57_CD8_Apoptosis_group2.pdf",p,width = 6,height = 3)

data <- data.frame(group1 = CD8_T_cells$group1,celltype = CD8_T_cells$celltype,Apoptosis_scores = CD8_T_cells$`Apoptosis scores`)
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

p <- p + stat_compare_means(method = "wilcox.test", label = "p.signif", comparisons = my_comparisons)
p
ggsave("z57_Apoptosis-group1-violin-celltype.pdf",plot = p,width = 10,height = 5)

data2 <- data.frame(group1 = CD8_T_cells$group2,celltype = CD8_T_cells$celltype,Apoptosis_scores = CD8_T_cells$`Apoptosis scores`)
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
ggsave("z57_Apoptosis-group2-celltype.pdf",plot = p,width = 10,height = 5)

data3 <- data.frame(group1 = CD8_T_cells$celltype,celltype = CD8_T_cells$celltype,Migration_scores = CD8_T_cells$`Apoptosis scores`)
data3$group1 <- factor(data3$group1, levels = c("CD8_Naive","CD8_Te", "CD8_Tem","CD8_Tem_exh"))
p <- ggviolin(data3, x = "group1", y = "Migration_scores", fill = "group1",
              palette = c("#AB3282", "#53A85F","#F1BB72","#D6E7A3"),  # 设置颜色
              ylab = "Average Expression", xlab = "Group",
              add = "boxplot", add.params = list(fill = "white"))  # 添加箱线图
p
ggsave("z57_Apoptosis-celltype.pdf",plot = p,width = 5,height = 4)

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
CD8_T_cells <- AddMetaData(CD8_T_cells, metadata = t(gsvaResult))

p <- FeaturePlot(CD8_T_cells, split.by = "group2", features = rownames(gsvaResult), reduction = "umap", ncol=2)
p
ggsave("z73_CD8_Cytokine_group2.pdf",p,width = 6,height = 3)

p <- FeaturePlot(CD8_T_cells, split.by = "group1", features = rownames(gsvaResult), reduction = "umap", ncol=2)
p
ggsave("z73_CD8_Cytokine_group1.pdf",p,width = 10,height = 3)

data <- data.frame(group1 = CD8_T_cells$group1,celltype = CD8_T_cells$celltype,Cytokine_scores = CD8_T_cells$`Cytokine scores`)
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

p <- p + stat_compare_means(method = "wilcox.test", label = "p.signif", comparisons = my_comparisons)
p
ggsave("z73_Cytokine-group1-violin-celltype.pdf",plot = p,width = 10,height = 5)

data2 <- data.frame(group1 = CD8_T_cells$group2,celltype = CD8_T_cells$celltype,Cytokine_scores = CD8_T_cells$`Cytokine scores`)
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
ggsave("z73_Cytokine-group2-celltype.pdf",plot = p,width = 10,height = 5)

data3 <- data.frame(group1 = CD8_T_cells$celltype,celltype = CD8_T_cells$celltype,Migration_scores = CD8_T_cells$`Cytokine scores`)
data3$group1 <- factor(data3$group1, levels = c("CD8_Naive","CD8_Te", "CD8_Tem","CD8_Tem_exh"))
p <- ggviolin(data3, x = "group1", y = "Migration_scores", fill = "group1",
              palette = c("#AB3282", "#53A85F","#F1BB72","#D6E7A3"),  # 设置颜色
              ylab = "Average Expression", xlab = "Group",
              add = "boxplot", add.params = list(fill = "white"))  # 添加箱线图
p
ggsave("z73_Cytokine-celltype.pdf",plot = p,width = 5,height = 4)

library(ComplexHeatmap)
library(pheatmap)
library(circlize)
exhaustion <- c("PDCD1","TIGIT","LAG3","CTLA4","HAVCR2","BTLA","CD80")
group1_means <- AggregateExpression(
  object = CD8_T_cells,
  features = exhaustion,
  group.by = "group1",
  slot = "data"
)
group1_means <- as.data.frame(group1_means$RNA)
df <- data.frame(colnames(group1_means))
colnames(df) <- 'class'
top_anno = HeatmapAnnotation(df = df,
                             border = T,
                             show_annotation_name = F,
                             gp = gpar(col = 'black'),
                             col = list(class = c('Health'="#00AFBB",
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
pdf("z32_exhaustion-heatmap-group1.pdf", width = 3, height = 2)
print(p)
dev.off()

group1_means_celltype <- AggregateExpression(
  object = CD8_T_cells,
  features = exhaustion,
  group.by = "celltype.group",
  slot = "data"
)
group1_means_celltype <- as.data.frame(group1_means_celltype$RNA)

CD8_Naive <- group1_means_celltype[c("CD8-Naive-Health","CD8-Naive-Mild","CD8-Naive-Severe")]
df <- data.frame(colnames(CD8_Naive))
colnames(df) <- 'class'
top_anno = HeatmapAnnotation(df = df,
                             border = T,
                             show_annotation_name = F,
                             gp = gpar(col = 'black'),
                             col = list(class = c('CD8-Naive-Health'="#00AFBB",
                                                  'CD8-Naive-Severe'="#FC4E07",
                                                  'CD8-Naive-Mild'="#E7B800")))#颜色设置
marker_exp <- t(scale(t(CD8_Naive),scale = T,center = F))
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
pdf("z32_exhaustion-heatmap-group1-CD8_Naive.pdf", width = 3.8, height = 2)
print(p)
dev.off()

CD8_Tem <- group1_means_celltype[c("CD8-Tem-Health","CD8-Tem-Mild","CD8-Tem-Severe")]
df <- data.frame(colnames(CD8_Tem))
colnames(df) <- 'class'
top_anno = HeatmapAnnotation(df = df,#细胞名/cluster
                             border = T,
                             show_annotation_name = F,
                             gp = gpar(col = 'black'),
                             col = list(class = c('CD8-Tem-Health'="#00AFBB",
                                                  'CD8-Tem-Severe'="#FC4E07",
                                                  'CD8-Tem-Mild'="#E7B800")))#颜色设置
marker_exp <- t(scale(t(CD8_Tem),scale = T,center = F))
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
pdf("z32_exhaustion-heatmap-group1-CD8_Tem.pdf", width = 3.7, height = 2)
print(p)
dev.off()

CD8_Tex <- group1_means_celltype[c("CD8-Tex-Health","CD8-Tex-Mild","CD8-Tex-Severe")]
df <- data.frame(colnames(CD8_Tex))
colnames(df) <- 'class'
top_anno = HeatmapAnnotation(df = df,
                             border = T,
                             show_annotation_name = F,
                             gp = gpar(col = 'black'),
                             col = list(class = c('CD8-Tex-Health'="#00AFBB",
                                                  'CD8-Tex-Severe'="#FC4E07",
                                                  'CD8-Tex-Mild'="#E7B800")))#颜色设置
marker_exp <- t(scale(t(CD8_Tex),scale = T,center = F))
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
pdf("z32_exhaustion-heatmap-group1-CD8_Tex.pdf", width = 4, height = 2)
print(p)
dev.off()

CD8_Te <- group1_means_celltype[c("CD8-Te-Health","CD8-Te-Mild","CD8-Te-Severe")]
df <- data.frame(colnames(CD8_Te))
colnames(df) <- 'class'
top_anno = HeatmapAnnotation(df = df,
                             border = T,
                             show_annotation_name = F,
                             gp = gpar(col = 'black'),
                             col = list(class = c('CD8-Te-Health'="#00AFBB",
                                                  'CD8-Te-Severe'="#FC4E07",
                                                  'CD8-Te-Mild'="#E7B800")))#颜色设置
marker_exp <- t(scale(t(CD8_Te),scale = T,center = F))
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
pdf("z32_exhaustion-heatmap-group1-CD8_Te.pdf", width = 3.6, height = 2)
print(p)
dev.off()


cytotoxicity <- c("PRF1",	"IFNG",	"GNLY",	"NKG7",	"GZMA",	"GZMK", "GZMM",	"GZMH",	"KLRK1",	"KLRB1",	"KLRD1",	"CTSW",	"CST7",	"FCGR3A",	"FGFBP2",	"ZEB2")
group1_means <- AggregateExpression(
  object = CD8_T_cells,
  features = cytotoxicity,
  group.by = "group1",
  slot = "data"
)
group1_means <- as.data.frame(group1_means$RNA)
df <- data.frame(colnames(group1_means))
colnames(df) <- 'class'
top_anno = HeatmapAnnotation(df = df,
                             border = T,
                             show_annotation_name = F,
                             gp = gpar(col = 'black'),
                             col = list(class = c( 'Health'="#00AFBB",
                                                   'Severe'="#FC4E07",
                                                   'Mild'="#E7B800")))#颜色设置
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
pdf("z41_cytotoxicity-heatmap-group1.pdf", width = 3, height = 4)
print(p)
dev.off()


group1_means_celltype <- AggregateExpression(
  object = CD8_T_cells,
  features = cytotoxicity,
  group.by = "celltype.group",
  slot = "data"
)
group1_means_celltype <- as.data.frame(group1_means_celltype$RNA)

CD8_Naive <- group1_means_celltype[c("CD8-Naive-Health","CD8-Naive-Mild","CD8-Naive-Severe")]
df <- data.frame(colnames(CD8_Naive))
colnames(df) <- 'class'
top_anno = HeatmapAnnotation(df = df,#细胞名/cluster
                             border = T,
                             show_annotation_name = F,
                             gp = gpar(col = 'black'),
                             col = list(class = c('CD8-Naive-Health'="#00AFBB",
                                                  'CD8-Naive-Severe'="#FC4E07",
                                                  'CD8-Naive-Mild'="#E7B800")))#颜色设置
marker_exp <- t(scale(t(CD8_Naive),scale = T,center = F))
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
pdf("z41_cytotoxicity-heatmap-group1-CD8_Naive.pdf", width = 3.8, height = 4)
print(p)
dev.off()

CD8_Tem <- group1_means_celltype[c("CD8-Tem-Health","CD8-Tem-Mild","CD8-Tem-Severe")]
df <- data.frame(colnames(CD8_Tem))
colnames(df) <- 'class'
top_anno = HeatmapAnnotation(df = df,#细胞名/cluster
                             border = T,
                             show_annotation_name = F,
                             gp = gpar(col = 'black'),
                             col = list(class = c('CD8-Tem-Health'="#00AFBB",
                                                  'CD8-Tem-Mild'="#E7B800",
                                                  'CD8-Tem-Severe'="#FC4E07")))#颜色设置
marker_exp <- t(scale(t(CD8_Tem),scale = T,center = F))
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
pdf("z41_cytotoxicity-heatmap-group1-CD8_Tem.pdf", width = 3.7, height = 4)
print(p)
dev.off()

CD8_Tex <- group1_means_celltype[c("CD8-Tex-Health","CD8-Tex-Mild","CD8-Tex-Severe")]
df <- data.frame(colnames(CD8_Tem_exh))
colnames(df) <- 'class'
top_anno = HeatmapAnnotation(df = df,
                             border = T,
                             show_annotation_name = F,
                             gp = gpar(col = 'black'),
                             col = list(class = c('CD8-Tex-Health'="#00AFBB",
                                                  'CD8-Tex-Mild'="#E7B800",
                                                  'CD8-Tex-Severe'="#FC4E07")))#颜色设置
marker_exp <- t(scale(t(CD8_Tex),scale = T,center = F))
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
pdf("z41_cytotoxicity-heatmap-group1-CD8_Tex.pdf", width = 4, height = 4)
print(p)
dev.off()

CD8_Te <- group1_means_celltype[c("CD8-Te-Health","CD8-Te-Mild","CD8-Te-Severe")]
df <- data.frame(colnames(CD8_Te))
colnames(df) <- 'class'
top_anno = HeatmapAnnotation(df = df,
                             border = T,
                             show_annotation_name = F,
                             gp = gpar(col = 'black'),
                             col = list(class = c('CD8-Te-Health'="#00AFBB",
                                                  'CD8-Te-Mild'="#E7B800",
                                                  'CD8-Te-Severe'="#FC4E07")))#颜色设置
marker_exp <- t(scale(t(CD8_Te),scale = T,center = F))
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
pdf("z41_cytotoxicity-heatmap-group1-CD8_Te.pdf", width = 3.6, height = 4)
print(p)
dev.off()

migration <- c("CCL20",	"CXCL13",	"CXCL10",	"XCL2",	"CCL22",	"CCL5",	"CCL2",	"CXCL8",	"CCL26",	"CCL3",	"XCL1",	"CXCL1",	"CCL17",	"CXCL2",	"CCL23",	"CXCL2",	"CXCL16",	"CXCR1",	"CCL25",	"CCL13",	"CCL14",	"CCL4",	"CX3CL1")
group1_means <- AggregateExpression(
  object = CD8_T_cells,
  features = migration,
  group.by = "group1",
  slot = "data"
)
group1_means <- as.data.frame(group1_means$RNA)
df <- data.frame(colnames(group1_means))
colnames(df) <- 'class'
top_anno = HeatmapAnnotation(df = df,
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
pdf("z49_migration-heatmap-group1.pdf", width = 3, height = 4)
print(p)
dev.off()

group1_means_celltype <- AggregateExpression(
  object = CD8_T_cells,
  features = migration,
  group.by = "celltype.group",
  slot = "data"
)
group1_means_celltype <- as.data.frame(group1_means_celltype$RNA)

CD8_Naive <- group1_means_celltype[c("CD8-Naive-Health","CD8-Naive-Mild","CD8-Naive-Severe")]
df <- data.frame(colnames(CD8_Naive))
colnames(df) <- 'class'
top_anno = HeatmapAnnotation(df = df,
                             border = T,
                             show_annotation_name = F,
                             gp = gpar(col = 'black'),
                             col = list(class = c('CD8-Naive-Health'="#00AFBB",
                                                  'CD8-Naive-Severe'="#FC4E07",
                                                  'CD8-Naive-Mild'="#E7B800")))#颜色设置
marker_exp <- t(scale(t(CD8_Naive),scale = T,center = F))
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
pdf("z49_migration-heatmap-group1-CD8_Naive.pdf", width = 3.8, height = 4)
print(p)
dev.off()

CD8_Tem <- group1_means_celltype[c("CD8-Tem-Health","CD8-Tem-Mild","CD8-Tem-Severe")]
df <- data.frame(colnames(CD8_Tem))
colnames(df) <- 'class'
top_anno = HeatmapAnnotation(df = df,
                             border = T,
                             show_annotation_name = F,
                             gp = gpar(col = 'black'),
                             col = list(class = c('CD8-Tem-Health'="#00AFBB",
                                                  'CD8-Tem-Mild'="#E7B800",
                                                  'CD8-Tem-Severe'="#FC4E07")))#颜色设置
marker_exp <- t(scale(t(CD8_Tem),scale = T,center = F))
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
pdf("z49_migration-heatmap-group1-CD8_Tem.pdf", width = 3.7, height = 4)
print(p)
dev.off()

CD8_Tex <- group1_means_celltype[c("CD8-Tex-Health","CD8-Tex-Mild","CD8-Tex-Severe")]
df <- data.frame(colnames(CD8_Tem_exh))
colnames(df) <- 'class'
top_anno = HeatmapAnnotation(df = df,
                             border = T,
                             show_annotation_name = F,
                             gp = gpar(col = 'black'),
                             col = list(class = c('CD8-Tex-Health'="#00AFBB",
                                                  'CD8-Tex-Mild'="#E7B800",
                                                  'CD8-Tex-Severe'="#FC4E07")))#颜色设置
marker_exp <- t(scale(t(CD8_Tex),scale = T,center = F))
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
pdf("z49_migration-heatmap-group1-CD8_Tex.pdf", width = 4, height = 4)
print(p)
dev.off()

CD8_Te <- group1_means_celltype[c("CD8-Te-Health","CD8-Te-Mild","CD8-Te-Severe")]
df <- data.frame(colnames(CD8_Te))
colnames(df) <- 'class'
top_anno = HeatmapAnnotation(df = df,
                             border = T,
                             show_annotation_name = F,
                             gp = gpar(col = 'black'),
                             col = list(class = c('CD8-Te-Health'="#00AFBB",
                                                  'CD8-Te-Mild'="#E7B800",
                                                  'CD8-Te-Severe'="#FC4E07")))#颜色设置
marker_exp <- t(scale(t(CD8_Te),scale = T,center = F))
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
pdf("z49_migration-heatmap-group1-CD8_Te.pdf", width = 3.6, height = 4)
print(p)
dev.off()

apoptosis <- c("TNFSF10",	"TNFSF10",	"TNFSF12",	"TNFSF15",	"FAS",	"FASLG",	"FADD", "TRADD",	"CASP8",	"IRF1",	"TP53",	"BCL2L11",	"XAF1", "CASP3")
group1_means <- AggregateExpression(
  object = CD8_T_cells,
  features = apoptosis,
  group.by = "group1",
  slot = "data"
)
group1_means <- as.data.frame(group1_means$RNA)
df <- data.frame(colnames(group1_means))
colnames(df) <- 'class'
top_anno = HeatmapAnnotation(df = df,
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
pdf("z57_apoptosis-heatmap-group1.pdf", width = 3, height = 3)
print(p)
dev.off()

group1_means_celltype <- AggregateExpression(
  object = CD8_T_cells,
  features = apoptosis,
  group.by = "celltype.group",
  slot = "data"
)
group1_means_celltype <- as.data.frame(group1_means_celltype$RNA)

CD8_Naive <- group1_means_celltype[c("CD8-Naive-Health","CD8-Naive-Mild","CD8-Naive-Severe")]
df <- data.frame(colnames(CD8_Naive))
colnames(df) <- 'class'
top_anno = HeatmapAnnotation(df = df,
                             border = T,
                             show_annotation_name = F,
                             gp = gpar(col = 'black'),
                             col = list(class = c('CD8-Naive-Health'="#00AFBB",
                                                  'CD8-Naive-Severe'="#FC4E07",
                                                  'CD8-Naive-Mild'="#E7B800")))#颜色设置
marker_exp <- t(scale(t(CD8_Naive),scale = T,center = F))
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
pdf("z57_apoptosis-heatmap-group1-CD8_Naive.pdf", width = 3.8, height = 3)
print(p)
dev.off()

CD8_Tem <- group1_means_celltype[c("CD8-Tem-Health","CD8-Tem-Mild","CD8-Tem-Severe")]
df <- data.frame(colnames(CD8_Tem))
colnames(df) <- 'class'
top_anno = HeatmapAnnotation(df = df,
                             border = T,
                             show_annotation_name = F,
                             gp = gpar(col = 'black'),
                             col = list(class = c('CD8-Tem-Health'="#00AFBB",
                                                  'CD8-Tem-Mild'="#E7B800",
                                                  'CD8-Tem-Severe'="#FC4E07")))#颜色设置
marker_exp <- t(scale(t(CD8_Tem),scale = T,center = F))
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
pdf("z57_apoptosis-heatmap-group1-CD8_Tem.pdf", width = 3.7, height = 3)
print(p)
dev.off()

CD8_Tex <- group1_means_celltype[c("CD8-Tex-Health","CD8-Tex-Mild","CD8-Tex-Severe")]
df <- data.frame(colnames(CD8_Tem_exh))
colnames(df) <- 'class'
top_anno = HeatmapAnnotation(df = df,
                             border = T,
                             show_annotation_name = F,
                             gp = gpar(col = 'black'),
                             col = list(class = c('CD8-Tex-Health'="#00AFBB",
                                                  'CD8-Tex-Mild'="#E7B800",
                                                  'CD8-Tex-Severe'="#FC4E07")))#颜色设置
marker_exp <- t(scale(t(CD8_Tex),scale = T,center = F))
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
pdf("z57_apoptosis-heatmap-group1-CD8_Tex.pdf", width = 4, height = 3)
print(p)
dev.off()

CD8_Te <- group1_means_celltype[c("CD8-Te-Health","CD8-Te-Mild","CD8-Te-Severe")]
df <- data.frame(colnames(CD8_Te))
colnames(df) <- 'class'
top_anno = HeatmapAnnotation(df = df,
                             border = T,
                             show_annotation_name = F,
                             gp = gpar(col = 'black'),
                             col = list(class = c('CD8-Te-Health'="#00AFBB",
                                                  'CD8-Te-Mild'="#E7B800",
                                                  'CD8-Te-Severe'="#FC4E07")))#颜色设置
marker_exp <- t(scale(t(CD8_Te),scale = T,center = F))
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
pdf("z57_apoptosis-heatmap-group1-CD8_Te.pdf", width = 3.6, height = 3)
print(p)
dev.off()

library(Seurat)
library(dplyr)
library(writexl)

gene_expression <- GetAssayData(CD8_T_cells, slot = "data", assay = "RNA")[exhaustion, ]

group_info <- CD8_T_cells$group1
gene_df <- data.frame(Group = group_info)
gene_df <- cbind(gene_df, t(gene_expression))

gene_summary_list <- list()
for (gene in exhaustion) {
  gene_df$Gene_Expression <- gene_df[[gene]]
  severe_df <- gene_df %>% filter(Group == "Severe")
  gene_summary <- severe_df %>%
    summarise(
      Group = "Severe",
      Gene = gene,
      Mean = mean(Gene_Expression, na.rm = TRUE),
      Variance = var(Gene_Expression, na.rm = TRUE)
    )%>%
    ungroup() 
  gene_summary_list[[gene]] <- gene_summary
}

all_gene_summary <- do.call(rbind, gene_summary_list)

write_xlsx(all_gene_summary, "Multiple_Genes_Severe_Summary.xlsx")

for (gene in exhaustion) {
  gene_df$Gene_Expression <- gene_df[[gene]]
  severe_df <- gene_df %>% filter(Group == "Mild")
  gene_summary <- severe_df %>%
    summarise(
      Group = "Mild",
      Gene = gene,
      Mean = mean(Gene_Expression, na.rm = TRUE),
      Variance = var(Gene_Expression, na.rm = TRUE)
    )%>%
    ungroup()
  gene_summary_list[[gene]] <- gene_summary
}

all_gene_summary2 <- do.call(rbind, gene_summary_list)

write_xlsx(all_gene_summary2, "Multiple_Genes_Mild_Summary.xlsx")

for (gene in exhaustion) {
  gene_df$Gene_Expression <- gene_df[[gene]]
  severe_df <- gene_df %>% filter(Group == "Health")
  gene_summary <- severe_df %>%
    summarise(
      Group = "Health",
      Gene = gene,
      Mean = mean(Gene_Expression, na.rm = TRUE),
      Variance = var(Gene_Expression, na.rm = TRUE)
    )%>%
    ungroup()
  gene_summary_list[[gene]] <- gene_summary
}

all_gene_summary3 <- do.call(rbind, gene_summary_list)

write_xlsx(all_gene_summary3, "Multiple_Genes_Health_Summary.xlsx")

gene_expression <- GetAssayData(CD8_T_cells, slot = "data", assay = "RNA")[cytotoxicity, ]

group_info <- CD8_T_cells$group1
gene_df <- data.frame(Group = group_info)
gene_df <- cbind(gene_df, t(gene_expression))

gene_summary_list <- list()
for (gene in cytotoxicity) {
  gene_df$Gene_Expression <- gene_df[[gene]]
  severe_df <- gene_df %>% filter(Group == "Severe")
  gene_summary <- severe_df %>%
    summarise(
      Group = "Severe",
      Gene = gene,
      Mean = mean(Gene_Expression, na.rm = TRUE),
      Variance = var(Gene_Expression, na.rm = TRUE)
    )%>%
    ungroup()
  gene_summary_list[[gene]] <- gene_summary
}

all_gene_summary <- do.call(rbind, gene_summary_list)

write_xlsx(all_gene_summary, "Multiple_Genes_Severe_Summary.xlsx")

for (gene in cytotoxicity) {
  gene_df$Gene_Expression <- gene_df[[gene]]
  severe_df <- gene_df %>% filter(Group == "Mild")
  gene_summary <- severe_df %>%
    summarise(
      Group = "Mild",
      Gene = gene,
      Mean = mean(Gene_Expression, na.rm = TRUE),
      Variance = var(Gene_Expression, na.rm = TRUE)
    )%>%
    ungroup()
  gene_summary_list[[gene]] <- gene_summary
}

all_gene_summary2 <- do.call(rbind, gene_summary_list)

write_xlsx(all_gene_summary2, "Multiple_Genes_Mild_Summary.xlsx")

IFNG <- FetchData(CD8_T_cells, vars = "IFNG")
IFNG_data <- data.frame(IFNG = IFNG, 
                        group = CD8_T_cells$group1, 
                        celltype = CD8_T_cells$celltype)
IFNG_data$group <- factor(IFNG_data$group, levels = c("Health", "Mild", "Severe"))

p <- ggviolin(IFNG_data, x = "group", y = "IFNG", fill = "group",
              palette = c("#00AFBB", "#E7B800", "#FC4E07"),  # 设置颜色
              ylab = "Average Expression", xlab = "Group",
              add = "boxplot", add.params = list(fill = "white"))  # 添加箱线图
p
p <- p + stat_compare_means(method = "wilcox.test", comparisons = my_comparisons, label = "p.signif")
p
ggsave("z34_IFNG-group1.pdf",plot = p,width = 4,height = 6)

p <- ggviolin(IFNG_data, x = "group", y = "IFNG", fill = "group",
              palette = c("#00AFBB", "#E7B800", "#FC4E07"),  # 设置颜色
              ylab = "Average Expression", xlab = "Group",
              add = "boxplot", add.params = list(fill = "white", width = 0.5, size = 1)) +  # 添加箱线图
  facet_wrap(~ celltype, scales = "free")  # 按celltype分面，并允许每个分面有不同的y轴范围
p
# 假设你想在每个细胞类型中比较group1的分组差异
p <- p + stat_compare_means(method = "wilcox.test", label = "p.signif", comparisons = my_comparisons)
p
ggsave("z34_IFNG-group1-celltype.pdf",plot = p,width = 12,height = 6)