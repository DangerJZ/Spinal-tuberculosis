library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(viridisLite)
library(stringr)

setwd("D:/1/")

fs = paste0("D:/1/",dir("D:/1/"))
fs

samples = dir("D:/1") %>% str_split_i(pattern = "_",i = 2) %>% unique();
samples

ctr = function(s){  ns = paste0("01_data/",s)  
if(!file.exists(ns))dir.create(ns,recursive = T)}
lapply(samples, ctr)

lapply(fs, function(s){  
  for(i in 1:length(samples)){   
    if(str_detect(s,samples[[i]])){      
      file.copy(s,paste0("01_data/",samples[[i]]))    
      }  
    }
  })


on = paste0("01_data/",dir("01_data/",recursive = T));
on

nn = str_remove(on,"GSM\\d+_sample\\d_");
nn
file.rename(on,nn)

rm(list = ls())
rdaf = "sce.all.Rdata"
if(!file.exists(rdaf)){
f = dir("01_data/")  
  scelist = list() 
  for(i in 1:length(f)){    
    pda <- Read10X(paste0("01_data/",f[[i]]))    
    scelist[[i]] <- CreateSeuratObject(counts = pda,                                        
                                       project = f[[i]],                                      
                                       min.cells = 3,                                      
                                       min.features = 200)    
    print(dim(scelist[[i]]))  
    }  
    sce.all = merge(scelist[[1]],scelist[-1])  
    sce.all = JoinLayers(sce.all)
    save(sce.all,file = rdaf)}
  load(rdaf)
  head(sce.all@meta.data)
  table(sce.all$orig.ident)

  sce.all[["percent.mt"]] <- PercentageFeatureSet(sce.all, pattern = "^MT-")
  sce.all[["percent.rp"]] <- PercentageFeatureSet(sce.all, pattern = "^RP[SL]")
  sce.all[["percent.hb"]] <- PercentageFeatureSet(sce.all, pattern = "^HB[^(P)]")

  head(sce.all@meta.data, 3)

VlnPlot(sce.all,features = c("nFeature_RNA","nCount_RNA","percent.mt"),
          ncol = 3,pt.size = 0, group.by = "orig.ident")

sce.all <- subset(sce.all, subset = nFeature_RNA > 200 & nFeature_RNA < 9000 & percent.mt < 10)

VlnPlot(sce.all,features = c("nFeature_RNA","nCount_RNA","percent.mt"),
        ncol = 3,pt.size = 0, group.by = "orig.ident")
table(sce.all$orig.ident)

all.genes=rownames(sce.all)

g2m_genes <- cc.genes$g2m.genes
g2m_genes <- CaseMatch(search = g2m_genes, match = rownames(sce.all))

s_genes <- cc.genes$s.genes
s_genes <- CaseMatch(search = s_genes, match = rownames(sce.all))

options(future.globals.maxSize = 60000 * 1024^3)
sce.all <- NormalizeData(sce.all) %>% FindVariableFeatures() %>% ScaleData(features = all.genes) %>% RunPCA()

sce.all <- CellCycleScoring(sce.all, g2m.features = g2m_genes, s.features = s_genes)

colnames(sce.all@meta.data)
table(sce.all$Phase)
DimPlot(sce.all, group.by = "Phase")

sce.all <- NormalizeData(sce.all, normalization.method = "LogNormalize", scale.factor = 10000)
sce.all <- FindVariableFeatures(sce.all, selection.method = "vst", nfeatures = 2000)
sce.all <- ScaleData(sce.all, vars.to.regress = c("S.Score", "G2M.Score"))
DimPlot(sce.all, group.by = "Phase")

table(Idents(sce.all))

top10 <- head(VariableFeatures(sce.all), 10)
pdf(file = "y1.pdf", width = 7, height = 6)
VariableFeaturePlot(object = sce.all)
dev.off()

pdf(file = "y2.pdf", width = 7, height = 6)
LabelPoints(plot= VariableFeaturePlot(object = sce.all), points = top10, repel = TRUE)
dev.off()

sce.all <- RunPCA(sce.all, verbose = F)
pdf(file = "y3.pdf", width = 7, height = 6)
DimPlot(object = sce.all, reduction = "pca")
dev.off()

pdf(file = "y4.pdf", width = 10, height = 9)
VizDimLoadings(object = sce.all, dims = 1:4, reduction = "pca", nfeatures = 20)
dev.off()

pdf(file = "y5.pdf", width = 10, height = 9)
DimHeatmap(object = sce.all, dims = 1:15, cells = 500, balanced = TRUE, nfeatures = 30, ncol = 2)
dev.off()

pdf(file = "y6.pdf", width = 10, height = 9)
ElbowPlot(sce.all, ndims = 50)
dev.off()

pct <- sce.all [["pca"]]@stdev / sum(sce.all [["pca"]]@stdev) * 100
pct
cumu <- cumsum(pct)
cumu
co1 = which(cumu > 90 & pct < 5)[1]
co2 = sort(which((pct[1:length(pct)-1]-pct[2:length(pct)]) > 0.1), decreasing = T)[1]+1
pcs = min(co1, co2)
pcs
pcs = 1:13

library(harmony)
sce.all <- RunHarmony(sce.all, group.by.vars="orig.ident", assay.use="SCT", max.iter.harmony = 20)
table(sce.all@meta.data$orig.ident)

library(clustree)
seq = seq(0.1,1,by=0.1)
sce.all <- FindNeighbors(sce.all, dims = pcs)
for (res in seq){
  sce.all = FindClusters(sce.all, resolution = res)
}

library(patchwork)
p8 = clustree(sce.all, prefix = "RNA_snn_res.")+coord_flip()
p = p8+plot_layout(widths = c(3,1))
p
ggsave("y7_T_RNA_snn_res.pdf", p, width = 30, height = 10)

for (i in c(10,11,12,13,14,15,16,17,18,19,20)) {
  sce.all <- FindNeighbors(sce.all,reduction = "harmony", dims = 1:i) %>% FindClusters(resolution = 0.1)
  sce.all <- RunUMAP(sce.all, reduction = "harmony", dims = 1:i)
  plot_i <- print(DimPlot(sce.all,reduction = "umap", label = TRUE, pt.size = 1, repel = TRUE, label.size = 5))
  plotname <- paste("plot_", i, sep = "")
  assign(plotname, plot_i)
  print(plot_i)
}
p = plot_10+plot_11+plot_12+plot_13+plot_14+plot_15+plot_16+plot_17+plot_18+plot_19+plot_20+plot_layout(ncol = 3)
ggsave(p,filename="dim.10-20.pdf",width = 15,height = 16)

for (i in seq) {
  sce.all <- FindNeighbors(sce.all,reduction = "harmony", dims = pcs) %>% FindClusters(resolution = i)
  sce.all <- RunUMAP(sce.all, reduction = "harmony", dims = pcs)
  plot_i <- print(DimPlot(sce.all,reduction = "umap", label = TRUE, pt.size = 1, repel = TRUE, label.size = 5))
  plotname <- paste("plot_", i, sep = "")
  assign(plotname, plot_i)
  print(plot_i)
}
p = plot_0.1+plot_0.2+plot_0.3+plot_0.4+plot_0.5+plot_0.6+plot_0.7+plot_0.8+plot_0.9+plot_1+plot_layout(ncol = 4)
ggsave(p,filename="dim.13,res.0.1-1.pdf",width = 20,height = 12)

Biocols = c("#AB3282","#53A85F","#F1BB72","#D6E7A3","#57C3F3","#E95C59","#E59CC4","#9FA3A8","#CCE0F5",
            "#BD956A","#8C549C","#E0D4CA","#5F3D69","#C5DEBA","#58A4C3","#E4C755","#F3B1A0","#23452F",
            "#E5D2DD","#AA9A59","#E63863","#E39A35","#C1E6F3","#6778AE","#91D0BE","#712820","#DCC1DD",
            "#CCC9E6","#625D9E","#68A180","#968175")
sce.all <- FindNeighbors(sce.all, reduction = "harmony", dims = pcs) %>% FindClusters(resolution = 0.1)
sce.all <- RunUMAP(sce.all, reduction = "harmony", dims = pcs) %>% RunTSNE(dims = pcs, reduction = "harmony")
DimPlot(sce.all,reduction = "umap",cols = Biocols)

pdf(file = "y8.pdf", width = 7, height = 6)
DimPlot(sce.all, reduction = "umap", label = T,cols = Biocols)
dev.off()

pdf(file = "y9.pdf", width = 7, height = 6)
DimPlot(sce.all, reduction = "umap", label = F,group.by = "orig.ident",cols = Biocols)
dev.off()

pdf(file = "y10.pdf", width = 7, height = 6)
DimPlot(sce.all, reduction = "tsne", label = T,cols = Biocols)
dev.off()

pdf(file = "y11.pdf", width = 7, height = 6)
DimPlot(sce.all, reduction = "tsne", label = F,group.by = "orig.ident",cols = Biocols)
dev.off()

sce.all.markers1 <- FindAllMarkers(sce.all,only.pos = TRUE,logfc.threshold = 0.25)
write.csv(sce.all.markers1,file="sce.all_markers1.RNA.csv")

marker_genes <- c("CD3D","CD3E","CD163","CD68","LYZ","CST3","FCGR3B","CSF3R","NT5E","THY1","IGHG1","MZB1","SDC1","EPCAM","KIT","CPA3","CD19","MS4A1","CD79A","CD79B","PECAM1","VWF","ADAMTS4","ACAN","VPREB1","IGLL1","BACH2")

p11 <- FeaturePlot(sce.all,features = marker_genes,ncol = 3)
p11
ggsave("y14.pdf",p11,width = 10,height = 24)

p12 <- DotPlot(sce.all,features = marker_genes) + RotatedAxis()
p12
ggsave("y15.pdf",p12,width = 8,height = 4)

p13 <- VlnPlot(sce.all,features = marker_genes,stack = T,flip = T) + NoLegend()
p13
ggsave("y16.pdf",p13,width = 8,height = 12)

library(SingleR)
library(celldex)
load("ref_Human_all.RData")
testdata=GetAssayData(object = sce.all@assays$RNA)
clusters <- sce.all@meta.data$seurat_clusters
cellpred <- SingleR(test = testdata, ref = ref_Human_all,clusters = clusters, assay.type.test = "logcounts",
                    labels = ref_Human_all@colData@listData[["label.main"]], assay.type.ref = "logcounts")
celltype = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels,stringsAsFactors = F)
sce.all@meta.data$SingleR="NA"
for (i in 1:nrow(celltype)) {
  sce.all@meta.data[which(sce.all$seurat_clusters == celltype$ClusterID[i]), 'SingleR'] <- celltype$celltype[i]
}
library(pheatmap)
p=plotScoreHeatmap(cellpred)
ggsave("B16_SingleR_Heatmap.pdf",p,width = 12,height = 6)

p1 <- DimPlot(sce.all, group.by = "SingleR", label = T, reduction = "tsne")
p2 <- DimPlot(sce.all, group.by = "SingleR", label = T, reduction = "umap")
p <- p1|p2
ggsave("B17_SingleR.pdf",p,width = 13,height = 5)

sce.all$celltype <- recode(sce.all@meta.data$seurat_clusters,
                                 "0" = "T_cells",
                                 "1" = "Monocytes",
                                 "2" = "Neutrophils",
                                 "3" = "Tissue_stem_cells",
                                 "4" = "Plasma",
                                 "5" = "Mast_cells",
                                 "6" = "B_cells",
                                 "7" = "Endothelial_cells",
                                 "8" = "Chondrocytes",
                                 "9" = "Pre_B_cells")

group_info <- data.table::fread("group.csv",header = TRUE)
head(group_info)
metadata <- FetchData(sce.all, "orig.ident")
metadata$cell_id <- rownames(metadata)
metadata <- left_join(x= metadata, y=group_info, by="orig.ident")
rownames(metadata) <- metadata$cell_id
sce.all <- AddMetaData(sce.all, metadata = metadata)
table(sce.all@meta.data$group1)
table(sce.all@meta.data$group2)

Idents(sce.all) <- "celltype"
table(Idents(sce.all))

p12 <- DotPlot(sce.all,features = marker_genes) + RotatedAxis()
p12
ggsave("y17.pdf",p12,width = 10,height = 4)

sce.all.markers3 <- FindAllMarkers(sce.all, only.pos = TRUE, logfc.threshold = 1, min.pct = 0.3)
write.csv(sce.all.markers3, file = "markers.celltype.RNA.csv")

top10 <- sce.all.markers3 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

markers = as.data.frame(top10[,"gene"])
sce.all <- ScaleData(sce.all, features = as.character(unique(markers$gene)))
p = DoHeatmap(sce.all,
              features = as.character(unique(markers$gene)),
              group.by = "celltype")
p
ggsave("B20.pdf", p, width = 20,height = 20,limitsize = FALSE)

allCells = names(Idents(sce.all))
allType = levels(Idents(sce.all))
choose_Cells = unlist(lapply(allType, function(x){
  cgCells = allCells[Idents(sce.all) == x]
  n = min(table(sce.all@meta.data$celltype))
  if (length(cgCells) >= n) {
    cg = sample(cgCells, n)
  } else {
    cg = sample(cgCells, length(cgCells))
  }
  cg
}))

cg_sce = sce.all[, allCells %in% choose_Cells]
table(Idents(cg_sce))

p = DoHeatmap(cg_sce,
              features = as.character(unique(markers$gene)),
              group.by = "celltype")
p
ggsave("B21.pdf", p, width = 15,height = 15)

cell.prop <- as.data.frame(prop.table(table(sce.all@meta.data$celltype,sce.all@meta.data$orig.ident)))
colnames(cell.prop)<-c("cluster", "group", "proportion")

p= ggplot(cell.prop, aes(group,proportion,fill=cluster))+
  geom_bar(stat = "identity",position = "fill")+
  ggtitle("")+
  theme_bw()+
  theme(axis.ticks.length = unit(0.5, "cm"))+
  guides(fill=guide_legend(title = NULL))
p
ggsave("B22.pdf", p, width = 7,height = 4)

cell.prop <- as.data.frame(prop.table(table(sce.all@meta.data$celltype,sce.all@meta.data$group1)))
colnames(cell.prop)<-c("cluster", "group", "proportion")

p= ggplot(cell.prop, aes(group,proportion,fill=cluster))+
  geom_bar(stat = "identity",position = "fill")+
  ggtitle("")+
  theme_bw()+
  theme(axis.ticks.length = unit(0.5, "cm"))+
  guides(fill=guide_legend(title = NULL))
p
ggsave("B24.pdf", p, width = 4,height = 4)

cell.prop <- as.data.frame(prop.table(table(sce.all@meta.data$celltype,sce.all@meta.data$group2)))
colnames(cell.prop)<-c("cluster", "group", "proportion")

p= ggplot(cell.prop, aes(group,proportion,fill=cluster))+
  geom_bar(stat = "identity",position = "fill")+
  ggtitle("")+
  theme_bw()+
  theme(axis.ticks.length = unit(0.5, "cm"))+
  guides(fill=guide_legend(title = NULL))
p
ggsave("B25.pdf", p, width = 4,height = 4)

p=DimPlot(sce.all, split.by = "group1", label = TRUE,cols = Biocols)
p
ggsave("A13.pdf",p,width = 10,height = 4)

p=DimPlot(sce.all, split.by = "group1", label = F,cols = Biocols)
p
ggsave("A13'.pdf",p,width = 10,height = 4)

p=DimPlot(sce.all, split.by = "group2", label = T,cols = Biocols)
p
ggsave("A14.pdf",p,width = 8,height = 4)

p=DimPlot(sce.all, split.by = "group2", label = F,cols = Biocols)
p
ggsave("A14'.pdf",p,width = 8,height = 4)


pdf(file = "A7.pdf", width = 6, height = 4)
DimPlot(sce.all, reduction = "umap", label = T,cols = Biocols)
dev.off()

pdf(file = "A7'.pdf", width = 6, height = 4)
DimPlot(sce.all, reduction = "umap", label = F,cols = Biocols)
dev.off()

pdf(file = "A8.pdf", width = 5, height = 4)
DimPlot(sce.all, reduction = "umap", label = F,group.by = "orig.ident",cols = Biocols)
dev.off()

pdf(file = "A9.pdf", width = 5, height = 4)
DimPlot(sce.all, reduction = "umap", label = F,group.by = "group1",cols = Biocols)
dev.off()

pdf(file = "A9'.pdf", width = 5.5, height = 4)
DimPlot(sce.all, reduction = "umap", label = F,group.by = "group2",cols = Biocols)
dev.off()

pdf(file = "A10.pdf", width = 6, height = 4)
DimPlot(sce.all, reduction = "tsne", label = T,cols = Biocols)
dev.off()

pdf(file = "A10'.pdf", width = 6, height = 4)
DimPlot(sce.all, reduction = "tsne", label = F,cols = Biocols)
dev.off()

pdf(file = "A11.pdf", width = 5, height = 4)
DimPlot(sce.all, reduction = "tsne", label = F,group.by = "orig.ident",cols = Biocols)
dev.off()

pdf(file = "A11.pdf", width = 5, height = 4)
DimPlot(sce.all, reduction = "tsne", label = F,group.by = "group1",cols = Biocols)
dev.off()

pdf(file = "A11'.pdf", width = 5.5, height = 4)
DimPlot(sce.all, reduction = "tsne", label = F,group.by = "group2",cols = Biocols)
dev.off()

ncol(sce.all)
names(sce.all@meta.data)
Idents(sce.all) <- 'orig.ident'
scRNA <- subset(sce.all,downsample = 2000)
Idents(scRNA) <- 'celltype'
ncol(scRNA)

library(CellChat)
cellchat <- createCellChat(object = GetAssayData(scRNA, assay = "RNA", slot = "data"), 
                           meta = scRNA@meta.data, 
                           group.by = "celltype")
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
pdf("I3_NetVisual_overview_all.pdf", width = 10, height = 6)
par(xpd=TRUE)
netVisual_circle(cellchat@net$count,vertex.weight = groupSize, weight.scale = T,
                 label.edge = F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight,vertex.weight = groupSize, weight.scale = T,
                 label.edge = F, title.name = "Interaction weights/strength")
dev.off()

pdf("I4_NetVisual_overview_split.pdf", width = 8, height = 6)
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

pdf("I5_NetVisual_pathways_circle.pdf", width = 8, height = 6)
for (pathways.show in mypathways) {
  par(xpd=TRUE)
  netVisual_aggregate(cellchat,signaling = pathways.show,layout = "circle")
  
}
dev.off()

pdf("I6_NetVisual_pathways_chord.pdf", width = 8, height = 6)
for (pathways.show in mypathways) {
  netVisual_aggregate(cellchat,signaling = pathways.show,layout = "chord")
}
dev.off()

pdf("I7_NetVisual_pathways_heatmap.pdf", width = 6, height = 4)
for (pathways.show in mypathways) {
  par(xpd=TRUE)
  p <- netVisual_heatmap(cellchat,signaling = pathways.show,color.heatmap = "Reds")
  plot(p)
}
dev.off()

pdf("I8_Pathways.pdf",width = 8, height = 6)
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
ggsave("I9_CCI_all.pdf",p,width = 22,height = 60,limitsize = F)

cellchat2 <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
pdf("I15_SignalingRole.pdf",width = 6,height = 4.5)
for (pathways.show in mypathways) {
  netAnalysis_signalingRole_network(cellchat2,signaling = pathways.show,
                                    width = 8,height = 2.5,font.size = 10)
}                                          
dev.off()