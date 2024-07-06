library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(DoubletFinder)
library(cowplot)
library(clustree)
library(ggraph)


setwd("/Volumes/MV")



## GATHERING DATA TOGETHER
WT1.data<- Read10X(data.dir = "/Volumes/2TB/WT1/")
C1039G.data<- Read10X(data.dir = "/Volumes/2TB/C1039G/")
WT2.data<- Read10X(data.dir = "/Volumes/2TB/WT2/")
mgD.data<- Read10X(data.dir = "/Volumes/2TB/mgD/")


WT1 <- CreateSeuratObject(counts = WT1.data, project = "WT1")
WT2 <- CreateSeuratObject(counts = WT2.data, project = "WT2")
C1039G <- CreateSeuratObject(counts = C1039G.data, project = "C1039G")
mgD <- CreateSeuratObject(counts = mgD.data, project = "mgD")


WT1$orig.ident<-'WT1'
WT2$orig.ident<-'WT2'
mgD$orig.ident<-'mgD'
C1039G$orig.ident<-'C1039G'



WT1<-RenameCells(WT1,add.cell.id = "WT1")
WT2<-RenameCells(WT2,add.cell.id = "WT2")
mgD<-RenameCells(mgD,add.cell.id = "mgD")
C1039G<-RenameCells(C1039G,add.cell.id = "C1039G")



##Merged all samples to MV
MV <- merge(WT1, y = c(WT2,C1039G,mgD), add.cell.ids = c("WT1","WT2","C1039G","mgD"), project = "Mitral_valve")
head(MV@meta.data,5)
tail(MV@meta.data,5)
saveRDS(MV,file = "/Volumes/2TB/MV.rds")



###QC
MV[["percent.mt"]] <- PercentageFeatureSet(MV, pattern = "^mt-")
VlnPlot(MV, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0.000001)
MV <- subset(MV, subset = nCount_RNA > 500 & nCount_RNA < 15000 & percent.mt <10)

###NORMALIZATION
MV2 <- NormalizeData(MV2, normalization.method = "LogNormalize", scale.factor = 10000)
MV2 <- FindVariableFeatures(MV2, selection.method = "vst", nfeatures = 2000)
MV2 <- ScaleData(MV2)
MV2 <- RunPCA(MV2, features = VariableFeatures(object = MV2))
ElbowPlot(MV2,ndims=50)


###CLUSTERING
MV2 <- FindNeighbors(MV2, dims = 1:30)
MV2 <- FindClusters(MV2, resolution = c(seq(0,1,0.1)))
###UMAP 
MV2 <- RunUMAP(MV2, dims = 1:30)
DimPlot(MV2, reduction = "umap",label  = TRUE)


# Create the DimPlot and specify the custom colors
p <- DimPlot(MV2, reduction = "umap", group.by = "orig.ident", pt.size = 0.0005) +
  scale_color_manual(values = c("#d6873b","#347852","#425785","#52a5c1"
  ))

# Print the plot
print(p)

# Close the PDF device
dev.off()











##FINDING ANS SAVING MARKERS
Idents(MV2)<-"seurat_clusters"
MV.markers <- FindAllMarkers(MV2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
head(MV.markers,10)
write.table(MV.markers,file ="MV_markers.tsv",sep = "\t")
saveRDS(MV.markers,file = "MV_markers.rds")






###top5/10 DEGs heatmap
MV.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

top10 <- MV.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
pdf(file = "heatmap_top10.pdf", width=12, height=8)
DoHeatmap(MV2, features = top10$gene,size = 3) + NoLegend() + theme(text = element_text(size = 8))
dev.off()

top5 <- MV.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
pdf(file = "heatmap_top5.pdf", width=8, height=4)
DoHeatmap(MV2, features = top5$gene,size = 3) + NoLegend() + theme(text = element_text(size = 8))
dev.off()

### DotPlot(MV2, features = top10_gene, 
top10_gene<-top10$gene
top10_gene<-unique(top10_gene)
pdf(file = "dotplot_top10.pdf", width=12, height=7)
p<-DotPlot(MV2, features = top10_gene, dot.scale = 3) + RotatedAxis() + FontSize(x.text = 8, y.text = 12) 
p+theme(axis.text.x = element_text(vjust=0.5,angle = 90))
dev.off()


top5_gene<-top5$gene
top5_gene<-unique(top5_gene)
pdf(file = "dotplot_top5.pdf", width=8, height=5)
p<-DotPlot(MV2, features = top5_gene, dot.scale = 3) + RotatedAxis() + FontSize(x.text = 8, y.text = 12) 
p+theme(axis.text.x = element_text(vjust=0.5,angle = 90))
dev.off()




