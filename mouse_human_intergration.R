

library(Seurat)
library(nichenetr)
library(dplyr)

dir = "/vast/palmer/scratch/mouse_human"
setwd(dir)
mouse <- readRDS(file.path("/vast/palmer/mouse_human", "mouse.rds"))

human <- readRDS(file.path("/vast/palmer/mouse_human", "human.rds"))

gene_trans = rownames(mouse) %>% convert_mouse_to_human_symbols() %>% as.data.frame()
gene_mouse <- as.data.frame(rownames(mouse))
gene_use <- cbind(gene_trans, gene_mouse)
gene_use <- na.omit(gene_use)
colnames(gene_use) <- c('human','mouse')
mouse_data_trans <- subset(mouse,features=gene_use$mouse)

df = data.frame(rownames(mouse_data_trans@assays$RNA@counts), gene_use$human)

rownames(mouse_data_trans@assays$RNA@counts) <- gene_use$human
rownames(mouse_data_trans@assays$RNA@data) <- gene_use$human

mat <- mouse_data_trans@assays$RNA@counts
unique_features <- unique(rownames(mat))
dup_features <- which(duplicated(rownames(mouse_data_trans)) | duplicated(rownames(mouse_data_trans), fromLast = TRUE))

mouse_data_trans1 <- subset(mouse_data_trans, features = -dup_features)

human$species = "human"
mouse_data_trans1$species = "mouse"

mouse_data_trans1$status = mouse_data_trans1$orig.ident
human$status = human$human

shared.features <- intersect(rownames(mouse_data_trans1), rownames(human))
mouse_data_trans2 <- subset(mouse_data_trans1, features = shared.features)
human2 <- subset(human, features = shared.features)

mouse_data_trans = NULL
human = NULL

gc()

ifnb.list <- list(mouse_data_trans2, human2)

#DefaultAssay()

# normalize and identify variable features for each dataset independently
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 1000)
})

#for (i in 1:2) { DefaultAssay(ifnb.list[[i]]) <- 'RNA' }

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = ifnb.list )

immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features= features )

#DefaultAssay(immune.combined) <- "integrated"
# this command creates an 'integrated' data assay
immune.combined <- IntegrateData(anchorset = immune.anchors)

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(immune.combined) <- "integrated"
saveRDS(immune.combined, file = "mouse_human.rds")
saveRDS(immune.anchors, file = "mouse_human_anchors.rds")

# Run the standard workflow for visualization and clustering
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)

immune.combined <- RunHarmony(immune.combined, "species")
immune.combined <- RunUMAP(immune.combined, reduction = "harmony", dims = 1:30)
immune.combined <- FindNeighbors(immune.combined, reduction = "harmony", dims = 1:30)
immune.combined <- FindClusters(immune.combined, resolution = c(0.2))

immune.combined$type <- factor(immune.combined$species, levels = c("mouse", "human"))

# Visualization
p1 <- DimPlot(immune.combined, reduction = "umap", group.by = "species", cols = c("#347852","#d6873b" ), pt.size = 0.0000000000000001 )
p2 <- DimPlot(immune.combined, reduction = "umap", label = FALSE, repel = TRUE, cols = c("#d6873b","#347852","#425785","#9f2b39","#52a5c1",
                                                                                         "#984ea3","#c65341","#4daf4a","#92b8da","#bc9dcc",
                                                                                         "#b5aa82","#de9d3d","#ca8399","#409079","#984ea3"))
pdf(file = "DimPlot_harmoney_celltype.pdf", width=13, height=5)                         
p1 + p2
dev.off()


pdf(file = "DimPlot4_harmoney.pdf", width=12, height=5) 
DimPlot(immune.combined, reduction = "umap", group.by = "species", split.by = "species", cols = c("#347852","#d6873b","#425785","#9f2b39","#52a5c1",
                                                                                                  "#984ea3","#c65341","#4daf4a","#92b8da","#bc9dcc",
                                                                                                  "#b5aa82","#de9d3d","#ca8399","#409079","#984ea3"))

dev.off()



pbmc.markers <- FindAllMarkers(immune.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)


go <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)
write.table(go, "DEG1.txt",sep="\t")
write.table(go$gene, "DEG_gene.txt",sep="\t")

top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
#write.table(top10$gene, file = "top10.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
DoHeatmap(immune.combined, features = as.character(top10$gene)) + scale_fill_gradientn(colors = c("red", "black","green")) + NoLegend()

pdf(file = paste("heatmap2.pdf",sep =""), width=10, height=12)
DoHeatmap(immune.combined, features = as.character(top10$gene)) + scale_fill_gradientn(colors = c("red", "black","green")) + NoLegend()
dev.off()

pdf(file = paste("GENE.pdf",sep =""), width=12, height=12)
FeaturePlot(immune.combined, features = c("COL1A1","COL3A1","COL5A1","COL8A1","ELN","FN1","ACTA2"))
dev.off()

saveRDS(immune.combined,"mouse_human2.rds")









