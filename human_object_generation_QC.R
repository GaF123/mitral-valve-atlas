library(Seurat)
library(clustree)
library(ggraph)
library(dplyr)
library(patchwork)
library(ggplot2)
library(DoubletFinder)
library(cowplot)

###load samples 8 noraml 

N1A.data<- Read10X(data.dir = "/Volumes/4TB/MV RNAseq/N1A/filtered_feature_bc_matrix/")
N1A <- CreateSeuratObject(counts = N1A.data, project = "N1A", min.cells = 3, min.features = 200)
N1A@meta.data$human<-"Normal"
N1A@meta.data$sample<-"N1A"

N1P.data<- Read10X(data.dir = "/Volumes/4TB/MV RNAseq/N1P/filtered_feature_bc_matrix/")
N1P <- CreateSeuratObject(counts = N1P.data, project = "N1P", min.cells = 3, min.features = 200)
N1P@meta.data$human<-"Normal"
N1P@meta.data$sample<-"N1P"


N2A.data<- Read10X(data.dir = "/Volumes/4TB/MV RNAseq/N2A/filtered_feature_bc_matrix/")
N2A <- CreateSeuratObject(counts = N2A.data, project = "N2A", min.cells = 3, min.features = 200)
N2A@meta.data$human<-"Normal"
N2A@meta.data$sample<-"N2A"



N2P.data<- Read10X(data.dir = "/Volumes/4TB/MV RNAseq/N2P/")
N2P <- CreateSeuratObject(counts = N2P.data, project = "N2P", min.cells = 3, min.features = 200)
N2P@meta.data$human<-"Normal"
N2P@meta.data$sample<-"N2P"


N3A.data<- Read10X(data.dir = "/Volumes/4TB/MV RNAseq/N3A/filtered_feature_bc_matrix/")
N3A <- CreateSeuratObject(counts = N3A.data, project = "N3A", min.cells = 3, min.features = 200)
N3A@meta.data$human<-"Normal"
N3A@meta.data$sample<-"N3A"



N3P.data<- Read10X(data.dir = "/Volumes/4TB/MV RNAseq/N3P/filtered_feature_bc_matrix/")
N3P <- CreateSeuratObject(counts = N3P.data, project = "N3P", min.cells = 3, min.features = 200)
N3P@meta.data$human<-"Normal"
N3P@meta.data$sample<-"N3P"


N4A.data<- Read10X(data.dir = "/Volumes/4TB/MV RNAseq/N4A/filtered_feature_bc_matrix/")
N4A <- CreateSeuratObject(counts = N4A.data, project = "N4A", min.cells = 3, min.features = 200)
N4A@meta.data$human<-"Normal"
N4A@meta.data$sample<-"N4A"



N4P.data<- Read10X(data.dir = "/Volumes/4TB/MV RNAseq/N4P/filtered_feature_bc_matrix/")
N4P <- CreateSeuratObject(counts = N4P.data, project = "N4P", min.cells = 3, min.features = 200)
N4P@meta.data$human<-"Normal"
N4P@meta.data$sample<-"N4P"



M1.data<- Read10X(data.dir = "/Volumes/4TB/MV RNAseq/M1/filtered_feature_bc_matrix/")
M1 <- CreateSeuratObject(counts = M1.data, project = "M1", min.cells = 3, min.features = 200)
M1@meta.data$human<-"Myxomatous"
M1@meta.data$sample<-"M1"


M2.data<- Read10X(data.dir = "/Volumes/4TB/MV RNAseq/M2/filtered_feature_bc_matrix/")
M2 <- CreateSeuratObject(counts = M2.data, project = "M2", min.cells = 3, min.features = 200)
M2@meta.data$human<-"Myxomatous"
M2@meta.data$sample<-"M2"


M3.data<- Read10X(data.dir = "/Volumes/4TB/MV RNAseq/M3/filtered_feature_bc_matrix/")
M3 <- CreateSeuratObject(counts = M3.data, project = "M3", min.cells = 3, min.features = 200)
M3@meta.data$human<-"Myxomatous"
M3@meta.data$sample<-"M3"


M4.data<- Read10X(data.dir = "/Volumes/4TB/MV RNAseq/M4/filtered_feature_bc_matrix/")
M4 <- CreateSeuratObject(counts = M4.data, project = "M4", min.cells = 3, min.features = 200)
M4@meta.data$human<-"Myxomatous"
M4@meta.data$sample<-"M4"


M5.data<- Read10X(data.dir = "/Volumes/4TB/MV RNAseq/M5/filtered_feature_bc_matrix/")
M5 <- CreateSeuratObject(counts = M5.data, project = "M5", min.cells = 3, min.features = 200)
M5@meta.data$human<-"Myxomatous"
M5@meta.data$sample<-"M5"

M6.data<- Read10X(data.dir = "/Volumes/4TB/MV RNAseq/M6/filtered_feature_bc_matrix/")
M6 <- CreateSeuratObject(counts = M6.data, project = "M6", min.cells = 3, min.features = 200)
M6@meta.data$human<-"Myxomatous"
M6@meta.data$sample<-"M6"



Marfan1A.data<- Read10X(data.dir = "/Volumes/4TB/MV RNAseq/Marfan1A/filtered_feature_bc_matrix/")
Marfan1A <- CreateSeuratObject(counts = Marfan1A.data, project = "Marfan1A", min.cells = 3, min.features = 200)
Marfan1A@meta.data$human<-"Marfan"
Marfan1A@meta.data$sample<-"Marfan1A"


Marfan1P.data<- Read10X(data.dir = "/Volumes/4TB/MV RNAseq/Marfan1P/filtered_feature_bc_matrix/")
Marfan1P <- CreateSeuratObject(counts = Marfan1P.data, project = "Marfan1P", min.cells = 3, min.features = 200)
Marfan1P@meta.data$human<-"Marfan"
Marfan1P@meta.data$sample<-"Marfan1P"



Marfan2A.data<- Read10X(data.dir = "/Volumes/4TB/MV RNAseq/Marfan2A/filtered_feature_bc_matrix/")
Marfan2A <- CreateSeuratObject(counts = Marfan2A.data, project = "Marfan2A", min.cells = 3, min.features = 200)
Marfan2A@meta.data$human<-"Marfan"
Marfan2A@meta.data$sample<-"Marfan2A"


Marfan2P.data<- Read10X(data.dir = "/Volumes/4TB/MV RNAseq/Marfan2P/filtered_feature_bc_matrix/")
Marfan2P <- CreateSeuratObject(counts = Marfan2P.data, project = "Marfan2P", min.cells = 3, min.features = 200)
Marfan2P@meta.data$human<-"Marfan"
Marfan2P@meta.data$sample<-"Marfan2P"



Rheumatic1A.data<- Read10X(data.dir = "/Volumes/4TB/MV RNAseq/Rheumatic1A/filtered_feature_bc_matrix/")
Rheumatic1A <- CreateSeuratObject(counts = Rheumatic1A.data, project = "Rheumatic1A", min.cells = 3, min.features = 200)
Rheumatic1A@meta.data$human<-"Rheumatic"
Rheumatic1A@meta.data$sample<-"Rheumatic1A"






Rheumatic2A.data<- Read10X(data.dir = "/Volumes/4TB/MV RNAseq/Rheumatic2A/filtered_feature_bc_matrix/")
Rheumatic2A <- CreateSeuratObject(counts = Rheumatic2A.data, project = "Rheumatic2A", min.cells = 3, min.features = 200)
Rheumatic2A@meta.data$human<-"Rheumatic"
Rheumatic2A@meta.data$sample<-"Rheumatic2A"


Rheumatic3A.data<- Read10X(data.dir = "/Volumes/4TB/MV RNAseq/Rheumatic3A/filtered_feature_bc_matrix/")
Rheumatic3A <- CreateSeuratObject(counts = Rheumatic3A.data, project = "Rheumatic3A", min.cells = 3, min.features = 200)
Rheumatic3A@meta.data$human<-"Rheumatic"
Rheumatic3A@meta.data$sample<-"Rheumatic3A"




##merge all samples，name it as MV，and preattach on each cells
MV <- merge(N1A, y = c(N1P,N2A,N2P,N3A,N3P,N4A,N4P,M1,M2,M3,M4,M5,M6,Marfan1A,Marfan1P,Marfan2A,Marfan2P,Rheumatic1A,Rheumatic2A,Rheumatic3A), add.cell.ids = c("N1A","N1P","N2A","N2P","N3A","N3P","N4A","N4P","M1","M2","M3","M4","M5","M6","Marfan1A","Marfan1P","Marfan2A","Marfan2P","Rheumatic1A","Rheumatic2A","Rheumatic3A"), project = "MV_all")



###QC
levels(MV)
MV[["percent.mt"]] <- PercentageFeatureSet(MV, pattern = "^MT-")


VlnPlot(MV, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

MV2 <- subset(MV, subset = nCount_RNA > 500 & nCount_RNA < 15000 & percent.mt <10)


###NORMALIZE
MV2_1 <- NormalizeData(MV2, normalization.method = "LogNormalize", scale.factor = 10000)

MV2_1 <- FindVariableFeatures(MV2_1, selection.method = "vst", nfeatures = 2000)




MV2_1 <- ScaleData(MV2_1)

MV2_1 <- RunPCA(MV2_1, features = VariableFeatures(object = MV2_1))


ElbowPlot(MV2_1,ndims=50)


###Harmony

MV = MV %>% RunHarmony("sample", plot_convergence = TRUE)


MV <- FindNeighbors(MV, reduction = "harmony", dims = 1:30)


MV <- FindClusters(MV, resolution = c(seq(0,1,0.1)))


MV <- RunUMAP(MV, reduction = "harmony", dims = 1:30)











###UMAP 
MV2_1 <- RunUMAP(MV2_1, dims = 1:30)



pdf(file = "UMAP.pdf", width=5.3, height=4)
DimPlot(MV2_1, reduction = "umap")
dev.off()





