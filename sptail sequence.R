library(Seurat)
library(ggplot2)
library(plyr)
library(gridExtra)
library(magrittr)
library(tidyr)
library(raster)
library(OpenImageR)
library(ggpubr)
library(grid)
library(wesanderson)
library(dplyr)



imported_raster=OpenImageR::readImage("S4.jpg")

my_data <- read.table(file = 'S4_stdata.tsv', sep = '\t', header = TRUE, stringsAsFactors=FALSE)

data_filtered <- my_data

#data_filtered <- my_data
##remove genes has less than 10 expression
count <- rowSums(data_filtered[,2:(ncol(data_filtered)-1)])
data_filtered_binary <- data_filtered[,2:ncol(data_filtered)] %>% mutate_all(as.logical)
gene_count <- rowSums(data_filtered_binary)


##UMI Count
region <- 500  #change the x axis maxium
test <- data_filtered %>% separate(X, c("A", "B"),  sep = "x")
df <- data.frame(number=1, c=count)
pdf(file = paste("UMI_prerun.pdf",sep =""), width=8.6, height=8.6)
ggplot(df,aes(x=c),color='blue', xlab="Gene") + 
  geom_histogram(aes(y=..density..),binwidth=region/20, color="black", fill="white",size=1)+ 
  geom_density(alpha=.2, fill="#FF6666",size=1,color ="red") +
  scale_x_continuous(name="UMI",limits = c(0,region)) + 
  scale_y_continuous(name="Density", expand = c(0, 0)) + 
  #xlim(0,4000) +
  #expand_limits(x = 0, y = 0) +
  theme(plot.title = element_text(hjust = 0.5, size = 25, face = "bold"),
        axis.text=element_text(colour="black",size=20),
        axis.title=element_text(colour="black",size=25,face="bold"),
        legend.text=element_text(colour="black",size=20),
        legend.title = element_text(colour="black", size=20, face="bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))
dev.off()

##Gene Count
df <- data.frame(number=1, c=gene_count)
region = 500 #change the x axis maxium
pdf(file = paste("Gene_prerun.pdf",sep =""), width=8.6, height=8.6)
ggplot(df,aes(x=c),color='blue', xlab="Gene") + 
  geom_histogram(aes(y=..density..),binwidth=region/20, color="black", fill="white",size=1)+ 
  geom_density(alpha=.2, fill="#FF6666",size=1,color ="red") +
  scale_x_continuous(name="Gene",limits = c(0,region)) + 
  scale_y_continuous(name="Density", expand = c(0, 0)) + 
#xlim(0,4000) +
  #expand_limits(x = 0, y = 0) +
  theme(plot.title = element_text(hjust = 0.5, size = 25, face = "bold"),
        axis.text=element_text(colour="black",size=20),
        axis.title=element_text(colour="black",size=25,face="bold"),
        legend.text=element_text(colour="black",size=20),
        legend.title = element_text(colour="black", size=20, face="bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))
dev.off()

#imported_raster=OpenImageR::readImage("ventricle.jpg")
g <- rasterGrob(imported_raster, width=unit(1,"npc"), height=unit(1,"npc"), interpolate = FALSE)
#Get the plot information so the image will fill the plot box, and draw it

#UMI heatmap
pdf(file = paste("UMI_heatmap_prerun.pdf",sep =""), width=8.6, height=8.6)
ggplot(test, aes(x = as.numeric(A), y = as.numeric(B), color=count)) +
  #scale_color_gradientn(colours = c("black", "green")) + 
  scale_color_gradientn(colours = c("blue","green", "red"),limits=c(0,1000),
                        oob = scales::squish) +
  ggtitle("UMI") +
  #annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
  guides(colour = guide_colourbar(barwidth = 1, barheight = 30)) +
  geom_point(shape = 15, size = 1.5)+
  expand_limits(x = 0, y = 0) +
  scale_x_continuous(name="X", limits = c(NA, NA), expand = expansion(mult = c(-0.013, -0.013))) +
  scale_y_reverse(name="Y", limits = c(NA, NA), expand = expansion(mult = c(-0.013, 0.008))) +
  coord_equal(xlim=c(0,101),ylim=c(101,1)) +
  theme(plot.title = element_text(hjust = 0.5, size = 25, face = "bold"),
        axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"),
        legend.text=element_text(size=20),
        legend.title = element_blank(),
        #legend.title = element_text(colour="black", size=15, face="bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())
dev.off()

#Gene heatmap
pdf(file = paste("Gene_heatmap_prerun.pdf",sep =""), width=8.6, height=8.6)
ggplot(test, aes(x = as.numeric(A), y = as.numeric(B), color=gene_count)) +
  #scale_color_gradientn(colours = c("black", "green")) + 
  scale_color_gradientn(colours = c("blue","green", "red"),limits=c(0,1000),
                        oob = scales::squish) +
  ggtitle("Gene") +
  #annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
  guides(colour = guide_colourbar(barwidth = 1, barheight = 30)) +
  geom_point(shape = 15, size = 1.5)+
  expand_limits(x = 0, y = 0) +
  scale_x_continuous(name="X", limits = c(NA, NA), expand = expansion(mult = c(-0.013, -0.013))) +
  scale_y_reverse(name="Y", limits = c(NA, NA), expand = expansion(mult = c(-0.013, 0.008))) +
  coord_equal(xlim=c(0,101),ylim=c(101,1)) +
  theme(plot.title = element_text(hjust = 0.5, size = 25, face = "bold"),
        axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"),
        legend.text=element_text(size=20),
        legend.title = element_blank(),
        #legend.title = element_text(colour="black", size=15, face="bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())
dev.off()


#use No1_generate_expression_matrix.R to generate expression_matrix.tsv
data1  <- my_data

#extract the coordinates of each pixel
temp1 <- data1 %>% separate(X, c("A", "B"),  sep = "x")

#repair the strips of the data, here is the col = 34
#the repair was done by averaging the two neighboring columns.
col1 = 34
col2 = 28
row = 15
row2 = 1
row3 = 50







temp1 = data.frame(X=data1$X, temp1)

temp1$A = NULL
temp1$B = NULL

location <- read.table("position.txt", sep =",", header = FALSE, dec =".", stringsAsFactors = F)
x <- as.character(location[1,])
x = x[-1]

data_filtered <- temp1[temp1$X %in% x,]
duplicated(data_filtered$X)
write.table(data_filtered, file = 'Filtered_matrix_correct.tsv', sep = '\t',col.names=TRUE, row.names = FALSE, quote = FALSE)



##read expression matrix
my_data <- read.table(file = 'Filtered_matrix_correct.tsv', sep = '\t', header = TRUE, stringsAsFactors=FALSE)
data_filtered <- my_data


#data_filtered <- my_data
##remove Proteins has less than 10 expression
count <- rowSums(data_filtered[,2:ncol(data_filtered)])
data_filtered_binary <- data_filtered[,2:ncol(data_filtered)] %>% mutate_all(as.logical)
gene_count <- rowSums(data_filtered_binary)


#log_count <- log(count)

##UMI Count
region <- 1000  #change the x axis maxium
test <- data_filtered %>% separate(X, c("A", "B"),  sep = "x")
df <- data.frame(number=1, c=count)
pdf(file = paste("UMI.pdf",sep =""), width=8.6, height=8.6)
ggplot(df,aes(x=c),color='blue', xlab="Gene") + 
  geom_histogram(aes(y=..density..),binwidth=region/20, color="black", fill="white",size=1)+ 
  geom_density(alpha=.2, fill="#FF6666",size=1,color ="red") +
  scale_x_continuous(name="UMI",limits = c(0,region)) + 
  scale_y_continuous(name="Density", expand = c(0, 0)) + 
  #xlim(0,4000) +
  #expand_limits(x = 0, y = 0) +
  theme(plot.title = element_text(hjust = 0.5, size = 25, face = "bold"),
        axis.text=element_text(colour="black",size=20),
        axis.title=element_text(colour="black",size=25,face="bold"),
        legend.text=element_text(colour="black",size=20),
        legend.title = element_text(colour="black", size=20, face="bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))
dev.off()

##Protein Count
df <- data.frame(number=1, c=gene_count)
region = 1000 #change the x axis maxium
pdf(file = paste("Gene.pdf",sep =""), width=8.6, height=8.6)
ggplot(df,aes(x=c),color='blue', xlab="Gene") + 
  geom_histogram(aes(y=..density..),binwidth=region/20, color="black", fill="white",size=1)+ 
  geom_density(alpha=.2, fill="#FF6666",size=1,color ="red") +
  scale_x_continuous(name="Gene",limits = c(0,region)) + 
  scale_y_continuous(name="Density", expand = c(0, 0)) + 
  #xlim(0,4000) +
  #expand_limits(x = 0, y = 0) +
  theme(plot.title = element_text(hjust = 0.5, size = 25, face = "bold"),
        axis.text=element_text(colour="black",size=20),
        axis.title=element_text(colour="black",size=25,face="bold"),
        legend.text=element_text(colour="black",size=20),
        legend.title = element_text(colour="black", size=20, face="bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))
dev.off()

#imported_raster=OpenImageR::readImage("ventricle.jpg")
g <- rasterGrob(imported_raster, width=unit(1,"npc"), height=unit(1,"npc"), interpolate = FALSE)
#Get the plot information so the image will fill the plot box, and draw it

#UMI heatmap
pdf(file = paste("UMI_heatmap.pdf",sep =""), width=8.6, height=8.6)
ggplot(test, aes(x = as.numeric(A), y = as.numeric(B), color=count)) +
  #scale_color_gradientn(colours = c("black", "green")) + 
  scale_color_gradientn(colours = c("blue","green", "red"),limits=c(0,1000),
                        oob = scales::squish) +
  ggtitle("UMI") +
  #annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
  guides(colour = guide_colourbar(barwidth = 1, barheight = 30)) +
  geom_point(shape = 15, size = 1.5)+
  expand_limits(x = 0, y = 0) +
  scale_x_continuous(name="X", limits = c(NA, NA), expand = expansion(mult = c(-0.013, -0.013))) +
  scale_y_reverse(name="Y", limits = c(NA, NA), expand = expansion(mult = c(-0.013, 0.008))) +
  coord_equal(xlim=c(0,101),ylim=c(101,1)) +
  theme(plot.title = element_text(hjust = 0.5, size = 25, face = "bold"),
        axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"),
        legend.text=element_text(size=20),
        legend.title = element_blank(),
        #legend.title = element_text(colour="black", size=15, face="bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())
dev.off()

#Protein heatmap
pdf(file = paste("Gene_heatmap.pdf",sep =""), width=8.6, height=8.6)
ggplot(test, aes(x = as.numeric(A), y = as.numeric(B), color=gene_count)) +
  #scale_color_gradientn(colours = c("black", "green")) + 
  scale_color_gradientn(colours = c("blue","green", "red"),limits=c(0,1000),
                        oob = scales::squish) +
  ggtitle("Gene") +
  #annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
  guides(colour = guide_colourbar(barwidth = 1, barheight = 30)) +
  geom_point(shape = 15, size = 1.5)+
  expand_limits(x = 0, y = 0) +
  scale_x_continuous(name="X", limits = c(NA, NA), expand = expansion(mult = c(-0.013, -0.013))) +
  scale_y_reverse(name="Y", limits = c(NA, NA), expand = expansion(mult = c(-0.013, 0.008))) +
  coord_equal(xlim=c(0,101),ylim=c(101,1)) +
  theme(plot.title = element_text(hjust = 0.5, size = 25, face = "bold"),
        axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"),
        legend.text=element_text(size=20),
        legend.title = element_blank(),
        #legend.title = element_text(colour="black", size=15, face="bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())
dev.off()



#change filename1 to name of txt file you want to load
data1 <- read.table("Filtered_matrix_correct.tsv", header = TRUE, sep = "\t", row.names = 1)

data3 <- data.frame(X = rownames(data1), data1)
temp1 <- data3 %>% separate(X, c("A", "B"),  sep = "x")


#clustering based on SCT
pbmc <- PercentageFeatureSet(pbmc, pattern = "^MT-", col.name = "percent.mt")
pbmc <- SCTransform(pbmc, vars.to.regress = "percent.mt", verbose = FALSE)
pbmc <- RunPCA(pbmc, verbose = FALSE)
pbmc <- RunUMAP(pbmc, dims = 1:10, verbose = FALSE)
pbmc <- FindNeighbors(pbmc, dims = 1:10, verbose = FALSE)
pbmc <- FindClusters(pbmc, resolution=0.8, verbose = FALSE)
DimPlot(pbmc, label = TRUE) + NoLegend()
ident <- Idents(pbmc)
df <- data.frame(ident[])
df1 <-data.frame(X =row.names(df), count= df$ident..)
test <- df1 %>% separate(X, c("A", "B"),  sep = "x")


g <- rasterGrob(imported_raster, width=unit(1,"npc"), height=unit(1,"npc"), interpolate = FALSE)
pdf(file = paste("clustering_SCT_3_dim20.pdf",sep =""), width=8.6, height=8.6)
ggplot(test, aes(x = as.numeric(A), y = as.numeric(B), color=count)) +
  #scale_color_gradientn(colours = c("black", "green")) + 
  #scale_color_gradientn(colours = c("blue","green", "red"),
  #                      oob = scales::squish) +
  ggtitle("UMAP") +
  annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
  geom_point(shape = 15, size = 1.5)+
  expand_limits(x = 0, y = 0) +
  scale_x_continuous(name="X", limits = c(NA, NA), expand = expansion(mult = c(-0.013, -0.013))) +
  scale_y_reverse(name="Y", limits = c(NA, NA), expand = expansion(mult = c(-0.013, 0.008))) +
  coord_equal(xlim=c(0,101),ylim=c(101,1)) +
  theme(plot.title = element_text(hjust = 0.8, size = 25, face = "bold"),
        axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"),
        legend.text=element_text(size=20),
        legend.title = element_blank(),
        #legend.title = element_text(colour="black", size=15, face="bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())
dev.off()


#clustering based on CLR
pbmc <- CreateSeuratObject(matrix1.data, min.cells = 10, project = sample1.name)
pbmc <- NormalizeData(pbmc, normalization.method = 'CLR', margin = 2) 
pbmc <- ScaleData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- RunPCA(pbmc, verbose = FALSE)
pbmc <- RunUMAP(pbmc, dims = 1:10, verbose = FALSE)

pbmc <- FindNeighbors(pbmc, dims = 1:10, verbose = FALSE)
pbmc <- FindClusters(pbmc, resolution=0.8, verbose = FALSE)
DimPlot(pbmc, label = TRUE) + NoLegend()

# pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0, logfc.threshold = 0.02)
# 
# pbmc <-ScaleData(object=pbmc, features = rownames(pbmc))
# 
# top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
# DoHeatmap(pbmc, features = top10$gene) + scale_fill_gradientn(colors = c("red", "black","green"))
# 
# 
# VlnPlot(pbmc, features = c("DCX", "NR2F2"))
# FeaturePlot(pbmc, features = c("DCX", "GAD1", "GAD2", "NR2F2", "PROX1", "SP8", "LHX6", "SOX6", 
#                                "SST","PVALB", "VIP", "LAMP5","PAX6", "SLC17A7","PROX1"),pt.size = 0.1)


ident <- Idents(pbmc)
df <- data.frame(ident[])
df1 <-data.frame(X =row.names(df), count= df$ident..)
test <- df1 %>% separate(X, c("A", "B"),  sep = "x")

pal = ArchR::ArchRPalettes$stallion2[order(as.integer(names(ArchR::ArchRPalettes$stallion2)))]



pdf(file = paste("clustering-CLR2000mincell100.pdf",sep =""), width=8.6, height=8.6)
ggplot(test, aes(x = as.numeric(A), y = as.numeric(B), color=count)) +
  #scale_color_gradientn(colours = c("black", "green")) + 
  #scale_color_gradientn(colours = c("blue","green", "red"),
  #                      oob = scales::squish) +
  #scale_color_manual(values=pal[0:6]) + 
  ggtitle("UMAP") +
  annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
  geom_point(shape = 15, size = 1.5)+
  expand_limits(x = 0, y = 0) +
  scale_x_continuous(name="X", limits = c(NA, NA), expand = expansion(mult = c(-0.013, -0.013))) +
  scale_y_reverse(name="Y", limits = c(NA, NA), expand = expansion(mult = c(-0.013, 0.008))) +
  coord_equal(xlim=c(0,101),ylim=c(101,1)) +
  theme(plot.title = element_text(hjust = 0.5, size = 25, face = "bold"),
        axis.text=element_text(size=20),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title = element_blank(),
        legend.text=element_text(size=20),
        legend.title = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())
dev.off()




