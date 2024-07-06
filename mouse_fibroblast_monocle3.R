library(monocle3)
library(Seurat)
library(SeuratWrappers)
library(ggplot2)

# Load your data
load("Fibroblast")


# Prepare the data
data <- GetAssayData(Fibroblast, assay = "RNA", layer = "counts")
cell_metadata <- Fibroblast@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data, cell_metadata = cell_metadata, gene_metadata = gene_annotation)

# Preprocess and reduce dimensions 
cds <- preprocess_cds(cds, num_dim = 15)
cds <- reduce_dimension(cds)

# Plot UMAP with celltype_all coloring
pdf("plot1.pdf", width = 8, height = 5)
plot_cells(cds, label_groups_by_cluster = FALSE, color_cells_by = "celltype_all")
dev.off()



pdf("plot2.pdf", width = 8, height = 6)
plot_cells(cds, genes = ECM_genes, label_cell_groups = FALSE, show_trajectory_graph = FALSE)
dev.off()

# Cluster the cells and plot by partition
cds <- cluster_cells(cds)
pdf("plot3.pdf", width = 8, height = 5)
plot_cells(cds, color_cells_by = "partition")
dev.off()

# Learn the trajectory graph
cds <- learn_graph(cds)
pdf("plot4.pdf", width = 8, height = 5)
plot_cells(cds, color_cells_by = "celltype_all", label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE)
dev.off()

pdf("plot5.pdf", width = 8, height = 5)
plot_cells(cds, color_cells_by = "partition", label_cell_groups = FALSE, label_leaves = TRUE, label_branch_points = TRUE, graph_label_size = 1.5)
dev.off()

# Order cells in pseudotime
cds <- order_cells(cds)

pdf("plot6.pdf", width = 8, height = 5)
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE, label_branch_points = FALSE, graph_label_size = 1.5)
dev.off()



# Helper function to identify the root principal points
get_earliest_principal_node <- function(cds, time_bin = "130-170") {
  cell_ids <- which(colData(cds)[, "embryo.time.bin"] == time_bin)
  
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids, ]))))]
  
  root_pr_nodes
}

# Order cells in pseudotime using the identified root principal nodes
root_pr_node <- get_earliest_principal_node(cds)
cds <- order_cells(cds, root_pr_nodes = root_pr_node)

# Plot cells in pseudotime
pdf("plot7.pdf", width = 8, height = 6)
plot_cells(cds, color_cells_by = "pseudotime", label_branch_points = FALSE, label_roots = TRUE, label_leaves = FALSE, show_trajectory_graph = FALSE)
dev.off()

pdf("plot8.pdf", width = 8, height = 6)
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE, label_branch_points = FALSE, graph_label_size = 1.5)
dev.off()






