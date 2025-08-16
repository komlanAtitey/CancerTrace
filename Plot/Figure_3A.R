
library(ggplot2)
library(dplyr)
library(Seurat)

# Load required data
load("tsne_reduced_data/tsne.lung.meta.rdata")
load("gene_marker_data/find.gene.markers.meta.rdata")

load("tsne_data/tsne.lung.meta.rdata")
load("lung_cancer/gene_marker_data/find.gene.markers.meta.rdata")

load("tsne_data/tsne.lung.meta.rdata")
load("gene_marker_data/find.gene.markers.meta.rdata")


# Target marker gene
target_gene <- "EPCAM"

# Ensure the gene exists in the dataset
if (!target_gene %in% rownames(tsne.lung.meta)) {
  stop(paste("Gene", target_gene, "not found in Seurat object."))
}

# Extract tSNE coordinates, cluster identity, and EPCAM expression
expr_data <- FetchData(tsne.lung.meta, vars = c("ident", target_gene, "tSNE_1", "tSNE_2"))

# Compute mean EPCAM per cluster
cluster_means <- expr_data %>%
  group_by(ident) %>%
  summarise(mean_expr = mean(.data[[target_gene]], na.rm = TRUE))

# Define high expression threshold (e.g., top 25% of clusters)
threshold <- quantile(cluster_means$mean_expr, 0.75)

# Select clusters with high EPCAM expression
high_epcam_clusters <- cluster_means %>%
  filter(mean_expr >= threshold) %>%
  pull(ident)

# Add flag for high EPCAM clusters
expr_data <- expr_data %>%
  mutate(high_epcam_cluster = ident %in% high_epcam_clusters)

# Create convex hull for each high EPCAM cluster
hull_data <- expr_data %>%
  filter(high_epcam_cluster) %>%
  group_by(ident) %>%
  slice(chull(tSNE_1, tSNE_2))

# Plot t-SNE with per-cluster EPCAM outlines
ggplot(expr_data, aes(x = tSNE_1, y = tSNE_2)) +
  geom_point(aes(color = .data[[target_gene]]), size = 1) +
  scale_color_viridis_c(option = "plasma", na.value = "lightgray") +
  #geom_polygon(
  #  data = hull_data,
  #  aes(x = tSNE_1, y = tSNE_2, fill = ident, group = ident),
  #  color = "red", alpha = 0.2, inherit.aes = FALSE
  #) +
  guides(fill = "none") +
  ggtitle("EPCAM Expression with High-Expression Cluster Outlines") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line=element_line(color = "black"),
    panel.border = element_blank()
  )