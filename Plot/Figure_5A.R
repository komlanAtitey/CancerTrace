
library(visNetwork)
library(igraph)
library(dplyr)
library(scales)

# A. Prepare edge list
all_edges <- do.call(rbind, top_influencers)
edge_list <- dplyr::select(
  all_edges,
  from = Non_Driver_Gene,
  to = Driver_Gene,
  weight = Influence_Score
)

# B. Create igraph object
g <- graph_from_data_frame(edge_list, directed = TRUE)

# C. Identify node types
non_driver_nodes <- unique(edge_list$from)
driver_nodes <- unique(edge_list$to)

# D. Assign colors
non_driver_color <- "#1f77b4"  # blue for all non-drivers
driver_palette <- hue_pal()(length(driver_nodes))
driver_colors <- setNames(driver_palette, driver_nodes)

# E. Assign group and color to each node
node_df <- data.frame(
  id = V(g)$name,
  label = V(g)$name,
  color = sapply(V(g)$name, function(node) {
    if (node %in% driver_nodes) driver_colors[[node]] else non_driver_color
  }),
  group = sapply(V(g)$name, function(node) {
    if (node %in% driver_nodes) node else "Non-Driver Genes"
  }),
  font.size = 0,#20 
  stringsAsFactors = FALSE
)

# F. Prepare edge data with width scaled
edge_df <- edge_list %>%
  dplyr::mutate(width = scales::rescale(weight, to = c(1, 10))) %>%  # scale edge width
  dplyr::rename(from = from, to = to) %>%
  dplyr::select(from, to, width)

# G. Build custom node legend: non-drivers and drivers
node_legend <- data.frame(
  label = c("Non-Driver Genes", driver_nodes),
  color = c(non_driver_color, driver_colors),
  shape = "dot",
  font.size = 0, #20
  stringsAsFactors = FALSE
)

# H. Build custom edge legend for influence
edge_legend <- data.frame(
  #label = c("Low influence →", "Highe →"),
  label = c(" →", " →"),
  arrows = "to",
  width = c(1, 10),
  color = "black",
  stringsAsFactors = FALSE
)

# I. Plot the interactive network with labeled legend
visNetwork(node_df, edge_df) %>%
  visEdges(arrows = "to", scaling = list(min = 1, max = 10)) %>%
  visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
  visLegend(
    addNodes = node_legend,
    addEdges = edge_legend,
    useGroups = FALSE,  # We manually specify nodes
    position = "right"
  ) %>%
  visLayout(randomSeed = 123)
