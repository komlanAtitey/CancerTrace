
source("github/cancertrace_algorithm_2.R")

# ------------------------------------------------------------
# Dependencies
# ------------------------------------------------------------
libs <- c("MASS", "tidyr", "dplyr", "caret", "viridis", "magrittr", "plyr", "forcats",
          "vars", "tidyverse", "reshape2", "zoo", "HDInterval", "mvtnorm", "Matrix", "ggplot2")
lapply(libs, library, character.only=TRUE)

# ------------------------------------------------------------
# load data
# ------------------------------------------------------------
load("epithelial.level.time1.rdata")
load("epithelial.level.time2.rdata")
load("epithelial.level.time3.rdata")
epithelial_level_1 <- data.frame(epithelial.level.time1)
epithelial_level_2 <- data.frame(epithelial.level.time2)
epithelial_level_3 <- data.frame(epithelial.level.time3)
epithelial.gene.level <- cbind(epithelial_level_1$level_1,
                               epithelial_level_2$level_2,
                               epithelial_level_3$level_3)
epithelial.gene.level <- data.frame(epithelial.gene.level)

# ------------------------------------------------------------
# run cancertrace_algorithm_2
# ------------------------------------------------------------
epithelial.gene <- lapply(1:dim(epithelial.gene.level)[1], 
                          function(w) {cancertrace_algorithm_2(epithelial.gene.level[w,])})

dr.coef <- sapply(1:dim(epithelial.gene.level)[1], function(m) {epithelial.gene[[m]]$driver.effect})
dr.coef[is.infinite(dr.coef)] <- 0
gene.id <- epithelial_level_1$id.time1
gene.dr <- cbind(gene.id, dr.coef)
gene.dr <- data.frame(gene.dr)
gene.dr <- gene.dr[order(gene.dr$dr.coef, decreasing = TRUE), ]
head(gene.dr)

# ------------------------------------------------------------
# visualization
# ------------------------------------------------------------
# --- Data prep (top 20 + ordering high→low) ---
df <- gene.dr %>%
  slice_head(n = 20) %>%
  mutate(
    dr.coef = as.numeric(dr.coef),
    gene.id  = fct_reorder(gene.id, dr.coef, .desc = TRUE)  # <<— high at top
  )

# list genes to mark with an "X"
known_drivers <- c("UHRF1","CD82","TRIM44","APC","CHEK1","TP53INP1")
df <- df %>% mutate(is_driver = gene.id %in% known_drivers)

mx  <- max(df$dr.coef, na.rm = TRUE)
pad <- mx * 0.08

# --- Plot ---
ggplot(df, aes(x = dr.coef, y = gene.id)) +
  # left strip of boxes
  geom_tile(aes(x = -pad*0.5, y = gene.id),
            width = pad*0.7, height = 0.9,
            fill = "white", color = "black") +
  # X marks for known drivers
  geom_point(data = subset(df, is_driver),
             aes(x = -pad*0.5, y = gene.id),
             shape = 4, size = 3, stroke = 1, color = "black") +
  # bars
  geom_col(fill = "darkred", width = 0.8) +
  # tiny right-end error bars
  {eps <- 0.001; geom_segment(aes(x = dr.coef - eps, xend = dr.coef + eps,
                                  y = gene.id,       yend = gene.id),
                              color = "red", linewidth = 0.6)} +
  scale_x_continuous(limits = c(-pad, mx * 1.05), expand = c(0, 0)) +
  labs(x = "dr.coef", y = NULL,
       title = "Top 20 genes (high → low) with known drivers marked by X") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 10),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank()
  )






