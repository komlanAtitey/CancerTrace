
# ------------------------------------------------------------
# Dependencies
# ------------------------------------------------------------
libs <- c("vars", "caret", "pROC", "dplyr", "tidyr", "tibble", "reshape2", "ggplot2")
lapply(libs, library, character.only=TRUE)

source("cancertrace_algorithm_3.R")
# ------------------------------------------------------------
# load data
# ------------------------------------------------------------
load("epithelial.level.time1.rdata")
load("epithelial.level.time2.rdata")
load("epithelial.level.time3.rdata")

# ------------------------------------------------------------
# run cancertrace_algorithm_3
# ------------------------------------------------------------
result <- cancertrace_algorithm_3(
  data_vector_1 = epithelial.level.time1,
  data_vector_2 = epithelial.level.time2,
  data_vector_3 = epithelial.level.time3,
  gene_vector   = epithelial.level.time1$id.time1,
  driver_genes  = c("TP53INP1", "CCNL2", "VPS37D", "ATP11AUN", "FBXO6")
)

# ------------------------------------------------------------
# Inspect results
# ------------------------------------------------------------
names(result$knockout_results$per_driver)   # All drivers processed
head(result$knockout_results$table)         # Combined knockout table
result$knockout_results$logp_orig           # Vector of -log10(p) original
result$knockout_results$logp_knock          # Vector of -log10(p) after knockout
result$auc_mean   

# ------------------------------------------------------------
# visualization
# ------------------------------------------------------------
# Reshape the data for plotting
df_plot <- reshape2::melt(result$knockout_results, id.vars = c("non_driver", "driver"),
                          measure.vars = c("logp_orig", "logp_knock"),
                          variable.name = "Condition",
                          value.name = "-log10(p)")
df_plot <- na.omit(df_plot)
# Rename conditions for plot clarity
df_plot$Condition <- factor(df_plot$Condition, levels = c("logp_orig", "logp_knock"),
                            labels = c("Original", "Knockout"))

# Plot bar chart comparing -log10(p-value) before and after knockout
ggplot(df_plot, aes(x = interaction(non_driver, driver), y = `-log10(p)`, fill = Condition)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(
    values = c(
      "Original" = "blue",
      "Knockout" = "red"
      # Add more if needed
    )
  ) +
  labs(
    x = "Non-driver → Driver Pair",
    y = "-log10(p-value)",
    title = "Effect of Non-driver Knockout on Granger Causality"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    panel.border = element_blank()
  )


