library(reshape2)
library(vars)      
library(ggplot2)   

# Reshape the data for plotting
df_plot <- reshape2::melt(results, id.vars = c("non_driver", "driver"),
                          measure.vars = c("logp_orig", "logp_knock"),
                          variable.name = "Condition",
                          value.name = "-log10(p)")

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

