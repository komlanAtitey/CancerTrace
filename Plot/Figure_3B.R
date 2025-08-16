library(ggplot2)
library(reshape2)
library(dplyr)

load("cells.fraction.time1.rdata")
load("cells.fraction.time2.rdata")
load("cells.fraction.time3.rdata")

# Combine the time point matrices into one tidy data frame
time_list <- list(
  time1 = cells.fraction.time1,
  time2 = cells.fraction.time2,
  time3 = cells.fraction.time3
)

cell.fraction.patient1 <- bind_rows(lapply(names(time_list), function(t) {
  df <- as.data.frame(time_list[[t]])
  df$CellType <- rownames(df)
  df$Time <- t
  df
}))

# Convert to long format for plotting
cell_long <- melt(cell.fraction.patient1, id.vars = c("CellType", "Time"),
                  variable.name = "Condition", value.name = "Fraction")

# Plot: one facet per cell type, grouped by time
ggplot(cell_long, aes(x = Time, y = Fraction, fill = Condition)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_wrap(~ CellType, scales = "free_y") +
  scale_fill_manual(values = c("darkgreen", "brown")) + # "darkolivegreen", "firebrick3"
  labs(y = "Fraction", x = "Time Point") +
  theme_classic() +
  theme(strip.text = element_text(size = 8))