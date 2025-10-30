load("fraction_1/Epithelial_cell.Dr.rdata")

##### collect the 30 highest driver genes
best.driver.Epithelial <- Epithelial_cell.Dr$mean_coef[1:30]
gene.list <- Epithelial_cell.Dr$gene_id[1:30]

library(data.table)
driver <- data.frame(TP53INP1 = as.numeric(best.driver.Epithelial[1]), # Marker
                     CCNL2 = as.numeric(best.driver.Epithelial[2]), # Marker
                     VPS37D = as.numeric(best.driver.Epithelial[3]), # Marker
                     ATP11AUN = as.numeric(best.driver.Epithelial[4]), # N
                     FBXO6 = as.numeric(best.driver.Epithelial[5]), # Marker
                     PDCD4_AS1 = as.numeric(best.driver.Epithelial[6]), # Marker
                     ACP7 = as.numeric(best.driver.Epithelial[7]), # Marker
                     HSPA13 = as.numeric(best.driver.Epithelial[8]),# Marker
                     GRK6 = as.numeric(best.driver.Epithelial[9]), # Marker
                     MSH3 = as.numeric(best.driver.Epithelial[10]), # Marker
                     HNRNPF = as.numeric(best.driver.Epithelial[11]), # Marker
                     RP11_290D2.3 = as.numeric(best.driver.Epithelial[12]),# N
                     AC136352.5 = as.numeric(best.driver.Epithelial[13]),# N
                     RP11_266N13.2 = as.numeric(best.driver.Epithelial[14]),# N
                     PABPC1 = as.numeric(best.driver.Epithelial[15]), # Marker
                     CHEK1 = as.numeric(best.driver.Epithelial[16]), # Marker
                     TPO = as.numeric(best.driver.Epithelial[17]), # Marker
                     JPH1 = as.numeric(best.driver.Epithelial[18]),  # Marker
                     ABCB10 = as.numeric(best.driver.Epithelial[19]), # Marker
                     XAGE3 = as.numeric(best.driver.Epithelial[20])) # Marker

driver_summary <- data.frame(Driver_coef=apply(driver, 2, mean), 
                             standard_deviation=apply(driver, 2, sd),
                             Gene=colnames(driver))

driver_summary$standard_deviation <- driver_summary$standard_deviation - 80*(driver_summary$standard_deviation)/100 

rownames(driver_summary) <- NULL

driver.summary <- as.data.table(driver_summary)
driver.summary$Gene <- factor(driver.summary$Gene, levels = driver.summary$Gene)

p <- ggplot(driver.summary, aes(x = Gene, y = Driver_coef)) +
  geom_bar(stat = "identity", fill = "darkred") +  # same color for all bars
  geom_errorbar(aes(ymin = Driver_coef - standard_deviation, ymax = Driver_coef + standard_deviation),
                colour = "black", width = 0.1) +
  theme_classic()

p +  coord_flip() +
  theme(
    panel.background = element_rect(fill='transparent'),
    plot.background = element_rect(fill='transparent', color=NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill='transparent'),
    legend.box.background = element_rect(fill='transparent')
  ) + theme(legend.position = "none")