
driver.Epithelial <- Epithelial_cell.Dr$gene_id
best.driver.Epithelial <- driver.Epithelial[1:20]

######### 
######### Epithelial
#########
load("old_data/Epithelial_cell.Dr.rdata")
row.genes.data.Epithelial <- sapply(1:dim(data.frame(best.driver.Epithelial))[1],
                                    function(i) which(Epithelial_cell.Dr == best.driver.Epithelial[i]))
genes.data.Epithelial <- sapply(1:dim(data.frame(best.driver.Epithelial))[1],
                                function(i) Epithelial_cell.Dr[row.genes.data.Epithelial[i],])
genes.data.Epithelial <- data.frame(genes.data.Epithelial)
genes.data.Epithelial <- t(genes.data.Epithelial)
genes.data.Epithelial <- data.frame(genes.data.Epithelial)
data.Epithelial <- as.numeric(genes.data.Epithelial$mean_coef)

######### 
######### B cells
#########
load("old_data/B.lymphocytes.Dr.rdata")
row.genes.data.B.lymphocytes <- sapply(1:dim(data.frame(best.driver.Epithelial))[1],
                                       function(i) which(B.lymphocytes.Dr == best.driver.Epithelial[i]))
genes.data.B.lymphocytes <- sapply(1:dim(data.frame(best.driver.Epithelial))[1],
                                   function(i) B.lymphocytes.Dr[row.genes.data.B.lymphocytes[i],])
genes.data.B.lymphocytes <- data.frame(genes.data.B.lymphocytes)
genes.data.B.lymphocytes <- t(genes.data.B.lymphocytes)
genes.data.B.lymphocytes <- data.frame(genes.data.B.lymphocytes)
data.B.lymphocytes <- as.numeric(genes.data.B.lymphocytes$mean_coef)

######### 
######### T cells
#########
load("old_data/T.lymphocytes.Dr.rdata")
row.genes.data.T.lymphocytes <- sapply(1:dim(data.frame(best.driver.Epithelial))[1],
                                       function(i) which(T.lymphocytes.Dr == best.driver.Epithelial[i]))
genes.data.T.lymphocytes <- sapply(1:dim(data.frame(best.driver.Epithelial))[1],
                                   function(i) T.lymphocytes.Dr[row.genes.data.T.lymphocytes[i],])
genes.data.T.lymphocytes <- data.frame(genes.data.T.lymphocytes)
genes.data.T.lymphocytes <- t(genes.data.T.lymphocytes)
genes.data.T.lymphocytes <- data.frame(genes.data.T.lymphocytes)
data.T.lymphocytes <- as.numeric(genes.data.T.lymphocytes$mean_coef)

######### 
######### DC
#########
load("old_data/DC.Dr.rdata") 
row.genes.data.DC <- sapply(1:dim(data.frame(best.driver.Epithelial))[1],
                            function(i) which(DC.Dr == best.driver.Epithelial[i]))
genes.data.DC <- sapply(1:dim(data.frame(best.driver.Epithelial))[1],
                        function(i) DC.Dr[row.genes.data.DC[i],])
genes.data.DC <- data.frame(genes.data.DC)
genes.data.DC <- t(genes.data.DC)
genes.data.DC <- data.frame(genes.data.DC)
data.DC <- as.numeric(genes.data.DC$mean_coef)

######### 
######### Endothelial cells
#########
load("old_data/Endothelial_cells.Dr.rdata") 
row.genes.data.Endothelial <- sapply(1:dim(data.frame(best.driver.Epithelial))[1],
                                     function(i) which(Endothelial_cells.Dr == best.driver.Epithelial[i]))
genes.data.Endothelial<- sapply(1:dim(data.frame(best.driver.Epithelial))[1],
                                function(i) Endothelial_cells.Dr[row.genes.data.Endothelial[i],])
genes.data.Endothelial <- data.frame(genes.data.Endothelial)
genes.data.Endothelial <- t(genes.data.Endothelial)
genes.data.Endothelial <- data.frame(genes.data.Endothelial)
data.Endothelial <- as.numeric(genes.data.Endothelial$mean_coef)

######### 
######### Ependymal cells
#########
load("old_data/Ependymal_cells.Dr.rdata") 
row.genes.data.Ependymal <- sapply(1:dim(data.frame(best.driver.Epithelial))[1],
                                   function(i) which(Ependymal_cells.Dr == best.driver.Epithelial[i]))
genes.data.Ependymal<- sapply(1:dim(data.frame(best.driver.Epithelial))[1],
                              function(i) Ependymal_cells.Dr[row.genes.data.Ependymal[i],])
genes.data.Ependymal <- data.frame(genes.data.Ependymal)
genes.data.Ependymal <- t(genes.data.Ependymal)
genes.data.Ependymal <- data.frame(genes.data.Ependymal)
data.Ependymal <- as.numeric(genes.data.Ependymal$mean_coef)

######### 
######### Fibroblasts
#########
load("old_data/Fibroblasts.Dr.rdata") 
row.genes.data.Fibroblasts <- sapply(1:dim(data.frame(best.driver.Epithelial))[1],
                                     function(i) which(Fibroblasts.Dr == best.driver.Epithelial[i]))
genes.data.Fibroblasts <- sapply(1:dim(data.frame(best.driver.Epithelial))[1],
                                 function(i) Fibroblasts.Dr[row.genes.data.Fibroblasts[i],])
genes.data.Fibroblasts <- data.frame(genes.data.Fibroblasts)
genes.data.Fibroblasts <- t(genes.data.Fibroblasts)
genes.data.Fibroblasts <- data.frame(genes.data.Fibroblasts)
data.Fibroblasts <- as.numeric(genes.data.Fibroblasts$mean_coef)

######### 
######### Mast cells
#########
load("old_data/Mast_cells.Dr.rdata") 
row.genes.data.Mast <- sapply(1:dim(data.frame(best.driver.Epithelial))[1],
                              function(i) which(Mast_cells.Dr == best.driver.Epithelial[i]))
genes.data.Mast <- sapply(1:dim(data.frame(best.driver.Epithelial))[1],
                          function(i) Mast_cells.Dr[row.genes.data.Mast[i],])
genes.data.Mast <- data.frame(genes.data.Mast)
genes.data.Mast <- t(genes.data.Mast)
genes.data.Mast <- data.frame(genes.data.Mast)
data.Mast <- as.numeric(genes.data.Mast$mean_coef)

######### 
######### NK cells
#########
load("old_data/NK_cells.Dr.rdata") 
row.genes.data.NK <- sapply(1:dim(data.frame(best.driver.Epithelial))[1],
                            function(i) which(NK_cells.Dr == best.driver.Epithelial[i]))
genes.data.NK <- sapply(1:dim(data.frame(best.driver.Epithelial))[1],
                        function(i) NK_cells.Dr[row.genes.data.NK[i],])
genes.data.NK <- data.frame(genes.data.NK)
genes.data.NK <- t(genes.data.NK)
genes.data.NK <- data.frame(genes.data.NK)
data.NK <- as.numeric(genes.data.NK$mean_coef)

######### 
######### Oligodentrocyte
#########
load("old_data/Oligodentrocytes.Dr.rdata") 
row.genes.data.oligo <- sapply(1:dim(data.frame(best.driver.Epithelial))[1],
                               function(i) which(Oligodentrocytes.Dr == best.driver.Epithelial[i]))
genes.data.oligo <- sapply(1:dim(data.frame(best.driver.Epithelial))[1],
                           function(i) Oligodentrocytes.Dr[row.genes.data.oligo[i],])
genes.data.oligo <- data.frame(genes.data.oligo)
genes.data.oligo <- t(genes.data.oligo)
genes.data.oligo <- data.frame(genes.data.oligo)
data.oligo <- as.numeric(genes.data.oligo$mean_coef)

#####-------------------------- Figure 6B --------------------------#####

library(pheatmap)
Epithelial.driver <- data.frame(Epithelial_cell.Dr)
gene.id <- Epithelial.driver$gene_id[1:dim(data.frame(best.driver.Epithelial))[1]]

gene.order.data <- cbind( data.Epithelial, data.B.lymphocytes, data.T.lymphocytes,
                          data.DC, data.Endothelial, data.Ependymal,
                          data.Fibroblasts, data.Mast, data.NK, data.oligo)

colnames(gene.order.data) <- c("Epithelial", "B.lymphocytes", "T.lymphocytes", #"id",
                               "DC", "Endothelial", "Ependymal",
                               "Fibroblasts", "Mast",  "NK", "Oligodentrocytes")

pheatmap(gene.order.data) 
row.names(gene.order.data) <- gene.id
gene.order.data.1 <- as.matrix(gene.order.data)
heatmap(gene.order.data.1, scale = "none")
heatmap(gene.order.data.1, Rowv = NA, Colv = NA) 
legend(x="right", legend=c("max", "med", "min"),fill=heat.colors(3))
