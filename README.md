# CancerTrace: Multi-Stage Single-Cell Analysis of Networked Cancer Evolution for Driver and Modulator Gene Identification 

CancerTrace is a computational framework designed to identify cancer driver genes and their upstream regulators from longitudinal single-cell RNA sequencing (scRNA-seq) datasets. By integrating patient-specific, time-resolved data, CancerTrace enables dynamic mapping of gene regulatory networks

## Features
- 🔎 Identify cancer-originating clusters from longitudinal scRNA-seq data  
- 🧬 Reconstruct time-resolved gene regulatory networks using **Transfer Entropy** + **Sparse Precision Matrix Inference**  
- 📊 Quantify causal influence of non-driver → driver genes with **Bayesian Logistic Modeling**  
- 🛠 Variational Bayesian inference ensures robust predictions  
- 📈 Rank both known and novel driver genes across multiple tumor stages  
- 🌐 Reveal upstream modulators and regulatory hierarchies in tumor progression  

## Dataset
- 9 longitudinal scRNA-seq datasets  
- 3 LUAD patients (Normal → Early → Mid → Late stages)  
- Patient-specific cell states preserved across time  

## Validation
- ✅ Cross-validation and ROC analysis for predictive accuracy  
- 🔄 In silico perturbation for functional relevance  
- 📚 More than half genes matches literature-reported oncogenes & tumor suppressors  


$~~$

## Steps in the CancerTrace Framework 

![](Figure/github.png)

$~~$


# Run CancerTrace

#---------------<br>
**R Packages**<br>
#---------------<br>
libs <- c("MASS", "tidyr", "dplyr", "caret", "viridis", "magrittr", "plyr", "vars", "tidyverse", "reshape2", "zoo")<br>
lapply(libs, library, character.only=TRUE)<br>

#---------------<br>
**load output data of algorithm 1**<br>
#---------------<br>
load("dataset/epithelial.level.time1.rdata")<br>
load("dataset/epithelial.level.time2.rdata")<br>
load("dataset/epithelial.level.time3.rdata")<br>

#---------------<br>
**load genes id**<br>
#---------------<br>
load("gene.id.rdata")<br>

#---------------<br>
**call function of algorithm 2**<br>
#---------------<br>
source("Algorithm/algorithm_2.R")<br>
epithelial.gene.level <- cbind(epithelial.level.time1$level_1, epithelial.level.time2$level_2, epithelial.level.time3$level_3)<br>
epithelial.gene.level <- data.frame(epithelial.gene.level)<br>

epithelial.gene <- lapply(1:dim(epithelial.gene.level)[1], function(w) {algorithm_2(epithelial.gene.level[w,])}) <br>
dr.coef <- sapply(1:dim(epithelial.gene.level)[1], function(m) {epithelial.gene[[m]]$driver.effect})<br>

gene.dr <- cbind(gene.id[1:5,], dr.coef)<br>
colnames(gene.dr) <- c("gene_id", "coef_dr")<br>
gene.dr <- data.frame(gene.dr)<br>
gene.dr <- gene.dr[order(gene.dr$coef_dr, decreasing = TRUE), ]<br>

#---------------<br>
**call function of algorithm 3**<br>
#---------------<br>
source("Algorithm/algorithm_3.R")<br>
num_data <- generate_evolved_matrix(<br>
  df1 = Epithelial.level.time1, col1 = "level_1",<br>
  df2 = Epithelial.level.time2, col2 = "level_2",<br>
  df3 = Epithelial.level.time3, col3 = "level_3",<br>
  noise_sd = 0.1,<br>
  seed = 42<br>
)<br>
genes <- Epithelial.level.time1$gene<br>
rownames(num_data) <- genes<br>
driver_genes <- c("CA12", "LINC01620", "FKRP", "GPX1", "RP11-146F11.1")<br>
non_driver_genes <- setdiff(rownames(num_data), driver_genes)<br>
CIS_matrix <- compute_CIS_matrix(num_data, non_driver_genes, driver_genes)<br> 
top_influencers <- get_top_influencers_per_driver(CIS_matrix, top_n = 5)<br>

#---------------<br>
**validation**<br>
#---------------<br>
likelihood_output <- compute_transformation_likelihood(CIS_matrix, driver_genes, rownames(num_data))<br> 
likelihood_df <- likelihood_output$model_df<br>
auc_mean <- evaluate_model_performance(likelihood_df)<br>

$~~$

## Running the CancerTrace Algorithms
Example run scripts are located in the Run/ folder. )<br>
They demonstrate how to execute both pipelines using example data and visualize key outputs. )<br>

#---------------<br>
**Run Algorithm 2 — Driver Scoring**<br>
#---------------<br>
Script: run/run_algorithm_2.R<br>
This script: <br>
1.	Loads stage-ordered epithelial datasets (time1, time2, time3). <br>
2.	Runs Algorithm 2 on each gene to estimate the driver effect coefficient (dr.coef). <br>
3.	Outputs a ranked gene table and a bar plot of the top 20 genes (known drivers marked with “X”). <br>
Output:
- gene.dr - data frame of genes × driver scores<br>
- Top 20 Drivers.png - visualization of ranked genes<br>

$~~$

#---------------<br>
**Run Algorithm 3 — Dynamic Causal Inference**<br>
#---------------<br>
Script: Run/run_algorithm_3.R<br>
This script: <br>
1.	Loads the same stage-ordered data. <br>
2.	Runs Algorithm 3, which constructs a time-densified expression matrix, infers causal gene-gene influences, and performs knockout Granger-causality analysis. <br>
3.	Produces a summary of top modulators per driver and a bar chart comparing −log10 (p-values) before vs after knockout. <br>
Output: <br>
- CIS_matrix — non-driver → driver causal scores<br>
- top_influencers — top modulators for each driver<br>
- knockout_results.table — combined knockout statistics<br>
- Knockout_Effect.png — plot of original vs knockout significance<br>






