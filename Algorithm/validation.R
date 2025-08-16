# ----------------------------
# Step 1: Simulate synthetic gene expression data
# ----------------------------
generate_evolved_matrix_with_causality <- function(df1, df2, df3, col1, col2, col3, noise_sd = 0.1, seed = 42, causal_weight = 0.8) {
  
  evolve_markov_noise <- function(v, noise_sd = 0.05, seed = NULL) {
    # v: numeric vector (initial state)
    # noise_sd: standard deviation of Gaussian noise added to transition matrix
    # seed: optional, for reproducibility
    
    if (!is.null(seed)) set.seed(seed)
    
    n <- length(v)
    
    # Create noisy transition matrix
    P <- diag(n) + matrix(rnorm(n^2, mean = 0, sd = noise_sd), nrow = n)
    P <- P / rowSums(P)  # Make it stochastic
    
    # Evolve the vector
    evolved_v <- as.vector(v %*% P)
    
    # Rescale evolved_v to match the input range
    input_min <- min(v)
    input_max <- max(v)
    
    evolved_min <- min(evolved_v)
    evolved_max <- max(evolved_v)
    
    if (evolved_max == evolved_min) {
      # Avoid divide-by-zero if constant vector
      evolved_v_scaled <- rep((input_max + input_min) / 2, n)
    } else {
      # Min-max rescaling
      evolved_v_scaled <- (evolved_v - evolved_min) / (evolved_max - evolved_min)
      evolved_v_scaled <- evolved_v_scaled * (input_max - input_min) + input_min
    }
    return(evolved_v_scaled)
  }
  
  # Convert to numeric
  data_time_1 <- as.numeric(df1[[col1]])
  data_time_2 <- as.numeric(df2[[col2]])
  data_time_3 <- as.numeric(df3[[col3]])
  
  # Apply Markov evolution
  data_time_1_1 <- evolve_markov_noise(data_time_1, noise_sd = noise_sd, seed = seed)
  data_time_1_2 <- evolve_markov_noise(data_time_1_1, noise_sd = noise_sd, seed = seed)
  
  data_time_2_1 <- evolve_markov_noise(data_time_2, noise_sd = noise_sd, seed = seed)
  data_time_2_2 <- evolve_markov_noise(data_time_2_1, noise_sd = noise_sd, seed = seed)
  
  data_time_3_1 <- evolve_markov_noise(data_time_3, noise_sd = noise_sd, seed = seed)
  
  # ✅ Inject causal signal from data_time_1_2 into data_time_3_2
  # This creates a synthetic edge: T3 (GeneX) ← T9 (GeneY)
  data_time_3_2 <- causal_weight * data_time_1_2 + (1 - causal_weight) * data_time_3_1 + rnorm(length(data_time_3_1), 0, 0.01)
  
  # Combine into a matrix
  num_data <- cbind(
    data_time_1, data_time_1_1, data_time_1_2,     # T1, T2, T3
    data_time_2, data_time_2_1, data_time_2_2,     # T4, T5, T6
    data_time_3, data_time_3_1, data_time_3_2      # T7, T8, T9
  )
  
  # Name columns
  colnames(num_data) <- paste0("T", 1:9)
  
  return(num_data)
}
causality_data <- generate_evolved_matrix(
  df1 = Epithelial_cells.level.time1, col1 = "level_1",
  df2 = Epithelial_cells.level.time2, col2 = "level_2",
  df3 = Epithelial_cells.level.time3, col3 = "level_3",
  noise_sd = 0.1,
  seed = 42
)
genes <- Epithelial_cells.level.time1$gene
rownames(causality_data) <- genes

# ----------------------------
# Step 2: Define function to compute Granger causality
# ----------------------------
compute_gc <- function(geneA, geneB, expr_data, max_lag = 2) {
  # Create a data frame with time-series for geneA and geneB
  df <- data.frame(t(expr_data[c(geneA, geneB), ]))
  colnames(df) <- c("A", "B")
  
  # Fit VAR model (A and B are time series)
  model <- tryCatch(
    VAR(df, p = max_lag, type = "const"),
    error = function(e) return(NULL)
  )
  if (is.null(model)) return(c(F = NA, p = NA))
  
  # Perform Granger causality test: Does A Granger-cause B?
  test <- tryCatch(
    causality(model, cause = "A"),
    error = function(e) return(NULL)
  )
  if (is.null(test) || is.null(test$Granger)) return(c(F = NA, p = NA))
  
  # Extract test statistics
  stat <- tryCatch(test$Granger$statistic, error = function(e) NA)
  pval <- tryCatch(test$Granger$p.value, error = function(e) NA)
  
  return(c(F = stat, p = pval))
}

# ----------------------------
# Step 3: Define non-driver and driver genes
# ----------------------------
# Non-driver genes: genes that may influence others (potential causes)
non_drivers <- c("FRMD4A", "NKX2-5", "NMNAT1", "SDCBP2-AS1", "RP11-449J10.1")

# Driver genes: target genes (potentially affected)
drivers <- c("RP11-146F11.1")

# ----------------------------
# Step 4: Run perturbation consistency analysis
# ----------------------------
results <- data.frame()

for (nd in non_drivers) {
  for (d in drivers) {
    # Compute Granger causality with original data
    stat_orig <- compute_gc(nd, d, causality_data)
    
    # Knock out the non-driver by replacing its expression with near-zero noise
    expr_knockout <- causality_data
    expr_knockout[nd, ] <- rnorm(ncol(expr_knockout), mean = 0, sd = 1e-4)
    
    # Re-compute Granger causality after knockout
    stat_knock <- compute_gc(nd, d, expr_knockout)
    
    # Collect results
    results <- rbind(results, data.frame(
      non_driver = nd,
      driver = d,
      F_orig = stat_orig[1],
      p_orig = stat_orig[2],
      F_knock = stat_knock[1],
      p_knock = stat_knock[2]
    ))
  }
}

# ----------------------------
# Step 5: Post-process for visualization
# ----------------------------
# Convert p-values to -log10(p-values) for easier visualization
results$logp_orig <- -log10(results$p_orig)
results$logp_knock <- -log10(results$p_knock)