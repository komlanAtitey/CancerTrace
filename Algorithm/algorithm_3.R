
# ------------------------------------------------------
# DATA REPRESENTATION: Extrapolate gene expression dynamics 
# ------------------------------------------------------

generate_evolved_matrix <- function(df1, df2, df3, col1, col2, col3, noise_sd = 0.1, seed = 42) {

  # Evolve via Markov-like Transition with Noise
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
  
  # Apply Markov evolution twice for each
  data_time_1_1 <- evolve_markov_noise(data_time_1, noise_sd = noise_sd, seed = seed)
  data_time_1_2 <- evolve_markov_noise(data_time_1_1, noise_sd = noise_sd, seed = seed)
  
  data_time_2_1 <- evolve_markov_noise(data_time_2, noise_sd = noise_sd, seed = seed)
  data_time_2_2 <- evolve_markov_noise(data_time_2_1, noise_sd = noise_sd, seed = seed)
  
  data_time_3_1 <- evolve_markov_noise(data_time_3, noise_sd = noise_sd, seed = seed)
  data_time_3_2 <- evolve_markov_noise(data_time_3_1, noise_sd = noise_sd, seed = seed)
  
  # Combine into a matrix
  num_data <- cbind(
    data_time_1, data_time_1_1, data_time_1_2,
    data_time_2, data_time_2_1, data_time_2_2,
    data_time_3, data_time_3_1, data_time_3_2
  )
  
  # Name columns
  colnames(num_data) <- paste0("T", 1:9)
  
  return(num_data)
}


# ------------------------------------------------------
# TIME-VARYING GENE INTERACTION NETWORKS (EDGE DEFINITION)
# ------------------------------------------------------

# Granger causality score function
compute_granger_score <- function(x, y, max_lag = 2) {
  df <- data.frame(x = x, y = y)
  model <- try(VAR(df, p = max_lag, type = "const"), silent = TRUE)
  if (inherits(model, "try-error")) return(0)
  test <- try(causality(model, cause = "x"), silent = TRUE)
  if (inherits(test, "try-error")) return(0)
  score <- -log10(test$Granger$p.value + 1e-10)
  return(score)
}

# ------------------------------------------------------
# CAUSAL INFLUENCE SCORE (CIS)
# ------------------------------------------------------

# Compute CIS matrix between non-drivers and drivers
compute_CIS_matrix <- function(expr_mat, non_drivers, drivers, max_lag = 2) {
  CIS_matrix <- matrix(0, nrow = length(non_drivers), ncol = length(drivers))
  rownames(CIS_matrix) <- non_drivers
  colnames(CIS_matrix) <- drivers
  
  for (i in seq_along(non_drivers)) {
    x <- expr_mat[non_drivers[i], ]
    for (j in seq_along(drivers)) {
      y <- expr_mat[drivers[j], ]
      CIS_matrix[i, j] <- compute_granger_score(x, y, max_lag)
    }
  }
  return(CIS_matrix)
}

# ------------------------------------------------------
# Identify top-k non-driver influencers per driver gene
# ------------------------------------------------------

get_top_influencers_per_driver <- function(CIS_matrix, top_n = number_gene_modulator) {
  results <- list()
  for (driver in colnames(CIS_matrix)) {
    top_df <- data.frame(
      Driver_Gene = driver,
      Non_Driver_Gene = rownames(CIS_matrix),
      Influence_Score = CIS_matrix[, driver]
    ) %>%
      arrange(desc(Influence_Score)) %>%
      slice_head(n = top_n)
    
    results[[driver]] <- top_df
  }
  return(results)
}

# ------------------------------------------------------
# TRANSFORMATION LIKELIHOOD FUNCTION
# ------------------------------------------------------

# Fit logistic regression to predict transformation likelihood
compute_transformation_likelihood <- function(CIS_matrix, driver_genes, all_genes) {
  influence_score <- rep(0, length(all_genes))
  names(influence_score) <- all_genes
  
  influence_score[rownames(CIS_matrix)] <- rowSums(CIS_matrix)
  label <- ifelse(all_genes %in% driver_genes, 1, 0)
  
  model_df <- data.frame(
    Gene = all_genes,
    Label = label,
    Influence = influence_score
  )
  
  fit <- glm(Label ~ Influence, data = model_df, family = "binomial")
  model_df$Transformation_Likelihood <- predict(fit, type = "response")
  return(list(model_df = model_df, model = fit))
}


# ------------------------------------------------------
# OPTIMIZATION FRAMEWORK
# ------------------------------------------------------

# Perform k-fold cross-validation and compute AUC
evaluate_model_performance <- function(model_df, k = 5) {
  library(pROC)    # For computing ROC curves and AUC
  library(caret)   # For creating stratified folds
  
  folds <- createFolds(model_df$Label, k = k, list = TRUE)  # Create k stratified folds based on the 'Label' column
  aucs <- c()  # Initialize an empty vector to store AUCs for each fold
  
  for (i in seq_along(folds)) {  # Loop over each fold
    test_idx <- folds[[i]]  # Get indices for the test set
    train_df <- model_df[-test_idx, ]  # Training data (excluding the test set)
    test_df <- model_df[test_idx, ]    # Testing data
    
    # Skip fold if test set has only one class (ROC requires both positive and negative cases)
    if (length(unique(test_df$Label)) < 2) {
      warning(paste("Skipping fold", i, "due to single-class test set"))
      next
    }
    
    fit <- glm(Label ~ Influence, data = train_df, family = "binomial")  # Fit logistic regression using Influence to predict Label
    preds <- predict(fit, newdata = test_df, type = "response")  # Predict probabilities on the test set
    
    roc_obj <- roc(test_df$Label, preds)  # Compute ROC curve
    aucs <- c(aucs, auc(roc_obj))  # Store AUC for the current fold
  }
  
  # If no valid folds were found (e.g., every test set had only one class), stop with an error
  if (length(aucs) == 0) {
    stop("No valid folds with both classes. Try increasing sample size or stratifying folds.")
  }
  
  return(mean(aucs))  # Return mean AUC across all valid folds
}

# ------------------------------------------------------
# MAIN PIPELINE EXECUTION
# ------------------------------------------------------

# Step 1-3: CIS Computation and Top Influencers
CIS_matrix <- compute_CIS_matrix(num_data, non_driver_genes, driver_genes) 
top_influencers <- get_top_influencers_per_driver(CIS_matrix, top_n = number_gene_modulator)

# Step 4: Likelihood estimation
likelihood_output <- compute_transformation_likelihood(CIS_matrix, driver_genes, rownames(num_data)) 
likelihood_df <- likelihood_output$model_df

# Step 5: Optimization evaluation
auc_mean <- evaluate_model_performance(likelihood_df)
cat("Mean AUC from cross-validation:", auc_mean, "\n")


