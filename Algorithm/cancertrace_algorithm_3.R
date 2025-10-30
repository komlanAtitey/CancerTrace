cancertrace_algorithm_3 <- function(
    data_vector_1,   # Early
    data_vector_2,   # Mid
    data_vector_3,   # Late
    gene_vector,     # gene names (vector or dataframe column)
    driver_genes,    # known/putative drivers
    # ---- evolution / simulation controls ----
    noise_sd       = 0.1,
    seed           = 42,
    inject_causality = TRUE,
    causal_weight  = 0.8,
    # ---- modeling controls ----
    max_lag        = 2,
    top_n          = 5,
    k_folds        = 5,
    knock_sd       = 1e-4
) {
  suppressPackageStartupMessages({
    require(vars)
    require(dplyr)
    require(tidyr)
    require(pROC)
    require(caret)
  })
  
  # --------------------------------------------------------------
  # Sanitize inputs
  # --------------------------------------------------------------
  sanitize_numeric_vector <- function(x, preferred_names = c("level_1","level_2","level_3","value")) {
    if (is.data.frame(x)) {
      have <- intersect(preferred_names, names(x))
      if (length(have) >= 1) {
        x <- x[[have[1]]]
      } else {
        num_cols <- names(x)[sapply(x, is.numeric)]
        if (length(num_cols) == 0) stop("No numeric columns found in provided data.frame.")
        x <- x[[num_cols[1]]]
      }
    }
    as.numeric(x)
  }
  
  sanitize_gene_vector <- function(g) {
    if (is.data.frame(g)) {
      likely <- intersect(c("gene","Gene","symbol","Symbol","id","ID"), names(g))
      if (length(likely) == 0) stop("gene_vector is a data.frame but no obvious gene column was found.")
      g <- g[[likely[1]]]
    }
    as.character(g)
  }
  
  # --------------------------------------------------------------
  # Evolve via Markov transition with optional causal injection
  # --------------------------------------------------------------
  evolve_markov_noise <- function(v, noise_sd = 0.05, seed = NULL) {
    v <- as.numeric(v)
    if (!is.null(seed)) set.seed(seed)
    n <- length(v)
    P <- diag(n) + matrix(rnorm(n^2, 0, noise_sd), nrow = n)
    P <- sweep(P, 1, rowSums(P), `/`)
    evolved_v <- as.vector(v %*% P)
    rng_in  <- range(v, na.rm = TRUE)
    rng_out <- range(evolved_v, na.rm = TRUE)
    if (!all(is.finite(rng_out)) || diff(rng_out) == 0)
      rep(mean(rng_in, na.rm = TRUE), n)
    else
      (evolved_v - rng_out[1]) / diff(rng_out) * diff(rng_in) + rng_in[1]
  }
  
  generate_evolved_matrix_with_causality <- function(v1, v2, v3, noise_sd = 0.1, seed = 42,
                                                     inject_causality = TRUE, causal_weight = 0.8) {
    v1_1 <- evolve_markov_noise(v1, noise_sd, seed)
    v1_2 <- evolve_markov_noise(v1_1, noise_sd, seed)
    v2_1 <- evolve_markov_noise(v2, noise_sd, seed)
    v2_2 <- evolve_markov_noise(v2_1, noise_sd, seed)
    v3_1 <- evolve_markov_noise(v3, noise_sd, seed)
    
    if (inject_causality) {
      set.seed(seed)
      eps <- rnorm(length(v3_1), 0, 0.01)
      v3_2 <- causal_weight * v1_2 + (1 - causal_weight) * v3_1 + eps
      rng <- range(v3, na.rm = TRUE)
      v3_2 <- pmin(pmax(v3_2, rng[1]), rng[2])
    } else {
      v3_2 <- evolve_markov_noise(v3_1, noise_sd, seed)
    }
    
    out <- cbind(v1, v1_1, v1_2, v2, v2_1, v2_2, v3, v3_1, v3_2)
    colnames(out) <- paste0("T", 1:9)
    out
  }
  
  # --------------------------------------------------------------
  #  Granger causality computation
  # --------------------------------------------------------------
  compute_granger_score <- function(x, y, max_lag = 2) {
    df <- data.frame(x = as.numeric(x), y = as.numeric(y))
    model <- try(VAR(df, p = max_lag, type = "const"), silent = TRUE)
    if (inherits(model, "try-error")) return(0)
    test <- try(causality(model, cause = "x"), silent = TRUE)
    if (inherits(test, "try-error")) return(0)
    pval <- tryCatch(test$Granger$p.value, error = function(e) 1)
    -log10(pval + 1e-10)
  }
  
  compute_CIS_matrix <- function(expr_mat, non_drivers, drivers, max_lag = 2) {
    CIS_matrix <- matrix(0, nrow = length(non_drivers), ncol = length(drivers),
                         dimnames = list(non_drivers, drivers))
    for (i in seq_along(non_drivers)) {
      x <- expr_mat[non_drivers[i], ]
      for (j in seq_along(drivers)) {
        y <- expr_mat[drivers[j], ]
        CIS_matrix[i, j] <- compute_granger_score(x, y, max_lag)
      }
    }
    CIS_matrix
  }
  
  # --------------------------------------------------------------
  #  Identify top influencers per driver
  # --------------------------------------------------------------
  get_top_influencers_per_driver <- function(CIS_matrix, top_n = 5) {
    res <- lapply(colnames(CIS_matrix), function(driver) {
      tibble::tibble(
        Driver_Gene = driver,
        Non_Driver_Gene = rownames(CIS_matrix),
        Influence_Score = CIS_matrix[, driver]
      ) |> arrange(desc(Influence_Score)) |> slice_head(n = top_n)
    })
    names(res) <- colnames(CIS_matrix)
    res
  }
  
  # --------------------------------------------------------------
  #  Likelihood + Model Performance
  # --------------------------------------------------------------
  compute_transformation_likelihood <- function(CIS_matrix, driver_genes, all_genes) {
    influence_score <- rep(0, length(all_genes)); names(influence_score) <- all_genes
    influence_score[rownames(CIS_matrix)] <- rowSums(CIS_matrix)
    label <- as.integer(all_genes %in% driver_genes)
    model_df <- data.frame(Gene = all_genes, Label = label, Influence = as.numeric(influence_score))
    fit <- glm(Label ~ Influence, data = model_df, family = "binomial")
    model_df$Transformation_Likelihood <- predict(fit, type = "response")
    list(model_df = model_df, model = fit)
  }
  
  evaluate_model_performance <- function(model_df, k = 5) {
    folds <- caret::createFolds(model_df$Label, k = k, list = TRUE)
    aucs <- c()
    for (i in seq_along(folds)) {
      test_idx <- folds[[i]]
      train_df <- model_df[-test_idx, ]
      test_df  <- model_df[test_idx, ]
      if (length(unique(test_df$Label)) < 2) next
      fit <- glm(Label ~ Influence, data = train_df, family = "binomial")
      preds <- predict(fit, newdata = test_df, type = "response")
      roc_obj <- pROC::roc(test_df$Label, preds, quiet = TRUE)
      aucs <- c(aucs, as.numeric(pROC::auc(roc_obj)))
    }
    list(aucs = aucs, auc_mean = mean(aucs))
  }
  
  # --------------------------------------------------------------
  #  Knockout GC module
  # --------------------------------------------------------------
  compute_gc <- function(non_driver_series, driver_series, max_lag = 2) {
    df <- data.frame(x = as.numeric(non_driver_series),
                     y = as.numeric(driver_series))
    fit <- try(VAR(df, p = max_lag, type = "const"), silent = TRUE)
    if (inherits(fit, "try-error")) return(c(NA_real_, NA_real_))
    tst <- try(causality(fit, cause = "x"), silent = TRUE)
    if (inherits(tst, "try-error")) return(c(NA_real_, NA_real_))
    Fval <- suppressWarnings(tryCatch(tst$Granger$statistic, error = function(e) NA_real_))
    pval <- suppressWarnings(tryCatch(tst$Granger$p.value,  error = function(e) NA_real_))
    c(Fval, pval)
  }
  
  run_knockout_gc <- function(expr_mat, nd_vec, d_vec, max_lag = 2, knock_sd = 1e-4) {
    if (is.null(rownames(expr_mat))) stop("expr_mat must have gene rownames.")
    nd_vec <- nd_vec[nd_vec %in% rownames(expr_mat)]
    d_vec  <- d_vec[d_vec %in% rownames(expr_mat)]
    if (!length(nd_vec) || !length(d_vec))
      stop("No valid non-driver or driver genes found in expression matrix.")
    out <- list(); idx <- 1
    for (nd in nd_vec) {
      for (d in d_vec) {
        stat_orig <- compute_gc(expr_mat[nd, ], expr_mat[d, ], max_lag)
        expr_knock <- expr_mat
        expr_knock[nd, ] <- rnorm(ncol(expr_mat), mean = 0, sd = knock_sd)
        stat_knock <- compute_gc(expr_knock[nd, ], expr_knock[d, ], max_lag)
        out[[idx]] <- data.frame(
          non_driver = nd,
          driver = d,
          F_orig = stat_orig[1],
          p_orig = stat_orig[2],
          F_knock = stat_knock[1],
          p_knock = stat_knock[2],
          logp_orig = -log10(as.numeric(stat_orig[2]) + 1e-12),
          logp_knock = -log10(as.numeric(stat_knock[2]) + 1e-12)
        ); idx <- idx + 1
      }
    }
    df <- do.call(rbind, out)
    list(logp_orig = df$logp_orig, logp_knock = df$logp_knock, table = df)
  }
  
  # --------------------------------------------------------------
  # Main execution
  # --------------------------------------------------------------
  v1 <- sanitize_numeric_vector(data_vector_1)
  v2 <- sanitize_numeric_vector(data_vector_2)
  v3 <- sanitize_numeric_vector(data_vector_3)
  g  <- sanitize_gene_vector(gene_vector)
  
  n <- length(g)
  if (!all(length(v1) == n, length(v2) == n, length(v3) == n))
    stop("Length mismatch between gene_vector and expression data vectors.")
  
  # Build time-augmented matrix
  causality_data <- generate_evolved_matrix_with_causality(
    v1, v2, v3,
    noise_sd = noise_sd, seed = seed,
    inject_causality = inject_causality, causal_weight = causal_weight)
  rownames(causality_data) <- g
  
  # Filter for valid driver/non-driver genes
  driver_genes <- intersect(driver_genes, rownames(causality_data))
  if (length(driver_genes) == 0)
    stop("None of the specified driver genes are present in expression matrix rownames.")
  non_driver_genes <- setdiff(rownames(causality_data), driver_genes)
  if (length(non_driver_genes) == 0)
    stop("No non-driver genes remain after filtering.")
  
  # Run CIS + downstream modules
  CIS_matrix <- compute_CIS_matrix(causality_data, non_driver_genes, driver_genes, max_lag)
  top_influencers <- get_top_influencers_per_driver(CIS_matrix, top_n)
  likelihood_output <- compute_transformation_likelihood(CIS_matrix, driver_genes, g)
  likelihood_df <- likelihood_output$model_df
  auc_result <- evaluate_model_performance(likelihood_df, k = k_folds)
  
  # --------------------------------------------------------------
  # Knockout analysis for ALL drivers
  # --------------------------------------------------------------
  driver_names <- names(top_influencers)
  knockout_results_list <- vector("list", length(driver_names))
  names(knockout_results_list) <- driver_names
  
  for (drv in driver_names) {
    nd_vec <- top_influencers[[drv]]$Non_Driver_Gene
    ko_res <- run_knockout_gc(
      expr_mat = causality_data,
      nd_vec   = nd_vec,
      d_vec    = drv,
      max_lag  = max_lag,
      knock_sd = knock_sd
    )
    knockout_results_list[[drv]] <- ko_res
  }
  
  knockout_table <- do.call(
    rbind,
    lapply(seq_along(knockout_results_list), function(i) {
      df <- knockout_results_list[[i]]$table
      df$Driver_Group <- names(knockout_results_list)[i]
      df
    })
  )
  
  knockout_results <- list(
    per_driver = lapply(knockout_results_list, `[[`, "table"),
    table      = knockout_table,
    logp_orig  = knockout_table$logp_orig,
    logp_knock = knockout_table$logp_knock
  )
  
  # --------------------------------------------------------------
  # Outputs
  # --------------------------------------------------------------
  list(
    causality_data    = causality_data,
    CIS_matrix        = CIS_matrix,
    top_influencers   = top_influencers,
    likelihood_output = likelihood_output,
    likelihood_df     = likelihood_df,
    auc_mean          = auc_result$auc_mean,
    knockout_results  = knockout_results
  )
}