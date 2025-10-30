cancertrace_algorithm_2 <- function(dataset) {
  
  # ------------------------------------------------------------
  # 0) Input sanitation: allow vector OR 1-row/1-col data.frame/matrix
  # ------------------------------------------------------------
  if (is.data.frame(dataset) || is.matrix(dataset)) {
    if (nrow(dataset) == 1L) dataset <- as.numeric(dataset[1, ])
    else if (ncol(dataset) == 1L) dataset <- as.numeric(dataset[, 1])
    else stop("`dataset` must be a numeric vector or a 1-row/1-col data.frame/matrix.")
  }
  dataset <- as.numeric(dataset)
  if (!length(dataset) || any(!is.finite(dataset))) {
    stop("`dataset` must be a finite numeric vector.")
  }
  
  # ------------------------------------------------------------
  # 1) Light augmentation (keeps your original idea)
  # ------------------------------------------------------------
  arma.data <- dataset * 10
  add_n <- max(0, 17 - length(arma.data))
  if (add_n > 0) {
    set.seed(123)
    rnd <- runif(add_n, min = min(arma.data), max = max(arma.data))
    arma.data <- c(arma.data, rnd)
  }
  
  # ------------------------------------------------------------
  # 2) Hyperparameters
  # ------------------------------------------------------------
  coef <- 1.15               # your original scalar scale
  N    <- 1000               # number of particles
  n    <- length(arma.data)  # time length
  H    <- 0.9                # Hurst exponent
  nu2  <- 1.78               # variance scale
  epsW <- 1e-16              # tiny epsilon for weights/numerics
  minVar <- 1e-12            # floor for variance
  
  # ------------------------------------------------------------
  # 3) Build fractional-Brownian-like covariance (Toeplitz), stabilize
  # ------------------------------------------------------------
  cor_increments <- sapply(1:n, function(j) {
    0.5 * (abs(j + 1)^(2 * H) - 2 * abs(j)^(2 * H) + abs(j - 1)^(2 * H))
  })
  Sigma <- nu2 * stats::toeplitz(cor_increments)
  Sigma <- as.matrix(Matrix::nearPD(Sigma, corr = FALSE)$mat)  # ensure PSD
  
  # ------------------------------------------------------------
  # 4) Draw particle trajectories: n x N
  # ------------------------------------------------------------
  noise.dist <- t(mvtnorm::rmvnorm(n = N, mean = rep(0, n), sigma = Sigma))  # n x N
  noise.dist[!is.finite(noise.dist)] <- 0
  
  # ------------------------------------------------------------
  # 5) Buffers
  # ------------------------------------------------------------
  obs.vec   <- arma.data
  particles <- matrix(0, nrow = n, ncol = N)
  variance  <- matrix(0, nrow = n, ncol = N)
  inv.var   <- matrix(0, nrow = n, ncol = N)
  sir.est   <- numeric(n)
  
  # Safe uniform weights
  w <- rep(1 / N, N)
  
  # ------------------------------------------------------------
  # 6) Scale prior from observations (kept from your idea)
  # ------------------------------------------------------------
  gam.fun        <- abs(obs.vec) / max(coef, epsW)
  arma.coef.var  <- exp(gam.fun)
  arma.coef      <- mean(arma.coef.var)
  
  # ------------------------------------------------------------
  # 7) SIR with robust guards
  # ------------------------------------------------------------
  effective_size <- function(wts) {
    wts <- as.numeric(wts)
    if (any(!is.finite(wts)) || sum(wts) <= 0) return(NA_real_)
    denom <- sum((wts / sum(wts))^2)
    if (!is.finite(denom) || denom <= 0) return(NA_real_)
    1 / denom
  }
  
  for (t in 1:n) {
    # Propagate particles (importance = transition)
    particles <- arma.coef * noise.dist
    
    # Observation variance at time t
    vt <- exp(particles[t, ])
    vt[!is.finite(vt)] <- minVar
    vt <- pmax(vt, minVar)
    variance[t, ] <- vt
    inv.var[t, ]  <- 1 / vt
    
    # Log-likelihood under N(0, vt) for y_t
    ll <- -0.5 * (log(2 * pi) + log(vt) + (obs.vec[t]^2) * inv.var[t, ])
    ll[!is.finite(ll)] <- -1e6  # floor
    
    # Update weights (log-sum-exp stabilized)
    lw <- log(pmax(w, epsW)) + ll
    m  <- max(lw)
    w  <- exp(lw - m)
    s  <- sum(w)
    if (!is.finite(s) || s <= 0) {
      w <- rep(1 / N, N)
    } else {
      w <- w / s
    }
    w[!is.finite(w)] <- 1 / N
    w <- w / sum(w)
    
    # SIR estimate at time t
    sir.est[t] <- sum(w * particles[t, ])
    
    # ESS and resample (guard NA)
    ess <- effective_size(w)
    if (!is.finite(ess)) ess <- N
    if (ess < 0.25 * N) {
      idx <- sample.int(N, size = N, replace = TRUE, prob = pmax(w, epsW))
      particles  <- particles[, idx, drop = FALSE]
      noise.dist <- noise.dist[, idx, drop = FALSE]
      w <- rep(1 / N, N)
    }
  }
  
  # ------------------------------------------------------------
  # 8) Post-processing
  # ------------------------------------------------------------
  sir.est[!is.finite(sir.est)] <- 0
  sir.est[sir.est < 0] <- 0
  
  density.function <- w
  density.function[!is.finite(density.function)] <- 0
  density.function <- density.function / max(sum(density.function), epsW)
  
  density.new <- density.function[density.function > exp(-20)]
  if (!length(density.new)) density.new <- density.function
  
  # ------------------------------------------------------------
  # 9) HDI score (credMass in [0,1])
  # ------------------------------------------------------------
  # Use weights as a proxy sample; if too degenerate, fallback safely
  hdi.range <- tryCatch(
    hdi(density.new, credMass = 0.95),
    error = function(e) c(min(density.new), max(density.new))
  )
  score.min <- min(hdi.range)
  score.max <- max(hdi.range)
  sel <- which(density.new >= score.min & density.new <= score.max)
  if (!length(sel)) sel <- seq_along(density.new)
  hdi.slice <- density.new[sel]
  
  sd_slice <- stats::sd(hdi.slice)
  sco.HDI <- if (!is.finite(sd_slice) || sd_slice == 0) {
    mean(hdi.slice)
  } else {
    mean((hdi.slice / sd_slice)^2)
  }
  
  # ------------------------------------------------------------
  # 10) Return
  # ------------------------------------------------------------
  list(
    sir.est       = sir.est,                       # filtered latent estimate per time
    obs.vec       = obs.vec,                       # observations used
    hdi.range     = c(score.min, score.max),       # HDI endpoints on final weight proxy
    driver.effect = sco.HDI,                       # HDI-based stability/strength score
    density.new   = density.new                    # cleaned final weights
  )
}
