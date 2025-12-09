# ==============================================================================
# MATRIX A_i CONSTRUCTION
# ==============================================================================

# ==============================================================================
# A_m: Monthly series (no aggregation)
# One row per observed month, one column per latent monthly value
# ==============================================================================

A_m <- function(x_i) {
  Tt  <- length(x_i)
  obs <- which(!is.na(x_i))          # months where the series is observed
  
  A <- matrix(0, nrow = length(obs), ncol = Tt)
  for (k in seq_along(obs)) {
    A[k, obs[k]] <- 1                # identity row (selection)
  }
  
  A
}


# ==============================================================================
# A_q_stock: Quarterly "stock" series
# Here: quarterly value = average of the 3 months in the quarter
#       (change to [0,0,1] if you prefer "last month only")
# ==============================================================================

A_q_stock <- function(x_i) {
  Tt  <- length(x_i)
  obs <- which(!is.na(x_i))          # months where the quarterly value is observed
  
  # keep only quarter-end months (m >= 3); first valid quarter ends at month 3
  obs <- obs[obs >= 3]
  n_q <- length(obs)
  
  A <- matrix(0, nrow = n_q, ncol = Tt)
  for (k in seq_along(obs)) {
    m <- obs[k]
    A[k, (m - 2):m] <- 1 / 3         # average of the 3 months
  }
  
  A
}


# ==============================================================================
# A_q_flow: Quarterly "flow" series
# quarterly value = sum of the 3 months in the quarter
# ==============================================================================

A_q_flow <- function(x_i) {
  Tt  <- length(x_i)
  obs <- which(!is.na(x_i))          # months where the quarterly value is observed
  
  # keep only quarter-end months (m >= 3)
  obs <- obs[obs >= 3]
  n_q <- length(obs)
  
  A <- matrix(0, nrow = n_q, ncol = Tt)
  for (k in seq_along(obs)) {
    m <- obs[k]
    A[k, (m - 2):m] <- 1             # sum over the 3 months
  }
  
  A
}

# ==============================================================================
# A_list: build A_i for all series
#
# X_na   : T x N data matrix with NAs
# N_m    : number of purely monthly series (first N_m columns)
# agg_q  : length N_q vector (N_q = N - N_m) with quarterly types:
#          "1" = stock, "2" = flow  (your coding)
# ==============================================================================

A_list <- function(X_na, N_q, agg_q) {
  
  Tt <- nrow(X_na)
  N  <- ncol(X_na)
  
  N_m <- N - N_q
  
  if (length(agg_q) != N_q) {
    stop("Length of 'agg_q' must be equal to number of quarterly series (N - N_m).")
  }
  
  A_list <- vector("list", N)   # one A_i per series
  
  # -------- MONTHLY SERIES (1,...,N_m) ----------
  for (i in 1:N_m) {
    x_i        <- X_na[, i]
    A_list[[i]] <- A_m(x_i)
  }
  
  # -------- QUARTERLY SERIES (N_m+1,...,N) ------
  for (j in 1:N_q) {
    idx <- N_m + j              # column index in X_na
    x_i <- X_na[, idx]
    
    type_j <- agg_q[j]
    
    if (type_j == "2") {
      # flow: sum over 3 months
      A_list[[idx]] <- A_q_flow(x_i)
    } else if (type_j == "1") {
      # stock: average over 3 months (or last month, if you change A_q_stock)
      A_list[[idx]] <- A_q_stock(x_i)
    } else {
      stop("Unknown quarterly type code in agg_q[", j, "]: ", type_j)
    }
  }
  
  A_list
}



# ==============================================================================
# EM ALGORITHM 
# ==============================================================================

# ==============================================================================
# 1. E-STEP: reconstruct completed X given (F_hat, Phi_hat, A_list)
# X_obs : T x N with NA in missing entries (fixed, raw data)
# X_old : T x N current completed panel (previous iteration)
# F_hat : T x r
# Phi_hat : N x r
# A_list : list of length N; A_list[[i]] is the aggregation/selection matrix for series i
# ==============================================================================

E_step <- function(X_old, X_obs, F_hat, Phi_hat, A_list) {
  
  Tt <- nrow(X_old)
  N  <- ncol(X_old)
  
  X_new <- matrix(0, nrow = Tt, ncol = N)
  
  for (i in seq_len(N)) {
    A_i <- A_list[[i]]          # rows = observed values, cols = T months
    
    # Factor loadings for series i (r x 1)
    phi_i <- matrix(Phi_hat[i, ], ncol = 1)
    
    # Predicted monthly path from factors (T x 1)
    pred_i <- F_hat %*% phi_i
    
    # If no observations at all for this series:
    if (nrow(A_i) == 0) {
      X_new[, i] <- pred_i
      next
    }
    
    # Indices of observed entries in the original data
    idx_obs_i <- which(!is.na(X_obs[, i]))
    
    # Observed values (monthly or quarterly) as a column vector
    x_obs_i <- matrix(X_obs[idx_obs_i, i], ncol = 1)
    
    # Check: length(x_obs_i) == nrow(A_i)  (should hold by construction)
    # stopifnot(length(idx_obs_i) == nrow(A_i))
    
    # Matrix to invert: A_i A_i'
    M <- A_i %*% t(A_i)
    
    # Constraint correction
    correction <- t(A_i) %*% solve(M) %*% (x_obs_i - A_i %*% pred_i)
    
    # Updated monthly series
    X_new[, i] <- pred_i + correction
  }
  
  X_new
}


# ==============================================================================
# 2. M-STEP: PCA estimation of factors and loadings
# X_old : T x N completed panel (no NA)
# r     : number of factors
# ==============================================================================

M_step <- function(X_old, r) {
  
  # PCA without re-centering or re-scaling (X_old already standardized)
  pca_res <- prcomp(X_old, center = FALSE, scale. = FALSE)
  
  # First r principal components
  F_hat   <- pca_res$x[, 1:r, drop = FALSE]         # T x r (factor scores)
  Phi_hat <- pca_res$rotation[, 1:r, drop = FALSE]  # N x r (loadings)
  
  list(F = F_hat, Phi = Phi_hat)
}



# ==============================================================================
# Initialization for EM (Hepenstrick & Marcellino)
# XP covariance + ER + XP factors → replace NAs with common component
# ==============================================================================

init_XP_ER <- function(X_std, Kmax = 15) {
  # 1. Covariance estimator
  cov_out <- all_purpose_covariance(X_std)
  
  # 2. Number of factors via ER
  ER_out <- select_num_factors_ER(cov_out$Sigma_tilde, Kmax = Kmax)
  r      <- ER_out$r
  
  # 3. Factors and common component (XP)
  fac_out <- estimate_factors_XP(X_std, r)
  C_hat   <- fac_out$C_hat
  
  # 4. Replace only missing entries with the common component
  X_init <- X_std
  na_idx <- is.na(X_std)
  X_init[na_idx] <- C_hat[na_idx]
  
  list(
    X_init = X_init,          # T x N, no NAs (for EM)
    r      = r,
    Sigma  = cov_out$Sigma_tilde,
    Lambda = fac_out$Lambda,
    F_hat  = fac_out$F_hat,
    C_hat  = C_hat,
    ER     = ER_out$ER
  )
}


# ==============================================================================
# 3. EM ALGORITHM LOOP
# X_init : T x N initial completed panel (no NA), e.g. from XP+ER
# X_obs  : T x N original data with NA (fixed)
# A_list : list of A_i matrices (one per series)
# r      : number of factors
# ==============================================================================


EM_algorithm <- function(X_init, X_obs, A_list, r,
                         max_iter = 1000, tol = 1e-4) {
  
  X_old <- as.matrix(X_init)
  diffs <- numeric(max_iter)   # store max changes
  
  for (iter in 1:max_iter) {
    
    # --- M-step: estimate factors and loadings ---
    mstep   <- M_step(X_old, r)
    F_hat   <- mstep$F          # T x r
    Phi_hat <- mstep$Phi        # N x r
    
    # --- E-step: reconstruct X given (F_hat, Phi_hat) ---
    X_new <- E_step(X_old, X_obs, F_hat, Phi_hat, A_list)
    
    # --- Convergence check (max change in X) ---
    diff_max    <- max(abs(X_new - X_old))
    diffs[iter] <- diff_max
    
    # Calcola anche la violazione dei vincoli A_i X = x_obs (diagnostica)
    viol_A <- 0
    
    for (i in seq_len(ncol(X_new))) {
      A_i <- A_list[[i]]
      if (nrow(A_i) == 0) next
      
      idx_obs_i <- which(!is.na(X_obs[, i]))
      x_obs_i   <- X_obs[idx_obs_i, i]
      
      # predizione aggregata dal pannello X_new
      x_hat_i <- A_i %*% X_new[, i, drop = FALSE]
      
      if (length(x_obs_i) != nrow(A_i)) {
        warning("EM: mismatch tra length(x_obs_i) e nrow(A_i) per serie ", i)
      } else {
        viol_A <- viol_A + sum((x_obs_i - x_hat_i)^2, na.rm = TRUE)
      }
    }
    
    cat("Iter ", iter, 
        ": diff_max = ", signif(diff_max, 3),
        "  viol_A = ", signif(viol_A, 3), "\n")
    
    if (diff_max < tol) {
      message("Converged at iteration ", iter,
              " with max change = ", signif(diff_max, 3))
      break
    }
    
    X_old <- X_new
  }
  
  # keep only the iterations actually used
  diffs <- diffs[1:iter]
  
  list(
    X_completed = X_new,
    F          = F_hat,
    Phi        = Phi_hat,
    iterations = iter,
    diffs      = diffs
  )
}
