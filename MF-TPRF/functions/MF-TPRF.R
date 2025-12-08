
# Input: 
#   X_lf   : matrice T x N dei predittori (dopo EM, eventualmente standardizzata)
#   y_q    : vettore T x 1 della target (stazionaria, non necessariamente centrata)
#   Lproxy : numero totale di proxy automatiche da costruire (inclusa y_q)
build_autoproxy_3prf <- function(X_lf, y_q, Lproxy) {
  
  T_obs <- nrow(X_lf)
  
  # Inizializza matrice delle proxy
  Z <- matrix(nrow = T_obs, ncol = 0)
  
  # L = 1: prima proxy è la target (target-proxy)
  Z <- cbind(Z, y_q)
  colnames(Z) <- "y_q"
  
  # Se voglio solo la target-proxy, esco
  if (Lproxy == 1) {
    return(Z)
  }
  
  # Caso Lproxy > 1: costruisco proxy automatiche
  for (L in 1:(Lproxy - 1)) {
    
    # -------------------------
    # Step 1: OLS di X_lf su [1, Z]
    # -------------------------
    # Z_tilde = [1, Z], dimensione T x (k+1)
    Z_tilde <- cbind(1, Z)
    
    # Coefficienti OLS: B_full ha dimensione N x (k+1)
    A      <- t(Z_tilde) %*% Z_tilde    
    B      <- t(Z_tilde) %*% X_lf
    B_full <- t(solve(A, B))             # N x (k+1)
    
    # Phi: solo le slope rispetto alle proxy (scarto la colonna dell'intercetta)
    Phi <- B_full[, -1, drop = FALSE]    # N x k, k = ncol(Z)
    
    # -------------------------
    # Step 2: fattore mirato
    # F_t = X_lf * Phi * (Phi' Phi)^{-1}
    # -------------------------
    Phi_cross <- t(Phi) %*% Phi                      # k x k
    F_hat     <- X_lf %*% Phi %*% solve(Phi_cross)   # T x k
    
    # -------------------------
    # Step 3: regressione y_t su F_{t-1} (one-step-ahead), con intercetta
    # -------------------------
    y_t   <- y_q[-1]                               # T-1
    F_lag <- F_hat[-nrow(F_hat), , drop = FALSE]   # T-1 x k
    
    fit   <- lm(y_t ~ F_lag)                       # intercetta inclusa di default
    y_hat <- c(NA, as.numeric(predict(fit)))       # riallineo alla lunghezza T
    
    # residui e^(L)
    e <- y_q - y_hat
    
    # Primo periodo non prevedibile: sostituisco l'NA
    e[is.na(e)] <- y_q[1]
    
    # Aggiungo la nuova proxy alla matrice Z
    Z <- cbind(Z, e)
    colnames(Z)[ncol(Z)] <- paste0("e", L)
  }
  
  # Output: matrice T x Lproxy (colonne: y_q, e1, e2, ..., e_{Lproxy-1})
  return(Z)
}






# ==============================================================================
# LAG U-MIDAS MF-3PRF: 
#===============================================================================

choose_UMIDAS_lag <- function(X_lf, X_hf, y_q, Lmax = 3, Lproxy = Lproxy) {
  
  X_lf <- as.matrix(X_lf)   # T_q x N (low frequency, dopo EM)
  X_hf <- as.matrix(X_hf)   # T_m x N (high frequency, mensile)
  y_q  <- as.numeric(y_q)   # T_q
  
  T_q <- length(y_q)
  
  # -------------------------
  # STEP 0: Build L-autoproxy (con intercette già gestite dentro)
  # -------------------------
  Z <- build_autoproxy_3prf(X_lf, y_q, Lproxy)   # T_q x Lproxy
  
  # -------------------------
  # STEP 1: First pass
  # xi,τ = α0_i + Z_τ α_i   (regressione con intercetta)
  # -------------------------
  # Z_tilde = [1, Z]  (T_q x (Lproxy + 1))
  Z_tilde <- cbind(1, Z)
  
  # OLS: B_full = (Z_tilde' Z_tilde)^(-1) Z_tilde' X_lf, dimensione N x (Lproxy+1)
  A <- t(Z_tilde) %*% Z_tilde
  B <- t(Z_tilde) %*% X_lf
  B_full <- t(qr.solve(A, B))       # N x (Lproxy + 1)
  
  # Phi_hat: prendo SOLO le slope rispetto alle proxy (scarto l'intercetta)
  Phi_hat <- B_full[, -1, drop = FALSE]   # N x Lproxy
  
  # -------------------------
  # STEP 2: Second pass: factor extraction
  # Ft = X_hf * Phi_hat * (Phi_hat' Phi_hat)^(-1)
  # -------------------------
  Phi_cross <- t(Phi_hat) %*% Phi_hat              # Lproxy x Lproxy
  
  Ft <- X_hf %*% Phi_hat %*% qr.solve(Phi_cross)   # T_m x Lproxy
  
  # -------------------------
  # STEP 3: split monthly -> quarterly (3 mesi per trimestre)
  # Assumiamo T_m = 3 * T_q
  # -------------------------
  F1 <- Ft[seq(1, 3 * T_q, by = 3), , drop = FALSE]  # mese 1 del trimestre
  F2 <- Ft[seq(2, 3 * T_q, by = 3), , drop = FALSE]  # mese 2
  F3 <- Ft[seq(3, 3 * T_q, by = 3), , drop = FALSE]  # mese 3
  
  # -------------------------
  # LOOP over L (lag MF-UMIDAS)
  # -------------------------
  results <- data.frame(L = integer(),
                        AIC = numeric(),
                        BIC = numeric())
  
  for (L in 1:Lmax) {
    
    Xreg <- NULL
    T_eff <- T_q - L      # numero effettivo di osservazioni utilizzabili
    
    # Per ogni lag 0,...,L-1 costruisco i regressori
    # target: y_dep = y_q[(L+1):T_q], lunghezza T_q - L = T_eff
    # per ogni ell, uso F(τ - ell) con τ = L+1,...,T_q
    for (ell in 1:L) {
      # ell0 = ell - 1 è il lag "vero"
      # indici corretti: (L+1-ell0):(T_q-ell0) con ell0 = ell-1
      idx <- (L + 2 - ell):(T_q + 1 - ell)   # lunghezza T_eff
      
      Xreg <- cbind(
        Xreg,
        F1[idx, , drop = FALSE],
        F2[idx, , drop = FALSE],
        F3[idx, , drop = FALSE]
      )
    }
    
    y_dep <- y_q[(L + 1):T_q]    # T_eff osservazioni
    
    # Regressione U-MIDAS con intercetta (di default in lm)
    fit <- lm(y_dep ~ Xreg)
    
    results <- rbind(
      results,
      data.frame(
        L   = L,
        AIC = AIC(fit),
        BIC = BIC(fit)
      )
    )
  }
  
  return(list(
    results = results,
    lag_AIC = results$L[which.min(results$AIC)],
    lag_BIC = results$L[which.min(results$BIC)]
  ))
}


# ==============================================================================
# MF-3PRF: 
#===============================================================================
#' PROXY
#' LOADINGS 
#' FATTORI
#' TARGET PREDICTION
#' add the high friquency predictors
# ==============================================================================

# X_hf <- X_em   
# X_lf <- X_em_agg 

MF_TPRF <- function(X_lf, X_hf, y_q, Lproxy = 1, L_midas = 1) {
  
  # --------------------------------------------------------
  # Preliminari
  # --------------------------------------------------------
  X_lf <- as.matrix(X_lf)   # trimestrali (T_q x N) - dopo EM + aggregazione
  X_hf <- as.matrix(X_hf)   # mensili (T_m x N)
  y_q  <- as.numeric(y_q)   # target trimestrale (T_q)
  T_q  <- length(y_q)
  T_m  <- nrow(X_hf)
  
  if (L_midas < 1) stop("L_midas deve essere >= 1")
  if (L_midas > T_q) stop("L_midas non può superare T_q")
  
  # --------------------------------------------------------
  # STEP 0 – Costruzione L-autoproxy (usa intercetta internamente)
  # --------------------------------------------------------
  Z <- build_autoproxy_3prf(X_lf, y_q, Lproxy)   # T_q x Lproxy
  
  # --------------------------------------------------------
  # STEP 1 – First Pass con intercetta
  #
  #  x_{i,τ} = α0_i + Z_τ α_i + u_{i,τ}
  #
  # Matricialmente:
  #  Beta_full = (Z_tilde' Z_tilde)^(-1) Z_tilde' X_lf
  #  con Z_tilde = [1 , Z]
  #  Beta_full ha dimensione (1+Lproxy) x N
  #  Prendiamo solo le slope α_i (righe 2..)
  # --------------------------------------------------------
  Z_tilde <- cbind(1, Z)                     # T_q x (1+Lproxy)
  XtX_1   <- t(Z_tilde) %*% Z_tilde          # (1+Lproxy) x (1+Lproxy)
  XtY_1   <- t(Z_tilde) %*% X_lf             # (1+Lproxy) x N
  
  Beta_full <- solve(XtX_1, XtY_1)           # (1+Lproxy) x N
  Beta_full <- t(Beta_full)                  # N x (1+Lproxy)
  
  Phi_hat <- Beta_full[, -1, drop = FALSE]   # N x Lproxy  (solo slope α_i)
  
  # --------------------------------------------------------
  # STEP 2 – Second Pass: fattori mensili
  #
  #  F_t = X_hf Φ_hat (Φ_hat' Φ_hat)^(-1)
  # --------------------------------------------------------
  XtX_2   <- t(Phi_hat) %*% Phi_hat          # Lproxy x Lproxy
  XtY_2   <- t(Phi_hat) %*% t(X_hf)          # Lproxy x T_m
  
  # Solviamo: (Phi'Phi) * M = Phi' X_hf'
  #  => M = (Phi'Phi)^(-1) Phi' X_hf'
  M_2     <- solve(XtX_2, XtY_2)             # Lproxy x T_m
  F_hat   <- t(M_2)                          # T_m x Lproxy
  
  # --------------------------------------------------------
  # STEP 3.1 – Fattori trimestrali per trimestri COMPLETI
  # --------------------------------------------------------
  if (T_m < 3 * T_q) {
    stop("Ci sono meno mesi di quelli necessari per coprire tutti i trimestri di y_q.")
  }
  
  # trimestri 1..T_q (completi) → usiamo i primi 3*T_q mesi
  F1 <- F_hat[seq(1, 3 * T_q, by = 3), , drop = FALSE]  # mese 1
  F2 <- F_hat[seq(2, 3 * T_q, by = 3), , drop = FALSE]  # mese 2
  F3 <- F_hat[seq(3, 3 * T_q, by = 3), , drop = FALSE]  # mese 3
  
  K <- ncol(F1)  # numero di fattori
  
  # mesi extra → trimestre T_q+1 (ragged edge)
  rem <- T_m - 3 * T_q               # 0,1,2 (o 3) mesi del trimestre successivo
  F_next1 <- if (rem >= 1) F_hat[3 * T_q + 1, , drop = FALSE] else NULL
  F_next2 <- if (rem >= 2) F_hat[3 * T_q + 2, , drop = FALSE] else NULL
  F_next3 <- if (rem >= 3) F_hat[3 * T_q + 3, , drop = FALSE] else NULL
  
  # --------------------------------------------------------
  # STEP 3.2 – U-MIDAS trimestrale con L_midas trimestri (corrente + lag)
  #
  # Modello:
  #  y_τ = β0 + sum_{ℓ=0}^{L_midas-1} [ β_{ℓ,1}' F1_{τ-ℓ} 
  #                                     + β_{ℓ,2}' F2_{τ-ℓ}
  #                                     + β_{ℓ,3}' F3_{τ-ℓ} ] + η_τ
  #
  # Stima su trimestri COMPLETI: τ = L_midas .. T_q
  # --------------------------------------------------------
  rows_list <- list()
  row_id    <- 1
  
  for (tau in L_midas:T_q) {
    row_vec <- c()
    for (ell_id in 1:L_midas) {
      ell   <- ell_id - 1       # ell = 0,1,...,L_midas-1
      lag_q <- tau - ell
      # ogni blocco: [F1_{lag_q}, F2_{lag_q}, F3_{lag_q}]
      row_vec <- c(
        row_vec,
        F1[lag_q, ],
        F2[lag_q, ],
        F3[lag_q, ]
      )
    }
    rows_list[[row_id]] <- row_vec
    row_id <- row_id + 1
  }
  
  Xreg <- do.call(rbind, rows_list)          # (T_q - L_midas + 1) x (3*K*L_midas)
  y_dep <- y_q[L_midas:T_q]                  # vettore dipendente
 
  # OLS con intercetta: y = β0 + Xreg β + errore
  X_tilde_3 <- cbind(1, Xreg)                # aggiungo colonna di 1
  XtX_3     <- t(X_tilde_3) %*% X_tilde_3    # (1+3*K*L_midas) x (1+3*K*L_midas)
  XtY_3     <- t(X_tilde_3) %*% y_dep        # (1+3*K*L_midas) x 1
  
  beta_hat  <- solve(XtX_3, XtY_3)           # (1+3*K*L_midas) x 1
  
  beta0     <- as.numeric(beta_hat[1])
  beta_vec  <- as.numeric(beta_hat[-1])
  
  # beta_mat: righe = ℓ_id = 1..L_midas (ell = 0..L_midas-1)
  # colonne per ogni lag: [β_{ℓ,1} (K), β_{ℓ,2} (K), β_{ℓ,3} (K)]
  beta_mat <- matrix(beta_vec,
                     nrow = L_midas,
                     ncol = 3 * K,
                     byrow = TRUE)
  
  # --------------------------------------------------------
  # STEP 4 – Nowcasting mensile
  #
  # y_nowcast: lunghezza T_m (un nowcast per ogni mese disponibile).
  #
  # Per ciascun trimestre τ = L_midas..T_q (storici):
  #   al mese m=1: usa solo F1_τ come "corrente" (ell=0), + lag
  #   al mese m=2: F1_τ + F2_τ
  #   al mese m=3: F1_τ + F2_τ + F3_τ
  #
  # Per il trimestre T_q+1 (se rem>=1):
  #   usa F_next1/F_next2/F_next3 allo stesso modo.
  # --------------------------------------------------------
  
  y_nowcast <- rep(NA_real_, T_m)
  
  # 4.a) Trimestri COMPLETI 1..T_q (backtest pseudo real time)
  for (tau in L_midas:T_q) {
    
    month_idx <- ((tau - 1) * 3 + 1):(tau * 3)  # posizioni dei 3 mesi di τ
    
    for (m in 1:3) {   # m = 1,2,3 (mese nel trimestre τ)
      
      contrib <- 0
      
      # somma sui lag ℓ_id=1..L_midas
      for (ell_id in 1:L_midas) {
        ell   <- ell_id - 1
        lag_q <- tau - ell
        if (lag_q < 1) next
        
        # per ogni mese mm=1,2,3
        for (mm in 1:3) {
          
          # Ragged edge "interno" al trimestre:
          # - per ell=0 (trimestre τ stesso), includo solo i mesi mm <= m
          # - per i lag (ell>=1), includo sempre mm=1,2,3
          if (ell_id == 1 && mm > m) next
          
          F_qm <- switch(
            mm,
            `1` = F1[lag_q, ],
            `2` = F2[lag_q, ],
            `3` = F3[lag_q, ]
          )
          
          start_col <- (mm - 1) * K + 1
          end_col   <- mm * K
          beta_block <- beta_mat[ell_id, start_col:end_col]
          
          contrib <- contrib + sum(F_qm * beta_block)
        }
      }
      
      y_nowcast[month_idx[m]] <- beta0 + contrib
    }
  }
  
  # 4.b) Trimestre CORRENTE T_q+1 (se ho mesi extra)
  if (rem > 0) {
    
    tau_curr   <- T_q + 1
    month_curr <- 3 * T_q + seq_len(rem)   # indici mesi disponibili del trimestre T_q+1
    
    for (m in 1:rem) {   # m=1,2 (o 3) mesi osservati nel trimestre corrente
      
      contrib <- 0
      
      for (ell_id in 1:L_midas) {
        ell   <- ell_id - 1
        lag_q <- tau_curr - ell
        
        # Per i lag (ell>=1) uso i trimestri COMPLETI in 1..T_q
        if (ell_id >= 2 && (lag_q < 1 || lag_q > T_q)) next
        
        for (mm in 1:3) {
          
          if (ell_id == 1) {
            # Trimestre corrente T_q+1: uso solo mesi <= m
            if (mm > m) next
            
            F_qm <- switch(
              mm,
              `1` = F_next1,
              `2` = F_next2,
              `3` = F_next3
            )
            if (is.null(F_qm)) next   # mese non ancora osservato
            
          } else {
            # Trimestri laggati: uso sempre i 3 mesi F1,F2,F3
            F_qm <- switch(
              mm,
              `1` = F1[lag_q, ],
              `2` = F2[lag_q, ],
              `3` = F3[lag_q, ]
            )
          }
          
          start_col <- (mm - 1) * K + 1
          end_col   <- mm * K
          beta_block <- beta_mat[ell_id, start_col:end_col]
          
          contrib <- contrib + sum(F_qm * beta_block)
        }
      }
      
      y_nowcast[month_curr[m]] <- beta0 + contrib
    }
  }
  
  return(list(
    Z         = Z,
    Phi_hat   = Phi_hat,
    F_hat     = F_hat,
    F1        = F1,
    F2        = F2,
    F3        = F3,
    beta0     = beta0,
    beta_mat  = beta_mat,
    y_nowcast = y_nowcast
  ))
}
