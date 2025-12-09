# ==============================================================================
# UTILITY TO ADD AN F-TEST IN THE FIRST STEP OF THEMF-3PRF: 
#===============================================================================

step1_select_Phi <- function(X_lf, Z,
                             alpha_level = 0.10,
                             robust_type = c("White", "NW"),
                             nw_lag = 1) {
  robust_type <- match.arg(robust_type)
  
  X_lf   <- as.matrix(X_lf)
  Z      <- as.data.frame(Z)   # preserva i nomi di colonna
  T_q    <- nrow(X_lf)
  N      <- ncol(X_lf)
  Lproxy <- ncol(Z)
  
  # Qui ci aspettiamo che Z abbia già i nomi: "y_q", "e1", "e2", ...
  if (is.null(colnames(Z))) {
    stop("Z deve avere nomi di colonna (ad es. 'y_q', 'e1', ...).")
  }
  
  Phi_hat_sel <- matrix(NA_real_, nrow = N, ncol = Lproxy)
  colnames(Phi_hat_sel) <- colnames(Z)  # stessi nomi delle proxy
  
  for (i in 1:N) {
    df_i <- data.frame(x = X_lf[, i], Z)  # x ~ y_q + e1 + e2 + ...
    mod_i <- lm(x ~ ., data = df_i)       # intercetta + tutte le proxy
    
    # matrice di covarianza robusta (White o Newey-West)
    if (robust_type == "White") {
      Vcov_i <- vcovHC(mod_i, type = "HC1")
    } else { # "NW"
      Vcov_i <- NeweyWest(mod_i, lag = nw_lag,
                          prewhite = FALSE, adjust = TRUE)
    }
    
    # slope (escludo l'intercetta)
    alpha_hat_i <- coef(mod_i)[-1]
    
    if (Lproxy == 1) {
      # ---- caso 1 proxy: t-test robusto ----
      test_i <- coeftest(mod_i, vcov. = Vcov_i)
      # la seconda riga è il coefficiente del primo proxy
      p_val  <- test_i[2, "Pr(>|t|)"]
      
      if (is.na(p_val) || p_val > alpha_level) {
        alpha_hat_i[] <- 0
      }
      
    } else {
      # ---- caso Lproxy > 1: F/Wald congiunto su TUTTE le proxy ----
      # es: c("y_q = 0", "e1 = 0", "e2 = 0", ...)
      hyp <- paste0(colnames(Z), " = 0")
      
      Ftest_i <- linearHypothesis(mod_i, hyp,
                                  vcov. = Vcov_i,
                                  test  = "F")
      p_val_F <- Ftest_i[2, "Pr(>F)"]
      
      if (is.na(p_val_F) || p_val_F > alpha_level) {
        # non significativo congiuntamente -> azzero tutta la serie
        alpha_hat_i[] <- 0
      }
      # (se un domani vuoi fare sparsity fine, qui guardi anche i singoli t-test)
    }
    
    Phi_hat_sel[i, ] <- alpha_hat_i
  }
  
  return(Phi_hat_sel)
}



# ==============================================================================
# NUMBER OF PROXIES IN THE MF-3PRF: 
#===============================================================================

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
# Fixing p_AR as a parameter
# considering the restriction after the Wald Test

choose_UMIDAS_lag <- function(X_lf, X_hf, y_q,
                              Lmax         = 5,
                              Lproxy       = 1, 
                              p_AR         = 1,      # ordine AR(y) scelto FUORI
                              Robust_F     = FALSE,
                              alpha        = 0.10,
                              robust_type  = c("White", "NW"),
                              nw_lag       = 1){
  
  robust_type <- match.arg(robust_type)
  
  # --------------------------------------------------------
  # Preliminari
  # --------------------------------------------------------
  X_lf <- as.matrix(X_lf)   # trimestrali (T_q x N) - dopo EM + aggregazione
  X_hf <- as.matrix(X_hf)   # mensili (T_m x N)
  y_q  <- as.numeric(y_q)   # target trimestrale (T_q)
  T_q  <- length(y_q)
  T_m  <- nrow(X_hf)
  
  if (Lmax < 1) stop("Lmax deve essere >= 1")
  if (Lmax >= T_q) stop("Lmax non può essere >= T_q")
  if (p_AR < 0) stop("p_AR deve essere >= 0")
  
  # --------------------------------------------------------
  # STEP 0 – Costruzione L-autoproxy (usa intercetta internamente)
  # --------------------------------------------------------
  Z <- build_autoproxy_3prf(X_lf, y_q, Lproxy)   # T_q x Lproxy
  
  # --------------------------------------------------------
  # STEP 1 – First Pass con intercetta
  #
  #  x_{i,τ} = α0_i + Z_τ α_i + u_{i,τ}
  # --------------------------------------------------------
  Z_tilde <- cbind(1, Z)                          # T_q x (1+Lproxy)
  XtX_1   <- t(Z_tilde) %*% Z_tilde               # (1+Lproxy) x (1+Lproxy)
  XtY_1   <- t(Z_tilde) %*% X_lf                  # (1+Lproxy) x N
  
  Beta_full <- solve(XtX_1, XtY_1)                # (1+Lproxy) x N
  Beta_full <- t(Beta_full)                       # N x (1+Lproxy)
  
  Phi_hat_full <- Beta_full[, -1, drop = FALSE]   # N x Lproxy (solo slope)
  
  # --------------------------------------------------------
  # STEP 1bis – First Pass con Robust F-test (opzionale)
  # --------------------------------------------------------
  if (Robust_F == FALSE) {
    Phi_hat <- Phi_hat_full
  } else {
    Phi_hat <- step1_select_Phi(X_lf        = X_lf,
                                Z           = Z,
                                alpha_level = alpha,
                                robust_type = robust_type,
                                nw_lag      = nw_lag)
  }
  
  # --------------------------------------------------------
  # STEP 2 – Second Pass: fattori mensili
  #
  #  F_t = X_hf Φ_hat (Φ_hat' Φ_hat)^(-1)
  # --------------------------------------------------------
  XtX_2   <- t(Phi_hat) %*% Phi_hat          # Lproxy x Lproxy
  XtY_2   <- t(Phi_hat) %*% t(X_hf)          # Lproxy x T_m
  
  # Uso qr.solve per sicurezza (Phi'Phi potrebbe essere quasi singolare)
  M_2     <- qr.solve(XtX_2, XtY_2)          # Lproxy x T_m
  F_hat   <- t(M_2)                          # T_m x Lproxy
  
  # --------------------------------------------------------
  # STEP 3.1 – Fattori trimestrali per trimestri COMPLETI
  # --------------------------------------------------------
  if (T_m < 3 * T_q) {
    stop("Ci sono meno mesi di quelli necessari per coprire tutti i trimestri di y_q.")
  }
  
  F1 <- F_hat[seq(1, 3 * T_q, by = 3), , drop = FALSE]  # mese 1
  F2 <- F_hat[seq(2, 3 * T_q, by = 3), , drop = FALSE]  # mese 2
  F3 <- F_hat[seq(3, 3 * T_q, by = 3), , drop = FALSE]  # mese 3
  
  # --------------------------------------------------------
  # LOOP over L (lag MF-UMIDAS) con AR(p_AR) su y
  # --------------------------------------------------------
  results <- data.frame(L = integer(),
                        AIC = numeric(),
                        BIC = numeric())
  
  for (L in 1:Lmax) {
    
    # punto di partenza in termini di trimestre:
    # serve avere disponibili sia i L lag dei fattori sia i p_AR lag di y
    start_tau <- max(L, p_AR) + 1
    T_eff     <- T_q - max(L, p_AR)   # numero osservazioni effettive
    
    # variabile dipendente: y_tau per tau = start_tau,...,T_q
    y_dep <- y_q[start_tau:T_q]
    
    # ---- blocco AR in y: lag trimestrali di y ----
    Y_lag <- NULL
    if (p_AR > 0) {
      for (j in 1:p_AR) {
        Y_lag <- cbind(
          Y_lag,
          y_q[(start_tau - j):(T_q - j)]
        )
      }
      colnames(Y_lag) <- paste0("y_lag", 1:p_AR)
    }
    
    # ---- blocco fattori U-MIDAS: lag di F1,F2,F3 ----
    Xreg_F <- NULL
    if (L > 0) {
      for (ell in 0:(L-1)) {
        idx <- (start_tau - ell):(T_q - ell)   # lunghezza T_eff
        
        Xreg_F <- cbind(
          Xreg_F,
          F1[idx, , drop = FALSE],
          F2[idx, , drop = FALSE],
          F3[idx, , drop = FALSE]
        )
      }
    }
    
    # Combino AR + fattori
    if (p_AR > 0 && !is.null(Xreg_F)) {
      Xreg <- cbind(Y_lag, Xreg_F)
    } else if (p_AR > 0 && is.null(Xreg_F)) {   # solo AR
      Xreg <- Y_lag
    } else {
      Xreg <- Xreg_F                            # solo fattori
    }
    
    # Regressione U-MIDAS con intercetta
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
# Robust_F   = TRUE
# alpha      = 0.05
# robust_type = "NW"
# nw_lag      = 1

MF_TPRF <- function(X_lf, X_hf, y_q,
                    Lproxy = 1,
                    L_midas = 1,
                    p_AR   = 1,          # nuovo: ordine AR scelto FUORI
                    Robust_F   = FALSE,
                    alpha      = 0.10,
                    robust_type = c("White", "NW"),
                    nw_lag      = 1) {
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
  Z_tilde <- cbind(1, Z)                          # T_q x (1+Lproxy)
  XtX_1   <- t(Z_tilde) %*% Z_tilde               # (1+Lproxy) x (1+Lproxy)
  XtY_1   <- t(Z_tilde) %*% X_lf                  # (1+Lproxy) x N
  
  Beta_full <- solve(XtX_1, XtY_1)                # (1+Lproxy) x N
  Beta_full <- t(Beta_full)                       # N x (1+Lproxy)
  
  Phi_hat_full <- Beta_full[, -1, drop = FALSE]   # N x Lproxy (solo slope)
  
  # --------------------------------------------------------
  # STEP 1bis – First Pass con intercetta considering a Robust F-test
  # --------------------------------------------------------
  
  if (Robust_F == FALSE) {
    Phi_hat <- Phi_hat_full
  } else {
    Phi_hat <- step1_select_Phi(X_lf        = X_lf,
                                Z           = Z,
                                alpha_level = alpha,
                                robust_type = robust_type,
                                nw_lag      = nw_lag)
  }
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
  # STEP 3.2 – U-MIDAS trimestrale con AR(p_AR) in y
  #
  # Modello:
  #  y_τ = β0
  #        + sum_{j=1}^{p_AR} ρ_j y_{τ-j}
  #        + sum_{ℓ=0}^{L_midas-1} [ β_{ℓ,1}' F1_{τ-ℓ} 
  #                                  + β_{ℓ,2}' F2_{τ-ℓ}
  #                                  + β_{ℓ,3}' F3_{τ-ℓ} ] + η_τ
  #
  # Stima su trimestri COMPLETI:
  #  τ = start_tau, ..., T_q
  #  con start_tau = max(L_midas, p_AR) + 1
  # --------------------------------------------------------
  
  if (p_AR < 0) stop("p_AR deve essere >= 0")
  
  start_tau <- max(L_midas, p_AR) + 1
  T_eff     <- T_q - max(L_midas, p_AR)
  
  # variabile dipendente: y_τ per τ = start_tau,...,T_q
  y_dep <- y_q[start_tau:T_q]   # lunghezza T_eff
  
  # ---- blocco AR in y ----
  Y_lag <- NULL
  if (p_AR > 0) {
    for (j in 1:p_AR) {
      Y_lag <- cbind(
        Y_lag,
        y_q[(start_tau - j):(T_q - j)]
      )
    }
    colnames(Y_lag) <- paste0("y_lag", 1:p_AR)
  }
  
  # ---- blocco fattori U-MIDAS (F1, F2, F3 con lag 0..L_midas-1) ----
  Xreg_F <- NULL
  for (ell_id in 1:L_midas) {
    ell   <- ell_id - 1                 # ell = 0,...,L_midas-1
    idx   <- (start_tau - ell):(T_q - ell)   # lunghezza T_eff
    
    Xreg_F <- cbind(
      Xreg_F,
      F1[idx, , drop = FALSE],
      F2[idx, , drop = FALSE],
      F3[idx, , drop = FALSE]
    )
  }
  
  # Combino AR + fattori
  if (p_AR > 0) {
    Xreg <- cbind(Y_lag, Xreg_F)
  } else {
    Xreg <- Xreg_F
  }
  
  # OLS con intercetta: y = β0 + Y_lag ρ + Xreg_F β + errore
  X_tilde_3 <- cbind(1, Xreg)
  XtX_3     <- t(X_tilde_3) %*% X_tilde_3
  XtY_3     <- t(X_tilde_3) %*% y_dep
  
  beta_hat  <- solve(XtX_3, XtY_3)
  
  beta0     <- as.numeric(beta_hat[1])
  
  # primi p_AR coefficienti dopo l'intercetta = AR in y
  if (p_AR > 0) {
    rho_hat  <- as.numeric(beta_hat[2:(1 + p_AR)])   # ρ_1,...,ρ_{p_AR}
    beta_vec <- as.numeric(beta_hat[(2 + p_AR):length(beta_hat)])
  } else {
    rho_hat  <- numeric(0)
    beta_vec <- as.numeric(beta_hat[-1])
  }
  
  # beta_mat: righe = ℓ_id = 1..L_midas (ell = 0..L_midas-1)
  # colonne per ogni lag: [β_{ℓ,1} (K), β_{ℓ,2} (K), β_{ℓ,3} (K)]
  beta_mat <- matrix(beta_vec,
                     nrow = L_midas,
                     ncol = 3 * K,
                     byrow = TRUE)
  
  # --------------------------------------------------------
  # STEP 4 – Nowcasting mensile con AR(p_AR) in y
  #
  # y_nowcast: lunghezza T_m (un nowcast per ogni mese disponibile).
  #
  # Per ciascun trimestre τ (storici) in cui il modello è definito:
  #   τ = start_tau, ..., T_q
  #   con start_tau = max(L_midas, p_AR) + 1
  #
  # Ad ogni trimestre τ:
  #   - parte AR:    AR_part(τ) = sum_{j=1}^{p_AR} ρ_j * y_{τ-j}
  #   - parte fattori: come nel codice originale (con ragged edge interno)
  #   - nowcast mensile m = 1,2,3: β0 + AR_part(τ) + contrib_fattori(τ,m)
  #
  # Per il trimestre T_q+1 (se rem>=1):
  #   - parte AR:    usa y_{T_q}, y_{T_q-1}, ..., y_{T_q+1-p_AR}
  #   - parte fattori: F_next1/F_next2/F_next3 + lag storici
  # --------------------------------------------------------
  
  y_nowcast <- rep(NA_real_, T_m)
  
  # punto di partenza coerente con la stima di STEP 3.2
  start_tau <- max(L_midas, p_AR) + 1
  
  # 4.a) Trimestri COMPLETI (backtest pseudo real time)
  if (start_tau <= T_q) {
    for (tau in start_tau:T_q) {
      
      # mesi del trimestre τ
      month_idx <- ((tau - 1) * 3 + 1):(tau * 3)
      
      # -------- parte AR in y per il trimestre τ --------
      if (p_AR > 0) {
        # y_{τ-1}, ..., y_{τ-p_AR}
        y_lags_tau <- sapply(1:p_AR, function(j) y_q[tau - j])
        AR_part_tau <- sum(rho_hat * y_lags_tau)
      } else {
        AR_part_tau <- 0
      }
      
      # -------- parte fattoriale mese per mese --------
      for (m in 1:3) {   # m = 1,2,3 (mese nel trimestre τ)
        
        contrib <- 0
        
        # somma sui lag ℓ_id = 1..L_midas
        for (ell_id in 1:L_midas) {
          ell   <- ell_id - 1
          lag_q <- tau - ell
          if (lag_q < 1) next
          
          # per ogni mese mm = 1,2,3
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
        
        # nowcast per il mese m del trimestre τ:
        # β0 + parte AR (solo y) + parte fattori
        y_nowcast[month_idx[m]] <- beta0 + AR_part_tau + contrib
      }
    }
  }
  
  # 4.b) Trimestre CORRENTE T_q+1 (se ho mesi extra)
  if (rem > 0) {
    
    tau_curr   <- T_q + 1
    month_curr <- 3 * T_q + seq_len(rem)   # indici mesi disponibili del trimestre T_q+1
    
    # -------- parte AR in y per il trimestre corrente --------
    if (p_AR > 0) {
      # y_{T_q+1-1}, ..., y_{T_q+1-p_AR} = y_{T_q}, y_{T_q-1}, ...
      y_lags_curr <- sapply(1:p_AR, function(j) y_q[tau_curr - j])
      AR_part_curr <- sum(rho_hat * y_lags_curr)
    } else {
      AR_part_curr <- 0
    }
    
    # -------- parte fattoriale con F_next1/F_next2/F_next3 --------
    for (m in 1:rem) {   # m = 1,2 (o 3) mesi osservati nel trimestre corrente
      
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
      
      # nowcast per il mese m del trimestre T_q+1
      y_nowcast[month_curr[m]] <- beta0 + AR_part_curr + contrib
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
