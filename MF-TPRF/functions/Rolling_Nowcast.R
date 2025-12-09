# ==============================================================================
# UNBALANCEDNESS
# ==============================================================================

# X_full <- X
# current_t <- 303

unbalancedness <- function(X_full, dates, Freq, Unb, current_t) {
  
  stopifnot(nrow(X_full) >= current_t)
  X_cut <- X_full[1:current_t, , drop = FALSE]
  
  for (j in seq_len(ncol(X_cut))) {
    
    delay_j <- Unb[j]          # ritardo di pubblicazione (in mesi)
    release_t <- current_t - delay_j
    
    if (release_t < 1) {
      # la variabile NON è mai stata pubblicata entro current_t
      X_cut[, j] <- NA
    } else {
      # la variabile è disponibile solo fino a release_t
      if (release_t + 1 <= current_t) {
        X_cut[(release_t + 1):current_t, j] <- NA
      }
    }
  }
  
  return(X_cut)
}


compute_m_tr <- function(date_t, dates_q) {
  
  # ultimo trimestre rilasciato PRIMA del mese corrente
  last_q_date <- max(dates_q[dates_q < date_t])
  
  # indice del trimestre
  last_q_idx <- which(dates_q == last_q_date)
  
  # mesi del trimestre successivo (target quarter)
  M1 <- dates_q[last_q_idx] %m+% months(1)
  M2 <- dates_q[last_q_idx] %m+% months(2)
  M3 <- dates_q[last_q_idx] %m+% months(3)
  
  if (date_t == M1) return(1)
  if (date_t == M2) return(2)
  if (date_t == M3) return(3)
  
  return(NA)
}


# ==============================================================================
# EXPANDING NOWCAST
# ==============================================================================

# X_full <- X
# dates <- dates_m
# tt <- 304 # [M1]
# tt <- 305 # [M2]
# tt <- 306 # [M3]

pseudo_realtime_TPRF_EM <- function(
    X_full,      # matrice mensile completa (T_m x N)
    y_q,         # vettore trimestrale GDP (T_q)
    params,
    dates,       # date mensili (lunghezza T_m)
    dates_q,     # date trimestrali (lunghezza T_q)
    Freq,        # vettore "M"/"Q" per ciascuna colonna di X_full
    Unb,         # ritardo di pubblicazione (in mesi) per ciascuna colonna di X_full
    Type,        # (se ti serve per altro, qui non lo uso)
    agg_m,       # oggetto per aggregazione mensile->trimestrale (per X_m)
    agg_q        # oggetto per aggregazione trimestrale (per X_q)
) {
  # --------------------------------------------
  # 0. Evaluation window
  # --------------------------------------------
  t_start <- which(dates == params$start_eval)
  t_end   <- which(dates == params$end_eval)
  
  if (length(t_start) == 0 | length(t_end) == 0) {
    stop("start_eval o end_eval non trovati in 'dates'.")
  }
  
  # --------------------------------------------
  # 0.b Estimation window per Lproxy, L_midas
  #     (tutte le date < start_eval)
  # --------------------------------------------
  t_est_end <- max(which(dates < params$start_eval))
  if (length(t_est_end) == 0 || t_est_end < 24) {
    stop("Estimation sample troppo corto o non definito.")
  }
  
  cat("\n>>> HYPER-PARAM SELECTION using data up to",
      as.character(dates[t_est_end]), "\n")
  
  ## --- Pre-step: costruisco dataset mensile/trimestrale "estimation" ---
  
  # Unbalancedness alla data t_est_end
  X_cut_est <- unbalancedness(
    X_full    = X_full,
    dates     = dates,
    Freq      = Freq,
    Unb       = Unb,
    current_t = t_est_end
  )
  
  # A_list per EM
  N_q  <- sum(Freq == "Q")
  A_est <- A_list(X_cut_est, N_q, agg_q)
  
  # Standardization
  out_std_est <- standardize_with_na(X_cut_est)
  X_std_est   <- out_std_est$X_std
  
  # Init EM (qui init_XP_ER seleziona r)
  init_est   <- init_XP_ER(X_std_est)
  X_init_est <- init_est$X_init
  r_est      <- init_est$r
  
  # EM
  EM_out_est <- EM_algorithm(X_init_est, X_std_est, A_est,
                             r = r_est, max_iter = 50, tol = 1e-4)
  X_em_est   <- EM_out_est$X_completed
  T_m_est    <- nrow(X_em_est)
  
  # Separa M / Q e aggrega a trimestrale
  N_m  <- sum(Freq == "M")
  X_m_em_est <- X_em_est[, Freq == "M", drop = FALSE]
  X_q_em_est <- X_em_est[, Freq == "Q", drop = FALSE]
  
  X_mq_em_est <- agg_mq(X_m_em_est, agg_m)
  X_qq_em_est <- agg_qq(X_q_em_est, agg_q)
  
  X_em_agg_est <- cbind(X_mq_em_est, X_qq_em_est)
  T_q_em_est   <- nrow(X_em_agg_est)
  
  # Quanti trimestri di PIL sono pubblicati a t_est_end?
  idx_pub_est <- which(dates_q < dates[t_est_end])
  if (length(idx_pub_est) < 2) {
    stop("Troppo pochi trimestri di PIL nell'estimation sample.")
  }
  T_q_current_est <- tail(idx_pub_est, 1)
  T_q_current_est <- min(T_q_current_est, T_q_em_est)
  
  y_q_est   <- y_q[1:T_q_current_est]
  X_lf_est  <- X_em_agg_est[1:T_q_current_est, , drop = FALSE]
  X_hf_est  <- X_em_est[1:T_m_est, , drop = FALSE]
  
  ## --- Scelta Lproxy su base mensile (estimation) ---
  
  # Costruisco y mensile espansa tipo (NA, NA, y_tau)
  y_m_expand_est <- as.numeric(kronecker(y_q_est, c(NA, NA, 1)))
  y_m_cut_est    <- y_m_expand_est[1:T_m_est]
  
  data_proxy_est <- cbind(X_hf_est, y_m_cut_est)
  Y_out_std_est  <- standardize_with_na(data_proxy_est)
  Y_std_est      <- Y_out_std_est$X_std
  
  cov_proxy_est  <- all_purpose_covariance(Y_std_est)
  ER_proxy_est   <- select_num_factors_ER(cov_proxy_est$Sigma_tilde)
  Lproxy_fix     <- ER_proxy_est$r
  
  cat("# [Hyper] Lproxy selected on estimation sample:", Lproxy_fix, "\n")
  
  ## --- Scelta L_midas su base trimestrale (estimation) ---
  
  lag_sel_est <- choose_UMIDAS_lag(
    X_lf       = X_lf_est,
    X_hf       = X_hf_est,
    y_q        = y_q_est,
    Lproxy     = Lproxy_fix,
    Lmax       = params$Lmax,
    Robust_F   = params$Robust_F,
    alpha      = params$alpha,
    robust_type = params$robust_type,
    nw_lag      = params$nw_lag,
    p_AR        = params$p_AR   # se l'hai aggiunto in choose_UMIDAS_lag
  )
  
  L_midas_fix <- lag_sel_est$lag_BIC
  cat("# [Hyper] L_midas selected on estimation sample (BIC):",
      L_midas_fix, "\n")
  
  # ------------------------------------------------
  # Ora passo alla evaluation window con hyper fissi
  # ------------------------------------------------
  
  now_M1 <- list()
  now_M2 <- list()
  now_M3 <- list()
  
  N      <- ncol(X_full)
  
  # --------------------------------------------
  # 1. EXPANDING NOWCAST LOOP sui mesi tt
  # --------------------------------------------
  for (tt in seq(t_start, t_end)) {
    
    date_t <- dates[tt]
    cat("\n>>> REAL-TIME at", as.character(date_t),
        "| Lproxy_fix =", Lproxy_fix,
        "| L_midas_fix =", L_midas_fix, "\n")
    
    # --------------------------
    # Step 1: Unbalanced cut
    # --------------------------
    X_cut <- unbalancedness(
      X_full    = X_full,
      dates     = dates,
      Freq      = Freq,
      Unb       = Unb,
      current_t = tt
    )
    
    # --------------------------
    # Step 2: Build A_list (per EM)
    # --------------------------
    A <- A_list(X_cut, N_q, agg_q)
    
    # --------------------------
    # Step 3: Standardization
    # --------------------------
    out_std <- standardize_with_na(X_cut)
    X_std   <- out_std$X_std
    
    # --------------------------
    # Step 4: Init EM
    # --------------------------
    init   <- init_XP_ER(X_std, Kmax = r_est)
    X_init <- init$X_init
    r      <- init$r      # numero di fattori (può variare nel tempo se vuoi)
    
    # --------------------------
    # Step 5: EM
    # --------------------------
    EM_out <- EM_algorithm(X_init, X_std, A, r, max_iter = 50, tol = 1e-4)
    X_em   <- EM_out$X_completed
    T_m_cur <- nrow(X_em)
    
    # --------------------------
    # Step 6: separa mensili / trimestrali e aggrega
    # --------------------------
    X_m_em <- X_em[, Freq == "M", drop = FALSE]
    X_q_em <- X_em[, Freq == "Q", drop = FALSE]
    
    X_mq_em <- agg_mq(X_m_em, agg_m)
    X_qq_em <- agg_qq(X_q_em, agg_q)
    
    X_em_agg <- cbind(X_mq_em, X_qq_em)
    T_q_em   <- nrow(X_em_agg)
    
    # --------------------------
    # Step 7: trimestri di PIL pubblicati a date_t
    # --------------------------
    idx_pub <- which(dates_q < date_t)
    if (length(idx_pub) < 2) {
      next
    }
    T_q_current <- tail(idx_pub, 1)
    T_q_current <- min(T_q_current, T_q_em)
    
    # --------------------------
    # Step 8: taglia y_q e X_lf/X_hf ai trimestri utilizzabili
    # --------------------------
    y_q_cut   <- y_q[1:T_q_current]
    X_lf_cut  <- X_em_agg[1:T_q_current, , drop = FALSE]
    X_hf_cut  <- X_em[1:T_m_cur, , drop = FALSE]
    
    if (length(y_q_cut) < 2) next
    
    # --------------------------
    # Step 9 & 10: ORA NON seleziono più Lproxy e L_midas
    #              Li uso fissi: Lproxy_fix, L_midas_fix
    # --------------------------
    
    Lproxy  <- Lproxy_fix
    L_midas <- L_midas_fix
    
    # --------------------------
    # Step 11: MF-TPRF su sotto-campione (con hyper fissi)
    # --------------------------
    MF_TPRF_RT_out <- MF_TPRF(
      X_lf      = X_lf_cut,
      X_hf      = X_hf_cut,
      y_q       = y_q_cut,
      Lproxy    = Lproxy,
      L_midas   = L_midas,
      p_AR      = params$p_AR,
      Robust_F  = params$Robust_F,
      alpha     = params$alpha,
      robust_type = params$robust_type,
      nw_lag      = params$nw_lag
    )
    
    y_rt_full <- MF_TPRF_RT_out$y_nowcast
    y_rt_last <- tail(y_rt_full, 1)
    
    # --------------------------
    # Step 12: attribuisco il nowcast a M1/M2/M3 del trimestre target
    # --------------------------
    m_tr <- compute_m_tr(date_t, dates_q)    # 1,2,3 oppure NA
    if (is.na(m_tr)) next
    
    key <- as.character(date_t)
    
    if (m_tr == 1) now_M1[[key]] <- y_rt_last
    if (m_tr == 2) now_M2[[key]] <- y_rt_last
    if (m_tr == 3) now_M3[[key]] <- y_rt_last
  }
  
  return(list(
    Lproxy_fix  = Lproxy_fix,
    L_midas_fix = L_midas_fix,
    M1 = now_M1,
    M2 = now_M2,
    M3 = now_M3
  ))
}






