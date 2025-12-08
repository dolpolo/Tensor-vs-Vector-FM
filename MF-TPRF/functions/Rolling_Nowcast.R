# ==============================================================================
# NOWCAST
# ==============================================================================

# X_full <- X
# y_q <- y_q
# dates <- dates_m
# Unb  <- Unb    # deve contenere le colonne: Frequency, Type, M1
# Freq <- Freq    # deve contenere le colonne: Frequency, Type, M1
# Lmax = Lmax
# tt <- 304 # [M1]
# tt <- 305 # [M2]
# tt <- 306 # [M3]
# L_midas <- 1
# dates_q <- all_countries$data[[country]]$DatesQ
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
  
  now_M1 <- list()
  now_M2 <- list()
  now_M3 <- list()
  
  N      <- ncol(X_full)
  N_m    <- sum(Freq == "M")
  N_q    <- sum(Freq == "Q")
  
  # --------------------------------------------
  # 1. ROLLING NOWCAST LOOP sui mesi tt
  # --------------------------------------------
  for (tt in seq(t_start, t_end)) {
    
    date_t <- dates[tt]
    cat("\n>>> REAL-TIME at", as.character(date_t), "\n")
    
    # --------------------------
    # Step 1: Unbalanced cut (applica ritardi di pubblicazione)
    # --------------------------
    X_cut <- unbalancedness(
      X_full   = X_full,
      dates    = dates,
      Freq     = Freq,
      Unb      = Unb,
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
    init   <- init_XP_ER(X_std)
    X_init <- init$X_init
    r      <- init$r      # numero di fattori
    
    # --------------------------
    # Step 5: EM
    # --------------------------
    EM_out <- EM_algorithm(X_init, X_std, A, r, max_iter = 50, tol = 1e-4)
    X_em   <- EM_out$X_completed          # T_m_current x N
    T_m_cur <- nrow(X_em)
    
    # --------------------------
    # Step 6: separa mensili / trimestrali e aggrega
    # --------------------------
    X_m_em <- X_em[, Freq == "M", drop = FALSE]
    X_q_em <- X_em[, Freq == "Q", drop = FALSE]
    
    X_mq_em <- agg_mq(X_m_em, agg_m)      # T_q_em x N_m
    X_qq_em <- agg_qq(X_q_em, agg_q)      # T_q_em x N_q
    
    X_em_agg <- cbind(X_mq_em, X_qq_em)   # T_q_em x (N_m+N_q)
    T_q_em   <- nrow(X_em_agg)            # numero di trimestri che posso costruire dai dati
    
    # --------------------------
    # Step 7: quanti trimestri di PIL sono pubblicati a date_t?
    #         (PIL ritardato di 1 mese)
    # --------------------------
    idx_pub <- which(dates_q < date_t)    # trimestri usciti PRIMA del mese corrente
    if (length(idx_pub) < 2) {
      # non ho almeno 2 trimestri per stimare
      next
    }
    T_q_current <- tail(idx_pub, 1)       # ultimo trimestre di PIL osservato
    
    # uso il minimo tra PIL disponibile e trimestri aggregati da X
    T_q_current <- min(T_q_current, T_q_em)
    
    # --------------------------
    # Step 8: taglia y_q e X_lf/X_hf ai trimestri utilizzabili
    # --------------------------
    y_q_cut   <- y_q[1:T_q_current]                       # T_q_current
    X_lf_cut  <- X_em_agg[1:T_q_current, , drop = FALSE]  # T_q_current x N
    X_hf_cut  <- X_em[1:T_m_cur, , drop = FALSE]          # tutti i mesi fino a tt
    
    if (length(y_q_cut) < 2) next
    
    # --------------------------
    # Step 9: selezione Lproxy via ER su DATASET MENSILE
    #         (non più su [X_lf_cut, y_q_cut] trimestrali)
    # --------------------------
    
    # 9.1. Costruisci y_m_cut a frequenza mensile
    # Se hai già una ricostruzione monthly di y (Xiong & Pelger), usala qui, es:
    #   y_m_full <- y_recon_monthly  # T_m x 1
    #   y_m_cut  <- y_m_full[1:T_m_cur]
    # Altrimenti, come fallback, usa una espansione tipo (NA, NA, y_tau):
    
    y_m_expand <- as.numeric(kronecker(y_q_cut, c(NA, NA, 1)))  # lunghezza = 3*T_q_current
    y_m_cut    <- y_m_expand[1:T_m_cur]                         # tronco ai mesi disponibili
    
    # 9.2. Dataset mensile per ER: [X_hf_cut, y_m_cut]
    data_proxy <- cbind(X_hf_cut, y_m_cut)   # T_m_cur x (N+1)
    
    Y_out_std     <- standardize_with_na(data_proxy)
    Y_std         <- Y_out_std$X_std
    cov_proxy_out <- all_purpose_covariance(Y_std)
    ER_proxy_out  <- select_num_factors_ER(cov_proxy_out$Sigma_tilde)
    Lproxy        <- ER_proxy_out$r
    
    message("# [", as.character(date_t), "] Lproxy (monthly ER) selected: ", Lproxy)
    
    # --------------------------
    # Step 10: selezione lag U-MIDAS (L_midas) sul sotto-campione trimestrale
    # --------------------------
    lag_sel <- choose_UMIDAS_lag(
      X_lf   = X_lf_cut,    # trimestrale (T_q_current x N)
      X_hf   = X_hf_cut,    # mensile (T_m_cur x N)
      y_q    = y_q_cut,     # PIL trimestrale
      Lproxy = Lproxy,
    )
    
    L_midas <- lag_sel$lag_BIC
    message("# [", as.character(date_t), "] U-MIDAS lag selected: ", L_midas)
    
    # --------------------------
    # Step 11: TPRF su sotto-campione
    # --------------------------
    MF_TPRF_RT_out <- MF_TPRF(
      X_lf    = X_lf_cut,
      X_hf    = X_hf_cut,
      y_q     = y_q_cut,
      Lproxy  = Lproxy,
      L_midas = L_midas
    )
    
    y_rt_full <- MF_TPRF_RT_out$y_nowcast    # vettore mensile (lunghezza T_m_cur)
    y_rt_last <- tail(y_rt_full, 1)          # nowcast del trimestre target alla data_t
    
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
    M1 = now_M1,
    M2 = now_M2,
    M3 = now_M3
  ))
}




