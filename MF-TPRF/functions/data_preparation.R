# ==============================================================================
# STANDARDIZATION HANDLING NA
# ==============================================================================

standardize_with_na <- function(X) {
  # X : matrix T x N with NA
  
  # Compute column means and sds (ignoring NAs)
  col_means <- apply(X, 2, mean, na.rm = TRUE)
  col_sds   <- apply(X, 2, sd,   na.rm = TRUE)
  
  # Replace zero or NA std devs with 1 (prevents division by zero)
  col_sds[col_sds == 0 | is.na(col_sds)] <- 1
  
  # Standardize
  X_centered <- sweep(X, 2, col_means, "-")
  X_std      <- sweep(X_centered, 2, col_sds, "/")
  
  return(list(
    X_std   = X_std,
    mean    = col_means,
    sd      = col_sds
  ))
}



# ==============================================================================
# COVID PERIOD MASKING
# ==============================================================================

mask_covid <- function(X, dates, class, start, end) {
  
  real <- tolower(class) %in% c("r", "real")
  per  <- dates >= start & dates <= end
  
  # LOGICAL mask (important!)
  M <- outer(per, real, FUN = "&")
  
  X[M] <- NA
  
  list(data = X, mask = M)
}




# ==============================================================================
# MONTHLY AND QUARTERLY VARIABLES TYPES OF AGGREGATION
# ==============================================================================
# Xm <- Ym
# Type_m <- leg_type[1:ncol(Ym)]
# j <- 1
agg_mq <- function(Xm, agg) {
  
  Xm <- as.matrix(Xm)
  T  <- nrow(Xm)
  N  <- ncol(Xm)
  
  # Numero di trimestri completi
  T_q <- floor(T / 3)
  
  Xq <- matrix(NA, T_q, N)
  
  for (j in 1:N) {
    
    for (tau in 1:T_q) {
      
      m3 <- 3 * tau
      m2 <- m3 - 1
      m1 <- m3 - 2
      
      if (agg[j] == 2) {
        # FLOW → SOMMA
        Xq[tau, j] <- Xm[m1, j] + Xm[m2, j] + Xm[m3, j]
        
      } else if (agg[j] == 1) {
        # STOCK → MEDIA
        Xq[tau, j] <- mean(c(Xm[m1, j], Xm[m2, j], Xm[m3, j]), na.rm = TRUE)
        
      } else {
        stop("Aggregation code must be 1 (stock) or 2 (flow).")
      }
    }
  }
  
  colnames(Xq) <- colnames(Xm)
  return(Xq)
}

agg_qq <- function(Xm, agg) {
  
  Xm <- as.matrix(Xm)
  T  <- nrow(Xm)
  N  <- ncol(Xm)
  
  # numero di trimestri completi
  T_q <- floor(T / 3)
  
  Xq <- matrix(NA, T_q, N)
  
  for (j in 1:N) {
    
    for (tau in 1:T_q) {
      
      m3 <- 3 * tau
      m2 <- m3 - 1
      m1 <- m3 - 2
      
      if (agg[j] == 2) {
        # FLOW → somma
        Xq[tau, j] <- Xm[m1, j] + Xm[m2, j] + Xm[m3, j]
        
      } else if (agg[j] == 1) {
        # STOCK → media
        Xq[tau, j] <- mean(c(Xm[m1, j], Xm[m2, j], Xm[m3, j]), na.rm = TRUE)
      }
    }
  }
  
  colnames(Xq) <- colnames(Xm)
  return(Xq)
}





# ==============================================================================
# VARIABLE SELECTION BASED ON GDP CORRELATION
# ==============================================================================
# path   <- path_data
# covid_mask = TRUE
# cc <- "FR"

# Selects variables most correlated with GDP for each country.
# Returns base (country-independent) names for monthly and quarterly datasets.

select_vars <- function(countries, params, path) {
  
  target <- params$target
  
  # dates limits
  start_lim <- params$start_est
  end_lim   <- params$end_eval
  
  sel_m <- list()
  sel_q <- list()
  
  for (cc in countries) {
    
    # ---- Load data ----
    dm  <- read_excel(file.path(path, paste0(cc, "dataM_LT.xlsx")))
    dq  <- read_excel(file.path(path, paste0(cc, "dataQ_LT.xlsx")))
    
    # ---- select the correct time span
    dm  <-  dm %>%
      filter(Time >= start_lim & Time <= end_lim)
    dq  <-  dq %>%
      filter(Time >= start_lim & Time <= end_lim)
    
    # ---- Remove date column ----
    Ym <- as.data.frame(dm[, -1])
    Yq <- as.data.frame(dq[, -1])
     
    # ---- Legend ----
    leg <- read_excel(file.path(path, paste0(cc, "data.xlsx")), sheet = "info")
    TR   <- floor(leg$TR)
    freq <- leg$Frequency
    agg  <- leg$Aggregation
    
    # ---- Identify target quarterly variable ----
    target_id <- grep(paste0("^", target, "_"), colnames(Yq))[1]
    names_q   <- colnames(Yq)
    
    # ---- Match quarterly dataset with legend ----
    idx_leg_in_q <- match(tolower(names_q), tolower(leg$Name))
    
    # ---- Identify quarterly series (all Q variables) ----
    freq_q <- freq[idx_leg_in_q]
    need_shift <- which(freq_q == "Q")
    
    # --------------------------------------------------------
    # SHIFT TRIMESTRALI
    # --------------------------------------------------------
    
    Yq_shifted_full <- as.data.frame(
      apply(Yq, 2, function(col){
        Y <- col                      # valori trimestrali
        kronecker(Y, c(NA, NA, 1))         # monthly expansion: NA,NA,Y
      })
    )
    
    # estrai SOLO i mesi finali dei trimestri (3, 6, 9, ...)
    rows_quarter_end <- seq(3, nrow(Yq_shifted_full), by = 3)
    Yq_shifted <- Yq_shifted_full[rows_quarter_end, , drop = FALSE]
    
    # --------------------------------------------------------
    # MONTHLY → QUARTERLY aggregation
    # --------------------------------------------------------
    
    Y_mq <- agg_mq(Ym, agg)   # dimensione Tq
    
    Tq <- nrow(Y_mq)
    
    # ---- Target quarterly series (aligned) ----
    y <- as.numeric(Yq_shifted[[target_id]])[1:Tq]
    
    # ---- Selection Method ----
    method <- params$sel_method
    
    #──────────────────────────────────────────────
    # METHOD 1: NONE (pick top n correlations)
    #──────────────────────────────────────────────
    if (method == "none") {
      
      # -- Monthly --
      cm <- abs(cor(Y_mq, y, use="pairwise.complete.obs"))
      cm <- setNames(as.numeric(cm), colnames(Y_mq))
      vars_m_base <- names(sort(cm, decreasing=TRUE))[1:params$n_m]
      
      # -- Quarterly (shifted) --
      cq <- abs(cor(Yq_shifted[1:Tq, -target_id, drop=FALSE], y, use="pairwise.complete.obs"))
      cq <- setNames(as.numeric(cq), colnames(Yq_shifted)[-target_id])
      vars_q_base <- names(sort(cq, decreasing=TRUE))[1:params$n_q]
      
      #──────────────────────────────────────────────
      # METHOD 2: CORRELATION THRESHOLD
      #──────────────────────────────────────────────
    } else if (method == "corr_threshold") {
      
      cm <- abs(cor(Y_mq, y, use="pairwise.complete.obs"))
      cm <- setNames(as.numeric(cm), colnames(Y_mq))
      vars_m_base <- names(cm[cm >= params$thr_m])
      
      cq <- abs(cor(Yq_shifted[1:Tq, -target_id, drop=FALSE], y, use="pairwise.complete.obs"))
      cq <- setNames(as.numeric(cq), colnames(Yq_shifted)[-target_id])
      vars_q_base <- names(cq[cq >= params$thr_q])
      
      #──────────────────────────────────────────────
      # METHOD 3: F-TEST
      #──────────────────────────────────────────────
    } else if (method == "F-Test") {
      
      pvals_m <- sapply(1:ncol(Y_mq), function(j) {
        f <- summary(lm(Y_mq[, j] ~ y))$fstatistic
        pf(f[1], f[2], f[3], lower.tail = FALSE)
      })
      names(pvals_m) <- colnames(Y_mq)
      vars_m_base <- names(pvals_m[pvals_m <= params$thr_F_test])
      
      Yq_nt <- Yq_shifted[1:Tq, -target_id, drop=FALSE]
      
      pvals_q <- sapply(1:ncol(Yq_nt), function(j) {
        f <- summary(lm(Yq_nt[, j] ~ y))$fstatistic
        pf(f[1], f[2], f[3], lower.tail = FALSE)
      })
      names(pvals_q) <- colnames(Yq_nt)
      vars_q_base <- names(pvals_q[pvals_q <= params$thr_F_test])
    }
    
    #──────────────────────────────────────────────
    # FINALIZE
    #──────────────────────────────────────────────
    
    vars_m_cc <- vars_m_base
    vars_q_cc <- vars_q_base
    
    target_cc <- colnames(Yq)[target_id]
    
    sel_m[[cc]] <- vars_m_cc
    sel_q[[cc]] <- c(target_cc, vars_q_cc)
  }
  
  return(list(
    m = sel_m,
    q = sel_q
  ))
}


# ==============================================================================
# PREPARE DATA FOR A SINGLE COUNTRY
# ==============================================================================
# sel_m <- vars_sel$m
# sel_q <- vars_sel$q
prepare_country_data <- function(cc, params, sel_m, sel_q, path,
                                 covid_mask    = TRUE,
                                 covid_mask_m  = covid_mask,  # mask mensili
                                 covid_mask_q  = covid_mask)  # mask trimestrali
{
  # ============================================================
  # 1. LOAD DATA (MONTHLY LT + QUARTERLY LT)
  # ============================================================
  dm  <- read_excel(file.path(path, paste0(cc, "dataM_LT.xlsx")))
  dq  <- read_excel(file.path(path, paste0(cc, "dataQ_LT.xlsx")))
  leg <- read_excel(file.path(path, paste0(cc, "data.xlsx")), sheet = "info")
  
  # ============================================================
  # 2. LEGEND PARSING
  # ============================================================
  var_col   <- intersect("Name", names(leg))[1]
  leg_names <- tolower(trimws(leg[[var_col]]))
  leg_base  <- sub(paste0("_", tolower(cc), "$"), "", leg_names)
  
  TR        <- floor(leg$TR)
  leg_class <- leg$Class
  leg_type  <- leg$Type
  leg_freq  <- leg$Frequency
  leg_unb   <- leg$M1
  freq      <- leg$Frequency
  agg       <- leg$Aggregation
  
  # ============================================================
  # 3. DATES LIMITS
  # ============================================================
  start_lim <- params$start_est
  end_lim   <- params$end_eval
  
  # ---- select the correct time span
  dm <- dm %>%
    filter(Time >= start_lim & Time <= end_lim)
  dq <- dq %>%
    filter(Time >= start_lim & Time <= end_lim)
  
  # ============================================================
  # 4. EXTRACT MONTHLY DATA
  # ============================================================
  date_m <- as.Date(dm[[1]])
  Ym_raw <- as.matrix(dm[, -1, drop = FALSE])
  Xm     <- Ym_raw[, sel_m[[cc]], drop = FALSE]
  
  # metadata mensili
  base_m  <- sub(paste0("_", cc, "$"), "", sel_m[[cc]])
  class_m <- leg_class[match(tolower(base_m), leg_base)]
  type_m  <- leg_type[match(tolower(base_m), leg_base)]
  freq_m  <- leg_freq[match(tolower(base_m), leg_base)]
  unb_m   <- leg_unb[match(tolower(base_m), leg_base)]
  agg_m   <- agg[match(tolower(base_m), leg_base)]
  
  # COVID mask monthly (controllata da covid_mask_m)
  if (covid_mask_m) {
    res_m <- mask_covid(Xm, date_m, class_m,
                        params$covid_start, params$covid_end)
    Xm     <- res_m$data
    mask_m <- res_m$mask
  } else {
    mask_m <- matrix(FALSE, nrow(Xm), ncol(Xm))
  }
  
  # ============================================================
  # 5. EXTRACT QUARTERLY DATA
  # ============================================================
  date_q <- as.Date(dq[[1]])
  Xq_all <- as.matrix(dq[, -1, drop = FALSE])
  
  # colonne trimestrali effettivamente presenti
  valid_cols <- sel_q[[cc]][sel_q[[cc]] %in% colnames(Xq_all)]
  
  Xq <- Xq_all[, valid_cols, drop = FALSE]
  Xq <- matrix(as.numeric(Xq), nrow = nrow(Xq_all),
               dimnames = list(NULL, valid_cols))
  
  # metadata trimestrali
  base_q  <- sub(paste0("_", cc, "$"), "", valid_cols)
  class_q <- leg_class[match(tolower(base_q), leg_base)]
  type_q  <- leg_type[match(tolower(base_q), leg_base)]
  freq_q  <- leg_freq[match(tolower(base_q), leg_base)]
  unb_q   <- leg_unb[match(tolower(base_q), leg_base)]
  agg_q   <- agg[match(tolower(base_q), leg_base)]
  
  # --------------------------------------------------------
  # 5B. IDENTIFICA LA TARGET TRA LE TRIMESTRALI
  #     (senza assumere che sia la prima colonna)
  # --------------------------------------------------------
  target_pattern <- paste0("^", tolower(params$target), "_")
  target_idx_q   <- grep(target_pattern, tolower(colnames(Xq)))
  
  # vettore classi usato SOLO per la maschera (target esclusa)
  class_q_mask <- class_q
  if (length(target_idx_q) > 0) {
    class_q_mask[target_idx_q] <- "target"  # qualcosa ≠ "r"/"real"
  }
  
  # date del terzo mese del trimestre
  date_q_third <- date_q %m+% months(2)
  
  # ============================================================
  # 6. CREATE FULL MONTHLY TIME AXIS
  # ============================================================
  date_full <- seq.Date(
    from = min(date_m),
    to   = max(date_m),
    by   = "month"
  )
  
  # ============================================================
  # 7. ALIGN MONTHLY
  # ============================================================
  Xm_full <- matrix(NA_real_, length(date_full), ncol(Xm))
  rownames(Xm_full) <- as.character(date_full)
  colnames(Xm_full) <- colnames(Xm)
  Xm_full[as.character(date_m), ] <- Xm
  
  mask_m_full <- matrix(FALSE, length(date_full), ncol(Xm))
  rownames(mask_m_full) <- as.character(date_full)
  colnames(mask_m_full) <- colnames(Xm)
  mask_m_full[as.character(date_m), ] <- mask_m
  
  # ============================================================
  # 8. ALIGN QUARTERLY (portate al 3° mese del trimestre)
  # ============================================================
  Xq_full <- matrix(NA_real_, length(date_full), ncol(Xq))
  rownames(Xq_full) <- as.character(date_full)
  colnames(Xq_full) <- colnames(Xq)
  
  # inserisci ciascun trimestre nel 3° mese
  for (t in seq_along(date_q_third)) {
    key <- as.character(date_q_third[t])
    if (key %in% rownames(Xq_full)) {
      Xq_full[key, ] <- Xq[t, ]
    }
  }
  
  # ---- COVID mask trimestrali (target esclusa tramite class_q_mask) ----
  if (covid_mask_q) {
    res_q <- mask_covid(Xq_full, date_full, class_q_mask,
                        params$covid_start, params$covid_end)
    Xq_full     <- res_q$data
    mask_q_full <- res_q$mask
  } else {
    mask_q_full <- matrix(FALSE, nrow(Xq_full), ncol(Xq_full))
    rownames(mask_q_full) <- rownames(Xq_full)
    colnames(mask_q_full) <- colnames(Xq_full)
  }
  
  # ============================================================
  # 9. APPLY FINAL TIME WINDOW
  # ============================================================
  idx <- which(date_full >= start_lim & date_full <= end_lim)
  
  Xm_out     <- Xm_full[idx, , drop = FALSE]
  Xq_out     <- Xq_full[idx, , drop = FALSE]
  mask_m_out <- mask_m_full[idx, , drop = FALSE]
  mask_q_out <- mask_q_full[idx, , drop = FALSE]
  dates_out  <- date_full[idx]
  
  dates_m_out <- date_m[date_m >= start_lim & date_m <= end_lim]
  
  # Date dei trimestri (3° mese) dentro la finestra
  dates_q_all <- date_q_third[date_q_third >= start_lim & date_q_third <= end_lim]
  dates_q_out <- dates_q_all
  
  # ============================================================
  # 10. IDENTIFY TARGET COLUMN (su tutte le serie)
  # ============================================================
  all_series <- c(sel_m[[cc]], sel_q[[cc]])
  
  target_pattern_all <- paste0("^", tolower(params$target), "_")
  target_col         <- grep(target_pattern_all, tolower(all_series))
  target_name        <- all_series[target_col]
  
  # ============================================================
  # 11. FINAL OUTPUT
  # ============================================================
  list(
    Data        = cbind(Xm_out, Xq_out),
    Dates       = dates_out,
    DatesM      = dates_m_out,
    DatesQ      = dates_q_out,
    Series      = all_series,
    nM          = length(sel_m[[cc]]),
    nQ          = length(sel_q[[cc]]),
    agg_m       = agg_m,
    agg_q       = agg_q,
    ClassM      = class_m,
    ClassQ      = class_q,
    TypeM       = type_m,
    TypeQ       = type_q,
    MaskM       = mask_m_out,
    MaskQ       = mask_q_out,
    freq_m      = freq_m,
    freq_q      = freq_q,
    unb_m       = unb_m,
    unb_q       = unb_q,
    idx_q       = ncol(Xm_out) + 1,
    target_col  = target_col,
    target_name = target_name
  )
}


# ==============================================================================
# WRAPPER FUNCTION: COMPLETE MULTI-COUNTRY DATA PREPARATION
# ==============================================================================


prepare_all_countries <- function(countries, params, path,
                                  covid_mask=TRUE, covid_mask_m, covid_mask_q) {
  
  # 1. variable selection country-by-country
  vars_sel <- select_vars(countries, params, path)
  
  # 2. prepare each country
  out <- lapply(countries, function(cc) {
    prepare_country_data(
      cc        = cc,
      params    = params,
      sel_m     = vars_sel$m,
      sel_q     = vars_sel$q,
      path      = path,
      covid_mask = covid_mask,
      covid_mask_m    = params$covid_mask_m,
      covid_mask_q    = params$covid_mask_q
    )
  })
  
  names(out) <- countries
  
  list(
    data   = out,
    sel    = vars_sel
  )
}

