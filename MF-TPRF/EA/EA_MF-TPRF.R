# ==============================================================================
# 0. SET WORKING DIRECTORY & PATHS
# ==============================================================================

path_main    <- "C:/Users/david/Desktop/Paper/Tensor-vs-Vector-FM/MF-TPRF"
setwd(path_main)

# Paths: data, functions, results, graphs
path_data    <- file.path(path_main, "EA/data")
path_func    <- file.path(path_main, "functions")
path_results <- file.path(path_main, "EA/results/outputs")
path_graph   <- file.path(path_main, "EA/results/graph")



# ==============================================================================
# 1. LOAD REQUIRED LIBRARIES
# ==============================================================================

## Data Handling
library(tidyverse)
library(lubridate)
library(abind)
library(dplyr)

## Time Series
library(tseries)
library(zoo)
library(MASS)

## I/O
library(readxl)

## Plotting
library(ggplot2)

## Resolve conflicts
library(conflicted)
conflict_prefer("select",  "dplyr", quiet = TRUE)
conflict_prefer("filter",  "dplyr", quiet = TRUE)



# ==============================================================================
# 2. SOURCE USER-DEFINED FUNCTIONS
# ==============================================================================

source(file.path(path_func, "data_preparation.R"))
source(file.path(path_func, "MF-TPRF.R"))
source(file.path(path_func, "EM_MF-TPRF.R"))
source(file.path(path_func, "Factor_Selection.R"))
source(file.path(path_func, "Rolling_Nowcast.R"))


# ==============================================================================
# 3.1 PARAMETERS & COUNTRY LIST & COUNTRY OF INTEREST SELECTION
# ==============================================================================

# Users's parameters
params <- list(
  start_est   = as.Date("2000-04-01"),   # questo dipende da come vogliamo trattare il primo trimestre nella crescita del PIL
  start_eval  = as.Date("2017-01-01"),
  end_eval    = as.Date("2025-10-01"),
  covid_start = as.Date("2020-03-01"),
  covid_end   = as.Date("2021-07-01"),
  covid_mask  = TRUE,    # boolean, not string
  covid_mask_m = TRUE,
  covid_mask_q = FALSE,
  target      = "GDP",   # quarterly target variable
  sel_method  = "none", # "none" | "corr_threshold" | "F-Test"
  
  # Variable selection parameters
  n_m        = 38,
  n_q        = 10,
  thr_m      = 0.10,
  thr_q      = 0.85,
  thr_F_test = 0.01
)

# Country List
countries <- c("DE", "FR", "IT", "ES")

# Selected country
country <- "FR"


# ==============================================================================
# 4. PREPARE DATA FOR ALL COUNTRIES AND EXCTARCT THE ONE OF INTEREST
# ==============================================================================
all_countries <- prepare_all_countries(
  countries       = countries,
  params          = params,
  path            = path_data,
  covid_mask      = params$covid_mask,
  covid_mask_m    = params$covid_mask_m,
  covid_mask_q    = params$covid_mask_q
)


# Extract selected country's data
data     <- all_countries$data[[country]]$Data
dates_m  <- all_countries$data[[country]]$Dates
dates_q  <- all_countries$data[[country]]$DatesQ
series   <- all_countries$data[[country]]$Series



# ==============================================================================
# 5. TARGET AND PREDICTORS IDENTIFICATION
# ==============================================================================

gdp_col <- all_countries$data[[country]]$target_col

y   <- as.matrix(data[, gdp_col, drop = FALSE])     # target

# The growth in the first quarter is 0 (NO EA before the 2000)
y_q <- as.matrix(y[!is.na(y)])                      # quarterly target

X   <- as.matrix(data[, -gdp_col, drop = FALSE])    # predictors


# ==============================================================================
# 6. METADATA EXTRACTION
# ==============================================================================

## Dimensions
N_m     <- all_countries$data[[country]]$nM
nQ_tot  <- all_countries$data[[country]]$nQ

# Quarterly series (all)
Q_series_all <- all_countries$data[[country]]$Series[(N_m + 1):(N_m + nQ_tot)]

# Remove target from quarterly list
q_cols <- which(!grepl(params$target, Q_series_all, ignore.case = TRUE))
N_q    <- length(q_cols)

# Total number of predictors
N <- N_m + N_q


agg_m   <- all_countries$data[[country]]$agg_m
agg_q   <- all_countries$data[[country]]$agg_q[q_cols]
agg     <- c(agg_m, agg_q)

Type_m  <- all_countries$data[[country]]$TypeM
Type_q  <- all_countries$data[[country]]$TypeQ[q_cols]
Type    <- c(Type_m, Type_q)

Freq_m  <- all_countries$data[[country]]$freq_m
Freq_q  <- all_countries$data[[country]]$freq_q[q_cols]
Freq    <- c(Freq_m, Freq_q)

Unb_m   <- all_countries$data[[country]]$unb_m
Unb_q   <- all_countries$data[[country]]$unb_q[q_cols]
Unb     <- c(Unb_m, Unb_q)

Class_m <- all_countries$data[[country]]$ClassM
Class_q <- all_countries$data[[country]]$ClassQ[q_cols]
Class   <- c(Class_m, Class_q)

# ==============================================================================
# 7. STANDARDIZATION OF PREDICTORS
# ==============================================================================

out_std <- standardize_with_na(X)
X_std   <- out_std$X_std



# ==============================================================================
# 8. EM ALGORITHM: RECONSTRUCT THE MATRIX OF PREDICTORS
# ==============================================================================

# Initialization using Xiong and Pelger "all Purposes Estimator"
init <- init_XP_ER(X_std)

X_init <- init$X_init
r      <- init$r

# for ES set just: 
# r <- 7

# Matrix A to map low frequency to high frequency variables
A <- A_list(X, N_q, agg_q)

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
    
    # Indici di osservazioni non-NA nella serie originale
    idx_obs_i <- which(!is.na(X_obs[, i]))
    
    # Valori osservati (mensili o trimestrali)
    x_obs_i <- matrix(X_obs[idx_obs_i, i], ncol = 1)
    
    # (Dovrebbe valere: length(idx_obs_i) == nrow(A_i))
    # se vuoi, puoi riattivare questo check:
    # if (length(idx_obs_i) != nrow(A_i)) {
    #   warning("E-step: mismatch tra length(x_obs_i) e nrow(A_i) per serie ", i)
    # }
    
    # Matrix to invert: A_i A_i'
    M <- A_i %*% t(A_i)
    
    # Diagnostica: condizionamento di M
    cond_M <- kappa(M)
    if (cond_M > 1e10) {
      warning("E-step: M molto mal condizionata per serie ", i,
              " (kappa = ", signif(cond_M, 3), ").")
    }
    
    # Constraint correction: pred_i + A_i' (A_i A_i')^{-1} (x_obs - A_i pred_i)
    correction <- t(A_i) %*% solve(M, (x_obs_i - A_i %*% pred_i))
    
    # Updated monthly series
    X_new[, i] <- pred_i + correction
  }
  
  X_new
}

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

# EM algorithm 
EM_out <- EM_algorithm(X_init, X_std, A, r, max_iter = 100, tol = 1e-4)

# Reconstruct the quarterly matrix of predictors
X_em   <- EM_out$X_completed
X_m_em <- X_em[,1:N_m]
X_q_em <- X_em[,(N_m+1):N]

X_mq_em <- agg_mq(X_m_em, agg_m)
X_qq_em <- agg_qq(X_q_em, agg_q)

X_em_agg <- cbind(X_mq_em, X_qq_em)


# ==============================================================================
# 8. FIND THE NUMBER OF FACTORS AND U-MIDAS LAGS IN THE WHOLE DATASET INCLUDING TARGET
# ==============================================================================

# Eigenvalue ratio by Ahn
Y_out_std     <- standardize_with_na(data)
Y_std         <- Y_out_std$X_std
cov_proxy_out <- all_purpose_covariance(Y_std)
ER_proxy_out  <- select_num_factors_ER(cov_proxy_out$Sigma_tilde)
Lproxy        <- ER_proxy_out$r

# Number of U-MIDAS LAG
lag_sel <- choose_UMIDAS_lag(
  X_lf   = X_em_agg,     # trimestrale (T_q x N)
  X_hf   = X_em,         # mensile (T_m x N)
  y_q    = y_q,          # GDP trimestrale (non centrato)
  Lproxy = Lproxy
)

L_midas <- lag_sel$lag_BIC


# ==============================================================================
# 9. MF-TPRF
# ==============================================================================

MF_TPRF_out <- MF_TPRF(
  X_lf    = X_em_agg,
  X_hf    = X_em,
  y_q     = y_q,
  Lproxy  = Lproxy,
  L_midas = L_midas
)


# ==============================================================================
# 14. SAVE ALL RESULTS TO COUNTRY-SPECIFIC FOLDER
# ==============================================================================

## 14.1 Create folder for the country
path_country <- file.path(path_results, country)
if (!dir.exists(path_country)) {
  dir.create(path_country, recursive = TRUE)
}

## 14.2 Build filename
file_out <- file.path(
  path_country,
  paste0(
    "MF_TPRF_",  country,
    "_sel-",     params$sel_method,
    "_Covid-",   params$covid_mask,
    "_Lproxy-",  Lproxy,
    "_L_midas-", L_midas,
    "_Nm-",      N_m,
    "_Nq-",      N_q,
    "_",         format(params$start_est, "%Y-%m"),
    "_to_",      format(params$end_eval, "%Y-%m"),
    ".rds"
  )
)

## 14.3 Save everything into one RDS
saveRDS(
  list(
    # Country info
    country = country,
    params  = params,
    
    # EM Results
    A_list     = A,
    X_init     = X_init,
    X_em       = X_em,
    F_EM       = EM_out$F,
    Phi_EM     = EM_out$Phi,
    EM_iter    = EM_out$iterations,
    EM_diffs   = EM_out$diffs,
    X_mq_em    = X_mq_em,
    X_qq_em    = X_qq_em,
    X_em_agg   = X_em_agg,
    
    # MF-3PRF Results
    UMIDAS_lag = lag_sel,
    MF_TPRF    = MF_TPRF_out
  ),
  file_out
)

cat("\n*** SAVED EM + MF-3PRF RESULTS TO ***\n", file_out, "\n")



# ==============================================================================
# 1. LOAD RESULTS (EM + MF-3PRF)
# ==============================================================================

path_country <- file.path(path_results, country)

file_results <- file.path(
  path_country,
  paste0(
    "MF_TPRF_",  country,
    "_sel-",     params$sel_method,
    "_Covid-",   params$covid_mask,
    "_Lproxy-",  Lproxy,
    "_L_midas-", L_midas,
    "_Nm-",      N_m,
    "_Nq-",      N_q,
    "_",         format(params$start_est, "%Y-%m"),
    "_to_",      format(params$end_eval, "%Y-%m"),
    ".rds"
  )
)

MF_TPRF_res <- readRDS(file_results)



# ==============================================================================
# 2. PARAMETERS AND TIME INDEXES
# ==============================================================================

# Prendo i parametri salvati nell'RDS (più robusto che usare l'oggetto in memoria)
params <- MF_TPRF_res$params        # <<<

covid_start <- params$covid_start
covid_end   <- params$covid_end

dates_m <- all_countries$data[[country]]$DatesM   # monthly dates
dates_q <- all_countries$data[[country]]$DatesQ   # quarterly dates

# Ground-truth quarterly target series
y_true_q <- y_q   # se y_q è già la serie target nel sample di valutazione
T_q_complete <- length(y_true_q)    # <<< numero di trimestri osservati


# ==============================================================================
# 3. EXTRACT MF-3PRF REAL-TIME OUTPUT
# ==============================================================================

y_now_full <- MF_TPRF_res$MF_TPRF$y_nowcast  # vettore mensile
T_m_full   <- length(y_now_full)

# Controllo che ci siano almeno 3*T_q mesi
if (T_m_full < 3 * T_q_complete) {
  stop("y_now_full ha meno mesi di 3 * T_q_complete.")
}

# Monthly dates corrispondenti a y_now_full
dates_m_full <- dates_m[1:T_m_full]

# ==============================================================================
# 4. BUILD DATA FRAMES FOR PLOTTING
# ==============================================================================

# Numero di mesi in-sample (trimestri completi)
n_in  <- 3 * T_q_complete
n_tot <- T_m_full

# Costruisce il vettore “period” che distingue in-sample da real-time
period_vec <- rep("in-sample", n_tot)       # default: tutti in-sample
if (n_tot > n_in) {
  period_vec[(n_in + 1):n_tot] <- "real-time"
}

df_now_full <- data.frame(
  date       = dates_m_full,
  y_now_full = y_now_full,
  period     = factor(period_vec, levels = c("in-sample", "real-time"))
)

df_quarterly <- data.frame(
  date   = dates_q,
  y_true = y_true_q
)


# ==============================================================================
# 5. PLOT: MONTHLY NOWCAST vs TRUE QUARTERLY GDP
# ==============================================================================

library(ggplot2)
library(dplyr)

plot_nowcast <- ggplot() +
  
  # --- Highlight COVID window
  annotate(
    "rect",
    xmin = covid_start, xmax = covid_end,
    ymin = -Inf, ymax = Inf,
    fill = "grey70", alpha = 0.15
  ) +
  
  # --- Monthly nowcast (line, con colore in base a in-sample / real-time)
  geom_line(
    data = df_now_full,
    aes(x = date, y = y_now_full, color = period),
    linewidth = 1.1
  ) +
  
  # --- Real-time last points (markers)
  geom_point(
    data = df_now_full %>% filter(period == "real-time"),
    aes(x = date, y = y_now_full),
    size = 2.5, color = "black"
  ) +
  
  # --- Quarterly true GDP
  geom_line(
    data = df_quarterly,
    aes(x = date, y = y_true, color = "Quarterly GDP"),
    linewidth = 1.2
  ) +
  
  scale_color_manual(
    values = c(
      "in-sample"     = "#1F77B4",  # blue
      "real-time"     = "#2ECC71",  # green
      "Quarterly GDP" = "#D62728"   # red
    ),
    name = ""
  ) +
  
  labs(
    title = paste0("MF-3PRF Monthly Nowcast (In-Sample + Real-Time) – ", country),
    x = "Date",
    y = "GDP Growth"
  ) +
  
  theme_minimal(base_size = 14)
# Se vuoi puoi tenere un ylim, ma io lo toglierei per non tagliare valori:
# + coord_cartesian(ylim = c(-0.05, 0.05))

print(plot_nowcast)

# ==============================================================================
# 6. PERFORMANCE METRICS — RMSFE (M1, M2, M3, IN-SAMPLE ONLY)
# ==============================================================================

# Nowcast solo sui trimestri per cui conosco il GDP (in-sample)
y_now_in <- y_now_full[1:n_in]    # primi 3*T_q_complete mesi

# Indexes per M1, M2, M3 dentro la parte in-sample
M1_idx <- seq(1, n_in, by = 3)
M2_idx <- seq(2, n_in, by = 3)
M3_idx <- seq(3, n_in, by = 3)

RMSFE_M1 <- sqrt(mean((y_true_q - y_now_in[M1_idx])^2, na.rm = TRUE))
RMSFE_M2 <- sqrt(mean((y_true_q - y_now_in[M2_idx])^2, na.rm = TRUE))
RMSFE_M3 <- sqrt(mean((y_true_q - y_now_in[M3_idx])^2, na.rm = TRUE))

cat("MF-3PRF RMSFE (Monthly, In-Sample) –", country, "\n",
    "--------------------------------------------\n",
    "RMSFE (M1) =", round(RMSFE_M1, 4), "\n",
    "RMSFE (M2) =", round(RMSFE_M2, 4), "\n",
    "RMSFE (M3) =", round(RMSFE_M3, 4), "\n\n")

# ==============================================================================
# 7. CONVERGENCE OF THE EM ALGORITHM
# ==============================================================================

# Se hai ricaricato dal file .rds:
em_diffs <- MF_TPRF_res$EM_diffs   # li avevi salvati come EM_diffs = EM_out$diffs

df_em <- data.frame(
  iter = seq_along(em_diffs),
  diff = em_diffs
)

plot_em <- ggplot(df_em, aes(x = iter, y = diff)) +
  geom_line() +
  geom_point() +
  labs(
    title = paste0("EM convergence – ", country),
    x = "Iteration",
    y = "Max change in X"
  ) +
  theme_minimal(base_size = 14)

print(plot_em)


# ==============================================================================
# 8. SAVE GRAPHS TO COUNTRY FOLDER
# ==============================================================================

path_graph_country <- file.path(path_graph, country)
if (!dir.exists(path_graph_country)) {
  dir.create(path_graph_country, recursive = TRUE)
}

# --- File nome per il nowcast ---
file_graph_now <- file.path(
  path_graph_country,
  paste0(
    "MF3PRF_Nowcast_", country,
    "_sel-",     params$sel_method,
    "_Covid-",   params$covid_mask,
    "_Lproxy-",  Lproxy,
    "_L_midas-", L_midas,
    "_Nm-",      N_m,
    "_Nq-",      N_q,
    ".png"
  )
)

# --- File nome per il plot EM ---
file_graph_em <- file.path(
  path_graph_country,
  paste0(
    "EM_Convergence_", country,
    "_sel-",     params$sel_method,
    "_Covid-",   params$covid_mask,
    "_Lproxy-",  Lproxy,
    "_L_midas-", L_midas,
    "_Nm-",      N_m,
    "_Nq-",      N_q,
    ".png"
  )
)

# Salva NOWCAST
ggsave(
  filename = file_graph_now,
  plot     = plot_nowcast,
  width    = 12,
  height   = 6,
  dpi      = 300
)

# Salva EM CONVERGENCE
ggsave(
  filename = file_graph_em,
  plot     = plot_em,
  width    = 8,
  height   = 5,
  dpi      = 300
)

cat("*** GRAPHS SAVED TO ***\n",
    file_graph_now, "\n",
    file_graph_em,  "\n")







# ==============================================================================
# 15. PSEUDO REAL-TIME FORECASTING EXERCISE (WITH RELEASE SCHEDULE)
# ==============================================================================

# Questo esercizio implementa il vero pseudo real-time MF-3PRF + EM,
# rispettando calendario dei rilasci, unbalancedness dati, revisione dei fattori
# e disponibilità delle variabili ad alta e bassa frequenza.

pseudo_realtime_raw <- pseudo_realtime_TPRF_EM(
  X_full  = X,         # full mixed-frequency predictor matrix (with NA)
  y_q     = y_q,       # true quarterly target
  params  = params,    # all model parameters (include start_eval, end_eval, Lmax_midas)
  dates   = dates_m,   # monthly dates (for release schedule)
  dates_q = dates_q,   # quarterly publication dates
  Freq    = Freq,      # vector: monthly/quarterly frequency classification
  Unb     = Unb,       # unbalancedness info (lag di pubblicazione in mesi)
  Type    = Type,      # flow/stock/custom ecc. (se serve a monte)
  agg_m   = agg_m,     # schema aggregazione mensile -> trimestrale
  agg_q   = agg_q      # schema aggregazione trimestrale
)

# piccola utility per rendere M1/M2/M3 in data frame puliti
list_to_df <- function(lst, tag) {
  if (length(lst) == 0L) {
    return(data.frame(
      date             = as.Date(character()),
      nowcast          = numeric(),
      month_in_quarter = character()
    ))
  }
  data.frame(
    date             = as.Date(names(lst)),
    nowcast          = as.numeric(unlist(lst)),
    month_in_quarter = tag,
    row.names        = NULL
  )
}


# ==============================================================================
# 16. SAVE REAL-TIME RESULTS (COUNTRY-SPECIFIC FOLDER)
# ==============================================================================

# Create directory for this country
path_country <- file.path(path_results, country)
if (!dir.exists(path_country)) {
  dir.create(path_country, recursive = TRUE)
}

# Construct output filename (molto parlante)
file_out <- file.path(
  path_country,
  paste0(
    "MF3PRF_pseudoRT_", country,
    "_sel-",        params$sel_method,
    "_Covid-",      params$covid_mask,
    "_Nm-",         sum(Freq == "M"),
    "_Nq-",         sum(Freq == "Q"),
    "_",            format(params$start_est, "%Y-%m"),
    "_eval_",       format(params$start_eval, "%Y-%m"),
    "_to_",         format(params$end_eval, "%Y-%m"),
    ".rds"
  )
)

# Remove existing file (fresh overwrite)
if (file.exists(file_out)) {
  file.remove(file_out)
}



df_M1 <- list_to_df(pseudo_realtime_raw$M1, "M1")
df_M2 <- list_to_df(pseudo_realtime_raw$M2, "M2")
df_M3 <- list_to_df(pseudo_realtime_raw$M3, "M3")

# tutti i nowcast insieme, ordinati per data e mese nel trimestre
df_pseudoRT_all <- rbind(df_M1, df_M2, df_M3)
if (nrow(df_pseudoRT_all) > 0L) {
  df_pseudoRT_all <- df_pseudoRT_all[order(df_pseudoRT_all$date,
                                           df_pseudoRT_all$month_in_quarter), ]
}

# Save results in modo super riconoscibile
saveRDS(
  list(
    country          = country,
    params           = params,
    dates_m          = dates_m,
    dates_q          = dates_q,
    y_q              = y_q,
    
    # output grezzo della funzione pseudo_realtime_TPRF_EM (liste M1/M2/M3)
    pseudo_realtime_raw = pseudo_realtime_raw,
    
    # versioni tidy dei nowcast
    pseudo_realtime_M1  = df_M1,
    pseudo_realtime_M2  = df_M2,
    pseudo_realtime_M3  = df_M3,
    pseudo_realtime_all = df_pseudoRT_all
  ),
  file = file_out
)

cat("\n*** PSEUDO REAL-TIME RESULTS SAVED TO ***\n", file_out, "\n")


# ==============================================================================
# 1. LOAD PSEUDO REAL-TIME RESULTS
# ==============================================================================

library(dplyr)
library(ggplot2)

path_country <- file.path(path_results, country)

file_realtime <- file.path(
  path_country,
  paste0(
    "MF3PRF_pseudoRT_", country,
    "_sel-",        params$sel_method,
    "_Covid-",      params$covid_mask,
    "_Nm-",         sum(Freq == "M"),
    "_Nq-",         sum(Freq == "Q"),
    "_",            format(params$start_est, "%Y-%m"),
    "_eval_",       format(params$start_eval, "%Y-%m"),
    "_to_",         format(params$end_eval, "%Y-%m"),
    ".rds"
  )
)

RT <- readRDS(file_realtime)

# uso sempre params/dates/y_q salvati nel file (più robusto)
params  <- RT$params
dates_q <- RT$dates_q
y_q     <- RT$y_q

# ==============================================================================
# 2. BUILD TIDY REAL-TIME NOWCAST DATAFRAME (M1–M2–M3)
# ==============================================================================

df_rt <- RT$pseudo_realtime_all %>%
  rename(
    type    = month_in_quarter,
    nowcast = nowcast,
    date    = date
  ) %>%
  arrange(date)

df_rt$type <- factor(df_rt$type, levels = c("M1", "M2", "M3"))

# ==============================================================================
# 3. BUILD TRUE QUARTERLY GDP SERIES (EVALUATION WINDOW)
# ==============================================================================

idx_eval <- which(dates_q >= params$start_eval & dates_q <= params$end_eval)

df_yq <- data.frame(
  date = dates_q[idx_eval],
  GDP  = y_q[idx_eval]
)

# ==============================================================================
# 4. VISUALIZATION — ROLLING REAL-TIME NOWCASTS (M1–M2–M3)
# ==============================================================================

plot_rt <- ggplot() +
  
  # Periodo COVID
  annotate(
    "rect",
    xmin = params$covid_start, xmax = params$covid_end,
    ymin = -Inf, ymax = Inf,
    fill = "grey80", alpha = 0.20
  ) +
  
  # Serie trimestrale GDP reale
  geom_line(
    data = df_yq,
    aes(x = date, y = GDP, color = "True GDP"),
    linewidth = 1.2
  ) +
  
  # Nowcast real-time M1/M2/M3
  geom_line(
    data = df_rt,
    aes(x = date, y = nowcast, color = type),
    linewidth = 1
  ) +
  geom_point(
    data = df_rt,
    aes(x = date, y = nowcast, color = type),
    size = 2
  ) +
  
  scale_color_manual(
    values = c(
      "True GDP" = "#D62728",   # rosso
      "M1"       = "#1F77B4",   # blu
      "M2"       = "#2ECC71",   # verde
      "M3"       = "#F1C40F"    # giallo
    ),
    name = "Series"
  ) +
  
  labs(
    title    = paste0("MF-3PRF Rolling Real-Time Nowcasts – ", country),
    subtitle = "M1: early • M2: mid-quarter • M3: end-quarter",
    x        = "Date",
    y        = "GDP Growth Rate"
  ) +
  
  theme_minimal(base_size = 15) +
  theme(
    legend.position  = "bottom",
    plot.title       = element_text(face = "bold"),
    plot.subtitle    = element_text(color = "gray30"),
    panel.grid.minor = element_blank()
  )+
  
  coord_cartesian(ylim = c(-0.05, 0.1))
# Se vuoi, puoi aggiungere un ylim:
# + coord_cartesian(ylim = c(-0.05, 0.05))

print(plot_rt)

# ==============================================================================
# 5. SAVE GRAPH TO COUNTRY FOLDER
# ==============================================================================

path_graph_country <- file.path(path_graph, country)
if (!dir.exists(path_graph_country)) {
  dir.create(path_graph_country, recursive = TRUE)
}

file_graph <- file.path(
  path_graph_country,
  paste0(
    "MF3PRF_RollingPseudoRT_M1M2M3_", country,
    "_sel-",      params$sel_method,
    "_Covid-",    params$covid_mask,
    "_LmaxMid-",  params$Lmax_midas,
    "_",          format(params$start_eval, "%Y-%m"),
    "_to_",       format(params$end_eval, "%Y-%m"),
    ".png"
  )
)

ggsave(
  filename = file_graph,
  plot     = plot_rt,
  width    = 12,
  height   = 6,
  dpi      = 300
)

cat("\n*** SAVED GRAPH TO ***\n", file_graph, "\n")














###################################################################


MF_TPRF <- function(X_lf, X_hf, y_q,
                    Lproxy       = 1,
                    L_midas      = 1,
                    p_y          = 0,      # numero di lags di y (AR(p))
                    signif_alpha = 0.05,   # livello di significatività
                    shrink_alpha = TRUE,   # azzera alphas non significative
                    use_robust   = TRUE    # usa varianze HC1 per i t-test
) {
  # --------------------------------------------------------
  # Preliminari
  # --------------------------------------------------------
  X_lf <- as.matrix(X_lf)   # trimestrali (T_q x N) - dopo EM + aggregazione
  X_hf <- as.matrix(X_hf)   # mensili (T_m x N)
  y_q  <- as.numeric(y_q)   # target trimestrale (T_q)
  T_q  <- length(y_q)
  T_m  <- nrow(X_hf)
  N    <- ncol(X_lf)
  
  if (L_midas < 1) stop("L_midas deve essere >= 1")
  if (L_midas > T_q) stop("L_midas non può superare T_q")
  if (p_y < 0) stop("p_y deve essere >= 0")
  
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
  #  Beta_full: (1+Lproxy) x N  -> trasposto in N x (1+Lproxy)
  # --------------------------------------------------------
  Z_tilde <- cbind(1, Z)                     # T_q x (1+Lproxy)
  XtX_1   <- t(Z_tilde) %*% Z_tilde          # (1+Lproxy) x (1+Lproxy)
  XtY_1   <- t(Z_tilde) %*% X_lf             # (1+Lproxy) x N
  
  Beta_full <- solve(XtX_1, XtY_1)           # (1+Lproxy) x N
  Beta_full <- t(Beta_full)                  # N x (1+Lproxy)
  
  # slope α_i (escludo intercept)
  Phi_hat <- Beta_full[, -1, drop = FALSE]   # N x Lproxy
  
  # --------------------------------------------------------
  # 1.b – Shrinkage delle α_i non significative (robust t-test)
  # --------------------------------------------------------
  if (shrink_alpha && Lproxy > 0) {
    
    # Provo a usare varianze robuste (HC1); se non ho 'sandwich', passo a classiche
    if (use_robust) {
      if (!requireNamespace("sandwich", quietly = TRUE)) {
        warning("Package 'sandwich' non trovato: passo a varianze classiche.")
        use_robust <- FALSE
      }
    }
    
    for (j in seq_len(N)) {   # per ogni variabile x_j:
      x_j <- X_lf[, j]
      
      # regressione con intercetta
      fit_j <- lm(x_j ~ Z)
      
      if (use_robust) {
        vcov_j <- sandwich::vcovHC(fit_j, type = "HC1")
      } else {
        vcov_j <- vcov(fit_j)
      }
      se_j <- sqrt(diag(vcov_j))
      
      beta_j <- coef(fit_j)
      t_j    <- beta_j / se_j
      
      # slope corrispondono alle componenti 2:(Lproxy+1)
      t_slopes <- t_j[-1]
      
      # soglia t critica (bilaterale)
      t_crit <- qt(1 - signif_alpha / 2, df = fit_j$df.residual)
      
      # indici (in 1..Lproxy) dei coefficienti non significativi
      non_sig <- which(abs(t_slopes) < t_crit)
      
      if (length(non_sig) > 0) {
        # azzero le corrispondenti α_jk
        Phi_hat[j, non_sig] <- 0
      }
    }
  }
  
  # --------------------------------------------------------
  # STEP 2 – Second Pass: fattori mensili
  #
  #  F_t = X_hf Φ_hat (Φ_hat' Φ_hat)^(-1)
  # --------------------------------------------------------
  XtX_2   <- t(Phi_hat) %*% Phi_hat          # Lproxy x Lproxy
  XtY_2   <- t(Phi_hat) %*% t(X_hf)          # Lproxy x T_m
  
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
  rem      <- T_m - 3 * T_q               # 0,1,2 (o 3) mesi del trimestre successivo
  F_next1  <- if (rem >= 1) F_hat[3 * T_q + 1, , drop = FALSE] else NULL
  F_next2  <- if (rem >= 2) F_hat[3 * T_q + 2, , drop = FALSE] else NULL
  F_next3  <- if (rem >= 3) F_hat[3 * T_q + 3, , drop = FALSE] else NULL
  
  # --------------------------------------------------------
  # STEP 3.2 – U-MIDAS + AR(p) con L_midas trimestri (corrente + lag)
  #
  # Modello di base (senza AR):
  #  y_τ = β0 + sum_{ℓ=0}^{L_midas-1} [ β_{ℓ,1}' F1_{τ-ℓ} 
  #                                     + β_{ℓ,2}' F2_{τ-ℓ}
  #                                     + β_{ℓ,3}' F3_{τ-ℓ} ] + η_τ
  #
  # Estendiamo con AR(p):
  #  y_τ = ... + sum_{j=1}^{p_y} γ_j y_{τ-j} + η_τ
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
      # blocco: [F1_{lag_q}, F2_{lag_q}, F3_{lag_q}]
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
  
  Xreg_raw <- do.call(rbind, rows_list)              # (T_q - L_midas + 1) x (3*K*L_midas)
  y_dep_raw <- y_q[L_midas:T_q]                      # stessa dimensione in righe
  
  n_raw <- length(y_dep_raw)
  
  # AR-lags di y (se p_y > 0)
  if (p_y > 0) {
    if (n_raw <= p_y) {
      stop("Non ci sono abbastanza osservazioni per stimare AR(", p_y, ") con L_midas = ", L_midas, ".", sep = "")
    }
    
    # y_{τ-j} per τ = L_midas..T_q e j=1..p_y
    y_lag_mat <- sapply(1:p_y, function(j) {
      y_q[(L_midas - j):(T_q - j)]
    })
    
    # allineo: scarto le prime p_y righe per avere tutte le lag definite
    keep_idx <- (p_y + 1):n_raw
    
    Xreg_eff   <- Xreg_raw[keep_idx, , drop = FALSE]
    y_dep_eff  <- y_dep_raw[keep_idx]
    y_lags_eff <- as.matrix(y_lag_mat[keep_idx, , drop = FALSE])
    
    X_tilde_3 <- cbind(1, Xreg_eff, y_lags_eff)
  } else {
    Xreg_eff   <- Xreg_raw
    y_dep_eff  <- y_dep_raw
    y_lags_eff <- NULL
    X_tilde_3  <- cbind(1, Xreg_eff)
  }
  
  XtX_3    <- t(X_tilde_3) %*% X_tilde_3
  XtY_3    <- t(X_tilde_3) %*% y_dep_eff
  
  beta_hat <- solve(XtX_3, XtY_3)
  beta_hat <- as.numeric(beta_hat)
  
  # split dei coefficienti: intercept, fattori, AR(y)
  n_factor_params <- 3 * K * L_midas
  idx_factor_end  <- 1 + n_factor_params
  
  beta0        <- beta_hat[1]
  beta_factors <- beta_hat[2:idx_factor_end]
  gamma_y      <- if (p_y > 0) beta_hat[(idx_factor_end + 1):length(beta_hat)] else numeric(0)
  
  beta_mat <- matrix(beta_factors,
                     nrow = L_midas,
                     ncol = 3 * K,
                     byrow = TRUE)
  
  # --------------------------------------------------------
  # STEP 4 – Nowcasting mensile con fattori + AR(p)
  #
  # y_nowcast: lunghezza T_m (un nowcast per ogni mese disponibile).
  #
  # Per ciascun trimestre τ = L_midas..T_q (storici):
  #   al mese m=1: usa solo F1_τ come "corrente" (ell=0), + lag dei fattori + AR(y)
  #   al mese m=2: F1_τ + F2_τ
  #   al mese m=3: F1_τ + F2_τ + F3_τ
  #
  # Per il trimestre T_q+1 (se rem>=1):
  #   usa F_next1/F_next2/F_next3 allo stesso modo, ed AR(y) con y_q reali
  # --------------------------------------------------------
  
  y_nowcast <- rep(NA_real_, T_m)
  
  # 4.a) Trimestri COMPLETI 1..T_q
  for (tau in L_midas:T_q) {
    
    month_idx <- ((tau - 1) * 3 + 1):(tau * 3)  # posizioni dei 3 mesi di τ
    
    for (m in 1:3) {   # m = 1,2,3
      
      contrib_f <- 0
      
      # contributo fattori
      for (ell_id in 1:L_midas) {
        ell   <- ell_id - 1
        lag_q <- tau - ell
        if (lag_q < 1) next
        
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
          
          contrib_f <- contrib_f + sum(F_qm * beta_block)
        }
      }
      
      # contributo AR(p) di y
      contrib_ar <- 0
      if (p_y > 0) {
        for (j in 1:p_y) {
          lag_idx <- tau - j
          if (lag_idx >= 1) {
            contrib_ar <- contrib_ar + gamma_y[j] * y_q[lag_idx]
          }
        }
      }
      
      y_nowcast[month_idx[m]] <- beta0 + contrib_f + contrib_ar
    }
  }
  
  # 4.b) Trimestre CORRENTE T_q+1 (mesi extra)
  if (rem > 0) {
    
    tau_curr   <- T_q + 1
    month_curr <- 3 * T_q + seq_len(rem)   # indici mesi disponibili del trimestre T_q+1
    
    for (m in 1:rem) {
      
      contrib_f <- 0
      
      for (ell_id in 1:L_midas) {
        ell   <- ell_id - 1
        lag_q <- tau_curr - ell
        
        # Per i lag (ell>=1) uso i trimestri COMPLETI 1..T_q
        if (ell_id >= 2 && (lag_q < 1 || lag_q > T_q)) next
        
        for (mm in 1:3) {
          
          if (ell_id == 1) {
            # trimestre corrente T_q+1: uso solo mesi <= m
            if (mm > m) next
            
            F_qm <- switch(
              mm,
              `1` = F_next1,
              `2` = F_next2,
              `3` = F_next3
            )
            if (is.null(F_qm)) next   # mese non ancora osservato
            
          } else {
            # lag sui trimestri passati 1..T_q
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
          
          contrib_f <- contrib_f + sum(F_qm * beta_block)
        }
      }
      
      # contributo AR(p) anche per il trimestre T_q+1
      contrib_ar <- 0
      if (p_y > 0) {
        for (j in 1:p_y) {
          lag_idx <- tau_curr - j
          if (lag_idx >= 1 && lag_idx <= T_q) {
            contrib_ar <- contrib_ar + gamma_y[j] * y_q[lag_idx]
          }
        }
      }
      
      y_nowcast[month_curr[m]] <- beta0 + contrib_f + contrib_ar
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
    gamma_y   = gamma_y,   # coefficienti AR di y
    p_y       = p_y,       # ordine AR
    y_nowcast = y_nowcast
  ))
}

MF_TPRF_RT_out <- MF_TPRF(
  X_lf    = X_em_agg,
  X_hf    = X_em,
  y_q     = y_q,
  Lproxy  = Lproxy,
  L_midas = L_midas,
  p_y     = 1,       # AR(1)
  shrink_alpha = TRUE,
  signif_alpha = 0.05,
  use_robust   = TRUE
)

library(ggplot2)
library(dplyr)

# Serie di nowcast mensili
y_now <- MF_TPRF_RT_out$y_nowcast
T_m   <- length(y_now)

df_now <- data.frame(
  date    = dates_m[1:T_m],
  nowcast = y_now
)

# Allineo il GDP trimestrale al 3° mese di ogni trimestre
idx_q <- seq(3, by = 3, length.out = length(y_q))  # 3,6,9,...
df_yq_plot <- data.frame(
  date = dates_m[idx_q],
  GDP  = y_q
)

plot_tprf <- ggplot() +
  geom_line(
    data = df_now,
    aes(x = date, y = nowcast, color = "Nowcast (monthly)"),
    linewidth = 1
  ) +
  geom_line(
    data = df_yq_plot,
    aes(x = date, y = GDP, color = "True GDP (quarterly)"),
    linewidth = 1.2
  ) +
  geom_point(
    data = df_yq_plot,
    aes(x = date, y = GDP, color = "True GDP (quarterly)"),
    size = 2
  ) +
  scale_color_manual(
    values = c(
      "Nowcast (monthly)"      = "#1F77B4",
      "True GDP (quarterly)"   = "#D62728"
    ),
    name = ""
  ) +
  labs(
    title = paste0("MF-TPRF – Nowcast mensile vs True GDP – ", country),
    x     = "Date",
    y     = "GDP"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    plot.title      = element_text(face = "bold")
  )

print(plot_tprf)




################

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
    agg_q,       # oggetto per aggregazione trimestrale (per X_q)
    p_y          = 1,        # ordine AR per y nel 3° step TPRF
    shrink_alpha = TRUE,     # shrink delle alphas non significative
    signif_alpha = 0.05,     # livello di test
    use_robust   = TRUE      # varianze HC1 per t-test
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
    init   <- init_XP_ER(X_std)
    X_init <- init$X_init
    r      <- init$r      # numero di fattori
    
    # --------------------------
    # Step 5: EM
    # --------------------------
    EM_out  <- EM_algorithm(X_init, X_std, A, r, max_iter = 100, tol = 1e-4)
    X_em    <- EM_out$X_completed          # T_m_current x N
    T_m_cur <- nrow(X_em)
    
    # --------------------------
    # Step 6: separa mensili / trimestrali e aggrega
    # --------------------------
    X_m_em <- X_em[, Freq == "M", drop = FALSE]
    X_q_em <- X_em[, Freq == "Q", drop = FALSE]
    
    X_mq_em <- agg_mq(X_m_em, agg_m)      # T_q_em x N_m
    X_qq_em <- agg_qq(X_q_em, agg_q)      # T_q_em x N_q
    
    X_em_agg <- cbind(X_mq_em, X_qq_em)   # T_q_em x (N_m+N_q)
    T_q_em   <- nrow(X_em_agg)            # numero di trimestri aggregati
    
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
    X_hf_cut  <- X_em[1:T_m_cur, , drop = FALSE]          # mesi 1..tt
    
    if (length(y_q_cut) < 2) next
    
    # --------------------------
    # Step 9: selezione Lproxy via ER su DATASET MENSILE
    # --------------------------
    y_m_expand <- as.numeric(kronecker(y_q_cut, c(NA, NA, 1)))  # 3*T_q_current
    y_m_cut    <- y_m_expand[1:T_m_cur]
    
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
    # Step 11: MF-TPRF su sotto-campione (con AR(p) + shrink_alpha)
    # --------------------------
    MF_TPRF_RT_out <- MF_TPRF(
      X_lf        = X_lf_cut,
      X_hf        = X_hf_cut,
      y_q         = y_q_cut,
      Lproxy      = Lproxy,
      L_midas     = L_midas,
      p_y         = p_y,
      shrink_alpha = shrink_alpha,
      signif_alpha = signif_alpha,
      use_robust   = use_robust
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

pseudo_realtime <- pseudo_realtime_TPRF_EM(
  X_full  = X,
  y_q     = y_q,
  params  = params,
  dates   = dates_m,
  dates_q = dates_q,
  Freq    = Freq,
  Unb     = Unb,
  Type    = Type,
  agg_m   = agg_m,
  agg_q   = agg_q,
  p_y          = 1,      # AR(1) come ti ha suggerito il prof
  shrink_alpha = TRUE,
  signif_alpha = 0.05,
  use_robust   = TRUE
)

# piccola utility per rendere M1/M2/M3 in data frame puliti
list_to_df <- function(lst, tag) {
  if (length(lst) == 0L) {
    return(data.frame(
      date             = as.Date(character()),
      nowcast          = numeric(),
      month_in_quarter = character()
    ))
  }
  data.frame(
    date             = as.Date(names(lst)),
    nowcast          = as.numeric(unlist(lst)),
    month_in_quarter = tag,
    row.names        = NULL
  )
}


# ==============================================================================
# 16. SAVE REAL-TIME RESULTS (COUNTRY-SPECIFIC FOLDER)
# ==============================================================================

# Create directory for this country
path_country <- file.path(path_results, country)
if (!dir.exists(path_country)) {
  dir.create(path_country, recursive = TRUE)
}

# Construct output filename (molto parlante)
file_out <- file.path(
  path_country,
  paste0(
    "MF3PRF_pseudoRT_", country,
    "p_y_ ", "1",      # AR(1) come ti ha suggerito il prof
    "shrink_alpha", "TRUE",
    "signif_alpha", "0.05",
    "use_robust", "TRUE",
    "_sel-",        params$sel_method,
    "_Covid-",      params$covid_mask,
    "_Nm-",         sum(Freq == "M"),
    "_Nq-",         sum(Freq == "Q"),
    "_",            format(params$start_est, "%Y-%m"),
    "_eval_",       format(params$start_eval, "%Y-%m"),
    "_to_",         format(params$end_eval, "%Y-%m"),
    ".rds"
  )
)

# Remove existing file (fresh overwrite)
if (file.exists(file_out)) {
  file.remove(file_out)
}

# Save results in modo super riconoscibile
saveRDS(
  list(
    country          = country,
    params           = params,
    dates_m          = dates_m,
    dates_q          = dates_q,
    y_q              = y_q,
    
    # output grezzo della funzione pseudo_realtime_TPRF_EM (liste M1/M2/M3)
    pseudo_realtime_raw = pseudo_realtime,
    
    # versioni tidy dei nowcast
    pseudo_realtime_M1  = df_M1,
    pseudo_realtime_M2  = df_M2,
    pseudo_realtime_M3  = df_M3,
    pseudo_realtime_all = df_pseudoRT_all
  ),
  file = file_out
)

cat("\n*** PSEUDO REAL-TIME RESULTS SAVED TO ***\n", file_out, "\n")

df_M1 <- list_to_df(pseudo_realtime_raw$M1, "M1")
df_M2 <- list_to_df(pseudo_realtime_raw$M2, "M2")
df_M3 <- list_to_df(pseudo_realtime_raw$M3, "M3")

# tutti i nowcast insieme, ordinati per data e mese nel trimestre
df_pseudoRT_all <- rbind(df_M1, df_M2, df_M3)
if (nrow(df_pseudoRT_all) > 0L) {
  df_pseudoRT_all <- df_pseudoRT_all[order(df_pseudoRT_all$date,
                                           df_pseudoRT_all$month_in_quarter), ]
}



# ==============================================================================
# 1. LOAD PSEUDO REAL-TIME RESULTS
# ==============================================================================

library(dplyr)
library(ggplot2)

path_country <- file.path(path_results, country)

file_realtime <- file.path(
  path_country,
  paste0(
    "MF3PRF_pseudoRT_", country,
    "p_y_ ", "1",      # AR(1) come ti ha suggerito il prof
    "shrink_alpha", "TRUE",
    "signif_alpha", "0.05",
    "use_robust", "TRUE",
    "_sel-",        params$sel_method,
    "_Covid-",      params$covid_mask,
    "_Nm-",         sum(Freq == "M"),
    "_Nq-",         sum(Freq == "Q"),
    "_",            format(params$start_est, "%Y-%m"),
    "_eval_",       format(params$start_eval, "%Y-%m"),
    "_to_",         format(params$end_eval, "%Y-%m"),
    ".rds"
  )
)

RT <- readRDS(file_realtime)

# uso sempre params/dates/y_q salvati nel file (più robusto)
params  <- RT$params
dates_q <- RT$dates_q
y_q     <- RT$y_q

# ==============================================================================
# 2. BUILD TIDY REAL-TIME NOWCAST DATAFRAME (M1–M2–M3)
# ==============================================================================

df_rt <- RT$pseudo_realtime_all %>%
  rename(
    type    = month_in_quarter,
    nowcast = nowcast,
    date    = date
  ) %>%
  arrange(date)

df_rt$type <- factor(df_rt$type, levels = c("M1", "M2", "M3"))

# ==============================================================================
# 3. BUILD TRUE QUARTERLY GDP SERIES (EVALUATION WINDOW)
# ==============================================================================

idx_eval <- which(dates_q >= params$start_eval & dates_q <= params$end_eval)

df_yq <- data.frame(
  date = dates_q[idx_eval],
  GDP  = y_q[idx_eval]
)

# ==============================================================================
# 4. VISUALIZATION — ROLLING REAL-TIME NOWCASTS (M1–M2–M3)
# ==============================================================================

plot_rt <- ggplot() +
  
  # Periodo COVID
  annotate(
    "rect",
    xmin = params$covid_start, xmax = params$covid_end,
    ymin = -Inf, ymax = Inf,
    fill = "grey80", alpha = 0.20
  ) +
  
  # Serie trimestrale GDP reale
  geom_line(
    data = df_yq,
    aes(x = date, y = GDP, color = "True GDP"),
    linewidth = 1.2
  ) +
  
  # Nowcast real-time M1/M2/M3
  geom_line(
    data = df_rt,
    aes(x = date, y = nowcast, color = type),
    linewidth = 1
  ) +
  geom_point(
    data = df_rt,
    aes(x = date, y = nowcast, color = type),
    size = 2
  ) +
  
  scale_color_manual(
    values = c(
      "True GDP" = "#D62728",   # rosso
      "M1"       = "#1F77B4",   # blu
      "M2"       = "#2ECC71",   # verde
      "M3"       = "#F1C40F"    # giallo
    ),
    name = "Series"
  ) +
  
  labs(
    title    = paste0("MF-3PRF Rolling Real-Time Nowcasts – ", country),
    subtitle = "M1: early • M2: mid-quarter • M3: end-quarter",
    x        = "Date",
    y        = "GDP Growth Rate"
  ) +
  
  theme_minimal(base_size = 15) +
  theme(
    legend.position  = "bottom",
    plot.title       = element_text(face = "bold"),
    plot.subtitle    = element_text(color = "gray30"),
    panel.grid.minor = element_blank()
  )+
  
  coord_cartesian(ylim = c(-0.05, 0.1))
# Se vuoi, puoi aggiungere un ylim:
# + coord_cartesian(ylim = c(-0.05, 0.05))

print(plot_rt)

# ==============================================================================
# 5. SAVE GRAPH TO COUNTRY FOLDER
# ==============================================================================

path_graph_country <- file.path(path_graph, country)
if (!dir.exists(path_graph_country)) {
  dir.create(path_graph_country, recursive = TRUE)
}

file_graph <- file.path(
  path_graph_country,
  paste0(
    "MF3PRF_RollingPseudoRT_M1M2M3_", country,
    "p_y_ ", "1",      # AR(1) come ti ha suggerito il prof
    "shrink_alpha", "TRUE",
    "signif_alpha", "0.05",
    "use_robust", "TRUE",
    "_sel-",      params$sel_method,
    "_Covid-",    params$covid_mask,
    "_LmaxMid-",  params$Lmax_midas,
    "_",          format(params$start_eval, "%Y-%m"),
    "_to_",       format(params$end_eval, "%Y-%m"),
    ".png"
  )
)

ggsave(
  filename = file_graph,
  plot     = plot_rt,
  width    = 12,
  height   = 6,
  dpi      = 300
)

cat("\n*** SAVED GRAPH TO ***\n", file_graph, "\n")






######################################


# UMIDAS with AR 

choose_UMIDAS_lag <- function(X_lf, X_hf, y_q,
                              Lmax   = 9,
                              Lproxy = Lproxy,
                              p_y    = 0) {
  
  X_lf <- as.matrix(X_lf)   # trimestrali (T_q x N, dopo EM + aggregazione)
  X_hf <- as.matrix(X_hf)   # mensili (T_m x N)
  y_q  <- as.numeric(y_q)   # T_q
  
  T_q_full <- length(y_q)
  T_m      <- nrow(X_hf)
  
  # numero di trimestri COMPLETI che posso costruire dai mesi disponibili
  T_q_eff <- min(T_q_full, floor(T_m / 3))
  if (T_q_eff <= 2) {
    stop("Troppo pochi trimestri completi per selezionare L_midas.")
  }
  
  # uso solo i trimestri completi
  y_q_eff  <- y_q[1:T_q_eff]
  X_lf_eff <- X_lf[1:T_q_eff, , drop = FALSE]
  
  # -------------------------
  # STEP 0: Build L-autoproxy (con intercette gestite dentro)
  # -------------------------
  Z <- build_autoproxy_3prf(X_lf_eff, y_q_eff, Lproxy)   # T_q_eff x Lproxy
  
  # -------------------------
  # STEP 1: First pass (x_i,τ ~ intercept + Z_τ)
  # -------------------------
  Z_tilde <- cbind(1, Z)                      # T_q_eff x (Lproxy+1)
  
  A <- t(Z_tilde) %*% Z_tilde                 # (Lproxy+1) x (Lproxy+1)
  B <- t(Z_tilde) %*% X_lf_eff                # (Lproxy+1) x N
  
  B_full <- t(qr.solve(A, B))                 # N x (Lproxy+1)
  Phi_hat <- B_full[, -1, drop = FALSE]       # N x Lproxy (solo slope)
  
  # -------------------------
  # STEP 2: Second pass – fattori mensili
  # Ft = X_hf * Phi_hat * (Phi_hat' Phi_hat)^(-1)
  # -------------------------
  Phi_cross <- t(Phi_hat) %*% Phi_hat         # Lproxy x Lproxy
  Ft        <- X_hf %*% Phi_hat %*% qr.solve(Phi_cross)   # T_m x Lproxy
  
  # Considero solo i mesi dei trimestri completi: primi 3*T_q_eff mesi
  T_m_eff <- 3 * T_q_eff
  Ft_eff  <- Ft[1:T_m_eff, , drop = FALSE]
  
  # -------------------------
  # STEP 3: split monthly -> quarterly (3 mesi per trimestre)
  # -------------------------
  F1 <- Ft_eff[seq(1, T_m_eff, by = 3), , drop = FALSE]  # mese 1
  F2 <- Ft_eff[seq(2, T_m_eff, by = 3), , drop = FALSE]  # mese 2
  F3 <- Ft_eff[seq(3, T_m_eff, by = 3), , drop = FALSE]  # mese 3
  
  K <- ncol(F1)  # numero di fattori
  
  # -------------------------
  # LOOP over L (lag MF-UMIDAS) con AR(p_y) opzionale
  # -------------------------
  Lmax_eff <- min(Lmax, T_q_eff - 1)   # non posso usare L >= T_q_eff
  out <- data.frame(L = integer(), AIC = numeric(), BIC = numeric())
  
  for (L in 1:Lmax_eff) {
    
    # Costruisco Xreg_raw in modo coerente con MF_TPRF:
    # tau va da L .. T_q_eff
    rows_list <- list()
    row_id    <- 1
    
    for (tau in L:T_q_eff) {
      row_vec <- c()
      for (ell_id in 1:L) {
        ell   <- ell_id - 1             # ell = 0,...,L-1
        lag_q <- tau - ell
        # blocco: [F1_{lag_q}, F2_{lag_q}, F3_{lag_q}]
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
    
    Xreg_raw  <- do.call(rbind, rows_list)   # (T_q_eff - L + 1) x (3*K*L)
    y_dep_raw <- y_q_eff[L:T_q_eff]          # stessa lunghezza
    
    n_raw <- length(y_dep_raw)
    if (n_raw <= p_y) {
      next   # troppo poche obs per stimare AR(p_y)
    }
    
    if (p_y > 0) {
      # costruiamo la matrice delle lag di y
      y_lag_mat <- sapply(1:p_y, function(j) {
        y_q_eff[(L - j):(T_q_eff - j)]
      })
      
      # allineamento: scarto le prime p_y righe
      keep_idx <- (p_y + 1):n_raw
      
      Xreg_eff   <- Xreg_raw[keep_idx, , drop = FALSE]
      y_dep_eff  <- y_dep_raw[keep_idx]
      y_lags_eff <- as.matrix(y_lag_mat[keep_idx, , drop = FALSE])
      
      # regressione con intercetta + fattori + AR(p_y)
      df_reg <- data.frame(
        y    = y_dep_eff,
        Xreg = I(Xreg_eff),
        ylag = I(y_lags_eff)
      )
      fit <- lm(y ~ Xreg + ylag, data = df_reg)
      
    } else {
      # senza AR(p_y)
      Xreg_eff  <- Xreg_raw
      y_dep_eff <- y_dep_raw
      
      df_reg <- data.frame(
        y    = y_dep_eff,
        Xreg = I(Xreg_eff)
      )
      fit <- lm(y ~ Xreg, data = df_reg)
    }
    
    out <- rbind(
      out,
      data.frame(
        L   = L,
        AIC = AIC(fit),
        BIC = BIC(fit)
      )
    )
  }
  
  if (nrow(out) == 0) {
    stop("Nessun L valido per la selezione di U-MIDAS (problemi di campione vs p_y).")
  }
  
  list(
    results = out,
    lag_AIC = out$L[which.min(out$AIC)],
    lag_BIC = out$L[which.min(out$BIC)]
  )
}

p_y <- 1

lag_sel <- choose_UMIDAS_lag(
  X_lf   = X_lf_cut,
  X_hf   = X_hf_cut,
  y_q    = y_q_cut,
  Lproxy = Lproxy,
  p_y    = p_y                  # stesso p_y che usi in MF_TPRF (es. 1)
)

L_midas <- lag_sel$lag_BIC