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

## test
library(sandwich)
library(lmtest)
library(car)

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
  start_est    = as.Date("2000-04-01"),    # Accounting from missingenss in growth rates wrt 1999
  start_eval   = as.Date("2017-01-01"),
  end_eval     = as.Date("2025-10-01"),
  covid_start  = as.Date("2020-03-01"),
  covid_end    = as.Date("2021-07-01"),
  covid_mask_m = TRUE,                    # Boolean: if empty == TRUE
  covid_mask_q = FALSE,                    # Boolean: if empty == TRUE
  target       = "GDP",                    # quarterly target variable
  sel_method   = "corr_threshold",         # "none" | "corr_threshold" | "F-Test"
  
  # Variable selection parameters
  n_m        = 38,
  n_q        = 10,
  thr_m      = 0.10,
  thr_q      = 0.85,
  thr_F_test = 0.01,
  
  # MF-TPRF parameters
  p_AR        = 1,                         # order AR(p) in U-MIDAS regression of y
  Lmax        = 3,                         # Max numbers of lags in the Factors in U-MIDAS
  Robust_F    = TRUE,                      # Robust F test in the first step
  alpha       = 0.1,
  robust_type = "NW",                      # White--> HC | NW(Newey West)--> HAC
  nw_lag      = 1                          # Lag in the error autocorrelation in step 1
)

# Country List
countries <- c("DE", "FR", "IT", "ES")

# Selected country
country <- "DE"


# ==============================================================================
# 4. PREPARE DATA FOR ALL COUNTRIES AND EXCTARCT THE ONE OF INTEREST
# ==============================================================================
all_countries <- prepare_all_countries(
  countries       = countries,
  params          = params,
  path            = path_data,
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

# Matrix A to map low frequency to high frequency variables
A <- A_list(X, N_q, agg_q)

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
  X_lf        = X_em_agg,
  X_hf        = X_em,
  y_q         = y_q,
  Lmax        = params$Lmax,
  Lproxy      = Lproxy,
  p_AR        = params$p_AR,
  Robust_F    = params$Robust_F,
  alpha       = params$alpha,
  robust_type = params$robust_type,
  nw_lag      = params$nw_lag
)

L_midas <- lag_sel$lag_BIC

# ==============================================================================
# 9. MF-TPRF
# ==============================================================================

MF_TPRF_out <- MF_TPRF(
  X_lf        = X_em_agg,
  X_hf        = X_em,
  y_q         = y_q,
  Lproxy      = Lproxy,
  L_midas     = L_midas,
  p_AR        = params$p_AR,
  Robust_F    = params$Robust_F,
  alpha       = params$alpha,
  robust_type = params$robust_type,
  nw_lag      = params$nw_lag
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
    "MF_TPRF_",      country,
    "_sel-",         params$sel_method,
    "_Lproxy-",      Lproxy,
    "_L_midas-",     L_midas,
    "_Nm-",          N_m,
    "_Nq-",          N_q,
    "_p-AR",         params$p_AR,
    "_Robust-F_",    params$Robust_F,
    "_alpha-",       params$alpha,
    "_robust-type_", params$robust_type,
    "_nw-lag_",      params$nw_lag,
    "_Covid_m-",     params$covid_mask_m,
    "_Covid_q-",     params$covid_mask_q,
    "_",             format(params$start_est, "%Y-%m"),
    "_to_",          format(params$end_eval, "%Y-%m"),
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
    "_sel-",         params$sel_method,
    "_Lproxy-",      pseudo_realtime_raw$Lproxy_fix,
    "_Lmidas-",      pseudo_realtime_raw$L_midas_fix,
    "_Nm-",          N_m,
    "_Nq-",          N_q,
    "_p-AR",         params$p_AR,
    "_Robust-F_",    params$Robust_F,
    "_alpha-",       params$alpha,
    "_robust-type_", params$robust_type,
    "_nw-lag_",      params$nw_lag,
    "_Covid_m-",     params$covid_mask_m,
    "_Covid_q-",     params$covid_mask_q,
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

cat("file_out =\n", file_out, "\n")
nchar(file_out)
dir.exists(path_country)

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
