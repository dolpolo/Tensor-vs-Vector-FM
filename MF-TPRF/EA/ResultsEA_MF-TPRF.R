
# ==============================================================================
# 1. LOAD RESULTS (EM + MF-3PRF)
# ==============================================================================

path_country <- file.path(path_results, country)

file_results <- file.path(
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

MF_TPRF_res <- readRDS(file_results)



# ==============================================================================
# 2. PARAMETERS AND TIME INDEXES
# ==============================================================================

# Prendo i parametri salvati nell'RDS (più robusto che usare l'oggetto in memoria)
params <- MF_TPRF_res$params        # <<<

start_eval  <- params$start_eval
end_eval    <- params$end_eval
covid_start <- params$covid_start
covid_end   <- params$covid_end

dates_m <- all_countries$data[[country]]$DatesM   # monthly dates
dates_q <- all_countries$data[[country]]$DatesQ   # quarterly dates

# Ground-truth quarterly target series
y_true_q <- y_q   # se y_q è già la serie target nel sample di valutazione
T_q_complete <- length(y_true_q)    # <<< numero di trimestri osservati


# PRE, COVID, POST (MONTHLY)
dates_m_PRE <- dates_m[dates_m >= start_eval & dates_m < covid_start]
dates_m_COVID <- dates_m[dates_m >= covid_start & dates_m <= covid_end]
dates_m_POST <- dates_m[dates_m > covid_end & dates_m <= end_eval]

# PRE, COVID, POST (QUARTERLY)
dates_q_PRE <- dates_q[dates_q >= start_eval & dates_q < covid_start]
dates_q_COVID <- dates_q[dates_q >= covid_start & dates_q <= covid_end]
dates_q_POST <- dates_q[dates_q > covid_end & dates_q <= end_eval]

# Logical indices per periodi sui TRIMESTRI
is_PRE   <- (dates_q >= start_eval  & dates_q <  covid_start)
is_COVID <- (dates_q >= covid_start & dates_q <= covid_end)
is_POST  <- (dates_q >  covid_end   & dates_q <= end_eval)
is_ALL   <- (dates_q >= start_eval & dates_q <= end_eval)


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
# 6. PERFORMANCE METRICS — RMSFE (M1, M2, M3) per periodi
# ==============================================================================

# Nowcast solo sui trimestri per cui conosco il GDP (in-sample)
y_now_in <- y_now_full[1:n_in]    # primi 3*T_q_complete mesi

# Indexes per M1, M2, M3 dentro la parte in-sample
M1_idx <- seq(1, n_in, by = 3)
M2_idx <- seq(2, n_in, by = 3)
M3_idx <- seq(3, n_in, by = 3)

stopifnot(length(M1_idx) == length(y_true_q))

# -----------------------------
# Helper per calcolare RMSFE in un periodo
# -----------------------------
rmsfe_period <- function(mask_q, y_true_q, y_now_in, M_idx) {
  if (!any(mask_q)) return(NA_real_)  # nessun trimestre nel periodo
  y_true_sub <- y_true_q[mask_q]
  y_now_sub  <- y_now_in[M_idx[mask_q]]
  sqrt(mean((y_true_sub - y_now_sub)^2, na.rm = TRUE))
}

# FULL SAMPLE (ALL)
RMSFE_M1_ALL <- rmsfe_period(is_ALL,   y_true_q, y_now_in, M1_idx)
RMSFE_M2_ALL <- rmsfe_period(is_ALL,   y_true_q, y_now_in, M2_idx)
RMSFE_M3_ALL <- rmsfe_period(is_ALL,   y_true_q, y_now_in, M3_idx)

# PRE-COVID
RMSFE_M1_PRE <- rmsfe_period(is_PRE,   y_true_q, y_now_in, M1_idx)
RMSFE_M2_PRE <- rmsfe_period(is_PRE,   y_true_q, y_now_in, M2_idx)
RMSFE_M3_PRE <- rmsfe_period(is_PRE,   y_true_q, y_now_in, M3_idx)

# COVID
RMSFE_M1_COV <- rmsfe_period(is_COVID, y_true_q, y_now_in, M1_idx)
RMSFE_M2_COV <- rmsfe_period(is_COVID, y_true_q, y_now_in, M2_idx)
RMSFE_M3_COV <- rmsfe_period(is_COVID, y_true_q, y_now_in, M3_idx)

# POST-COVID
RMSFE_M1_POST <- rmsfe_period(is_POST, y_true_q, y_now_in, M1_idx)
RMSFE_M2_POST <- rmsfe_period(is_POST, y_true_q, y_now_in, M2_idx)
RMSFE_M3_POST <- rmsfe_period(is_POST, y_true_q, y_now_in, M3_idx)

cat("MF-3PRF RMSFE by period –", country, "\n",
    "--------------------------------------------\n",
    "FULL sample: \n",
    "  M1 =", round(RMSFE_M1_ALL, 4),
    "  M2 =", round(RMSFE_M2_ALL, 4),
    "  M3 =", round(RMSFE_M3_ALL, 4), "\n",
    "PRE-COVID:\n",
    "  M1 =", round(RMSFE_M1_PRE, 4),
    "  M2 =", round(RMSFE_M2_PRE, 4),
    "  M3 =", round(RMSFE_M3_PRE, 4), "\n",
    "COVID:\n",
    "  M1 =", round(RMSFE_M1_COV, 4),
    "  M2 =", round(RMSFE_M2_COV, 4),
    "  M3 =", round(RMSFE_M3_COV, 4), "\n",
    "POST-COVID:\n",
    "  M1 =", round(RMSFE_M1_POST, 4),
    "  M2 =", round(RMSFE_M2_POST, 4),
    "  M3 =", round(RMSFE_M3_POST, 4), "\n\n")


# ==============================================================================
# 7. TABELLA LATEX RMSFE (M1, M2, M3) per periodo
# ==============================================================================

rmsfe_mat <- rbind(
  "Full sample" = c(RMSFE_M1_ALL, RMSFE_M2_ALL, RMSFE_M3_ALL),
  "Pre-COVID"   = c(RMSFE_M1_PRE, RMSFE_M2_PRE, RMSFE_M3_PRE),
  "COVID period"= c(RMSFE_M1_COV, RMSFE_M2_COV, RMSFE_M3_COV),
  "Post-COVID"  = c(RMSFE_M1_POST, RMSFE_M2_POST, RMSFE_M3_POST)
)

colnames(rmsfe_mat) <- c("M1", "M2", "M3")

# Costruzione stringa LaTeX
latex_tab <- paste0(
  "\\begin{table}[htbp]\n",
  "\\centering\n",
  sprintf("\\caption{MF-3PRF RMSFE by period -- %s}\n", country),
  sprintf("\\label{tab:RMSFE_MF3PRF_%s}\n", country),
  "\\begin{tabular}{lccc}\n",
  "\\toprule\n",
  "Period & M1 & M2 & M3 \\\\\n",
  "\\midrule\n",
  paste(
    sprintf("%s & %.4f & %.4f & %.4f \\\\",
            rownames(rmsfe_mat),
            rmsfe_mat[, 1], rmsfe_mat[, 2], rmsfe_mat[, 3]),
    collapse = "\n"
  ),
  "\n\\bottomrule\n",
  "\\end{tabular}\n",
  "\\end{table}\n"
)

cat(latex_tab)

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
    "_sel-",           params$sel_method,
    "_Lproxy-",        Lproxy,
    "_L_midas-",       L_midas,
    "_Nm-",            N_m,
    "_Nq-",            N_q,
    "_p-AR",           params$p_AR,
    "_Robust-F_",      params$Robust_F,
    "_alpha-",         params$alpha,
    "_robust-type_",   params$robust_type,
    "_nw-lag_",        params$nw_lag,
    "_Covid_m-",       params$covid_mask_m,
    "_Covid_q-",       params$covid_mask_q,
    "_",               format(params$start_est, "%Y-%m"),
    "_to_",            format(params$end_eval, "%Y-%m"),
    ".png"
  )
)

# --- File nome per il plot EM ---
file_graph_em <- file.path(
  path_graph_country,
  paste0(
    "EM_Convergence_", country,
    "_sel-",           params$sel_method,
    "_Lproxy-",        Lproxy,
    "_L_midas-",       L_midas,
    "_Nm-",            N_m,
    "_Nq-",            N_q,
    "_p-AR",           params$p_AR,
    "_Robust-F_",      params$Robust_F,
    "_alpha-",         params$alpha,
    "_robust-type_",   params$robust_type,
    "_nw-lag_",        params$nw_lag,
    "_Covid_m-",       params$covid_mask_m,
    "_Covid_q-",       params$covid_mask_q,
    "_",               format(params$start_est, "%Y-%m"),
    "_to_",            format(params$end_eval, "%Y-%m"),
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




################################################################################
################################################################################
################################################################################


# ==============================================================================
# 1. LOAD PSEUDO REAL-TIME RESULTS
# ==============================================================================

library(dplyr)
library(ggplot2)

path_country <- file.path(path_results, country)

# file_realtime <- file.path(
#  path_country,
#  "MF3PRF_pseudoRT_ES_sel-corr_threshold_Covid-TRUE_Nm-30_Nq-9_2000-04_eval_2017-01_to_2025-10.rds"
#)

file_realtime <- file.path(
  path_country,
  paste0(
    "MF3PRF_pseudoRT_", country,
    "_sel-",         params$sel_method,
    "_Lproxy-",      pseudo_realtime_raw$Lproxy_fix,     # !!! BISOGNA TROVARE UN MODO PER CHIAMARLO DAL PRINCIPIO
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
  
  coord_cartesian(ylim = c(-0.05, 0.05))
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



# ==============================================================================
# 6. PERFORMANCE METRICS — RMSFE (M1, M2, M3) per periodi (rolling pseudo-RT)
# ==============================================================================

library(dplyr)
library(tidyr)
library(lubridate)

# 6.1 Aggiungo il periodo (Pre / COVID / Post) e l'ID di trimestre alla serie trimestrale vera
df_yq <- df_yq %>%
  mutate(
    period = dplyr::case_when(
      date <  params$covid_start                            ~ "Pre-COVID",
      date >= params$covid_start & date <= params$covid_end ~ "COVID period",
      date >  params$covid_end                              ~ "Post-COVID",
      TRUE                                                  ~ NA_character_
    ),
    quarter_id = paste0(year(date), "Q", quarter(date))
  )

# 6.2 Aggiungo lo stesso ID di trimestre ai nowcast rolling M1/M2/M3
df_rt <- df_rt %>%
  mutate(
    quarter_id = paste0(year(date), "Q", quarter(date))
  )

# 6.3 Join per trimestre: per ogni quarter_id ho M1, M2, M3 + GDP + periodo
df_eval <- df_rt %>%
  inner_join(
    df_yq %>% select(quarter_id, GDP, period),
    by = "quarter_id"
  ) %>%
  arrange(date)

# 6.4 RMSFE per periodo e tipo (M1/M2/M3)
rmsfe_by_period <- df_eval %>%
  group_by(period, type) %>%
  summarise(
    RMSFE = sqrt(mean((nowcast - GDP)^2, na.rm = TRUE)),
    .groups = "drop"
  )

# 6.5 RMSFE FULL sample (tutti i trimestri nella finestra di valutazione)
rmsfe_full <- df_eval %>%
  group_by(type) %>%
  summarise(
    RMSFE = sqrt(mean((nowcast - GDP)^2, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(period = "Full sample")

# 6.6 Metto insieme full + periodi e sistemo fattori e ordine
rmsfe_all <- bind_rows(rmsfe_full, rmsfe_by_period) %>%
  filter(period %in% c("Full sample", "Pre-COVID", "COVID period", "Post-COVID")) %>%
  mutate(
    period = factor(period,
                    levels = c("Full sample", "Pre-COVID", "COVID period", "Post-COVID")),
    type   = factor(type, levels = c("M1", "M2", "M3"))
  ) %>%
  arrange(period, type)

# 6.7 Wide: righe = periodi, colonne = M1/M2/M3
rmsfe_mat_df <- rmsfe_all %>%
  tidyr::pivot_wider(names_from = type, values_from = RMSFE) %>%
  arrange(period)

rmsfe_mat <- as.matrix(rmsfe_mat_df[, c("M1", "M2", "M3")])
rownames(rmsfe_mat) <- rmsfe_mat_df$period

# ==============================================================================
# 7. TABELLA LATEX RMSFE (M1, M2, M3) per periodo (rolling pseudo-RT)
# ==============================================================================

latex_tab <- paste0(
  "\\begin{table}[htbp]\n",
  "\\centering\n",
  sprintf("\\caption{MF-3PRF Rolling RMSFE by period -- %s}\n", country),
  sprintf("\\label{tab:RMSFE_MF3PRF_Rolling_%s}\n", country),
  "\\begin{tabular}{lccc}\n",
  "\\toprule\n",
  "Period & M1 & M2 & M3 \\\\\n",
  "\\midrule\n",
  paste(
    sprintf("%s & %.4f & %.4f & %.4f \\\\",
            rownames(rmsfe_mat),
            rmsfe_mat[, "M1"], rmsfe_mat[, "M2"], rmsfe_mat[, "M3"]),
    collapse = "\n"
  ),
  "\n\\bottomrule\n",
  "\\end{tabular}\n",
  "\\end{table}\n"
)

cat(latex_tab)
