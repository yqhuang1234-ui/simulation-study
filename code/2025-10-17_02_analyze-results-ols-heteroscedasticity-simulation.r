## ================================
## 0. Load needed packages and data
## ================================
# Helper function: install if not already available
load_or_install <- function(pkg){
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

# List of needed packages
packages <- c("dplyr", "ggplot2", "knitr", "kableExtra", "webshot2","magick")

# Load or install each
lapply(packages, load_or_install)

true_beta <- 0.4
# Load simulation results
path <- "~/Library/CloudStorage/Dropbox/School/CU/fall 2025/BIOS 6618 adv biostatistical method/midterm simulation/simulation-study"
fits_file <- '/data/2025-10-18_centered-log-linear-fits_ols-heteroscedasticity-simulation_n100_reps1000_seed1234_c0-10-15-20.rds'
fits <- readRDS(paste0(path, fits_file))
# Extract parameter settings from filename for labeling outputs
param_data <- sub(".*_(n[0-9].*).rds$", "\\1", basename(fits_file))
message("Analyzing results for parameter settings: ", param_data)
# "n200_reps1000_seed1234_c0-0.5-1-2-4-8"

# ISO 8601 date (date only):
iso_date <- format(Sys.Date(), "%Y-%m-%d")
save_results_path <- paste0(path,"/results/",iso_date)

if (!dir.exists(save_results_path)) {
  dir.create(save_results_path, recursive = TRUE)
}

## ================================
## 1. Summarize results with MC-SEs
## ================================
fits_summ <- fits %>%
  mutate(
    covered   = as.integer(true_beta >= ci_lower & true_beta <= ci_upper),
    rejected  = as.integer(p_value <= 0.05),
    est_err   = beta1_hat - true_beta
  ) %>%
  group_by(c_param) %>%
  summarise(
    R = n(),  # number of replications per c
    
    # Bias and its MC-SE
    bias         = mean(est_err, na.rm = TRUE),
    mcse_bias    = sd(est_err, na.rm = TRUE) / sqrt(R),
    
    # SD of beta1_hat across reps, with delta-method MC-SE: sd / sqrt(2*(R-1))
    sd_beta1_hat = sd(beta1_hat, na.rm = TRUE),
    mcse_sd_beta1 = sd_beta1_hat / sqrt(2 * pmax(R - 1, 1)),
    
    # Mean (reported) SE and its MC-SE
    se_beta1_hat = mean(se, na.rm = TRUE),
    mcse_se_mean = sd(se, na.rm = TRUE) / sqrt(R),
    
    # Coverage and its MC-SE
    ci_coverage  = mean(covered, na.rm = TRUE),
    mcse_cov     = sqrt(ci_coverage * (1 - ci_coverage) / R),
    
    # Reject rate (Type I error if true_beta==0; otherwise power) + MC-SE
    reject_rate  = mean(rejected, na.rm = TRUE),
    mcse_reject  = sqrt(reject_rate * (1 - reject_rate) / R),
    
    .groups = "drop"
  )
fits_summ

## ================================
## 2. Create figure with metrics and 95% MC-SE error bars 
## ================================
# Create table with formatted metrics including MC-SEs
fits_table <- fits_summ %>%
  mutate(
    bias          = sprintf("%.4f (±%.4f)", bias, mcse_bias),
    sd_beta1_hat  = sprintf("%.4f (±%.4f)", sd_beta1_hat, mcse_sd_beta1),
    se_beta1_hat  = sprintf("%.4f (±%.4f)", se_beta1_hat, mcse_se_mean),
    ci_coverage   = sprintf("%.3f (±%.3f)", ci_coverage, mcse_cov),
    reject_rate   = sprintf("%.3f (±%.3f)", reject_rate, mcse_reject)
  ) %>%
  select(c_param, bias, sd_beta1_hat, se_beta1_hat, ci_coverage, reject_rate)

# Create kable with new column labels (without touching fits_summ)
tbl <- fits_table %>%
  kable("html",
        caption = "Simulation Results with Monte Carlo Standard Errors",
  col.names = c("c", "Bias", "SD(β̂1)", "Mean SE", "Coverage", "Power")) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"),
                full_width = FALSE, position = "center")


# Reshape data to long format for plotting
metrics_long <- bind_rows(
  fits_summ %>% transmute(c_param, metric = "Bias",           est = bias,           mcse = mcse_bias),
  fits_summ %>% transmute(c_param, metric = "SD(β̂1)",        est = sd_beta1_hat,   mcse = mcse_sd_beta1),
  fits_summ %>% transmute(c_param, metric = "Mean SE",        est = se_beta1_hat,   mcse = mcse_se_mean),
  fits_summ %>% transmute(c_param, metric = "Coverage",       est = ci_coverage,    mcse = mcse_cov),
  fits_summ %>% transmute(c_param, metric = "Power",   est = reject_rate,    mcse = mcse_reject)
) %>%
  mutate(
    lo = est - 1.96*mcse,
    hi = est + 1.96*mcse,
    ref = case_when(
      metric == "Bias"         ~ 0,
      metric == "Coverage"     ~ 0.95,
      TRUE ~ NA_real_
    ),
    metric = factor(metric, levels=c("Bias","SD(β̂1)","Mean SE","Coverage","Power")
  ))

p <- ggplot(metrics_long, aes(x = c_param, y = est)) +
  geom_hline(aes(yintercept = ref), linetype = "dashed", na.rm = TRUE) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0, linewidth = 0.5) +
  geom_point(size = 2) +
  facet_wrap(~ metric, scales = "free_y", ncol = 2) +
  labs(
    title = "Metrics with 95% Monte Carlo Error Bars",
    x = "Heteroscedasticity parameter (c)",
    y = "Estimate (point ± 1.96 × MC-SE)"
  ) +
  theme_minimal(base_size = 12)

print(p)

## ================================
## 3. Save results table and figure
## ================================
table_save <- paste0(save_results_path,"/01_centered-log-linear-table-simulation-results-mcse_",param_data,'.png')
save_kable(tbl, table_save , zoom = 4, vwidth = 1000, vheight = 1200)
plot_save <- paste0(save_results_path,"/02_centered-log-linear-plots-simulation-results-mcse_",param_data,'.png')
ggsave(plot_save, plot = p, width = 8, height = 6, dpi = 500)
message("Saved results table to: ", table_save,
        "\nSaved results plot to: ", plot_save)
