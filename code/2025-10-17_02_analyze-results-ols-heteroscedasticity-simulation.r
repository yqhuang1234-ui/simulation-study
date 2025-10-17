library(dplyr)
path <- "~/Library/CloudStorage/Dropbox/School/CU/fall 2025/BIOS 6618 adv biostatistical method/midterm simulation/simulation-study/data"
fits <- readRDS(paste0(path, '/2025-10-17_fits_ols-heteroscedasticity-simulation_n200_reps1000_seed1234_c0-0.5-1-2-4-8.rds'))
tail(fits)
true_beta <- 0.4
library(dplyr)

fits.grouped <- fits %>%
  mutate(
    covered = as.integer(true_beta >= ci_lower & true_beta <= ci_upper),
    rejected = as.integer(p_value <= 0.05),
    ci_width = ci_upper - ci_lower,
    est_err  = beta1_hat - true_beta
  ) %>%
  group_by(c_param) %>%
  summarise(
    R = n(),  # number of replications
    
    # Point estimate quality
    bias         = mean(est_err, na.rm = TRUE),
    mcse_bias    = sd(est_err, na.rm = TRUE) / sqrt(R),
    
    sd_beta1_hat = sd(beta1_hat, na.rm = TRUE),
    
    # Reported SE
    se_beta1_hat = mean(se, na.rm = TRUE),
    mcse_se_mean = sd(se, na.rm = TRUE) / sqrt(R),
    
    # CI performance
    ci_coverage  = mean(covered, na.rm = TRUE),
    mcse_cov     = sqrt(ci_coverage * (1 - ci_coverage) / R),
    
    mean_ci_lower = mean(ci_lower, na.rm = TRUE),
    mean_ci_upper = mean(ci_upper, na.rm = TRUE),
    mean_ci_width = mean(ci_width, na.rm = TRUE),
    mcse_ci_width = sd(ci_width, na.rm = TRUE) / sqrt(R),
    
    # Testing behavior
    reject_rate  = mean(rejected, na.rm = TRUE),
    mcse_reject  = sqrt(reject_rate * (1 - reject_rate) / R),
    
    .groups = "drop"
  )

fits.grouped
#------------------------
#Plot results
#------------------------
# install.packages(c("dplyr","ggplot2","tidyr"))
library(dplyr)
library(tidyr)
library(ggplot2)

# 1) Build per-c summary with MC-SEs
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

library(dplyr)
library(ggplot2)

# fits_summ must already exist with these columns:
# c_param, bias, mcse_bias, sd_beta1_hat, mcse_sd_beta1, 
# se_beta1_hat, mcse_se_mean, ci_coverage, mcse_cov, 
# reject_rate, mcse_reject

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
      metric == "Power" ~ 0.05,   
      TRUE ~ NA_real_
    ),
    metric = factor(metric, levels = c("Bias","SD(β̂1)","Mean SE","Coverage","Power/Type I"))
  )

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

# Optional if c spans orders of magnitude:
# p <- p + scale_x_continuous(trans = "log10", breaks = sort(unique(metrics_long$c_param)))

print(p)


