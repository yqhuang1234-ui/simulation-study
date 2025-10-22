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
packages <- c("dplyr", "ggplot2",  "kableExtra", "webshot2","magick","knitr","tidyverse","tidyr","ggh4x")
# Load or install each
lapply(packages, load_or_install)


# Load simulation results
path <- "~/Library/CloudStorage/Dropbox/School/CU/fall 2025/BIOS 6618 adv biostatistical method/midterm simulation/simulation-study/"
fits_file <- '/data/2025-10-22_optimized-seed-best-param-centered-log-linear-fits_ols-heteroscedasticity-simulation_n100_reps1000_seed1234_c0-0.3-1-6.rds'
simulated_data_file <- '/data/2025-10-22_optimized-seed-best-param-centered-log-linear-datasets_ols-heteroscedasticity-simulation_n100_reps1000_seed1234_c0-0.3-1-6.rds'
fits <- readRDS(paste0(path, fits_file))
simulated_data <- readRDS(paste0(path, simulated_data_file))
# Extract parameter settings from filename for labeling outputs
param_data <- sub(".*_(n[0-9].*).rds$", "\\1", basename(fits_file))
message("Analyzing results for parameter settings: ", param_data)
true_beta <- unique(fits$beta1_true)
c_params <- c(0,0.3,1,6)
# ISO 8601 date (date only):
iso_date <- format(Sys.Date(), "%Y-%m-%d")
save_results_path <- paste0(path,"/results/",iso_date)

if (!dir.exists(save_results_path)) {
  dir.create(save_results_path, recursive = TRUE)
}

## ================================
## 1. Summarize results
## ================================
# Summarize datasets to check variances
average_error_variance <- simulated_data %>%
  group_by(c_param) %>%
  summarise(
    var_y = round(var(y),0),                        # variance of Y
    mean_eps_var = round(mean(eps_var),0),          # average
    var_y_unnorm = round(var(y_unnorm),0),          
    mean_eps_var_unnorm = round(mean(eps_var_unnorm), 0),
    .groups = "drop"
  )

# Summarize simulation results with MC-SEs
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
    
    # Relative error of model SE vs empirical SD
    rel_err_se    = (se_beta1_hat - sd_beta1_hat) / sd_beta1_hat,
    
    # --- MCSE for rel_err_se via delta method ---
    # Notation: M = mean(se), S = sd(beta1_hat)
    # g(M, S) = M/S - 1; grad g = (1/S, -M/S^2)
    varM          = var(se, na.rm = TRUE) / R,
    varS          = (sd_beta1_hat^2) / (2 * pmax(R - 1, 1)),
    mcse_rel_err  = sqrt( (1 / sd_beta1_hat)^2 * varM +
                            (se_beta1_hat^2 / sd_beta1_hat^4) * varS ),
    # Coverage and its MC-SE
    ci_coverage  = mean(covered, na.rm = TRUE),
    mcse_cov     = sqrt(ci_coverage * (1 - ci_coverage) / R),
    
    
    
    .groups = "drop"
  )

## ================================
## 2. Create figure with metrics and 95% MC-SE error bars 
## ================================
# create variance table
variance_table <- average_error_variance %>%
  kbl(
    booktabs = TRUE,
    digits = 0,
    col.names = c("c", "Var(Y) normalized", "Mean error variance normalized", 
                  "Var(Y) unnormalized", "Mean error variance unnormalized")
  ) %>%
  kable_styling(latex_options = c("HOLD_position", "striped", "scale_down"))
variance_table
#-----------------------------------------------------
# Create table with formatted metrics including MC-SEs
## 1) Find rows using numeric values from fits_summ
row_min_relerr <- which.min(fits_summ$rel_err_se)      # rel_err_se is in proportion (not %)
row_min_cov    <- which.min(fits_summ$ci_coverage)     # proportion (0–1)
row_max_empse  <- which.max(fits_summ$sd_beta1_hat)

## 2) Build the display table (your code)
fits_table <- fits_summ %>%
  mutate(
    bias          = sprintf("%.3f (±%.3f)", bias, mcse_bias),
    se_beta1_hat  = sprintf("%.3f (±%.3f)", se_beta1_hat, mcse_se_mean),
    sd_beta1_hat  = sprintf("%.3f (±%.3f)", sd_beta1_hat, mcse_sd_beta1),
    rel_err_se    = sprintf("%.1f (±%.1f)", rel_err_se*100, mcse_rel_err*100),
    ci_coverage   = sprintf("%.1f (±%.1f)", ci_coverage*100, mcse_cov*100)
  ) %>%
  select(c_param, bias, se_beta1_hat, sd_beta1_hat, rel_err_se, ci_coverage)

## 3) Bold the target cells
fits_table_fmt <- fits_table %>%
  mutate(
    rel_err_se   = ifelse(row_number() == row_min_relerr,
                          cell_spec(rel_err_se, "html", bold = TRUE),
                          rel_err_se),
    ci_coverage  = ifelse(row_number() == row_min_cov,
                          cell_spec(ci_coverage, "html", bold = TRUE),
                          ci_coverage),
    sd_beta1_hat = ifelse(row_number() == row_max_empse,
                          cell_spec(sd_beta1_hat, "html", bold = TRUE),
                          sd_beta1_hat)
  )

## 4) Render (escape = FALSE to allow HTML styling)
tbl <- fits_table_fmt %>%
  kable("html",
        escape = FALSE,
        col.names = c("c", "Bias", "Model SE", "Empirical SE",
                      "Relative error in model SE(%)", "CI Coverage(%)")) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"),
                full_width = FALSE, position = "center")

tbl

#-----------------------------------------------------
# Create figure with metrics and 95% MC-SE confidence intervals
# Long data + factor levels + ref lines
metrics_long <- bind_rows(
  fits_summ %>% transmute(c_param, metric = "Bias",                 est = bias,             mcse = mcse_bias),
  fits_summ %>% transmute(c_param, metric = "Empirical SE",         est = sd_beta1_hat,     mcse = mcse_sd_beta1),
  fits_summ %>% transmute(c_param, metric = "Model SE",             est = se_beta1_hat,     mcse = mcse_se_mean),
  fits_summ %>% transmute(c_param, metric = "Relative Error in Model SE (%)",est = round(rel_err_se*100,1),   mcse = mcse_rel_err*100),
  fits_summ %>% transmute(c_param, metric = "Coverage(%)",          est = ci_coverage*100,  mcse = mcse_cov*100)
) %>%
  mutate(
    lo = est - 1.96 * mcse,
    hi = est + 1.96 * mcse,
    ref = case_when(
      metric == "Bias"        ~ 0,
      metric == "Coverage(%)" ~ 95,
      TRUE ~ NA_real_
    ),
    metric  = factor(metric, levels = c("Bias","Empirical SE","Model SE","Relative Error in Model SE (%)","Coverage(%)")),
    c_factor = factor(c_param, levels = c_params) 
  )

# Make sure c_param is treated as factor with correct order
metrics_long$c_param <- factor(
  metrics_long$c_param,
  levels = sort(unique(metrics_long$c_param))
)

# Labels only for c=0 and c=6
label_data <- metrics_long %>% 
  filter(c_param %in% c("0", "6")) %>%
  mutate(
    metric_chr = as.character(metric),
    label_val = case_when(
      grepl("%", metric_chr) ~ sprintf("%.1f%%", est),   # 1 decimal + %
      metric_chr == "Bias"   ~ sprintf("%.4f", est),     # 4 decimals for Bias
      TRUE                   ~ sprintf("%.3f", est)      # otherwise 3 decimals
    ),
    hjust_val = ifelse(c_param == "0", -0.1, 1.1),  
    vjust_val = -0.5
  )


p <- ggplot(metrics_long, aes(x = c_param, y = est)) +
  geom_hline(aes(yintercept = ref), linetype = "dashed", na.rm = TRUE, color = "grey50") +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0, linewidth = 0.5, color = "black") +
  geom_point(size = 2, color = "black") +
  geom_text(
    data = label_data,
    aes(label = label_val, hjust = hjust_val, vjust = vjust_val),
    size = 3, color = "black"
  ) +
  ggh4x::facet_wrap2(~ metric, scales = "free_y", ncol = 2, axes = "all") +
  # facet-specific y-limit for Model SE
  ggh4x::facetted_pos_scales(
    y = list(
      metric == "Model SE" ~ scale_y_continuous(limits = c(0.02, 0.072)),
      metric == "Bias" ~ scale_y_continuous(limits = c(-0.02, 0.02))
    )
  ) +
  labs(
    x = "Heteroscedasticity parameter (c)",
    y = "Estimate (point ± 1.96 × MC-SE)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.border = element_rect(color = "lightgray", fill = NA, linewidth = 0.5),
    panel.background = element_rect(fill = "white", color = NA),     # inside each panel
    plot.background  = element_rect(fill = "gray95", color = NA),
    strip.text = element_text(face = "bold", size = 10),            # bold facet titles
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

p


#-----------------------------------------------
# create residuals vs fitted plot
# 3) Residual vs y_hat subplots (one panel per c_param)
my_labels <- c(
  "0" = "c = 0 (Homoscedastic)",
  "0.4" = "c = 0.4 (Mild Heteroscedasticity)",
  "1" = "c = 1 (Moderate Heteroscedasticity)",
  "6" = "c = 6 (Strong Heteroscedasticity)"
)

ggplot(simulated_data, aes(y_hat, resid)) +
  geom_point(alpha = 0.3, size = 0.4) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_smooth(method = "lm", se = FALSE, formula = y ~ x, color = "#33A1E0") +
  facet_wrap(
    ~ c_param,
    ncol = 2,
    labeller = labeller(c_param = as_labeller(my_labels)),  
    scales = "fixed"
  ) +
  labs(x = expression(hat(y)), y = "Residual") +   
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),            # no grid lines
    strip.background = element_blank(),      # no gray background for facet labels
    strip.text = element_text(face = "bold"),
    panel.border = element_rect(color = "lightgray", fill = NA, linewidth = 0.5) # light gray box
  )

# Calculate sd of residuals for low and high fitted values
simulated_data %>%
  group_by(c_param) %>%
  summarise(sd_low  = sd(resid[y_hat < median(y_hat)]),
            sd_high = sd(resid[y_hat > median(y_hat)]),
            ratio   = sd_high / sd_low)


## ================================
## 3. Save results table and figure
## ================================

table_save <- paste0(save_results_path,"/01_optimized-seed-best-centered-log-linear-table-simulation-results-mcse_",param_data,'.png')
save_kable(tbl, table_save , zoom = 10, vwidth = 327, vheight = 206)
plot_save <- paste0(save_results_path,"/02_optimized-seed-best-centered-log-linear-plots-simulation-results-mcse_",param_data,'.png')
ggsave(plot_save, plot = p, width = 6, height = 6, dpi = 1000)
message("Saved results table to: ", table_save,
        "\nSaved results plot to: ", plot_save)
average_error_variance
