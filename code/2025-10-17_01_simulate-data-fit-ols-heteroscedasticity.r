## ================================
## 0. Setup
## ================================
seed <- 1234
sample_size <- 100
n_reps <- 1000
# true parameters of the linear model for data generation
beta0 <- 14
beta1 <- 0.4
# x is uniform distribution
x_min <- 1
x_max <- 6 
# mean of the error term
eps_mean <- 0
# heteroscedasticity parameters to simulate
# when c=0, variance is 1 and it is homoscedastic baseline.
c_params <- c(0,  2, 5, 8, 10,15)
# store parameters in a list for easy passing to functions
params <- list(
  n = sample_size,
  x_min = x_min,
  x_max = x_max,
  eps_mean = eps_mean,
  beta0 = beta0,
  beta1 = beta1
)
# path to save results
save_path <- "./data/"

## ================================
## 1. Function
## ================================
get_unif_moments <- function(x_min, x_max){
  expected_x2 <- (x_max^3 - x_min^3) / (3 * (x_max - x_min))
  return(expected_x2)
}
#------------------------------------------------------------
#' Compute error variance
#'
#' Given predictor values x, parameter c, and reference mean_x and var_x 
#' (the mean and variance of the predictor distribution), compute the
#' conditional variance of the error term.
#'
#' @param x Numeric vector of predictor values.
#' @param c Numeric scalar, heteroscedasticity parameter.
#' @param mean_x Mean of the predictor distribution (e.g., uniform mean).
#' @param var_x Variance of the predictor distribution (e.g., uniform variance).
#'
#' @return Numeric vector of conditional error variances for each x.
get_eps_var <- function(x, c, expected_x2) {
  g <- x^2
  eps_var <- (1 + c * g) / (1 + c * expected_x2)
  return(eps_var)
}
#-----------------------------------------------------------
#' Simulate linear regression data with controllable heteroscedasticity
#' 
#' This function generates synthetic data for a linear regression model
#' @param c Numeric. Heteroscedasticity parameter controlling the degree of
#'   variance inflation as a function of x.
#' @param x Numeric vector of length n. Predictor values used directly in
#'   the simulation. Supplying the same x across different values of c
#'   allows apples-to-apples comparisons by holding the predictor
#'   distribution constant.
#' @param params A list containing:
#'   - n: sample size
#'   - x_min, x_max: range of x (uniform)
#'   - eps_mean: mean of the error term
#'   - beta0, beta1: regression coefficients
#' @param verbose Logical (default = FALSE). If TRUE, the function prints a
#'   short status message ("Data simulated successfully") along with a
#'   diagnostic table of parameters and summary statistics
#'
#' @return A data.frame with columns:
#'   - x: predictor values
#'   - y: response values
#'   - c_param: heteroscedasticity parameter
#'   - eps: simulated error terms
#'   - eps_var: conditional variance of each error term
#'
#' @details
#' Why normalization is applied:
#'   The variance function is scaled so that the *average error variance across x*
#'   remains equal to 1 for all values of c. Without this normalization, larger c
#'   would both increase heteroscedasticity *and* inflate the overall noise level,
#'   confounding the comparison between homoscedastic and heteroscedastic settings.
#'   By holding the average variance constant, the simulation isolates the impact
#'   of heteroscedasticity on inference rather than on total noise.
#'
#' Why store extra parameters:
#'   Returning summary statistics (mean of eps_var, variance of y, etc.)
#'   provides a quick check that the simulation is working as expected,
#'   and documents the conditions under which the data were generated.
get_data <- function(c, x, params, verbose=FALSE) {
  n <- params$n
  x_min <- params$x_min
  x_max <- params$x_max
  eps_mean <- params$eps_mean
  beta0 <- params$beta0
  beta1 <- params$beta1
  # Pre-compute mean and variance of X ~ Uniform(x_min, x_max)
  mean_x <- (x_min + x_max) / 2
  var_x  <- (x_max - x_min)^2 / 12
  
  # Compute conditional error variance for each x
  expected_x2 <- get_unif_moments(x_min, x_max)
  eps_var <- get_eps_var(x, c, expected_x2)
  
  # Draw error terms
  eps <- rnorm(n, mean = eps_mean, sd = sqrt(eps_var))
  
  # Generate response variable
  y <- beta0 + beta1 * x + eps
  
  # Return dataset
  data <- data.frame(
    x = x,
    y = y,
    c_param = c,   # renamed for clarity
    eps = eps,
    eps_var = eps_var
  )
  params = data.frame(
      c = c, n = n, x_min = x_min, x_max = x_max, x_bar=mean(x), sample_var_x=var(x),
      beta0 = beta0, beta1 = beta1, eps_mean = eps_mean, mean_eps_var = mean(eps_var), var_y=var(y)
    )
  if (verbose){
    message("Data simulated successfully")
    print(params)}
  return(data)
}
#-----------------------------------------------------------
#' Fit OLS and extract inference metrics for the slope coefficient
#'
#' This function fits a simple linear regression model of the form
#'   \eqn{y = \beta_0 + \beta_1 x + \varepsilon} using ordinary least squares.
#' It returns the estimated slope, its standard error, 95% confidence interval,
#' and p-value based on classical OLS assumptions.
#'
#' @param dat A data.frame containing the variables:
#'   - x: numeric predictor
#'   - y: numeric response
#'
#' @return A data.frame with one row containing:
#'   - beta1_hat: estimated slope coefficient for x
#'   - se: standard error of the slope estimate
#'   - ci_lower: lower bound of the 95% confidence interval
#'   - ci_upper: upper bound of the 95% confidence interval
#'   - p_value: p-value testing H0: beta1 = 0
#'   - n: sample size used in the regression
#'   - r_squared: model R-squared
#'
#' @details
#' Why this is useful:
#'   These metrics allow assessment of bias, efficiency, coverage, and test size
#'   across different simulation settings. In particular:
#'   - beta1_hat can be compared to beta1_true to check bias.
#'   - se reflects estimation precision.
#'   - ci_lower and ci_upper allow CI coverage evaluation.
#'   - p_value is used to assess power (when beta1_true ≠ 0).
#'
#' Confidence intervals are computed using R’s built-in \code{confint()} function,
#' which relies on the t distribution with residual degrees of freedom.
fit_lm <- function(dat, verbose=FALSE) {
  stopifnot(all(c("x","y") %in% names(dat)))
  
  fit <- lm(y ~ x, data = dat)
  n <- nobs(fit)
  r2 <- summary(fit)$r.squared
  beta1_hat <- coef(fit)[["x"]]
  
  summ <- summary(fit)$coefficients
  se <- summ["x", "Std. Error"]
  pval <- summ["x", "Pr(>|t|)"]
  
  # 95% CI for slope directly from confint()
  ci <- confint(fit, "x", level = 0.95)
  ci_lower <- ci[1]
  ci_upper <- ci[2]
  
  result <- data.frame(
    beta1_hat = beta1_hat,
    se = se,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    p_value = pval,
    n = n,
    r_squared = r2,
    c_param = unique(dat$c_param)
  )
  if (verbose){
    message("model run successfully")
    print(result)
  }
  return(result)
}
#-----------------------------------------------------------
#' Simulate datasets and model fits across multiple c values
#'
#' This function loops over a vector of heteroscedasticity parameters (`c_values`)
#' and repeated replications (`n_reps`). For each replication, a vector of predictors
#' `x` is drawn once from Uniform(x_min, x_max) and reused across all `c` values. 
#' For each combination of `rep` and `c`, data are simulated via `get_data()` and
#' models are fitted with `fit_lm()`. Results are collected into two stacked data frames.
#'
#' A progress bar is shown in the console, and a one-line status message is updated
#' during the run to indicate current replication, c index, and iteration count.
#'
#' @param c_values Numeric vector. Values of heteroscedasticity parameter `c` to simulate.
#' @param params List. 
#'   parameters required by `get_data()`.
#' @param n_reps Integer. Number of replications to run for each `c` value.
#' @param seed Optional integer. If supplied, sets random seed for reproducibility.
#'
#' @return A list with two elements:
#'   \item{fits}{Data frame of stacked model fit results}
#'   \item{datasets}{Data frame of stacked simulated datasets}
#'
simulate_over_c <- function(c_values, params, n_reps, seed = NULL) {
  x_min <- params$x_min; x_max <- params$x_max; n <- params$n
  if (!is.null(seed)) set.seed(seed)

  fits_all     <- NULL
  datasets_all <- NULL

  total_iters <- n_reps * length(c_values)

  # progress bar + reserve a second line for status text
  pb <- utils::txtProgressBar(min = 0, max = total_iters, style = 3)
  on.exit(close(pb), add = TRUE)
  cat("\n")

  iter <- 0L

  for (r in seq_len(n_reps)) {
    # for each rep, generate x once. x is same across all c but different across reps
    x_rep <- runif(n, min = x_min, max = x_max)
    # for each c, generate data using same x_rep
    for (j in seq_along(c_values)) {
      c_val <- c_values[j]

      dat <- get_data(c = c_val, x = x_rep, params = params)
      fit <- fit_lm(dat)

      # tag
      dat$rep <- r; dat$seed <- seed
      fit$rep <- r; fit$seed <- seed

      # append
      datasets_all <- rbind(datasets_all, dat)
      fits_all     <- rbind(fits_all, fit)

      iter <- iter + 1L

      # throttle progress bar updates to reduce overhead
      if (iter %% 10L == 0L || iter == total_iters) {
        utils::setTxtProgressBar(pb, iter)
      }

      # live status line (overwrites same line)
      cat(sprintf(
        "\rrep %d/%d | c %d/%d (c = %s) | iter %d/%d",
        r, n_reps, j, length(c_values), as.character(c_val), iter, total_iters
      ))
      flush.console()
    }
  }

  # finish with a fresh line
  cat("\n")

  list(fits = fits_all, datasets = datasets_all)
}

#-----------------------------------------------------------
#' Sanity check for simulation results
#'
#' This function verifies that the simulation results produced by `simulate_over_c()`
#' behave as expected. It checks variation of `x`, `eps`, and `y` across replications
#' and consistency of `x` across different `c` values within each replication.
#'
#' Specifically, it performs three groups of checks:
#' \enumerate{
#'   \item \strong{Within each c}: verifies that `x`, `eps`, and `y` vary across
#'         replications (`rep`).
#'   \item \strong{Across c within each rep}: verifies that the vector of `x`
#'         values is identical for all `c` values within the same replication.
#'   \item \strong{Across c overall}: verifies that the variance of `eps` and
#'         the variance of `y` differ across `c` values (expected under heteroscedasticity).
#' }
#'
#' @param res A list returned by `simulate_over_c()`. Must contain an element
#'   `datasets`, which is a data frame with columns `x`, `eps`, `y`, `c`, and `rep`.
#'
#' @return Prints a summary of the checks to the console.
#'
#' @details
#' - If all checks return `TRUE`, it indicates that the simulation logic is functioning
#' - If any check returns `FALSE`, it may indicate a problem in the simulation logic
#'   (e.g., seed placement, incorrect reuse of `x`, or missing heteroscedasticity).
#'
sanity_check <- function(res) {
  dat <- res$datasets
  
  # --- A) varies across reps within each c ---
  vary_by_rep <- function(vec) {
    # different group means across reps ⇒ variation
    all(tapply(vec, list(dat$c, dat$rep), mean) |>
          apply(1, function(v) length(unique(v)) > 1))
  }
  x_varies_within_c   <- vary_by_rep(dat$x)
  eps_varies_within_c <- vary_by_rep(dat$eps)
  y_varies_within_c   <- vary_by_rep(dat$y)
  
  # --- B) for each rep, x is identical across all c (vector-wise) ---
  x_identical_across_c_per_rep <- {
    by_rep <- split(dat, dat$rep)
    all(vapply(by_rep, function(df_rep) {
      xs <- split(df_rep$x, df_rep$c)
      ref <- xs[[1]]
      # compare element-wise identity to the first c
      all(vapply(xs[-1], function(v) identical(as.numeric(v), as.numeric(ref)), logical(1)))
    }, logical(1)))
  }
  
  
  # --- C) across c: eps/y variance differ (expected if heteroscedastic by c) ---
  eps_var_diff_across_c <- length(unique(tapply(dat$eps, dat$c, var))) > 1
  y_var_diff_across_c   <- length(unique(tapply(dat$y,   dat$c, var))) > 1
  
  cat("Sanity check:\n")
  cat("  x varies across reps within each c:            ", x_varies_within_c,   "\n")
  cat("  eps varies across reps within each c:          ", eps_varies_within_c, "\n")
  cat("  y varies across reps within each c:            ", y_varies_within_c,   "\n")
  cat("  x identical across c for each rep (elementwise): ", x_identical_across_c_per_rep, "\n")
  cat("  eps variance differs across c:                 ", eps_var_diff_across_c, "\n")
  cat("  y variance differs across c:                   ", y_var_diff_across_c,   "\n")
}



## ================================
## 2. Data Generation & Simulation Loop
## ================================
res <- simulate_over_c(c_values = c_params, params = params, n_reps = n_reps, seed=seed)

## ================================
## 3. Sanity Check  
## ================================
sanity_check(res)

## ================================
## 4. Save Results
## ================================
if (!dir.exists(save_path)) {
  dir.create(save_path, recursive = TRUE)
}

# ISO 8601 date (date only):
iso_date <- format(Sys.Date(), "%Y-%m-%d")

# describe the run
experiment_tag <- "ols-heteroscedasticity-simulation"

# pull key params
n        <- unique(res$fits$n)
n_reps   <- max(res$fits$rep)
seed_val_input <- unique(res$fits$seed)
seed_val <- if (is.null(seed_val_input)) "na" else as.character(seed_val_input)
# get unique c values from fits
c_vals <- sort(unique(res$fits$c_param))  
# collapse into compact string

# compact c-values for filenames: 0, 0.5, 1, 2 -> "c0-0.5-1-2"
c_slug <- paste0("c", paste(c_vals, collapse = "-"))

# build a consistent suffix with key params
param_suffix <- sprintf("n%d_reps%d_seed%s_%s", n, n_reps, seed_val, c_slug)

# final filenames
fits_rds      <- file.path(save_path, sprintf("%s_level-dependent-fits_%s_%s.rds",      iso_date, experiment_tag, param_suffix))
datasets_rds  <- file.path(save_path, sprintf("%s_level-dependent-datasets_%s_%s.rds",  iso_date, experiment_tag, param_suffix))

# save (RDS for fidelity + CSV for sharing)
saveRDS(res$fits,     fits_rds)
saveRDS(res$datasets, datasets_rds)
# confirmation message
message("Results saved successfully:\n",
        "  - ", fits_rds, "\n",
        "  - ", datasets_rds, "\n")
