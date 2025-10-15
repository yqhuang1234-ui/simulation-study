## ================================
## 0. Setup
## ================================
seed <- 1234
sample_size <- 200
#---regression model parameters---
beta0 <- 14
beta1 <- 0.4
#---distribution of x--- 
# x is uniform distribution
x_min <- 1
x_max <- 6 
#---distribution of error term---
# error term is normally distributed N(0,1) for homoscedastic baseline and depend on x for heteroscedastic scenarios.
# variance of error term is detailed in the below function section.
eps_mean <- 0
# when c=0, variance is 1 and it is homoscedastic baseline.
c_params <- c(0, 0.5, 1, 2, 4, 8)
#used as input for the function below
params <- list(
  n = sample_size,
  x_min = x_min,
  x_max = x_max,
  eps_mean = eps_mean,
  beta0 = beta0,
  beta1 = beta1
)


## ================================
## 1. Function
## ================================
get_unif_moments <- function(x_min, x_max){
  mean_x <- (x_max + x_min) / 2
  var_x <- (x_max - x_min)^2 / 12
  return(list(mean=mean_x, var=var_x))
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
get_eps_var <- function(x, c, mean_x, var_x) {
  eps_var <- (1 + c * (x - mean_x)^2) / (1 + var_x * c)
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
#' @param verbose Logical (default = TRUE). If TRUE, the function prints a
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
get_data <- function(c, x, params, verbose=TRUE) {
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
  x_moments <- get_unif_moments(x_min, x_max)
  mean_x <- x_moments$mean
  var_x <- x_moments$var
  eps_var <- get_eps_var(x, c, mean_x, var_x)
  
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
      c = c, n = n, x_min = x_min, x_max = x_max, mean_x=mean_x, var_x=var_x, x_bar=mean(x), sample_var_x=var(x),
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
fit_lm <- function(dat) {
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
  
  return(data.frame(
    beta1_hat = beta1_hat,
    se = se,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    p_value = pval,
    n = n,
    r_squared = r2,
    row.names = NULL
  ))
}

## ================================
## 2. Data Generation
## ================================
set.seed(seed)
# Generate predictor values
x <- runif(sample_size, x_min, x_max)
set.seed(seed)
data_homo <- get_data(c=0, x=x, params = params)
set.seed(seed)
data_hete1 <- get_data(c=8, x=x, params = params)


## ================================
## 3. Simulation Loop
## ================================

## ================================
## 4. Summarize Results in Table
## ================================
