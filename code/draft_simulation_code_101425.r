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
params <- list(
  n = sample_size,
  x_min = x_min,
  x_max = x_max,
  eps_mean = eps_mean,
  beta0 = beta0,
  beta1 = beta1
)

# Access them like a dictionary
params$n
params$beta1


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
#' Generate simulated dataset with errors
#'
#' Draws predictors from a uniform distribution, computes error 
#' variance, simulates error terms, and generates responses from a linear model.
#'
#' @param c Numeric scalar, heteroscedasticity parameter.
#' @param n Sample size (number of observations).
#' @param x_min Lower bound of predictor distribution (uniform).
#' @param x_max Upper bound of predictor distribution (uniform).
#' @param eps_mean Mean of the error term.
#' @param beta0 Intercept parameter.
#' @param beta1 Slope parameter.
#'
#' @return Data frame with columns:
#'   - x: predictor values
#'   - y: response values
#'   - c_param: heteroscedasticity parameter
#'   - eps: error terms
#'   - eps_var: conditional error variance for each x
get_data <- function(c, paramsn, x_min, x_max, eps_mean, beta0, beta1) {
  n <- params$n
  x_min <- params$x_min
  x_max <- params$x_max
  eps_mean <- params$eps_mean
  beta0 <- 
  # Pre-compute mean and variance of X ~ Uniform(x_min, x_max)
  mean_x <- (x_min + x_max) / 2
  var_x  <- (x_max - x_min)^2 / 12
  
  # Generate predictor values
  x <- runif(n, x_min, x_max)
  
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
      c = c, n = n, x_min = x_min, x_max = x_max, mean_x=mean_x, var_x=var_x,
      eps_mean = eps_mean, mean_eps_var = mean(eps_var), var_y=var(y), beta0 = beta0, beta1 = beta1
    )
  print("data is simulated successfully")
  print(params)
  return(data)
}
#-----------------------------------------------------------

## ================================
## 2. Data Generation
## ================================
set.seed(1234)
data_homo <- get_data(c=0, n=sample_size, x_min=x_min, x_max=x_max, eps_mean=eps_mean, beta0=beta0, beta1=beta1)
set.seed(1234)
data_hete1 <- get_data(c=1, n=sample_size, x_min=x_min, x_max=x_max, eps_mean=eps_mean, beta0=beta0, beta1=beta1)
## ================================
## 3. Simulation Loop
## ================================

## ================================
## 4. Summarize Results in Table
## ================================
