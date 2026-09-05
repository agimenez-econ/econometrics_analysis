####################################################
# FINITE-SAMPLE PROPERTIES OF FGLS
# Monte Carlo Simulation
# Seed: 2026
####################################################


# --------------------------------------------------
# 1. Setup
# --------------------------------------------------

rm(list = ls())

# Reproducibility
set.seed(2026)


# --------------------------------------------------
# 2. Model parameters
# --------------------------------------------------

# True parameters
beta0 <- -3
beta1 <- 0.8

# Monte Carlo replications
R <- 5000

# Sample sizes
sample_sizes <- c(5, 10, 30, 100, 200, 500)

# Significance levels
crit_5 <- qnorm(0.975)
crit_1 <- qnorm(0.995)

# Variance structure
omega_values <- c(4, 9, 16, 25, 36)


# --------------------------------------------------
# 3. FGLS + Cholesky simulation function
# --------------------------------------------------

simulate_fgls <- function(n, beta1_true, R = 5000) {
  
  # Number of observations per group
  N <- n / 5
  
  # Variance structure:
  # Omega ⊗ I_N
  var_vec <- rep(omega_values, each = N)
  
  # Inverse of Omega
  Omega_inv <- diag(1 / var_vec)
  
  # Cholesky factor
  # Omega = P P'
  P <- diag(sqrt(var_vec))
  
  # Storage for FGLS
  beta0_fgls <- numeric(R)
  beta1_fgls <- numeric(R)
  t_fgls <- numeric(R)
  
  # Storage for Cholesky/MCC
  beta0_chol <- numeric(R)
  beta1_chol <- numeric(R)
  
  
  # ------------------------------------------------
  # Monte Carlo loop
  # ------------------------------------------------
  
  for (r in 1:R) {
    
    # ----------------------------------------------
    # Generate sample
    # ----------------------------------------------
    
    # Regressor
    x <- runif(n, 1, 50)
    
    # Design matrix
    X <- cbind(1, x)
    
    # Heteroskedastic errors
    u <- rnorm(
      n,
      mean = 0,
      sd = sqrt(var_vec)
    )
    
    # Dependent variable
    y <- beta0 + beta1_true * x + u
    
    
    # ----------------------------------------------
    # FGLS estimator
    # ----------------------------------------------
    
    XtOmegaInvX <- t(X) %*% Omega_inv %*% X
    
    beta_fgls <- solve(
      XtOmegaInvX
    ) %*%
      t(X) %*% Omega_inv %*% y
    
    # Variance of FGLS estimator
    V_fgls <- solve(XtOmegaInvX)
    
    # Standard error of beta1
    se_fgls <- sqrt(diag(V_fgls))
    
    # Test statistic:
    # H0: beta1 = 0.8
    t_fgls[r] <- (
      beta_fgls[2] - beta1
    ) / se_fgls[2]
    
    # Store FGLS estimates
    beta0_fgls[r] <- beta_fgls[1]
    beta1_fgls[r] <- beta_fgls[2]
    
    
    # ----------------------------------------------
    # Cholesky transformation
    # ----------------------------------------------
    
    # Since Omega = P P',
    # transform the model as:
    #
    # y* = P^(-1) y
    # X* = P^(-1) X
    
    y_star <- y / sqrt(var_vec)
    X_star <- X / sqrt(var_vec)
    
    
    # ----------------------------------------------
    # MCC on transformed model
    # ----------------------------------------------
    
    beta_chol <- solve(
      t(X_star) %*% X_star
    ) %*%
      t(X_star) %*%
      y_star
    
    # Store Cholesky/MCC estimates
    beta0_chol[r] <- beta_chol[1]
    beta1_chol[r] <- beta_chol[2]
  }
  
  
  # ------------------------------------------------
  # Return results
  # ------------------------------------------------
  
  return(
    list(
      
      # FGLS
      beta0_fgls = beta0_fgls,
      beta1_fgls = beta1_fgls,
      t_fgls = t_fgls,
      
      # Cholesky/MCC
      beta0_chol = beta0_chol,
      beta1_chol = beta1_chol
    )
  )
}


# --------------------------------------------------
# 4. Run Monte Carlo simulations
# --------------------------------------------------

fgls_results <- list()


for (n in sample_sizes) {
  
  cat("\n====================================\n")
  cat("Simulating FGLS, n =", n, "\n")
  cat("====================================\n")
  
  
  # ----------------------------------------------
  # H0: beta1 = 0.8
  # ----------------------------------------------
  
  cat("H0: beta1 = 0.8\n")
  
  res_H0 <- simulate_fgls(
    n = n,
    beta1_true = 0.8,
    R = R
  )
  
  
  # ----------------------------------------------
  # Alternative: beta1 = 0
  # ----------------------------------------------
  
  cat("Alternative: beta1 = 0\n")
  
  res_0 <- simulate_fgls(
    n = n,
    beta1_true = 0,
    R = R
  )
  
  
  # ----------------------------------------------
  # Alternative: beta1 = 0.4
  # ----------------------------------------------
  
  cat("Alternative: beta1 = 0.4\n")
  
  res_04 <- simulate_fgls(
    n = n,
    beta1_true = 0.4,
    R = R
  )
  
  
  # ----------------------------------------------
  # Test size
  # ----------------------------------------------
  
  size_5 <- mean(
    abs(res_H0$t_fgls) > crit_5
  )
  
  size_1 <- mean(
    abs(res_H0$t_fgls) > crit_1
  )
  
  
  # ----------------------------------------------
  # Power: beta1 = 0
  # ----------------------------------------------
  
  power0_5 <- mean(
    abs(res_0$t_fgls) > crit_5
  )
  
  power0_1 <- mean(
    abs(res_0$t_fgls) > crit_1
  )
  
  
  # ----------------------------------------------
  # Power: beta1 = 0.4
  # ----------------------------------------------
  
  power04_5 <- mean(
    abs(res_04$t_fgls) > crit_5
  )
  
  power04_1 <- mean(
    abs(res_04$t_fgls) > crit_1
  )
  
  
  # ----------------------------------------------
  # FGLS descriptive statistics
  # ----------------------------------------------
  
  beta0_mean <- mean(res_H0$beta0_fgls)
  beta0_median <- median(res_H0$beta0_fgls)
  beta0_sd <- sd(res_H0$beta0_fgls)
  
  beta1_mean <- mean(res_H0$beta1_fgls)
  beta1_median <- median(res_H0$beta1_fgls)
  beta1_sd <- sd(res_H0$beta1_fgls)
  
  
  # ----------------------------------------------
  # Cholesky/MCC descriptive statistics
  # ----------------------------------------------
  
  beta0_chol_mean <- mean(res_H0$beta0_chol)
  beta0_chol_median <- median(res_H0$beta0_chol)
  beta0_chol_sd <- sd(res_H0$beta0_chol)
  
  beta1_chol_mean <- mean(res_H0$beta1_chol)
  beta1_chol_median <- median(res_H0$beta1_chol)
  beta1_chol_sd <- sd(res_H0$beta1_chol)
  
  
  # ----------------------------------------------
  # FGLS vs. Cholesky
  # ----------------------------------------------
  
  diff_beta0 <- (
    res_H0$beta0_fgls -
      res_H0$beta0_chol
  )
  
  diff_beta1 <- (
    res_H0$beta1_fgls -
      res_H0$beta1_chol
  )
  
  
  # Maximum absolute difference
  max_diff_beta0 <- max(
    abs(diff_beta0)
  )
  
  max_diff_beta1 <- max(
    abs(diff_beta1)
  )
  
  
  # Mean absolute difference
  mean_abs_diff_beta0 <- mean(
    abs(diff_beta0)
  )
  
  mean_abs_diff_beta1 <- mean(
    abs(diff_beta1)
  )
  
  
  # ----------------------------------------------
  # Store results
  # ----------------------------------------------
  
  fgls_results[[paste0("n", n)]] <- list(
    
    n = n,
    
    # Test
    size_1 = size_1,
    size_5 = size_5,
    
    power0_1 = power0_1,
    power0_5 = power0_5,
    
    power04_1 = power04_1,
    power04_5 = power04_5,
    
    # FGLS
    beta0_mean = beta0_mean,
    beta0_median = beta0_median,
    beta0_sd = beta0_sd,
    
    beta1_mean = beta1_mean,
    beta1_median = beta1_median,
    beta1_sd = beta1_sd,
    
    # Cholesky/MCC
    beta0_chol_mean = beta0_chol_mean,
    beta0_chol_median = beta0_chol_median,
    beta0_chol_sd = beta0_chol_sd,
    
    beta1_chol_mean = beta1_chol_mean,
    beta1_chol_median = beta1_chol_median,
    beta1_chol_sd = beta1_chol_sd,
    
    # Differences
    max_diff_beta0 = max_diff_beta0,
    max_diff_beta1 = max_diff_beta1,
    
    mean_abs_diff_beta0 = mean_abs_diff_beta0,
    mean_abs_diff_beta1 = mean_abs_diff_beta1,
    
    # Simulation vectors
    beta0_simulation = res_H0$beta0_fgls,
    beta1_simulation = res_H0$beta1_fgls
  )
}


# --------------------------------------------------
# 5. Create final FGLS results table
# --------------------------------------------------

fgls_table <- do.call(
  rbind,
  lapply(
    fgls_results,
    function(x) {
      
      data.frame(
        
        n = x$n,
        
        size_1 = x$size_1,
        size_5 = x$size_5,
        
        power0_1 = x$power0_1,
        power0_5 = x$power0_5,
        
        power04_1 = x$power04_1,
        power04_5 = x$power04_5,
        
        beta0_mean = x$beta0_mean,
        beta0_median = x$beta0_median,
        beta0_sd = x$beta0_sd,
        
        beta1_mean = x$beta1_mean,
        beta1_median = x$beta1_median,
        beta1_sd = x$beta1_sd
      )
    }
  )
)

rownames(fgls_table) <- NULL


# --------------------------------------------------
# 6. Create FGLS vs. Cholesky table
# --------------------------------------------------

fgls_cholesky_table <- do.call(
  rbind,
  lapply(
    fgls_results,
    function(x) {
      
      data.frame(
        
        n = x$n,
        
        beta0_mean_FGLS = x$beta0_mean,
        beta0_mean_MCC = x$beta0_chol_mean,
        
        beta0_median_FGLS = x$beta0_median,
        beta0_median_MCC = x$beta0_chol_median,
        
        beta0_sd_FGLS = x$beta0_sd,
        beta0_sd_MCC = x$beta0_chol_sd,
        
        beta1_mean_FGLS = x$beta1_mean,
        beta1_mean_MCC = x$beta1_chol_mean,
        
        beta1_median_FGLS = x$beta1_median,
        beta1_median_MCC = x$beta1_chol_median,
        
        beta1_sd_FGLS = x$beta1_sd,
        beta1_sd_MCC = x$beta1_chol_sd,
        
        mean_abs_diff_beta0 = x$mean_abs_diff_beta0,
        mean_abs_diff_beta1 = x$mean_abs_diff_beta1,
        
        max_diff_beta0 = x$max_diff_beta0,
        max_diff_beta1 = x$max_diff_beta1
      )
    }
  )
)

rownames(fgls_cholesky_table) <- NULL


# --------------------------------------------------
# 7. Display final results
# --------------------------------------------------

cat("\n\n====================================\n")
cat("FINAL FGLS RESULTS\n")
cat("====================================\n\n")

print(fgls_table)


cat("\n\n====================================\n")
cat("FGLS VS. CHOLESKY / MCC RESULTS\n")
cat("====================================\n\n")

print(fgls_cholesky_table)

# --------------------------------------------------
# 8. Visualizations
# --------------------------------------------------

library(ggplot2)


# --------------------------------------------------
# 8.1. Empirical test size
# --------------------------------------------------

size_plot_data <- data.frame(
  n = factor(
    rep(fgls_table$n, 2),
    levels = sample_sizes
  ),
  
  size = c(
    fgls_table$size_1,
    fgls_table$size_5
  ),
  
  significance = rep(
    c("1%", "5%"),
    each = nrow(fgls_table)
  )
)

p_size <- ggplot(
  size_plot_data,
  aes(
    x = n,
    y = size,
    group = significance,
    linetype = significance
  )
) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.5) +
  geom_hline(
    yintercept = 0.01,
    linetype = "dashed"
  ) +
  geom_hline(
    yintercept = 0.05,
    linetype = "dashed"
  ) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1)
  ) +
  labs(
    title = "Empirical Test Size",
    subtitle = "Nominal significance levels: 1% and 5%",
    x = "Sample size (5N)",
    y = "Empirical rejection rate",
    linetype = "Nominal level"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

print(p_size)


# --------------------------------------------------
# 8.2. Test power
# --------------------------------------------------

power_plot_data <- data.frame(
  
  n = factor(
    rep(fgls_table$n, 4),
    levels = sample_sizes
  ),
  
  power = c(
    fgls_table$power0_1,
    fgls_table$power0_5,
    fgls_table$power04_1,
    fgls_table$power04_5
  ),
  
  alternative = c(
    rep("beta1 = 0", nrow(fgls_table)),
    rep("beta1 = 0", nrow(fgls_table)),
    rep("beta1 = 0.4", nrow(fgls_table)),
    rep("beta1 = 0.4", nrow(fgls_table))
  ),
  
  significance = rep(
    c("1%", "5%", "1%", "5%"),
    each = nrow(fgls_table)
  )
)

p_power <- ggplot(
  power_plot_data,
  aes(
    x = n,
    y = power,
    group = interaction(alternative, significance),
    linetype = significance
  )
) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.5) +
  scale_y_continuous(
    limits = c(0, 1),
    labels = scales::percent_format(accuracy = 1)
  ) +
  facet_wrap(
    ~ alternative
  ) +
  labs(
    title = "Test Power",
    subtitle = "Power under beta1 = 0 and beta1 = 0.4",
    x = "Sample size (5N)",
    y = "Empirical power",
    linetype = "Significance level"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

print(p_power)


# --------------------------------------------------
# 8.3. Sampling variability of beta1
# --------------------------------------------------

p_sd <- ggplot(
  fgls_table,
  aes(
    x = factor(n, levels = sample_sizes),
    y = beta1_sd
  )
) +
  geom_line(
    aes(group = 1),
    linewidth = 0.9
  ) +
  geom_point(size = 2.5) +
  labs(
    title = "Sampling Variability of beta1",
    subtitle = "Standard deviation of the FGLS estimator",
    x = "Sample size (5N)",
    y = "SD of beta1_hat"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold")
  )

print(p_sd)


# --------------------------------------------------
# 8.4. FGLS vs. Cholesky / MCC
# --------------------------------------------------

comparison_plot_data <- data.frame(
  
  n = factor(
    rep(fgls_cholesky_table$n, 4),
    levels = sample_sizes
  ),
  
  estimate = c(
    fgls_cholesky_table$beta1_mean_FGLS,
    fgls_cholesky_table$beta1_mean_MCC,
    fgls_cholesky_table$beta1_median_FGLS,
    fgls_cholesky_table$beta1_median_MCC
  ),
  
  estimator = c(
    rep(
      "FGLS",
      nrow(fgls_cholesky_table)
    ),
    rep(
      "Cholesky / MCC",
      nrow(fgls_cholesky_table)
    ),
    rep(
      "FGLS",
      nrow(fgls_cholesky_table)
    ),
    rep(
      "Cholesky / MCC",
      nrow(fgls_cholesky_table)
    )
  ),
  
  statistic = c(
    rep(
      "Mean",
      nrow(fgls_cholesky_table) * 2
    ),
    rep(
      "Median",
      nrow(fgls_cholesky_table) * 2
    )
  )
)

p_comparison <- ggplot(
  comparison_plot_data,
  aes(
    x = n,
    y = estimate,
    group = estimator,
    linetype = estimator
  )
) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.5) +
  geom_hline(
    yintercept = beta1,
    linetype = "dashed"
  ) +
  facet_wrap(
    ~ statistic
  ) +
  labs(
    title = "FGLS vs. Cholesky / MCC",
    subtitle = "Estimates of beta1 across sample sizes",
    x = "Sample size (5N)",
    y = "Estimate of beta1",
    linetype = "Estimator"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

print(p_comparison)

# --------------------------------------------------
# 8.5. Save visualizations
# --------------------------------------------------

# Create figures folder if it does not exist
if (!dir.exists("figures")) {
  dir.create("figures")
}

# Save empirical test size
ggsave(
  filename = "figures/fgls_test_size.png",
  plot = p_size,
  width = 8,
  height = 5,
  dpi = 300,
  bg = "white"
)

# Save test power
ggsave(
  filename = "figures/fgls_power.png",
  plot = p_power,
  width = 8,
  height = 5,
  dpi = 300,
  bg = "white"
)

# Save sampling variability
ggsave(
  filename = "figures/fgls_sampling_variability.png",
  plot = p_sd,
  width = 8,
  height = 5,
  dpi = 300,
  bg = "white"
)

# Save FGLS vs. Cholesky / MCC
ggsave(
  filename = "figures/fgls_vs_cholesky.png",
  plot = p_comparison,
  width = 8,
  height = 5,
  dpi = 300,
  bg = "white"
)

# --------------------------------------------------
# 9. Save numerical results
# --------------------------------------------------

# Create outputs folder if it does not exist
if (!dir.exists("outputs")) {
  dir.create("outputs")
}


# --------------------------------------------------
# 9.1. Save FGLS results
# --------------------------------------------------

write.csv(
  fgls_table,
  file = "outputs/fgls_results.csv",
  row.names = FALSE
)


# --------------------------------------------------
# 9.2. Save FGLS vs. Cholesky / MCC results
# --------------------------------------------------

write.csv(
  fgls_cholesky_table,
  file = "outputs/fgls_cholesky_comparison.csv",
  row.names = FALSE
)


# --------------------------------------------------
# 9.3. Confirmation
# --------------------------------------------------

cat("\n====================================\n")
cat("RESULTS SAVED SUCCESSFULLY\n")
cat("====================================\n")

cat("\nFiles created:\n")
cat("- outputs/fgls_results.csv\n")
cat("- outputs/fgls_cholesky_comparison.csv\n")
cat("- figures/fgls_test_size.png\n")
cat("- figures/fgls_power.png\n")
cat("- figures/fgls_sampling_variability.png\n")
cat("- figures/fgls_vs_cholesky.png\n")

####################################################
# FINITE-SAMPLE PROPERTIES OF WHITE'S TEST
# AND HETEROSKEDASTICITY-ROBUST VARIANCES
# Monte Carlo Simulation
# Seed: 2026
####################################################

# --------------------------------------------------
# 1. Setup
# --------------------------------------------------

rm(list = ls())

set.seed(2026)

library(sandwich)


# --------------------------------------------------
# 2. Model parameters
# --------------------------------------------------

# True parameters
beta0 <- 1
beta1 <- 1
beta2 <- 1

# Monte Carlo replications
R <- 5000

# Sample sizes
sample_sizes <- c(20, 60, 100, 200, 400, 600)


# --------------------------------------------------
# 3. Base regressors
# --------------------------------------------------

# Base sample: n = 20

x1_20 <- seq(
  -1.1,
  1.1,
  length.out = 20
)

x2_20 <- qnorm(
  runif(20)
)


# --------------------------------------------------
# 4. White test function
# --------------------------------------------------

white_test <- function(y, x1, x2) {
  
  # Original regression
  model <- lm(
    y ~ x1 + x2
  )
  
  # OLS residuals
  u_hat <- resid(model)
  
  # Auxiliary regression
  auxiliary_model <- lm(
    u_hat^2 ~
      x1 +
      x2 +
      I(x1^2) +
      I(x2^2) +
      I(x1 * x2)
  )
  
  # White LM statistic
  LM <- length(y) * summary(
    auxiliary_model
  )$r.squared
  
  # p-value
  p_value <- 1 - pchisq(
    LM,
    df = 5
  )
  
  return(p_value)
}

# --------------------------------------------------
# 5. Monte Carlo simulation for White's test
# --------------------------------------------------

simulate_white <- function(n, design, R = 5000) {
  
  # Number of replications of the base sample
  k <- n / 20
  
  # Construct regressors
  x1 <- rep(x1_20, k)
  x2 <- rep(x2_20, k)
  
  # Storage for rejection decisions
  rejection <- matrix(
    0,
    nrow = R,
    ncol = 3
  )
  
  # Monte Carlo loop
  for (r in 1:R) {
    
    # ----------------------------------------------
    # Design 0: homoskedasticity and normality
    # ----------------------------------------------
    
    if (design == 0) {
      
      nu <- rep(1, n)
      
      u <- rnorm(n)
    }
    
    
    # ----------------------------------------------
    # Design 1: heteroskedasticity and normality
    # ----------------------------------------------
    
    if (design == 1) {
      
      nu <- exp(
        0.25 * x1 +
          0.25 * x2
      )
      
      u <- rnorm(n)
    }
    
    
    # ----------------------------------------------
    # Design 2: heteroskedasticity and non-normality
    # ----------------------------------------------
    
    if (design == 2) {
      
      nu <- exp(
        0.25 * x1 +
          0.25 * x2
      )
      
      u <- rt(
        n,
        df = 5
      )
    }
    
    
    # ----------------------------------------------
    # Generate dependent variable
    # ----------------------------------------------
    
    y <- beta0 +
      beta1 * x1 +
      beta2 * x2 +
      sqrt(nu) * u
    
    
    # ----------------------------------------------
    # White test
    # ----------------------------------------------
    
    p_value <- white_test(
      y,
      x1,
      x2
    )
    
    
    # ----------------------------------------------
    # Rejection decisions
    # ----------------------------------------------
    
    rejection[r, ] <- c(
      p_value < 0.01,
      p_value < 0.05,
      p_value < 0.10
    )
  }
  
  # Empirical rejection rates
  colMeans(rejection)
}

# --------------------------------------------------
# 6. Run Monte Carlo simulations
# --------------------------------------------------

white_results <- list()

for (n in sample_sizes) {
  
  cat("\n====================================\n")
  cat("Simulating White's test, n =", n, "\n")
  cat("====================================\n")
  
  
  # Design 0: test size
  cat("Design 0: homoskedasticity\n")
  
  result_0 <- simulate_white(
    n = n,
    design = 0,
    R = R
  )
  
  
  # Design 1: power
  cat("Design 1: heteroskedasticity + normality\n")
  
  result_1 <- simulate_white(
    n = n,
    design = 1,
    R = R
  )
  
  
  # Design 2: power
  cat("Design 2: heteroskedasticity + non-normality\n")
  
  result_2 <- simulate_white(
    n = n,
    design = 2,
    R = R
  )
  
  
  # Store results
  white_results[[paste0("n", n)]] <- list(
    
    n = n,
    
    design_0 = result_0,
    design_1 = result_1,
    design_2 = result_2
  )
}

# --------------------------------------------------
# 7. Final results table
# --------------------------------------------------

white_table <- data.frame(
  
  n = sample_sizes,
  
  # Design 0: test size
  size_1 = sapply(
    white_results,
    function(x) x$design_0[1]
  ),
  
  size_5 = sapply(
    white_results,
    function(x) x$design_0[2]
  ),
  
  size_10 = sapply(
    white_results,
    function(x) x$design_0[3]
  ),
  
  # Design 1: power
  power_design1_1 = sapply(
    white_results,
    function(x) x$design_1[1]
  ),
  
  power_design1_5 = sapply(
    white_results,
    function(x) x$design_1[2]
  ),
  
  power_design1_10 = sapply(
    white_results,
    function(x) x$design_1[3]
  ),
  
  # Design 2: power
  power_design2_1 = sapply(
    white_results,
    function(x) x$design_2[1]
  ),
  
  power_design2_5 = sapply(
    white_results,
    function(x) x$design_2[2]
  ),
  
  power_design2_10 = sapply(
    white_results,
    function(x) x$design_2[3]
  )
)

print(white_table)

# --------------------------------------------------
# 8. Save results
# --------------------------------------------------

if (!dir.exists("outputs")) {
  dir.create("outputs")
}

write.csv(
  white_table,
  "outputs/white_test_results.csv",
  row.names = FALSE
)

# --------------------------------------------------
# 9. Visualizations: White's test
# --------------------------------------------------

library(ggplot2)


# --------------------------------------------------
# 9.1. Empirical test size
# --------------------------------------------------

size_white_data <- data.frame(
  
  n = factor(
    rep(white_table$n, 3),
    levels = sample_sizes
  ),
  
  rejection_rate = c(
    white_table$size_1,
    white_table$size_5,
    white_table$size_10
  ),
  
  significance = rep(
    c("1%", "5%", "10%"),
    each = nrow(white_table)
  )
)


p_white_size <- ggplot(
  size_white_data,
  aes(
    x = n,
    y = rejection_rate,
    group = significance,
    linetype = significance
  )
) +
  
  geom_line(linewidth = 0.9) +
  
  geom_point(size = 2.5) +
  
  # Nominal significance levels
  geom_hline(
    yintercept = 0.01,
    linetype = "dashed",
    linewidth = 0.6
  ) +
  
  geom_hline(
    yintercept = 0.05,
    linetype = "dashed",
    linewidth = 0.6
  ) +
  
  geom_hline(
    yintercept = 0.10,
    linetype = "dashed",
    linewidth = 0.6
  ) +
  
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, 0.12)
  ) +
  
  labs(
    title = "Empirical Size of White's Test",
    subtitle = "Design 0: homoskedastic errors",
    x = "Sample size (n)",
    y = "Empirical rejection rate",
    linetype = "Nominal level"
  ) +
  
  theme_minimal() +
  
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )


print(p_white_size)


# --------------------------------------------------
# 9.2. Empirical test power
# --------------------------------------------------

power_white_data <- data.frame(
  
  n = factor(
    rep(white_table$n, 6),
    levels = sample_sizes
  ),
  
  power = c(
    white_table$power_design1_1,
    white_table$power_design1_5,
    white_table$power_design1_10,
    white_table$power_design2_1,
    white_table$power_design2_5,
    white_table$power_design2_10
  ),
  
  design = c(
    rep(
      "Design 1: Heteroskedasticity + Normality",
      nrow(white_table) * 3
    ),
    rep(
      "Design 2: Heteroskedasticity + Non-normality",
      nrow(white_table) * 3
    )
  ),
  
  significance = rep(
    c("1%", "5%", "10%"),
    times = 2,
    each = nrow(white_table)
  )
)


p_white_power <- ggplot(
  power_white_data,
  aes(
    x = n,
    y = power,
    group = significance,
    linetype = significance
  )
) +
  
  geom_line(linewidth = 0.9) +
  
  geom_point(size = 2.5) +
  
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, 1)
  ) +
  
  facet_wrap(
    ~ design
  ) +
  
  labs(
    title = "Empirical Power of White's Test",
    subtitle = "Power under alternative heteroskedasticity designs",
    x = "Sample size (n)",
    y = "Empirical power",
    linetype = "Significance level"
  ) +
  
  theme_minimal() +
  
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )


print(p_white_power)


# --------------------------------------------------
# 9.3. Save visualizations
# --------------------------------------------------

if (!dir.exists("figures")) {
  dir.create("figures")
}


ggsave(
  filename = "figures/white_test_size.png",
  plot = p_white_size,
  width = 8,
  height = 5,
  dpi = 300,
  bg = "white"
)


ggsave(
  filename = "figures/white_test_power.png",
  plot = p_white_power,
  width = 8,
  height = 5,
  dpi = 300,
  bg = "white"
)


cat("\n====================================\n")
cat("WHITE'S TEST FIGURES SAVED\n")
cat("====================================\n")
cat("- figures/white_test_size.png\n")
cat("- figures/white_test_power.png\n")

# --------------------------------------------------
# 10. White robust variance simulation
# --------------------------------------------------

simulate_white_variance <- function(
    n,
    design,
    R = 5000,
    use_true_errors = FALSE
) {
  
  # Number of repetitions of the base sample
  k <- n / 20
  
  # Construct regressors
  x1 <- rep(x1_20, k)
  x2 <- rep(x2_20, k)
  
  
  # Storage for coefficient estimates
  beta_hat <- matrix(
    0,
    nrow = R,
    ncol = 3
  )
  
  # Storage for White covariance matrices
  vcov_white <- array(
    0,
    dim = c(3, 3, R)
  )
  
  
  # ------------------------------------------------
  # Monte Carlo loop
  # ------------------------------------------------
  
  for (r in 1:R) {
    
    # ----------------------------------------------
    # Design 1
    # ----------------------------------------------
    
    if (design == 1) {
      
      nu <- exp(
        0.25 * x1 +
          0.25 * x2
      )
      
      eps <- rnorm(n)
    }
    
    
    # ----------------------------------------------
    # Design 2
    # ----------------------------------------------
    
    if (design == 2) {
      
      nu <- exp(
        0.25 * x1 +
          0.25 * x2
      )
      
      eps <- rt(
        n,
        df = 5
      )
    }
    
    
    # ----------------------------------------------
    # True error
    # ----------------------------------------------
    
    u <- sqrt(nu) * eps
    
    
    # ----------------------------------------------
    # Dependent variable
    # ----------------------------------------------
    
    y <- beta0 +
      beta1 * x1 +
      beta2 * x2 +
      u
    
    
    # ----------------------------------------------
    # OLS estimation
    # ----------------------------------------------
    
    model <- lm(
      y ~ x1 + x2
    )
    
    
    # Store coefficient estimates
    beta_hat[r, ] <- coef(model)
    
    
    # ----------------------------------------------
    # White covariance matrix
    # ----------------------------------------------
    
    if (use_true_errors == FALSE) {
      
      # White HC0 using estimated residuals
      vcov_white[, , r] <- vcovHC(
        model,
        type = "HC0"
      )
      
    } else {
      
      # White covariance using true errors
      X <- model.matrix(model)
      
      Omega_true <- diag(
        u^2
      )
      
      XtX_inv <- solve(
        t(X) %*% X
      )
      
      vcov_white[, , r] <-
        XtX_inv %*%
        t(X) %*%
        Omega_true %*%
        X %*%
        XtX_inv
    }
  }
  
  
  # ------------------------------------------------
  # True Monte Carlo variance
  # ------------------------------------------------
  
  V_true <- cov(beta_hat)
  
  
  # ------------------------------------------------
  # Average White covariance matrix
  # ------------------------------------------------
  
  V_white <- apply(
    vcov_white,
    c(1, 2),
    mean
  )
  
  
  # ------------------------------------------------
  # Relative bias of variances
  # ------------------------------------------------
  
  relative_bias <- (
    diag(V_white) -
      diag(V_true)
  ) / diag(V_true)
  
  
  # ------------------------------------------------
  # Total absolute relative bias
  # ------------------------------------------------
  
  total_bias <- sum(
    abs(relative_bias)
  )
  
  
  # ------------------------------------------------
  # Return results
  # ------------------------------------------------
  
  data.frame(
    b0 = relative_bias[1],
    b1 = relative_bias[2],
    b2 = relative_bias[3],
    total = total_bias
  )
}

# --------------------------------------------------
# 11. White variance using estimated residuals
# --------------------------------------------------

white_variance_results <- list()

for (n in sample_sizes) {
  
  cat("\n====================================\n")
  cat("White variance simulation, n =", n, "\n")
  cat("====================================\n")
  
  
  # Design 1
  cat("Design 1\n")
  
  result_d1 <- simulate_white_variance(
    n = n,
    design = 1,
    R = R,
    use_true_errors = FALSE
  )
  
  
  # Design 2
  cat("Design 2\n")
  
  result_d2 <- simulate_white_variance(
    n = n,
    design = 2,
    R = R,
    use_true_errors = FALSE
  )
  
  
  white_variance_results[[paste0("n", n)]] <- list(
    design1 = result_d1,
    design2 = result_d2
  )
}

# --------------------------------------------------
# 12. Results table: estimated residuals
# --------------------------------------------------

white_variance_table <- data.frame(
  
  n = sample_sizes,
  
  
  # Design 1
  b0_d1 = sapply(
    white_variance_results,
    function(x) x$design1$b0
  ),
  
  b1_d1 = sapply(
    white_variance_results,
    function(x) x$design1$b1
  ),
  
  b2_d1 = sapply(
    white_variance_results,
    function(x) x$design1$b2
  ),
  
  total_d1 = sapply(
    white_variance_results,
    function(x) x$design1$total
  ),
  
  
  # Design 2
  b0_d2 = sapply(
    white_variance_results,
    function(x) x$design2$b0
  ),
  
  b1_d2 = sapply(
    white_variance_results,
    function(x) x$design2$b1
  ),
  
  b2_d2 = sapply(
    white_variance_results,
    function(x) x$design2$b2
  ),
  
  total_d2 = sapply(
    white_variance_results,
    function(x) x$design2$total
  )
)

print(white_variance_table)

# --------------------------------------------------
# 13. White variance using true errors
# --------------------------------------------------

white_true_results <- list()

for (n in sample_sizes) {
  
  cat("\n====================================\n")
  cat("White true-error simulation, n =", n, "\n")
  cat("====================================\n")
  
  
  # Design 1
  cat("Design 1\n")
  
  result_d1 <- simulate_white_variance(
    n = n,
    design = 1,
    R = R,
    use_true_errors = TRUE
  )
  
  
  # Design 2
  cat("Design 2\n")
  
  result_d2 <- simulate_white_variance(
    n = n,
    design = 2,
    R = R,
    use_true_errors = TRUE
  )
  
  
  white_true_results[[paste0("n", n)]] <- list(
    design1 = result_d1,
    design2 = result_d2
  )
}

# --------------------------------------------------
# 14. Results table: true errors
# --------------------------------------------------

white_true_table <- data.frame(
  
  n = sample_sizes,
  
  
  # Design 1
  b0_d1 = sapply(
    white_true_results,
    function(x) x$design1$b0
  ),
  
  b1_d1 = sapply(
    white_true_results,
    function(x) x$design1$b1
  ),
  
  b2_d1 = sapply(
    white_true_results,
    function(x) x$design1$b2
  ),
  
  total_d1 = sapply(
    white_true_results,
    function(x) x$design1$total
  ),
  
  
  # Design 2
  b0_d2 = sapply(
    white_true_results,
    function(x) x$design2$b0
  ),
  
  b1_d2 = sapply(
    white_true_results,
    function(x) x$design2$b1
  ),
  
  b2_d2 = sapply(
    white_true_results,
    function(x) x$design2$b2
  ),
  
  total_d2 = sapply(
    white_true_results,
    function(x) x$design2$total
  )
)

print(white_true_table)

# --------------------------------------------------
# 15. Save White variance results
# --------------------------------------------------

if (!dir.exists("outputs")) {
  dir.create("outputs")
}


write.csv(
  white_variance_table,
  "outputs/white_variance_residuals.csv",
  row.names = FALSE
)


write.csv(
  white_true_table,
  "outputs/white_variance_true_errors.csv",
  row.names = FALSE
)


cat("\n====================================\n")
cat("WHITE VARIANCE RESULTS SAVED\n")
cat("====================================\n")

cat("- outputs/white_variance_residuals.csv\n")
cat("- outputs/white_variance_true_errors.csv\n")


# ==================================================
# 16. VISUALIZATIONS - PART 2.2
# ==================================================

library(ggplot2)
library(tidyr)
library(dplyr)


# --------------------------------------------------
# Create folder for Part 2.2 figures
# --------------------------------------------------

if (!dir.exists("figures_22")) {
  dir.create("figures_22")
}


# ==================================================
# FIGURE 1
# Total relative bias - estimated residuals
# ==================================================

plot_total_residuals <- white_variance_table %>%
  
  select(
    n,
    total_d1,
    total_d2
  ) %>%
  
  pivot_longer(
    cols = c(total_d1, total_d2),
    names_to = "design",
    values_to = "total_bias"
  ) %>%
  
  mutate(
    design = case_when(
      design == "total_d1" ~ "Design 1",
      design == "total_d2" ~ "Design 2"
    )
  ) %>%
  
  ggplot(
    aes(
      x = factor(n, levels = sample_sizes),
      y = total_bias,
      linetype = design,
      shape = design
    )
  ) +
  
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  
  geom_line(
    linewidth = 1
  ) +
  
  geom_point(
    size = 3
  ) +
  
  labs(
    title = "Total Relative Bias of Variance Estimates",
    subtitle = "White covariance matrix using estimated residuals",
    x = "Sample size (n)",
    y = "Total relative bias",
    linetype = "Design",
    shape = "Design"
  ) +
  
  theme_minimal(base_size = 13) +
  
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )


print(plot_total_residuals)


ggsave(
  filename = "figures_22/figure_1_total_bias_residuals.png",
  plot = plot_total_residuals,
  width = 10,
  height = 6,
  dpi = 300,
  bg = "white"
)


# ==================================================
# FIGURE 2
# Individual relative bias - estimated residuals
# ==================================================

plot_individual_residuals <- white_variance_table %>%
  
  select(
    n,
    b0_d1,
    b1_d1,
    b2_d1,
    b0_d2,
    b1_d2,
    b2_d2
  ) %>%
  
  pivot_longer(
    cols = -n,
    names_to = "variable",
    values_to = "bias"
  ) %>%
  
  mutate(
    design = ifelse(
      grepl("_d1$", variable),
      "Design 1",
      "Design 2"
    ),
    
    parameter = case_when(
      grepl("^b0", variable) ~ "Beta0",
      grepl("^b1", variable) ~ "Beta1",
      grepl("^b2", variable) ~ "Beta2"
    )
  ) %>%
  
  ggplot(
    aes(
      x = factor(n, levels = sample_sizes),
      y = bias,
      linetype = parameter,
      shape = parameter
    )
  ) +
  
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  
  geom_line(
    linewidth = 0.9
  ) +
  
  geom_point(
    size = 2.5
  ) +
  
  facet_wrap(
    ~ design
  ) +
  
  labs(
    title = "Relative Bias of Individual Variance Estimates",
    subtitle = "White covariance matrix using estimated residuals",
    x = "Sample size (n)",
    y = "Relative bias",
    linetype = "Parameter",
    shape = "Parameter"
  ) +
  
  theme_minimal(base_size = 13) +
  
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )


print(plot_individual_residuals)


ggsave(
  filename = "figures_22/figure_2_individual_bias_residuals.png",
  plot = plot_individual_residuals,
  width = 10,
  height = 6,
  dpi = 300,
  bg = "white"
)


# ==================================================
# FIGURE 3
# Total relative bias - true errors
# ==================================================

plot_total_true <- white_true_table %>%
  
  select(
    n,
    total_d1,
    total_d2
  ) %>%
  
  pivot_longer(
    cols = c(total_d1, total_d2),
    names_to = "design",
    values_to = "total_bias"
  ) %>%
  
  mutate(
    design = case_when(
      design == "total_d1" ~ "Design 1",
      design == "total_d2" ~ "Design 2"
    )
  ) %>%
  
  ggplot(
    aes(
      x = factor(n, levels = sample_sizes),
      y = total_bias,
      linetype = design,
      shape = design
    )
  ) +
  
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  
  geom_line(
    linewidth = 1
  ) +
  
  geom_point(
    size = 3
  ) +
  
  labs(
    title = "Total Relative Bias of Variance Estimates",
    subtitle = "White covariance matrix using true errors",
    x = "Sample size (n)",
    y = "Total relative bias",
    linetype = "Design",
    shape = "Design"
  ) +
  
  theme_minimal(base_size = 13) +
  
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )


print(plot_total_true)


ggsave(
  filename = "figures_22/figure_3_total_bias_true_errors.png",
  plot = plot_total_true,
  width = 10,
  height = 6,
  dpi = 300,
  bg = "white"
)


# ==================================================
# FIGURE 4
# Individual relative bias - true errors
# ==================================================

plot_individual_true <- white_true_table %>%
  
  select(
    n,
    b0_d1,
    b1_d1,
    b2_d1,
    b0_d2,
    b1_d2,
    b2_d2
  ) %>%
  
  pivot_longer(
    cols = -n,
    names_to = "variable",
    values_to = "bias"
  ) %>%
  
  mutate(
    design = ifelse(
      grepl("_d1$", variable),
      "Design 1",
      "Design 2"
    ),
    
    parameter = case_when(
      grepl("^b0", variable) ~ "Beta0",
      grepl("^b1", variable) ~ "Beta1",
      grepl("^b2", variable) ~ "Beta2"
    )
  ) %>%
  
  ggplot(
    aes(
      x = factor(n, levels = sample_sizes),
      y = bias,
      linetype = parameter,
      shape = parameter
    )
  ) +
  
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  
  geom_line(
    linewidth = 0.9
  ) +
  
  geom_point(
    size = 2.5
  ) +
  
  facet_wrap(
    ~ design
  ) +
  
  labs(
    title = "Relative Bias of Individual Variance Estimates",
    subtitle = "White covariance matrix using true errors",
    x = "Sample size (n)",
    y = "Relative bias",
    linetype = "Parameter",
    shape = "Parameter"
  ) +
  
  theme_minimal(base_size = 13) +
  
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )


print(plot_individual_true)


ggsave(
  filename = "figures_22/figure_4_individual_bias_true_errors.png",
  plot = plot_individual_true,
  width = 10,
  height = 6,
  dpi = 300,
  bg = "white"
)


# --------------------------------------------------
# Confirmation
# --------------------------------------------------

cat("\n====================================\n")
cat("PART 2.2 FIGURES SAVED\n")
cat("====================================\n\n")

cat("1. figure_1_total_bias_residuals.png\n")
cat("2. figure_2_individual_bias_residuals.png\n")
cat("3. figure_3_total_bias_true_errors.png\n")
cat("4. figure_4_individual_bias_true_errors.png\n")