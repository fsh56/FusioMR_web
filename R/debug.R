# Comprehensive debugging function
debug_gibbs_results <- function(res, niter, burnin_prop) {
  cat("=== DEBUGGING GIBBS SAMPLER RESULTS ===\n")

  # Check raw results structure
  cat("\n1. Raw Results Structure:\n")
  cat("pst_tk dimensions:", dim(res$pst_tk), "\n")
  cat("beta_1_tk length:", length(res$beta_1_tk), "\n")
  cat("beta_2_tk length:", length(res$beta_2_tk), "\n")
  cat("eta_1_tk length:", length(res$eta_1_tk), "\n")
  cat("eta_2_tk length:", length(res$eta_2_tk), "\n")
  cat("alpha_1_tk length:", length(res$alpha_1_tk), "\n")
  cat("alpha_2_tk length:", length(res$alpha_2_tk), "\n")

  # Check for NA/NaN in raw results
  cat("\n2. NA/NaN Check in Raw Results:\n")
  cat("pst_tk has NA:", any(is.na(res$pst_tk)), "\n")
  cat("pst_tk has NaN:", any(is.nan(res$pst_tk)), "\n")
  cat("beta_1_tk has NA:", any(is.na(res$beta_1_tk)), "\n")
  cat("beta_2_tk has NA:", any(is.na(res$beta_2_tk)), "\n")
  cat("alpha_1_tk has NA:", any(is.na(res$alpha_1_tk)), "\n")
  cat("alpha_2_tk has NA:", any(is.na(res$alpha_2_tk)), "\n")

  # Check value ranges
  cat("\n3. Value Ranges in Raw Results:\n")
  cat("pst_tk range:", range(res$pst_tk, na.rm = TRUE), "\n")
  cat("beta_1_tk range:", range(res$beta_1_tk, na.rm = TRUE), "\n")
  cat("beta_2_tk range:", range(res$beta_2_tk, na.rm = TRUE), "\n")
  cat("alpha_1_tk range:", range(res$alpha_1_tk, na.rm = TRUE), "\n")
  cat("alpha_2_tk range:", range(res$alpha_2_tk, na.rm = TRUE), "\n")

  # Check first few values
  cat("\n4. First 5 values:\n")
  cat("pst_tk[1:5, ]:\n")
  print(res$pst_tk[1:5, ])
  cat("beta_1_tk[1:5]:", res$beta_1_tk[1:5], "\n")
  cat("beta_2_tk[1:5]:", res$beta_2_tk[1:5], "\n")

  # Check post-burnin results
  cat("\n5. Post-burnin Analysis:\n")
  burnin <- floor(niter * burnin_prop)
  cat("Burnin period:", burnin, "iterations\n")
  cat("Post-burnin samples:", niter - burnin, "iterations\n")

  post_res <- get_post_burnin_res(res, niter, burnin_prop)

  cat("Post-burnin pst_tk dimensions:", dim(post_res$pst_tk), "\n")
  cat("Post-burnin beta_1_tk length:", length(post_res$beta_1_tk), "\n")

  # Check probability columns
  if (is.matrix(post_res$pst_tk)) {
    cat("\n6. Probability Analysis:\n")
    p00_post <- post_res$pst_tk[, 1]
    p01_post <- post_res$pst_tk[, 2]
    p10_post <- post_res$pst_tk[, 3]
    p11_post <- post_res$pst_tk[, 4]

    cat("p00 range:", range(p00_post, na.rm = TRUE), "\n")
    cat("p01 range:", range(p01_post, na.rm = TRUE), "\n")
    cat("p10 range:", range(p10_post, na.rm = TRUE), "\n")
    cat("p11 range:", range(p11_post, na.rm = TRUE), "\n")

    # Check if probabilities sum to 1
    prob_sums <- rowSums(post_res$pst_tk)
    cat("Probability sums range:", range(prob_sums, na.rm = TRUE), "\n")
    cat("Mean probability sum:", mean(prob_sums, na.rm = TRUE), "\n")

    qq1 <- mean(p10_post + p11_post, na.rm = TRUE)
    qq2 <- mean(p01_post + p11_post, na.rm = TRUE)
    cat("qq1 (marginal prob 1):", qq1, "\n")
    cat("qq2 (marginal prob 2):", qq2, "\n")
  }

  # Check parameter combinations
  cat("\n7. Parameter Combination Analysis:\n")
  beta1_alpha1 <- post_res$beta_1_tk + post_res$alpha_1_tk
  beta2_alpha2 <- post_res$beta_2_tk + post_res$alpha_2_tk

  cat("beta_1_tk has NA:", any(is.na(post_res$beta_1_tk)), "\n")
  cat("alpha_1_tk has NA:", any(is.na(post_res$alpha_1_tk)), "\n")
  cat("beta_1 + alpha_1 has NA:", any(is.na(beta1_alpha1)), "\n")
  cat("beta_1 + alpha_1 range:", range(beta1_alpha1, na.rm = TRUE), "\n")

  cat("beta_2_tk has NA:", any(is.na(post_res$beta_2_tk)), "\n")
  cat("alpha_2_tk has NA:", any(is.na(post_res$alpha_2_tk)), "\n")
  cat("beta_2 + alpha_2 has NA:", any(is.na(beta2_alpha2)), "\n")
  cat("beta_2 + alpha_2 range:", range(beta2_alpha2, na.rm = TRUE), "\n")

  cat("\n=== END DEBUGGING ===\n")
}

# Enhanced label flip function with debugging
label_flip_joint_debug <- function(res, eta_true_1, eta_true_2, niter, burnin_prop) {
  cat("\n=== LABEL FLIP JOINT DEBUG ===\n")

  post_res <- get_post_burnin_res(res, niter, burnin_prop)

  # Check dimensions
  if (!is.matrix(post_res$pst_tk)) {
    stop("pst_tk is not a matrix after burnin removal")
  }

  if (ncol(post_res$pst_tk) != 4) {
    stop("pst_tk should have 4 columns, has ", ncol(post_res$pst_tk))
  }

  p00_post <- post_res$pst_tk[, 1]
  p01_post <- post_res$pst_tk[, 2]
  p10_post <- post_res$pst_tk[, 3]
  p11_post <- post_res$pst_tk[, 4]

  qq1 <- mean(p10_post + p11_post, na.rm = TRUE)
  qq2 <- mean(p01_post + p11_post, na.rm = TRUE)

  cat("qq1 (should trigger label flip if > 0.5):", qq1, "\n")
  cat("qq2 (should trigger label flip if > 0.5):", qq2, "\n")

  # Beta 1 estimation with debugging
  cat("\n--- Beta 1 Estimation ---\n")
  if (is.na(qq1)) {
    cat("ERROR: qq1 is NA\n")
    return(NULL)
  }

  if (qq1 > 0.5) {
    cat("Using flipped beta_1 (beta_1 + alpha_1)\n")
    beta1_values <- post_res$beta_1_tk + post_res$alpha_1_tk
  } else {
    cat("Using original beta_1\n")
    beta1_values <- post_res$beta_1_tk
  }

  cat("Beta1 values has NA:", any(is.na(beta1_values)), "\n")
  cat("Beta1 values range:", range(beta1_values, na.rm = TRUE), "\n")

  b1_mean <- mean(beta1_values, na.rm = TRUE)
  b1_sd <- sd(beta1_values, na.rm = TRUE)
  bci1 <- quantile(beta1_values, probs = c(0.025, 0.975), na.rm = TRUE)

  cat("b1_mean:", b1_mean, "\n")
  cat("b1_sd:", b1_sd, "\n")
  cat("bci1:", bci1, "\n")

  # Beta 2 estimation with debugging
  cat("\n--- Beta 2 Estimation ---\n")
  if (is.na(qq2)) {
    cat("ERROR: qq2 is NA\n")
    return(NULL)
  }

  if (qq2 > 0.5) {
    cat("Using flipped beta_2 (beta_2 + alpha_2)\n")
    beta2_values <- post_res$beta_2_tk + post_res$alpha_2_tk
  } else {
    cat("Using original beta_2\n")
    beta2_values <- post_res$beta_2_tk
  }

  cat("Beta2 values has NA:", any(is.na(beta2_values)), "\n")
  cat("Beta2 values range:", range(beta2_values, na.rm = TRUE), "\n")

  b2_mean <- mean(beta2_values, na.rm = TRUE)
  b2_sd <- sd(beta2_values, na.rm = TRUE)
  bci2 <- quantile(beta2_values, probs = c(0.025, 0.975), na.rm = TRUE)

  cat("b2_mean:", b2_mean, "\n")
  cat("b2_sd:", b2_sd, "\n")
  cat("bci2:", bci2, "\n")

  # P-value calculation with debugging
  cat("\n--- P-value Calculation ---\n")
  if (is.na(b1_mean) || is.na(b1_sd) || b1_sd == 0) {
    cat("ERROR: Cannot calculate p-value for beta1. Mean:", b1_mean, "SD:", b1_sd, "\n")
    est_pval1 <- NA
  } else {
    est_pval1 <- 2 * (1 - stats::pnorm(abs(b1_mean / b1_sd)))
    cat("est_pval1:", est_pval1, "\n")
  }

  if (is.na(b2_mean) || is.na(b2_sd) || b2_sd == 0) {
    cat("ERROR: Cannot calculate p-value for beta2. Mean:", b2_mean, "SD:", b2_sd, "\n")
    est_pval2 <- NA
  } else {
    est_pval2 <- 2 * (1 - stats::pnorm(abs(b2_mean / b2_sd)))
    cat("est_pval2:", est_pval2, "\n")
  }

  cat("=== END LABEL FLIP DEBUG ===\n")

  return(list(
    beta_est1 = b1_mean,
    beta_se1 = b1_sd,
    beta_est2 = b2_mean,
    beta_se2 = b2_sd,
    beta_pval1 = est_pval1,
    beta_pval2 = est_pval2,
    ci_emp1 = bci1,
    ci_emp2 = bci2
  ))
}

# Updated get_post_burnin_res function
get_post_burnin_res <- function(res, niter, burnin_prop) {
  burnin <- floor(niter * burnin_prop)

  post_res <- list(
    pst_tk = res$pst_tk[(burnin+1):niter, , drop = FALSE],
    beta_1_tk = res$beta_1_tk[(burnin+1):niter],
    beta_2_tk = res$beta_2_tk[(burnin+1):niter],
    eta_1_tk = res$eta_1_tk[(burnin+1):niter],
    eta_2_tk = res$eta_2_tk[(burnin+1):niter],
    alpha_1_tk = res$alpha_1_tk[(burnin+1):niter],
    alpha_2_tk = res$alpha_2_tk[(burnin+1):niter]
  )
  return(post_res)
}
