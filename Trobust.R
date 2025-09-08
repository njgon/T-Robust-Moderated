huber_estimation <- function(x, c = 1.345, tol = 1e-6, max_iter = 100) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  n <- length(x)
  mu <- median(x, na.rm = TRUE)
  mad <- median(abs(x - mu), na.rm = TRUE)
  sigma <- ifelse(is.finite(mad) && mad > 0, mad * 1.4826, sd(x, na.rm = TRUE))
  if (!is.finite(sigma) || sigma <= 0) sigma <- 1.0
  
  delta <- 2 * pnorm(c) - 1
  for (i in 1:max_iter) {
    r <- (x - mu) / sigma
    w <- ifelse(abs(r) <= c, 1, c / pmax(abs(r), .Machine$double.eps))
    mu_new <- sum(w * x, na.rm = TRUE) / sum(w, na.rm = TRUE)
    term_scale <- sum(w * (x - mu_new)^2, na.rm = TRUE)
    sigma_new <- sqrt(term_scale / (n * delta))
    if (!is.finite(mu_new) || !is.finite(sigma_new) || sigma_new <= 0) break
    if (abs(mu - mu_new) < tol && abs(sigma - sigma_new) < tol) { mu <- mu_new; sigma <- sigma_new; break }
    mu <- mu_new; sigma <- sigma_new
  }
  list(mu = mu, sigma = sigma)
}

robust_t_test <- function(x, y, c = 1.345) {
  n1 <- length(x)
  n2 <- length(y)
  
  est_x <- huber_estimation(x, c)
  est_y <- huber_estimation(y, c)
  
  mu_diff <- est_x$mu - est_y$mu
  
  var_x <- est_x$sigma^2
  var_y <- est_y$sigma^2
  
  se_rob <- sqrt(var_x/n1 + var_y/n2)
  t_rob <- mu_diff / se_rob
  
  num_df <- (var_x/n1 + var_y/n2)^2
  den_df <- (var_x/n1)^2/(n1-1) + (var_y/n2)^2/(n2-1)
  df_rob <- num_df / den_df
  
  p_value <- 2 * pt(-abs(t_rob), df_rob)
  
  result <- list(
    statistic = c(t_rob = t_rob),       
    parameter = c(df_rob = df_rob),    
    p.value = p_value,
    estimate = c(robust_mean_x = est_x$mu, robust_mean_y = est_y$mu),
    method = "Two-Sample Robust t-Test (Huber)",
    data.name = paste(deparse(substitute(x)), "and", deparse(substitute(y)))
  )
  
  class(result) <- "htest"
  return(result)
}

robust_huber_moderated <- function(counts, group, c_huber = 1.345, robust_prior = TRUE) {
  G <- nrow(counts)
  n1 <- sum(group == levels(group)[1])
  n2 <- sum(group == levels(group)[2])
  mu_diff <- numeric(G)
  v_diff <- numeric(G)
  for (g in seq_len(G)) {
    x <- log2(counts[g, group==levels(group)[1]] + 1)
    y <- log2(counts[g, group==levels(group)[2]] + 1)
    estx <- huber_estimation_safe(x, c = c_huber)
    esty <- huber_estimation_safe(y, c = c_huber)
    mu_diff[g] <- esty$mu - estx$mu
    v_diff[g] <- (estx$sigma^2) / n1 + (esty$sigma^2) / n2
  }
  # squeezeVar espera vector df por gen
  df_per_gene <- rep(n1 + n2 - 2, G)
  squeezed <- limma::squeezeVar(var = v_diff, df = df_per_gene, robust = robust_prior)
  df_total <- df_per_gene + squeezed$df.prior
  t_mod <- mu_diff / sqrt(squeezed$var.post)
  pvals <- 2 * pt(-abs(t_mod), df = df_total)
  pvals
}
