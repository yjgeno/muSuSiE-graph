## joint inference
# dta_1 and dta_2 are n x p data set
# scale_x : scale the data
# intercept: calculate the mean of Y
# sigma02_int is initialization for signal prior variance
# sigma2_int is initialization for error variance
# prior_vec is prior for common part and for single part
# itermax is the maximum iteration
# L is the largest number of parents
# tol is the threshold for ELBO
# sigma0_low_bd is the threshold for select effect l
# residual_variance_lowerbound is the lower bound for sigma2

## load variable selection function
source("twodatasets/sum_single_effect_two_graph.R")
joint_graph_fun_two <- function(dta_1, dta_2, scale_x = FALSE, intercept = TRUE,
                                sigma02_int = NULL, sigma2_int = NULL, prior_vec = NULL,
                                itermax = 100, L = 10, tol = 1e-4, sigma0_low_bd = 1e-8,
                                residual_variance_lowerbound = NULL) {
  ## Initialization
  p <- ncol(dta_1)
  if (p != ncol(dta_2)) stop("The number of features should be same!")
  n1 <- nrow(dta_1)
  n2 <- nrow(dta_2)
  ## define prior vector
  if (is.null(prior_vec)) {
    prior_vec <- c(1 / (2 * p^1.5), 1 / (p^2))
  }
  ## save matrix
  # probability of the edge exists
  alpha_res_1 <- matrix(0, nrow = p, ncol = p)
  alpha_res_2 <- matrix(0, nrow = p, ncol = p)
  # posterior mean of the edge
  A_res_1 <- matrix(0, nrow = p, ncol = p)
  A_res_2 <- matrix(0, nrow = p, ncol = p)
  # predicted value
  sigma2_vec <- rep(NA, p)
  # check intercept
  if (intercept) {
    mean_1 <- mean(dta_1[, 1])
    mean_2 <- mean(dta_2[, 1])
  } else {
    mean_1 <- 0
    mean_2 <- 0
  }
  # begin iteration
  for (iter_p in seq_len(p)) {
    X_1 <- dta_1[, -iter_p]
    Y_1 <- dta_1[, iter_p]
    X_2 <- dta_2[, -iter_p]
    Y_2 <- dta_2[, iter_p]
    ## variable selection
    res <- sum_single_effect_two(
      X_1 = X_1, Y_1 = Y_1, X_2 = X_2, Y_2 = Y_2,
      scale_x = scale_x, intercept = intercept, sigma02_int = sigma02_int,
      sigma2_int = sigma2_int, prior_vec = prior_vec, L = L,
      itermax = itermax, tol = tol, sigma0_low_bd = sigma0_low_bd,
      residual_variance_lowerbound = residual_variance_lowerbound
    )
    # save the matrix we want
    alpha_res_1[-iter_p, iter_p] <- res$alpha_1
    alpha_res_2[-iter_p, iter_p] <- res$alpha_2
    A_res_1[-iter_p, iter_p] <- res$post_mean1
    A_res_2[-iter_p, iter_p] <- res$post_mean2
    # calculate the likelihood
    sigma2_vec[iter_p] <- res$sigma2
  }
  ## return results
  return(list(
    alpha_res_1 = alpha_res_1, alpha_res_2 = alpha_res_2, A_res_1 = A_res_1, A_res_2 = A_res_2,
    sigma2_vec = sigma2_vec
  ))
}