## joint inference
# dta_list are n x p data set
# scale_x : scale the data
# intercept: calculate the mean of Y
# sigma02_int is initialization for signal prior variance
# sigma2_int is initialization for error variance
# prior_vec : prior for different models
# itermax is the maximum iteration
# L_max is the maximum number of parents
# tol is the threshold for ELBO
# sigma0_low_bd is the threshold for select effect l
# residual_variance_lowerbound is the lower bound for sigma2

source("mutipledatasets/sum_single_effect_multi_graph.R")
joint_graph_multi <- function(dta_list, scale_x = FALSE, intercept = TRUE,
                              sigma02_int = NULL, sigma2_int = NULL, prior_vec = NULL,
                              com_mat = NULL, com_list = NULL,
                              itermax = 100, L_max = 10, tol = 1e-4, sigma0_low_bd = 1e-8,
                              residual_variance_lowerbound = NULL) {
  ## Initialization
  K <- length(dta_list)
  n_group <- 2^K - 1
  p <- rep(NA, K)
  n <- rep(NA, K)
  for (iter_K in seq_len(K)) {
    p[iter_K] <- ncol(dta_list[[iter_K]])
    n[iter_K] <- nrow(dta_list[[iter_K]])
  }
  if (length(unique(p)) > 1) {
    stop("The number of features should be same!")
  } else {
    p <- p[1]
  }
  ## define prior vector
  if (is.null(prior_vec)) {
    prior_vec <- rep(1 / (n_group * p^1.5), n_group)
  }
  # combinatorics matrix
  if (is.null(com_mat)) {
    com_list <- list()
    com_mat <- matrix(c(0, 1), ncol = 1)
    for (iter in 2:K) {
      com_mat_copy <- com_mat
      com_mat <- cbind(1, com_mat)
      com_mat_copy <- cbind(0, com_mat_copy)
      com_mat <- rbind(com_mat_copy, com_mat)
    }
    com_mat <- com_mat[-1, ]

    for (iter_com in seq_len(n_group)) {
      com_list[[iter_com]] <- which(com_mat[iter_com, ] == 1)
    }
  }
  ## save list
  alpha_list <- list()
  A_list <- list()
  mean_list <- list()
  Y_list <- list()
  for (iter_K in seq_len(K)) {
    alpha_list[[iter_K]] <- matrix(0, nrow = p, ncol = p)
    A_list[[iter_K]] <- matrix(0, nrow = p, ncol = p)
    Y_list[[iter_K]] <- dta_list[[iter_K]][, 1]
    if (intercept) {
      mean_list[[iter_K]] <- mean(Y_list[[iter_K]])
    } else {
      mean_list[[iter_K]] <- 0
    }
  }
  sigma2_vec <- rep(NA, p)
  # begin iteration
  for (iter_p in seq_len(p)) {
    dta_vs_list <- list()
    for (iter_K in seq_len(K)) {
      dta_vs_list[[iter_K]] <- list()
      dta_vs_list[[iter_K]]$X <- dta_list[[iter_K]][, -iter_p]
      dta_vs_list[[iter_K]]$Y <- dta_list[[iter_K]][, iter_p]
    }
    ## variable selection
    res_vs <- sum_single_effect_multi(dta_vs_list,
      scale_x = scale_x, intercept = intercept,
      sigma02_int = sigma02_int, sigma2_int = sigma2_int,
      prior_vec = prior_vec, L = min(iter_p, L_max), itermax = itermax,
      tol = tol, sigma0_low_bd = sigma0_low_bd, com_mat = com_mat, com_list = com_list,
      residual_variance_lowerbound = residual_variance_lowerbound
    )
    ## save needed list
    sigma2_vec[iter_p] <- res_vs$sigma2
    for (iter_K in seq_len(K)) {
      alpha_list[[iter_K]][-iter_p, iter_p] <- res_vs$res[[iter_K]]$alpha
      A_list[[iter_K]][-iter_p, iter_p] <- res_vs$res[[iter_K]]$post_mean
    }
  }

  ## return results
  return(list(
    alpha_list = alpha_list, A_list = A_list, sigma2_vec = sigma2_vec
  ))
}