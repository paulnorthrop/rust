#' Internal rust functions
#'
#' Internal rust functions
#' @details
#' These functions are not intended to be called by the user.
#' @name rust-internal
#' @keywords internal
NULL

# =========================== box_cox ===========================

#' @keywords internal
#' @rdname rust-internal
box_cox <- function (x, lambda = 1, gm = 1, lambda_tol = 1e-6) {
  #
  # Computes the Box-Cox transformation of a vector.
  #
  # Args:
  #   x          : A numeric vector. (Positive) values to be Box-Cox
  #                transformed.
  #   lambda     : A numeric scalar.  Transformation parameter.
  #   gm         : A numeric scalar.  Optional scaling parameter.
  #   lambda_tol : A numeric scalar.  For abs(lambda) < lambda_tol use
  #                a Taylor series expansion.
  #
  # Returns:
  #   A numeric vector.  The transformed value
  #     (x^lambda - 1) / (lambda * gm ^ (lambda - 1))
  #
  if (abs(lambda) > lambda_tol) {
    retval <- (x ^ lambda - 1) / lambda / gm ^ (lambda - 1)
  } else {
    i <- 0:3
    retval <- sum(log(x) ^ (i + 1) * lambda ^ i / factorial(i + 1))
    retval <- retval / gm ^ (lambda - 1)
  }
  retval
}

# =========================== box_cox_vec ===========================

# Version of box_cox vectorized for lambda and gm.

#' @keywords internal
#' @rdname rust-internal
box_cox_vec <- Vectorize(box_cox, vectorize.args = c("lambda", "gm"))

# =========================== optim_box_cox ===========================

#' @keywords internal
#' @rdname rust-internal
optim_box_cox <- function(x, w, lambda_range = c(-3,3), start = NULL,
                          which_lam = 1:ncol(x)) {
  #
  # Finds the optimal value of the Box-Cox transformation parameter lambda,
  # based values of a probability density function.
  #
  # Args:
  #   x            : A numeric matrix. Values at which the density is evaluated.
  #                  each row contains a combination of the ncol(x) variables.
  #                  column numbers in which_lam must contain positive values.
  #   w            : A numeric vector. Density values corresponding to each
  #                  row of x (up to proportionality).
  #   lambda_range : A numeric vector (of length 2).  Range of lambda values
  #                  over which to search.
  #   start        : A numeric vector.  Optional starting value for lambda.
  #   which_lam    : A numeric vector.  Indicates which variables to Box-Cox
  #                  transform.
  #
  # Returns: a list containing
  #   lambda       : A numeric vector.  The optimal value of lambda.
  #   gm           : A numeric vector.  Geometric mean of x weighted by w.
  #
  n_var <- ncol(x)
  neg_loglik <- function(lambda_in, x, gm, w, which_lam) {
    lambda <- rep(1, n_var)
    lambda[which_lam] <- lambda_in
    for (j in 1:n_var) {
      x[, j] <- box_cox(x = x[, j], lambda = lambda[j], gm = gm[j])
    }
    (nrow(x) / 2) * log(det(stats::cov.wt(x, w, method = "ML")$cov))
  }
  w_n <- w / mean(w)
  n_lam <- length(which_lam)
  #
  if (any(x[, which_lam] <= 0)) {
    stop("All values must be > 0")
  }
  gm_w <- rep(1, n_var)
  x_mat <- x[, which_lam, drop = FALSE]
  gm_w[which_lam] <- apply(x_mat, 2, function(x) exp(mean(w_n * log(x))))
  #
  skew_wt <- function(x, w) {
    vp <- function(p) mean(w^p)
    v1 <- vp(1)
    xbar <- mean(w * x) / v1
    mp <- function(p) mean(w * (x - xbar) ^ p) / v1
    mp(3) / mp(2)^(3 / 2)
  }
  if (is.null(start)) {
    start <- 1 - apply(x, 2, skew_wt, w = w) / 2
  }
  if (n_lam == 1L) {
    method <- "Brent"
  } else {
    method <- "L-BFGS-B"
  }
  lower <- lambda_range[1]
  upper <- lambda_range[2]
  start <- start[which_lam]
  ret <- stats::optim(start, neg_loglik, method = method, x = x, gm = gm_w,
                      w = w_n, lower = lower, upper = upper,
                      which_lam = which_lam)
  lambda <- rep(1L, n_var)
  lambda[which_lam] <- ret$par
  ret$par <- lambda
  if (ret$convergence != 0) {
    warning(paste("Convergence failure: return code =", ret$convergence))
  }
  list(lambda = ret$par, gm = gm_w)
}

# =========================== n_grid_fun ===========================

#' @keywords internal
#' @rdname rust-internal
n_grid_fn <- function(d) {
  # Sets a default values value of n_grid.
  #
  # Args:
  #   d          : A numeric scalar.  Dimension of the target density.
  #
  # Returns:
  #   A numeric scalar.  The value of n_grid.
  #
  return(ceiling(2501 ^ (1 / d)))
}

# =========================== init_ps_calc ===========================

#' @keywords internal
#' @rdname rust-internal
init_psi_calc <- function(phi_to_psi, phi, lambda, gm, w, which_lam){
  # Estimates the mode and standard deviation of the Box-Cox transformed
  # target density.
  #
  # Args:
  #   phi_to_psi : A function.  The Box-Cox transformation function.
  #   phi        : A numeric matrix.  n_grid by d matrix of values of phi
  #   lambda     : A numeric vector.  The value of the Box-Cox transformation
  #                parameter lambda.
  #   gm         : A numeric vector.  The value of the Box-Cox scale
  #                parameter gm.
  #   w          : A numeric vector. The values of the target density at each
  #                value of the variables in phi.
  #   which_lam  : A numeric vector.  Indicates which variables to Box-Cox
  #                transform.
  #
  # Returns: A list containing
  #   init_psi : estimate of the mode of the Box-Cox transformed target
  #              density.
  #   sd_psi   : estimate of the standard deviation of the Box-Cox transformed
  #              target density.
  #
  psi <- phi_to_psi(phi)
  d <- ncol(phi)
  log_jac <- matrix(0, nrow(phi), d)
  for (j in which_lam) {
    log_jac[, j] <- (lambda[j] - 1) * log(phi[, j] / gm[j])
  }
  log_jac <- rowSums(log_jac)
  w_psi <- w * exp(-log_jac)
  temp_cov <- stats::cov.wt(psi, w_psi, method="ML")
  init_psi <- as.numeric(psi[which.max(w_psi), ])
  sd_psi <- sqrt(diag(temp_cov$cov))
  list(init_psi = init_psi, sd_psi = sd_psi)
}
