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

# =========================== log_gpd_mdi_prior ===========================

#' @keywords internal
#' @rdname rust-internal
log_gpd_mdi_prior <- function(pars) {
  # Generalized Pareto MDI prior log-density
  #
  # Calculates the log-density of a Maximal Data Information (MDI) prior
  # for generalized Pareto parameters.  The prior density is truncated,
  # (to xi >= -1), to produce a posterior density that is proper.
  #
  # Args:
  #   pars : A numeric vector containing the values of the generalized Pareto
  #          parameters sigma and xi.
  #
  # Returns:
  #   A numeric scalar. The value of the prior log-density.
  #
  if (pars[1] <= 0 | pars[2] < -1) {
    return(-Inf)
  }
  return(-log(pars[1]) - pars[2] - 1)
}

# =========================== gpd_loglik ===========================

#' @keywords internal
#' @rdname rust-internal
gpd_loglik <- function(pars, gpd_data, m, xm, sum_gp) {
  # Generalized Pareto log-likelihood
  #
  # Calculates the log-likelihood for a random sample from a Generalized Pareto.
  #
  # Args:
  #   pars     : A numeric vector containing the values of the generalized
  #              Pareto parameters sigma and xi.
  #   gpd_data : A numeric vector containing positive sample values.
  #   m        : A numeric scalar.  The sample size: the length of gpd_data.
  #   xm       : A numeric scalar. The sample maximum.
  #   sum_gp   : The sum of the sample values.
  #
  # Returns:
  #   A numeric scalar. The value of the log-likelihood.
  #
  if (pars[1] <= 0 | pars[2] <= -pars[1] / xm)
    return(-Inf)
  sdat <- gpd_data / pars[1]
  zz <- 1 + pars[2] * sdat
  if (abs(pars[2]) > 1e-6) {
    val <- -m * log(pars[1]) - (1 + 1 / pars[2]) * sum(log(zz))
  } else {
    i <- 1:4
    g_fun <- function(x) {
      t1 <- x ^ i
      t2 <- (i * x - i - 1)
      sum((-1) ^ i * t1 * t2 * pars[2] ^ i / i / (i + 1))
    }
    val <- -m * log(pars[1]) - sum_gp / pars[1] - sum(sapply(sdat, g_fun))
  }
  return(val)
}

# =========================== gpd_mle ===========================

#' @keywords internal
#' @rdname rust-internal
gpd_mle <- function(gpd_data) {
  # Maximum likelihood estimation for the generalized Pareto distribution
  #
  # Performs maximum likelihood estimation for the generalized Pareto
  # distribution.  Uses the function \code{gpdmle} associated with
  # Grimshaw (1993), which returns MLEs of sigma and k = - \xi.
  #
  # Grimshaw, S. D. (1993) Computing Maximum Likelihood Estimates
  #   for the Generalized Pareto Distribution.  Technometrics, 35(2), 185-191.
  #   and Computing (1991) 1, 129-133. http://dx.doi.org/10.1007/BF01889987.
  #
  # Args:
  #   gpd_data : A numeric vector containing positive values, assumed to be a
  #             random sample from a generalized Pareto distribution.
  #
  # Returns:
  #   A list with components
  #     mle  : A numeric vector.  MLEs of GP parameters sigma and xi.
  #     nllh : A numeric scalar.  The negated log-likelihood at the MLE.
  #
  # Call Grimshaw (1993) function, note: k is -xi, a is sigma
  #  pjn <- grimshaw_gpdmle(gpd_data)
  temp <- list()
  # Use revdbayes function if revdbayes is available, otherwise use fallback
  if (requireNamespace("revdbayes", quietly = TRUE)) {
    pjn <- revdbayes:::grimshaw_gp_mle(gpd_data)
    temp$mle <- c(pjn$a, -pjn$k)  # mle for (sigma,xi)
  } else {
    temp <- fallback_gp_mle(init = c(mean(gpd_data), 0), gpd_data = gpd_data,
                            m = length(gpd_data), xm = max(gpd_data),
                            sum_gp = sum(gpd_data))
  }
  sc <- rep(temp$mle[1], length(gpd_data))
  xi <- temp$mle[2]
  temp$nllh <- sum(log(sc)) + sum(log(1 + xi * gpd_data / sc) * (1 / xi + 1))
  return(temp)
}

# =========================== fallback_gp_mle ===========================

#' @keywords internal
#' @rdname rust-internal
fallback_gp_mle <- function(init, ...){
  x <- stats::optim(init, gpd_loglik, ..., control = list(fnscale = -1),
                    hessian = FALSE)
  temp <- list()
  temp$mle <- x$par
  temp$nllh <- -x$value
  return(temp)
}

# =========================== gpd_obs_info ===========================

#' @keywords internal
#' @rdname rust-internal
gpd_obs_info <- function(gpd_pars, y) {
  # Observed information for the generalized Pareto distribution
  #
  # Calculates the observed information matrix for a random sample \code{y}
  # from the generalized Pareto distribution, i.e. the negated Hessian matrix of
  # the generalized Pareto log-likelihood, evaluated at \code{gpd_pars}.
  #
  # Args:
  #   gpd_pars : A numeric vector. Parameters sigma and xi of the
  #   generalized Pareto distribution.
  #   y       : A numeric vector. A sample of positive data values.
  #
  # Returns:
  #   A 2 by 2 numeric matrix.  The observed information matrix.
  #
  s <- gpd_pars[1]
  x <- gpd_pars[2]
  i <- matrix(NA, 2, 2)
  i[1,1] <- -sum((1 - (1 + x) * y * (2 * s + x * y) / (s + x * y) ^ 2) / s ^ 2)
  i[1,2] <- i[2,1] <- -sum(y * (1 - y / s) / (1 + x * y / s) ^ 2 / s ^ 2)
  i[2,2] <- sum(2 * log(1 + x * y / s) / x ^ 3 - 2 * y / (s + x * y) / x ^ 2 -
                  (1 + 1 / x) * y ^ 2 / (s + x * y) ^ 2)
  return(i)
}
