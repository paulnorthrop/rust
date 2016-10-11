# =========================== gpd_sums_stats ===========================

#' Generalized Pareto summary statistics
#'
#' Calculates summary statistics involved in the Generalized Pareto
#' log-likelihood.
#'
#' @param gpd_data A numeric vector containing positive values.
#' @return A list with components
#'     \item{gpd_data}{A numeric vector. The input vector with any missings
#'     removed.}
#'     \item{m}{A numeric scalar. The sample size, i.e. the number of
#'     non-missing values.}
#'     \item{xm}{A numeric scalar. The sample maximum}
#'     \item{sum_gp}{A numeric scalar. The sum of the non-missing sample
#'     values.}
#' @examples
#' \dontrun{
#' # Sample data from a GP(sigma, xi) distribution
#' gpd_data <- rgpd(m = 100, xi = 0, sigma = 1)
#' # Calculate summary statistics for use in the log-likelihood
#' ss <- gpd_sum_stats(gpd_data)
#' }
#' @seealso \code{\link{rgpd}} for simulation from a generalized Pareto
#'   distribution.
gpd_sum_stats <- function(gpd_data) {
  ss <- list()
  nas <- is.na(gpd_data)
  if (any(nas))
    warning("Missing values have been removed from the data")
  gpd_data <- gpd_data[!nas]
  ss$gpd_data <- gpd_data
  ss$m <- length(gpd_data)                  # sample size
  ss$xm <- max(gpd_data)                    # maximum threshold excess
  ss$sum_gp <- sum(gpd_data)                # sum of threshold excesses
  return(ss)
}

# =========================== gpd_logpost ===========================

#' Generalized Pareto posterior log-density
#'
#' Calculates the generalized Pareto posterior log-density based on a particular
#' prior for the generalized Pareto parameters, a Maximal Data Information
#' (MDI) prior truncated to xi >= -1 in order to produce a posterior
#' density that is proper.
#'
#' @param pars A numeric vector containing the values of the generalized Pareto
#'   parameters sigma and xi.
#' @param ss A numeric list. Summary statistics to be passed to the generalized
#'   Pareto log-likelihood.  Calculated using \code{gpd_sum_stats}
#' @return A numeric scalar. The value of the log-likelihood.
#' @seealso \code{\link{gpd_sum_stats}} to calculate summary statistics for
#'   use in \code{gpd_loglik}.
#' @seealso \code{\link{rgpd}} for simulation from a generalized Pareto
#' @references Northrop, P. J. and Attalides, N. (2016) Posterior propriety in
#' Bayesian extreme value analyses using reference priors. Statistica Sinica,
#' 26(2), 721-743, \url{http://dx.doi.org/10.5705/ss.2014.034}.
#' @examples
#' \dontrun{
#' # Sample data from a GP(sigma, xi) distribution
#' gpd_data <- rgpd(m = 100, xi = 0, sigma = 1)
#' # Calculate summary statistics for use in the log-likelihood
#' ss <- gpd_sum_stats(gpd_data)
#' # Calculate the generalized Pareto log-posterior
#' gpd_logpost(pars = c(1, 0), ss = ss)
#' }
gpd_logpost <- function(pars, ss) {
  loglik <- do.call(gpd_loglik, c(list(pars = pars), ss))
  logprior <- log_gpd_mdi_prior(pars = pars)
  return(loglik + logprior)
}

# =========================== rgpd ===========================

#' Generalized Pareto simulation
#'
#' Simulates a sample of size \code{m} from a generalized Pareto distribution.
#'
#' @param m A numeric scalar.  The size of sample required.
#' @param sigma A numeric scalar.  The generalized Pareto scale parameter.
#' @param xi A numeric scalar.  The generalized Pareto shape parameter.
#' @return A numeric vector.  A generalized Pareto sample of size \code{m}.
#' @examples
#' \dontrun{
#' # Sample data from a GP(sigma, xi) distribution
#' gpd_data <- rgpd(m = 100, xi = 0, sigma = 1)
#' }
rgpd <- function (m = 1, sigma = 1, xi = 0) {
  if (min(sigma) <= 0) {
    stop("sigma must be positive")
  }
  if (length(xi) != 1) {
    stop("xi must be scalar")
  }
  if (xi==0) {
    return(sigma * stats::rexp(m))
  } else {
    return(sigma * (stats::runif(m) ^ (-xi) - 1) / xi)
  }
}

# =========================== gpd_init ===========================

#' Initial estimates for Generalized Pareto parameters
#'
#' Calculates initial estimates and estimated standard errors (SEs) for the
#' generalized Pareto parameters sigma and xi based on an
#' assumed random sample from this distribution.  Also, calculates
#' initial estimates and estimated standard errors for
#' phi1 = sigma and phi2 = xi + sigma / xm.
#'
#' @param gpd_data A numeric vector containing positive sample values.
#' @param m A numeric scalar.  The sample size, i.e. the length of gpd_data.
#' @param xm A numeric scalar. The sample maximum.
#' @param sum_gp A numeric scalar. The sum of the sample values.
#' @param xi_eq_zero A logical scalar.  If TRUE assume that the shape
#'   parameter xi = 0.
#' @param init_ests A numeric vector.  Initial estimate of
#'   theta = (sigma, xi).  If supplied \code{gpd_init()} just
#'   returns the corresponding initial estimate of phi = (phi1, phi2).
#' @details The main aim is to calculate an admissible estimate of theta at,
#'   i.e. one at which the log-likelihood is finite (necessary for the
#'   posterior log-density to be finite) at the estimate, and associated
#'   estimated SEs. These are converted into estimates and SEs for phi.  The
#'   latter can be used to set values of \code{min_phi} and \code{max_phi}
#'   for input to \code{find_lambda} and \code{find_lambda_one_d}.
#'
#'   In the default setting (\code{xi_eq_zero = FALSE} and
#'   \code{init_ests = NULL}) the methods tried are Maximum Likelihood
#'   Estimation (MLE) (Grimshaw, 1993), Probability-Weighted Moments (PWM)
#'   (Hosking and Wallis, 1987) and Linear Combinations of Ratios of Spacings
#'   (LRS) (Riess and Thomas, 2007, page 134) in that order.
#'
#'   For xi < -1 the likelihood is unbounded, MLE may fail when xi is not
#'   greater than -0.5 and the observed Fisher information for (sigma, xi) has
#'   finite variance only if xi > -0.25.  We use the ML estimate provided that
#'   the estimate of xi returned from \code{gpd_mle} is greater than -1. We only
#'   use the SE if the MLE of xi is greater than -0.25.
#'
#'   If either the MLE or the SE are not OK then we try PWM.  We use the PWM
#'   estimate only if is admissible, and the MLE was not OK.  We use the PWM SE,
#'   but this will be \code{c(NA, NA)} is the PWM estimate of xi is > 1/2.  If
#'   the estimate is still not OK then we try LRS.  As a last resort, which
#'   will tend to occur only when xi is strongly negative, we set xi = -1 and
#'   estimate sigma conditional on this.
#' @return If \code{init_ests} is not supplied by the user, a list is returned
#'   with components
#'     \item{init}{A numeric vector. Initial estimates of sigma
#'      and xi.}
#'     \item{se}{A numeric vector. Estimated standard errors of
#'      sigma and xi.}
#'     \item{init_phi}{A numeric vector. Initial estimates of
#'      phi1 = sigma and phi2 = xi + sigma / xm,
#'      where xm is the maximum of \code{gpd_data}.}
#'     \item{se_phi}{A numeric vector. Estimated standard errors of
#'      phi1 and phi2.}
#'   If \code{init_ests} is supplied then only the numeric vector
#'   \code{init_phi} is returned.
#' @references Grimshaw, S. D. (1993) Computing Maximum Likelihood Estimates
#'   for the Generalized Pareto Distribution.  Technometrics, 35(2), 185-191.
#'   and Computing (1991) 1, 129-133. \url{http://dx.doi.org/10.1007/BF01889987}.
#' @references Hosking, J. R. M. and Wallis, J. R. (1987) Parameter and Quantile
#'   Estimation for the Generalized Pareto Distribution. Technometrics, 29(3),
#'   339-349. \url{http://dx.doi.org/10.2307/1269343}.
#' @references Reiss, R.-D., Thomas, M. (2007) Statistical Analysis of Extreme Values
#'   with Applications to Insurance, Finance, Hydrology and Other Fields.Birkhauser.
#'   \url{http://dx.doi.org/10.1007/978-3-7643-7399-3}.
#' @seealso \code{\link{gpd_sum_stats}} to calculate summary statistics for
#'   use in \code{gpd_loglik}.
#' @seealso \code{\link{rgpd}} for simulation from a generalized Pareto
#' @seealso \code{\link{find_lambda_one_d}} to produce (somewhat) automatically
#'   a list for the argument \code{lambda} of \code{ru} for the
#'   \code{d} = 1 case.
#' @seealso \code{\link{find_lambda}} to produce (somewhat) automatically
#'   a list for the argument \code{lambda} of \code{ru} for any value of
#'   \code{d}.
#' @examples
#' \dontrun{
#' # Sample data from a GP(sigma, xi) distribution
#' gpd_data <- rgpd(m = 100, xi = 0, sigma = 1)
#' # Calculate summary statistics for use in the log-likelihood
#' ss <- gpd_sum_stats(gpd_data)
#' # Calculate initial estimates
#' do.call(gpd_init, ss)
#' }
gpd_init <- function(gpd_data, m, xm, sum_gp = NULL, xi_eq_zero = FALSE,
                    init_ests = NULL) {
  #
  theta_to_phi <- function(theta) c(theta[1], theta[2] + theta[1] / xm)
  #
  if (!is.null(init_ests)) {
    return(theta_to_phi(init_ests))
  }
  if (xi_eq_zero) {
    s_hat <- mean(gpd_data)
    init <- c(s_hat, 0)
    init_phi <- theta_to_phi(init)
    cov_mtx <- matrix(c(2*s_hat^2, -s_hat, -s_hat, 1), 2, 2) / m
    se <- sqrt(diag(cov_mtx))
    mat <- matrix(c(1, 0, 1 / xm, 1), 2, 2, byrow = TRUE)
    var_phi <- mat %*% cov_mtx %*% t(mat)
    se_phi <- sqrt(diag(var_phi))
    return(list(init = init, se = se, init_phi = init_phi, se_phi = se_phi))
  }
  #
  ests_ok <- ses_ok <- FALSE
  # First try MLE
  temp <- gpd_mle(gpd_data)
  init <- temp$mle
  # Check that the MLE is OK
  if (!is.na(init[1]) & init[2] > -1 & init[2] > -init[1]/xm &
      !is.infinite(temp$nllh)){
    ests_ok <- TRUE
    # Check whether or not we should use the SE
    if (init[2] > -0.25){
      cov_mtx <- solve(gpd_obs_info(init, gpd_data))
      se <- sqrt(diag(cov_mtx))
      mat <- matrix(c(1, 0, 1 / xm, 1), 2, 2, byrow = TRUE)
      var_phi <- mat %*% cov_mtx %*% t(mat)
      if (all(diag(var_phi) > 0)){
        se_phi <- sqrt(diag(var_phi))
        ses_ok <- TRUE
      }
    }
  }
  if (!ests_ok | !ses_ok){
    # Try PWM
    pwm <- gpd_pwm(gpd_data)
    se <- pwm$se
    mat <- matrix(c(1, 0, 1 / xm, 1), 2, 2, byrow = TRUE)
    var_phi <- mat %*% pwm$cov %*% t(mat)
    se_phi <- sqrt(diag(var_phi))
    # Note: se and se_phi will be NA if pwm$est > 1/2
    check <- gpd_loglik(pars = pwm$est, gpd_data = gpd_data, m = m, xm = xm,
                       sum_gp = sum_gp)
    # If MLE wasn't OK and PWM estimate is OK then use PWM estimate
    if (!ests_ok & init[2] > -1 & !is.infinite(check)) {
      init <- pwm$est
      ests_ok <- TRUE
    }
  }
  # If estimate is not OK then try LRS
  if (!ests_ok){
    init <- gpd_lrs(gpd_data)
    check <- gpd_loglik(pars = init, gpd_data = gpd_data, m = m, xm = xm,
                       sum_gp = sum_gp)
    if (init[2] > -1 & !is.infinite(check)) {
      ests_ok <- TRUE
    }
  }
  # If we get here without ests_ok = TRUE then the posterior mode is
  # probably close to xi = -1.  We want to avoid xi < -1 so we set
  # xi = -1 and use the (bias-adjusted) mle of sigma conditional on xi = -1.
  if (!ests_ok) {
    init <- c(xm * (m + 1) / m, -1)
  }
  init_phi <- theta_to_phi(init)
  return(list(init = init, se = se, init_phi = init_phi, se_phi = se_phi))
}

# =========================== log_gpd_mdi_prior ===========================

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

gpd_loglik <- function(pars, gpd_data, m, xm, sum_gp) {
  # Generalized Pareto log-likelihood
  #
  # Calculates the log-likelihood for a random sample from a Generalized Pareto.
  #
  # Args:
  #   pars    : A numeric vector containing the values of the generalized
  #             Pareto parameters sigma and xi.
  #   gpd_data : A numeric vector containing positive sample values.
  #   m       : A numeric scalar.  The sample size, i.e. the length of gpd_data.
  #   xm      : A numeric scalar. The sample maximum.
  #   sum_gp  : The sum of the sample values.
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
  pjn <- grimshaw_gpdmle(gpd_data)
  temp <- list()
  temp$mle <- c(pjn$a, -pjn$k)  # mle for (sigma,xi)
  sc <- rep(temp$mle[1], length(gpd_data))
  xi <- temp$mle[2]
  temp$nllh <- sum(log(sc)) + sum(log(1 + xi * gpd_data / sc) * (1 / xi + 1))
  return(temp)
}

# =========================== gpd_pwm ===========================

gpd_pwm <- function(gpd_data, u = 0) {
  # Probability weighted moments estimation for the generalized Pareto
  # distribution.
  #
  # Args:
  #   gpd_data : A numeric vector containing positive values, assumed to be a
  #             random sample from a generalized Pareto distribution.
  #   u       : A numeric scalar.  A threshold.  The GP distribution is
  #             fitted to the excesses of u.
  # Returns:
  #   A list with components
  #     est  : A numeric vector.  PWM estimates of GP parameters sigma and xi.
  #     se   : A numeric vector.  Estimated standard errors of sigma and xi.
  #    cov   : A numeric matrix.  Estimate covariance matrix of the the PWM
  #            estimators of sigma and xi.
  #
  n <- length(gpd_data)
  exceedances <- gpd_data[gpd_data > u]
  excess <- exceedances - u
  nu <- length(excess)
  xbar <- mean(excess)
  a0 <- xbar
  gamma <- -0.35
  delta <- 0
  pvec <- ((1:nu) + gamma)/(nu + delta)
  a1 <- mean(sort(excess) * (1 - pvec))
  xi <- 2 - a0/(a0 - 2 * a1)
  sigma <- (2 * a0 * a1)/(a0 - 2 * a1)
  pwm <- c(sigma, xi)
  names(pwm) = c("sigma","xi")
  denom <- nu * (1 - 2 * xi) * (3 - 2 * xi)
  if (xi > 0.5) {
    denom <- NA
    warning("Asymptotic Standard Errors not available for PWM when xi>0.5.")
  }
  one <- (7 - 18 * xi + 11 * xi^2 - 2 * xi^3) * sigma ^ 2
  two <- (1 - xi) * (1 - xi + 2 * xi^2) * (2 - xi) ^ 2
  cov <-  - sigma * (2 - xi) * (2 - 6 * xi + 7 * xi ^ 2 - 2 * xi ^ 3)
  pwm_varcov <- matrix(c(one, cov, cov, two), 2)/denom
  colnames(pwm_varcov) <- c("sigma","xi")
  rownames(pwm_varcov) <- c("sigma","xi")
  se <- sqrt(diag(pwm_varcov))
  return(list(est = pwm, se = se, cov = pwm_varcov))
}

# =========================== gpd_lrs ===========================

gpd_lrs <- function(x) {
  # LRS estimation for the generalized Pareto distribution.
  #
  # Args:
  #   x : A numeric vector containing positive values, assumed to be a
  #       random sample from a generalized Pareto distribution.
  #
  # Returns:
  #   A numeric vector.  Estimates of parameters sigma and xi.
  #
  n <- length(x)                             # sample size
  x <- sort(x)                               # put data in ascending order
  q0 <- 1 / (n + 1)
  q2 <- n / (n+1)                            # for sample minimum and maximum
  a <- sqrt((1 - q2) / (1 - q0))
  q1 <- 1 - a *(1 - q0)                      # `middle' quantile
  n0 <- 1
  n1 <- round((n + 1) * q1)
  n2 <- n                                    # corresponding order statistics
  ns <- c(n0,n1,n2); qs <- c(q0,q1,q2)
  xs <- x[ns]
  r_hat <- (xs[3] - xs[2]) / (xs[2] - xs[1])
  xi_hat <- -log(r_hat) / log(a)
  sigma_hat <- xi_hat * xs[3] / ((1 - q2) ^ (-xi_hat) - 1)
  return(c(sigma_hat, xi_hat))
}

# =========================== gpd_obs_info ===========================

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

# =========================== grimshaw_gpmle ===========================

grimshaw_gpdmle <- function(x) {
  #  Argument for function:
  #
  #  x     the sample values from the GPD
  #
  #
  #
  #  Returned from the function:
  #
  #  k       the mle of k
  #  a       the mle of a
  #
  n<-length(x)
  xq<-sort(x)
  xbar<-mean(x)
  sumx2<-sum(x^2)/n
  x1<-xq[1]
  xn<-xq[n]
  #
  #  Find the local maxima/minima of the likelihood.
  #
  #
  #  Initialize epsilon as the accuracy criterion
  #
  epsilon<-10^(-6)/xbar
  #
  #  The local maxima/minima must be found numerically by
  #  finding the zero(s) of h().
  #
  #  Algorithm for finding the zero(s) of h().
  #
  #
  #  Any roots that exist must be within the interval (lobnd,hibnd).
  #
  lobnd<-2*(x1-xbar)/x1^2
  if(lobnd>=0){
    lobnd<- -epsilon
  }
  hibnd<-(1/xn)-epsilon
  if(hibnd<=0){
    hibnd<-epsilon
  }
  #
  #  If h''(0)>0, look for one negative and one positive zero of h().
  #  If h''(0)<0, look for two negative and two positive zeros of h().
  #
  secderiv<-sumx2-2*xbar^2  #{ Evaluate h''(0). }
  if(secderiv>0){
    #
    #
    #  Look for one negative and one positive zero of h().
    #
    #
    thzeros<-cbind(c(0,0),c(0,0))
    nzeros<-2
    #
    #  Begin with the initial value at lobnd.
    #
    hlo<-(1+sum(log(1-lobnd*x))/n)*(sum(1/(1-lobnd*x))/n)-1
    if(hlo<0){
      thlo<-lobnd       #{  Orient the search so h(thlo)<0  }
      thhi<- -epsilon
    }
    else{
      thlo<- -epsilon
      thhi<-lobnd
    }
    thzero<-lobnd    #{  Initial value for modified Newton-Raphson is lobnd. }
    dxold<-abs(thhi-thlo)
    dx<-dxold
    temp1<-sum(log(1-thzero*x))/n
    temp2<-sum(1/(1-thzero*x))/n
    temp3<-sum(1/(1-thzero*x)^2)/n
    h<-(1+temp1)*(temp2)-1
    hprime<-(temp3-temp2^2-temp1*(temp2-temp3))/thzero

    #
    #  Newton-Raphson Algorithm to find the zero of the function h()
    #  for a given initial starting point.
    j<-1
    maxiter<-100  #{Maximum number of mod. Newton-Raphson iterations}
    while(j<=maxiter){
      #
      #  Determine whether it is better to use Bisection (if N-R is
      #  out of range or not decreasing fast enough) or Newton-Raphson.
      #
      c1<-(((thzero-thhi)*hprime-h)*((thzero-thlo)*hprime-h)>=0)
      c2<-(abs(2*h)>abs(dxold*hprime))
      if(c1+c2>=1){
        dxold<-dx
        dx<-(thhi-thlo)/2
        thzero<-thlo+dx
        if(thlo==thzero){       #{Change in root is negligible}
          j<-1000
        }
      }
      else{
        dxold<-dx
        dx<-h/hprime
        temp<-thzero
        thzero<-thzero-dx
        if(temp==thzero){       #{Change in root is negligible}
          j<-1001
        }
      }
      #
      #  Determine if convergence criterion is met.
      #
      if(abs(dx)<epsilon*abs(thlo+thhi)/2){
        j<-999
      }
      temp1<-sum(log(1-thzero*x))/n
      temp2<-sum(1/(1-thzero*x))/n
      temp3<-sum(1/(1-thzero*x)^2)/n
      h<-(1+temp1)*(temp2)-1
      hprime<-(temp3-temp2^2-temp1*(temp2-temp3))/thzero
      if(h<0){            #{Maintain the bracket on the root}
        thlo<-thzero
      }
      else{
        thhi<-thzero
      }
      j<-j+1
    }

    if(j>maxiter+1){
      thzeros[1,]<-cbind(thzero,j)
    }

    #
    #  Begin with the initial value at hibnd.
    #
    hlo<-(1+sum(log(1-epsilon*x))/n)*(sum(1/(1-epsilon*x))/n)-1
    if(hlo<0){
      thlo<-epsilon       #{  Orient the search so h(thlo)<0  }
      thhi<-hibnd
    }
    else{
      thlo<-hibnd
      thhi<-epsilon
    }
    thzero<-hibnd    #{  Initial value for modified Newton-Raphson is hibnd. }
    dxold<-abs(thhi-thlo)
    dx<-dxold
    temp1<-sum(log(1-thzero*x))/n
    temp2<-sum(1/(1-thzero*x))/n
    temp3<-sum(1/(1-thzero*x)^2)/n
    h<-(1+temp1)*(temp2)-1
    hprime<-(temp3-temp2^2-temp1*(temp2-temp3))/thzero

    #
    #  Newton-Raphson Algorithm to find the zero of the function h()
    #  for a given initial starting point.
    j<-1
    maxiter<-100  #{Maximum number of mod. Newton-Raphson iterations}
    while(j<=maxiter){
      #
      #  Determine whether it is better to use Bisection (if N-R is
      #  out of range or not decreasing fast enough) or Newton-Raphson.
      #
      c1<-(((thzero-thhi)*hprime-h)*((thzero-thlo)*hprime-h)>=0)
      c2<-(abs(2*h)>abs(dxold*hprime))
      if(c1+c2>=1){
        dxold<-dx
        dx<-(thhi-thlo)/2
        thzero<-thlo+dx
        if(thlo==thzero){       #{Change in root is negligible}
          j<-1000
        }
      }
      else{
        dxold<-dx
        dx<-h/hprime
        temp<-thzero
        thzero<-thzero-dx
        if(temp==thzero){       #{Change in root is negligible}
          j<-1001
        }
      }
      #
      #  Determine if convergence criterion is met.
      #
      if(abs(dx)<epsilon*abs(thlo+thhi)/2){
        j<-999
      }
      temp1<-sum(log(1-thzero*x))/n
      temp2<-sum(1/(1-thzero*x))/n
      temp3<-sum(1/(1-thzero*x)^2)/n
      h<-(1+temp1)*(temp2)-1
      hprime<-(temp3-temp2^2-temp1*(temp2-temp3))/thzero
      if(h<0){            #{Maintain the bracket on the root}
        thlo<-thzero
      }
      else{
        thhi<-thzero
      }
      j<-j+1
    }


    if(j>maxiter+1){
      thzeros[2,]=cbind(thzero,j)
    }
  }

  else{
    #
    #
    #  Look for two negative and two positive zeros of h().
    #
    #
    thzeros<-matrix(rep(0,8),ncol=2)
    nzeros<-4
    #
    #  Begin with the initial value at lobnd.
    #
    hlo<-(1+sum(log(1-lobnd*x))/n)*(sum(1/(1-lobnd*x))/n)-1
    if(hlo<0){
      thlo<-lobnd       #{  Orient the search so h(thlo)<0  }
      thhi<- -epsilon
    }
    else{
      thlo<- -epsilon
      thhi<-lobnd
    }
    thzero<-lobnd    #{  Initial value for modified Newton-Raphson is lobnd. }
    dxold<-abs(thhi-thlo)
    dx<-dxold
    temp1<-sum(log(1-thzero*x))/n
    temp2<-sum(1/(1-thzero*x))/n
    temp3<-sum(1/(1-thzero*x)^2)/n
    h<-(1+temp1)*(temp2)-1
    hprime<-(temp3-temp2^2-temp1*(temp2-temp3))/thzero

    #
    #  Newton-Raphson Algorithm to find the zero of the function h()
    #  for a given initial starting point.
    j<-1
    maxiter<-100  #{Maximum number of mod. Newton-Raphson iterations}
    while(j<=maxiter){
      #
      #  Determine whether it is better to use Bisection (if N-R is
      #  out of range or not decreasing fast enough) or Newton-Raphson.
      #
      c1<-(((thzero-thhi)*hprime-h)*((thzero-thlo)*hprime-h)>=0)
      c2<-(abs(2*h)>abs(dxold*hprime))
      if(c1+c2>=1){
        dxold<-dx
        dx<-(thhi-thlo)/2
        thzero<-thlo+dx
        if(thlo==thzero){       #{Change in root is negligible}
          j<-1000
        }
      }
      else{
        dxold<-dx
        dx<-h/hprime
        temp<-thzero
        thzero<-thzero-dx
        if(temp==thzero){       #{Change in root is negligible}
          j<-1001
        }
      }
      #
      #  Determine if convergence criterion is met.
      #
      if(abs(dx)<epsilon*abs(thlo+thhi)/2){
        j<-999
      }
      temp1<-sum(log(1-thzero*x))/n
      temp2<-sum(1/(1-thzero*x))/n
      temp3<-sum(1/(1-thzero*x)^2)/n
      h<-(1+temp1)*(temp2)-1
      hprime<-(temp3-temp2^2-temp1*(temp2-temp3))/thzero
      if(h<0){            #{Maintain the bracket on the root}
        thlo<-thzero
      }
      else{
        thhi<-thzero
      }
      j<-j+1
    }

    if(j>maxiter+1){
      thzeros[1,]<-cbind(thzero,j)
    }
    #
    #  Look at the derivative to determine where the second root lies.
    #   If h'(0)>0, second root lies between thzero and -epsilon.
    #   If h'(0)<0, second root lies between lobnd and thzero.
    #
    temp1<-sum(log(1-thzero*x))/n
    temp2<-sum(1/(1-thzero*x))/n
    temp3<-sum(1/(1-thzero*x)^2)/n
    hprime<-(temp3-temp2^2-temp1*(temp2-temp3))/thzero
    if(hprime>0){
      #
      #  h'(0)>0, so the second zero lies between thzero and -epsilon.
      #
      #
      #  Establish Initial Values.
      #
      thlo<-thzero
      thhi<- -epsilon
      thzero<-thhi
      dx<-thlo-thhi

      j<-1
      maxiter<-100  #{Maximum number of bisection iterations}
      while(j<=maxiter){
        dx<-.5*dx
        thmid<-thzero+dx
        hmid<-(1+sum(log(1-thmid*x))/n)*(sum(1/(1-thmid*x))/n)-1
        if(hmid<0){
          thzero<-thmid
        }
        if(hmid==0){   #{Zero of h() has been found}
          j<-999
        }
        if(abs(dx)<epsilon*abs(thlo+thhi)/2){
          j<-999
        }
        j<-j+1
      }

      if(j>maxiter+1){
        thzeros[2,]<-cbind(thzero,j)
      }
    }
    else{
      #
      #  h'(0)<0, so the second zero lies between lobnd and thzero.
      #
      #
      #  Establish Initial Values.
      #
      thlo<-lobnd
      thhi<-thzero
      thzero<-thlo
      dx<-thhi-thlo

      j<-1
      maxiter<-100  #{Maximum number of bisection iterations}
      while(j<=maxiter){
        dx<-.5*dx
        thmid<-thzero+dx
        hmid<-(1+sum(log(1-thmid*x))/n)*(sum(1/(1-thmid*x))/n)-1
        if(hmid<0){
          thzero<-thmid
        }
        if(hmid==0){   #{Zero of h() has been found}
          j<-999
        }
        if(abs(dx)<epsilon*abs(thlo+thhi)/2){
          j<-999
        }
        j<-j+1
      }

      if(j>maxiter+1){
        thzeros[2,]<-cbind(thzero,j)
      }
    }
    #
    #  Begin with the initial value at hibnd.
    #
    hlo<-(1+sum(log(1-epsilon*x))/n)*(sum(1/(1-epsilon*x))/n)-1
    if(hlo<0){
      thlo<-epsilon       #{  Orient the search so h(thlo)<0  }
      thhi<-hibnd
    }
    else{
      thlo<-hibnd
      thhi<-epsilon
    }
    thzero<-hibnd    #{  Initial value for modified Newton-Raphson is hibnd. }
    dxold<-abs(thhi-thlo)
    dx<-dxold
    temp1<-sum(log(1-thzero*x))/n
    temp2<-sum(1/(1-thzero*x))/n
    temp3<-sum(1/(1-thzero*x)^2)/n
    h<-(1+temp1)*(temp2)-1
    hprime<-(temp3-temp2^2-temp1*(temp2-temp3))/thzero

    #
    #  Newton-Raphson Algorithm to find the zero of the function h()
    #  for a given initial starting point.
    j<-1
    maxiter<-100  #{Maximum number of mod. Newton-Raphson iterations}
    while(j<=maxiter){
      #
      #  Determine whether it is better to use Bisection (if N-R is
      #  out of range or not decreasing fast enough) or Newton-Raphson.
      #
      c1<-(((thzero-thhi)*hprime-h)*((thzero-thlo)*hprime-h)>=0)
      c2<-(abs(2*h)>abs(dxold*hprime))
      if(c1+c2>=1){
        dxold<-dx
        dx<-(thhi-thlo)/2
        thzero<-thlo+dx
        if(thlo==thzero){       #{Change in root is negligible}
          j<-1000
        }
      }
      else{
        dxold<-dx
        dx<-h/hprime
        temp<-thzero
        thzero<-thzero-dx
        if(temp==thzero){       #{Change in root is negligible}
          j<-1001
        }
      }
      #
      #  Determine if convergence criterion is met.
      #
      if(abs(dx)<epsilon*abs(thlo+thhi)/2){
        j<-999
      }
      temp1<-sum(log(1-thzero*x))/n
      temp2<-sum(1/(1-thzero*x))/n
      temp3<-sum(1/(1-thzero*x)^2)/n
      h<-(1+temp1)*(temp2)-1
      hprime<-(temp3-temp2^2-temp1*(temp2-temp3))/thzero
      if(h<0){            #{Maintain the bracket on the root}
        thlo<-thzero
      }
      else{
        thhi<-thzero
      }
      j<-j+1
    }

    if(j>maxiter+1){
      thzeros[3,]<-cbind(thzero,j)
    }
    #
    #  Look at the derivative to determine where the second root lies.
    #   If h'(0)>0, second root lies between thzero and hibnd.
    #   If h'(0)<0, second root lies between epsilon and thzero.
    #
    temp1<-sum(log(1-thzero*x))/n
    temp2<-sum(1/(1-thzero*x))/n
    temp3<-sum(1/(1-thzero*x)^2)/n
    hprime<-(temp3-temp2^2-temp1*(temp2-temp3))/thzero
    if(hprime>0){
      #
      #  h'(0)>0, so the second zero lies between thzero and hibnd.
      #
      #
      #  Establish Initial Values.
      #
      thlo<-thzero
      thhi<-hibnd
      thzero<-thhi
      dx<-thlo-thhi

      j<-1
      maxiter<-100  #{Maximum number of bisection iterations}
      while(j<=maxiter){
        dx<-.5*dx
        thmid<-thzero+dx
        hmid<-(1+sum(log(1-thmid*x))/n)*(sum(1/(1-thmid*x))/n)-1
        if(hmid<0){
          thzero<-thmid
        }
        if(hmid==0){   #{Zero of h() has been found}
          j<-999
        }
        if(abs(dx)<epsilon*abs(thlo+thhi)/2){
          j<-999
        }
        j<-j+1
      }

      if(j>maxiter+1){
        thzeros[4,]<-cbind(thzero,j)
      }
    }
    else{
      #
      #  h'(0)<0, so the second zero lies between epsilon and thzero.
      #
      #
      #  Establish Initial Values.
      #
      thlo<-epsilon
      thhi<-thzero
      thzero<-thlo
      dx<-thhi-thlo

      j<-1
      maxiter<-100  #{Maximum number of bisection iterations}
      while(j<=maxiter){
        dx<-.5*dx
        thmid<-thzero+dx
        hmid<-(1+sum(log(1-thmid*x))/n)*(sum(1/(1-thmid*x))/n)-1
        if(hmid<0){
          thzero<-thmid
        }
        if(hmid==0){   #{Zero of h() has been found}
          j<-999
        }
        if(abs(dx)<epsilon*abs(thlo+thhi)/2){
          j<-999
        }
        j<-j+1
      }

      if(j>maxiter+1){
        thzeros[4,]<-cbind(thzero,j)
      }
    }
  }
  #
  #  Of the candidate zero(s) of h(), determine whether they correspond
  #  to a local maximum or minimum of the log-likelihood.
  #
  #  Eliminate any non-convergent roots}
  thetas<-thzeros[thzeros[,2]>maxiter+1,]
  nzeros<-nrow(thetas)
  proll<-rep(0,nzeros)
  mles<-matrix(rep(0,4*nzeros),ncol=4)
  i<-1
  while(i<=nzeros){
    temp1<-sum(log(1-thetas[i,1]*x))
    mles[i,1]<- -temp1/n
    mles[i,2]<-mles[i,1]/thetas[i,1]
    mles[i,3]<- -n*log(mles[i,2])+(1/mles[i,1]-1)*temp1
    mles[i,4]<-999
    i<-i+1
  }
  ind<-1:length(mles[,4])
  ind<-ind[mles[,4]==999]
  if(sum(ind)==0){   #{ Check to see if there are any local maxima. }
    nomle<-0          #{ If not, there is no mle. }
  }
  else{
    nomle<-1
  }
  if(nomle!=0){
    mles<-mles[ind,]
    nmles<-nrow(mles)
    #
    #  Add the boundary value where k=1 to the candidates for the
    #  maximum of the log-likelihood.
    #
    mles<-rbind(mles,c(1,xn,-n*log(xn),999))
    nmles<-nmles+1
    #
    #  Choose of the candidate mles whichever has the largest
    #  log-likelihood.
    #
    maxlogl<-max(mles[,3])
    ind<-order(mles[,3])
    ind<-ind[nmles]
    k<-mles[ind,1]
    #  label(k,'GPD mle of k')
    a<-mles[ind,2]
    #  label(a,'GPD mle of a')
  }
  else{
    #
    #  No Maximum Likelihood Estimators were found.
    #
    k<-NA
    a<-NA
  }
  return(list(k=k,a=a))
}
