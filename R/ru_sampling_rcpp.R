# =========================== ru_rucpp ===========================

#' Generalized ratio-of-uniforms sampling using C++ via Rcpp
#'
#' Uses the generalized ratio-of-uniforms method to simulate from a
#' distribution with log-density \eqn{log f} (up to an additive constant).
#' \eqn{f} must be bounded, perhaps after a transformation of variable.
#'
#' @param logf A pointer to the compiled C++ function returning the log
#'   of the target density \eqn{f}.
#' @param ... Further arguments to be passed to \code{logf} and related
#'   functions.
#' @param n A numeric scalar.  Number of simulated values required.
#' @param d A numeric scalar. Dimension of f.
#' @param init A numeric vector. Initial estimates of the mode of \code{logf}.
#'   If \code{trans="BC"} or \code{trans = "user"} this is \emph{after} Box-Cox
#'   transformation or user-defined transformation, but BEFORE any rotation
#'   of axes.
#' @param trans A character scalar. "none" for no transformation, "BC" for
#'   Box-Cox transformation, "user" for a user-defined transformation.
#'   If \code{trans = "user"} then the transformation should be specified
#'   using \code{phi_to_theta} and \code{log_j} and \code{user_args} may be
#'   used to pass arguments to \code{phi_to_theta} and \code{log_j}.
#' @param phi_to_theta A function returning (inverse) of the transformation
#'   from theta to phi used to ensure positivity of phi prior to Box-Cox
#'   transformation.  The argument is phi and the returned value is theta.
#' @param log_j A function returning the log of the Jacobian of the
#'  transformation from theta to phi, i.e. based on derivatives of phi with
#'  respect to theta. Takes theta as its argument.
#' @param user_args A list of numeric components. If \code{trans = ``user''}
#'   then \code{user_args} is a list providing arguments to the user-supplied
#'   functions \code{phi_to_theta} and \code{log_j}.  If \code{phi_to_theta}
#'   is undefined at the input value then the  function should return NA.
#' @param lambda Either
#' \itemize{
#'   \item {A numeric vector.  Box-Cox transformaton parameters, or}
#'   \item {A list with components}
#'   \describe{
#'     \item{lambda}{A numeric vector.  Box-Cox parameters (required).}
#'     \item{gm}{A numeric vector.  Box-cox scaling parameters (optional).
#'       If supplied this overrides any \code{gm} supplied by the individual
#'       \code{gm} argument described below.}
#'     \item{init_psi}{A numeric vector.  Initial estimate of mode \emph{after}
#'       Box-Cox transformation (optional).}
#'     \item{sd_psi}{A numeric vector.  Estimates of the marginal standard
#'       deviations of the Box-Cox transformed variables (optional).}
#'     \item{phi_to_theta}{as above (optional).}
#'     \item{log_j}{as above (optional).}
#'   }
#'   This list may be created using \link{find_lambda_one_d} (for \code{d} = 1)
#'   or \link{find_lambda} (for any \code{d}).
#' }
#' @param lambda_tol A numeric scalar.  Any values in lambda that are less
#'  than lambda_tol in magnitude are set to zero.
#' @param gm A numeric vector. Box-cox scaling parameters (optional). If
#'   \code{lambda$gm} is supplied in input list \code{lambda} then
#'   \code{lambda$gm} is used, not \code{gm}.
#' @param rotate A logical scalar. If TRUE (\code{d} > 1 only) use Choleski
#'   rotation.  If d = 1 and \code{rotate} = TRUE then rotate will be set to
#'   FALSE with a warning.
#' @param lower,upper Numeric vectors.  Lower/upper bounds on the arguments of
#'   the function \emph{after} any transformation from theta to phi implied by
#'   the inverse of \code{phi_to_theta}. If \code{rotate = FALSE} these
#'   are used in the optimizations used to construct the bounding box.  If
#'   \code{trans = "BC"} components of lower that are negative are set to zero
#'   without warning and the bounds implied after the Box-Cox transformation
#'   are calculated inside \code{ru}.  If \code{rotate = TRUE} all
#'   optimizations are unconstrained.
#' @param r A numeric scalar.  Parameter of generalized ratio-of-uniforms.
#' @param ep A numeric scalar.  Controls initial estimates for optimizations
#'   to find the b-bounding box parameters.  The default (\code{ep}=0)
#'   corresponds to starting at the mode of \code{logf} small positive values
#'   of \code{ep} move the constrained variable slightly away from the mode in
#'   the correct direction.  If \code{ep} is negative its absolute value is
#'   used, with no warning given.
#' @param a_algor,b_algor Character scalars.  Either "nlminb" or "optim".
#'   Respective optimization algorithms used to find a(r) and (bi-(r), bi+(r)).
#' @param a_method,b_method Character scalars.  Respective methods used by
#'   \code{optim} to find a(r) and (bi-(r), bi+(r)).  Only used if \code{optim}
#'   if the chosen algorithm.  If \code{d} = 1 then a_method and b_method are
#'   set to "Brent" without warning.
#' @param a_control,b_control  Lists of control arguments to \code{optim} or
#'   \code{nlminb} to find a(r) and (bi-(r), bi+(r)) respectively.
#' @param var_names A character vector.  Names to give to the column(s) of
#'   the simulated values.
#' @details If \code{trans = "none"} and \code{rotate = FALSE} then \code{rou}
#'   implements the (multivariate) generalized ratio of uniforms method
#'   described in Wakefield, Gelfand and Smith (1991) using a target
#'   density whose mode is relocated to the origin (`mode relocation') in the
#'   hope of increasing efficiency.
#'
#'   If \code{trans = "BC"} then marginal Box-Cox transformations of each of
#'   the \code{d} variables is performed, with parameters supplied in
#'   \code{lambda}.  The function \code{phi_to_theta} may be used, if
#'   necessary, to ensure positivity of the variables prior to Box-Cox
#'   transformation.
#'
#'   If \code{trans = "user"} then the function \code{phi_to_theta} enables
#'   the user to specify their own transformation.
#'
#'   In all cases the mode of the target function is relocated to the origin
#'   \emph{after} any user-supplied transformation and/or Box-Cox
#'   transformation.
#'
#'   If \code{d} is greater than one and \code{rotate = TRUE} then a rotation
#'   of the variable axes is performed \emph{after} mode relocation.  The
#'   rotation is based on the Choleski decomposition (see \link{chol}) of the
#'   estimated Hessian (computed using \link{optimHess} of the negated
#'   log-density after any user-supplied transformation or Box-Cox
#'   transformation.  If any of the eigenvalues of the estimated Hessian are
#'   non-positive (which may indicate that the estimated mode of \code{logf}
#'   is close to a variable boundary) then \code{rotate} is set to \code{FALSE}
#'   with a warning.  A warning is also given if this happens when
#'   \code{d} = 1.
#'
#'   The default value of the tuning parameter \code{r} is 1/2, which is
#'   likely to be close to optimal in many cases, particularly if
#'   \code{trans = "BC"}.
#'
#' See \code{vignette("rust-vignette", package = "rust")} for full details.
#'
#' @return An object of class "ru" is a list containing the following
#'   components:
#'     \item{sim_vals}{An \code{n} by \code{d} matrix of simulated values.}
#'     \item{box}{A (2 * \code{d} + 1) by \code{d} + 2 matrix of
#'       ratio-of-uniforms bounding box information, with row names indicating
#'       the box parameter.  The columns contain
#'       \describe{
#'         \item{column 1}{values of box parameters.}
#'         \item{columns 2 to (2+\code{d}-1)}{values of variables at which
#'          these box parameters are obtained.}
#'         \item{column 2+\code{d}}{convergence indicators.}
#'       }
#'       Scaling of f within \code{ru} and relocation of the
#'       mode to the origin means that the first row of \code{box} will always
#'       be \code{c(1, rep(0, d))}.
#'     }
#'     \item{pa}{A numeric scalar.  An estimate of the probability of
#'       acceptance.}
#'     \item{d}{A numeric scalar.  The dimension of \code{logf}.}
#'     \item{logf}{A function. \code{logf} function supplied by the user.}
#'     \item{logf_rho}{A function. The target function actually used in the
#'       ratio-of-uniforms algorithm.}
#'     \item{sim_vals_rho}{An \code{n} by \code{d} matrix of values simulated
#'       from the function used in the ratio-of-uniforms algorithm.}
#'     \item{logf_args}{A list of further arguments to \code{logf}.}
#'     \item{f_mode}{The estimated mode of the target density f, after any
#'       Box-Cox transformation and/or user supplied transformation, but before
#'       mode relocation.}
#' @references Wakefield, J. C., Gelfand, A. E. and Smith, A. F. M. (1991)
#'  Efficient generation of random variates via the ratio-of-uniforms method.
#'  Statistics and Computing (1991) 1, 129-133.
#'  \url{http://dx.doi.org/10.1007/BF01889987}.
#' @examples
#' Rcpp::sourceCpp("src/user_fns.cpp")
#' n <- 1000
#'
#' # Standard normal density ===================
#' x <- ru_rcpp(logf = ptr_N01, d = 1, n = 1000, init = 0.1)
#'
#' # Gamma (alpha, 1) density ===================
#' x <- ru_rcpp(logf = ptr_gam, alpha = alpha, d = 1, n = 1000,
#'   lower = 0, init = alpha)
#'
#' # two-dimensional normal with positive association ===================
#' x <- ru_rcpp(logf = ptr_bvn, d = 2, n = 1000, init = c(0.1, 0.1))
#'
#' # GP posterior density ===================
#' # Sample data from a GP(sigma, xi) distribution
#' gpd_data <- rgpd(m = 100, xi = -0.5, sigma = 1)
#' # Calculate summary statistics for use in the log-likelihood
#' ss <- gpd_sum_stats(gpd_data)
#' # Calculate an initial estimate
#' init <- c(mean(gpd_data), 0)
#'
#' for_ru_rcpp <- c(list(logf = ptr_gp, init = init, d = 2, n = n,
#'   lower = c(0, -Inf)), ss)
#' x <- do.call(ru_rcpp, for_ru_rcpp)
#'
#' library(microbenchmark)
#' microbenchmark(
#'   ru = ru(logf = gpd_logpost, ss = ss, d = 2, n = n, init = init,
#'     lower = c(0, -Inf)),
#'   ru_rcpp = do.call(ru_rcpp, for_ru_rcpp)
#' )
#'
#' @seealso \code{\link{summary.ru}} for summaries of the simulated values
#'   and properties of the ratio-of-uniforms algorithm.
#' @seealso \code{\link{plot.ru}} for a diagnostic plot (for \code{d} = 1
#'   and \code{d} = 2 only).
#' @seealso \code{\link{find_lambda_one_d}} to produce (somewhat) automatically
#'   a list for the argument \code{lambda} of \code{ru} for the
#'   \code{d} = 1 case.
#' @seealso \code{\link{find_lambda}} to produce (somewhat) automatically
#'   a list for the argument \code{lambda} of \code{ru} for any value of
#'   \code{d}.
#' @seealso \code{\link[stats]{optim}} for choices of the arguments
#'   \code{a_method}, \code{b_method}, \code{a_control} and \code{b_control}.
#' @seealso \code{\link[stats]{nlminb}} for choices of the arguments
#'   \code{a_control} and \code{b_control}.
#' @seealso \code{\link[stats]{optimHess}} for Hessian estimation.
#' @seealso \code{\link[base]{chol}} for the Choleski decomposition.
#'
#' @export
ru_rcpp <- function(logf, ..., n = 1, d = 1, init = NULL,
               trans = c("none", "BC", "user"),  phi_to_theta = NULL,
               log_j = NULL, user_args = list(), lambda = rep(1L, d),
               lambda_tol = 1e-6, gm = NULL,
               rotate = ifelse(d == 1, FALSE, TRUE),
               lower = rep(-Inf, d),
               upper = rep(Inf, d), r = 1 / 2, ep = 0L,
               a_algor = if (d == 1) "nlminb" else "optim",
               b_algor = c("nlminb", "optim"),
               a_method = c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN",
                            "Brent"),
               b_method = c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN",
                            "Brent"),
               a_control = list(), b_control = list(), var_names = NULL) {
  #
  # Extract list of parameters for logf
  #
  pars <- list(...)
  #
  # Check that the values of key arguments are suitable
  if (r < 0) {
    stop("r must be non-negative")
  }
  #
  a_algor <- match.arg(a_algor, c("optim", "nlminb"))
  a_method <- match.arg(a_method)
  b_algor <- match.arg(b_algor)
  b_method <- match.arg(b_method)
  if (any(upper <= lower)) {
    stop("upper must be greater than lower, componentwise.")
  }
  #
  trans <- match.arg(trans)
  # If Box-Cox scale parameter is not supplied (directly) set it to 1,
  # at least for the moment
  if (is.null(gm)) {
    gm <- rep(1, d)
  }
  # Set up Box-Cox transformation parameters (if necessary)
  if (trans == "BC") {
    lambda_type <- "numeric"
    # If lambda is a list then extract information from it
    if (is.list(lambda)) {
      lambda_type <- "list"
      if (is.null(lambda$lambda)) {
        stop("The list lambda must contain the object lambda$lambda")
      }
      if (!is.null(lambda$gm)) {
        gm <- lambda$gm
      }
      if (!is.null(lambda$init_psi)) {
        init <- lambda$init_psi
      }
      if (a_algor == "optim" & is.null(a_control$parscale)) {
        a_control <- c(a_control, list(parscale = lambda$sd_psi))
      }
      if (!is.null(lambda$phi_to_theta)) {
        phi_to_theta <- lambda$phi_to_theta
      }
      if (!is.null(lambda$log_j)) {
        log_j <- lambda$log_j
      }
      lambda <- lambda$lambda
    }
    # Set to zero values of lambda that are close to zero.
    lambda <- ifelse(abs(lambda) < lambda_tol, 0L, lambda)
    # Check that lambda is a vector of the correct length.
    if (!is.vector(lambda)) {
      stop("lambda must be a numeric vector")
    }
    if (!(length(lambda) %in% c(1, d))) {
      if (lambda_type == "numeric") {
        stop("lambda must be a numeric vector of length d")
      }
      if (lambda_type == "list") {
        stop("lambda$lambda must be a numeric vector of length d")
      }
    }
    if (length(lambda) == 1) {
      lambda <- rep(lambda, d)
    }
    # If rotate = FALSE, adjust lower and upper for the value of lambda
    if (!rotate) {
      # Check that all components of upper are positive
      if (any(upper <= 0)) {
        stop("when trans = ``BC'' all elements of upper must be positive")
      }
      # If any components of lower or upper are negative then set them to zero.
      lower <- pmax(0, lower)
      lower <- ifelse(lambda == 0, gm * log(lower),
                      (lower^lambda - 1) / (lambda * gm ^ (lambda -1)))
      upper <- ifelse(lambda == 0, gm * log(upper),
                      (upper^lambda - 1) / (lambda * gm ^ (lambda -1)))
    }
  }
  # If rotate = TRUE then don't use impose any (finite) bounds
  if (rotate) {
    lower <- rep(-Inf, d)
    upper <- rep(Inf, d)
  }
  # Check that the optimization algorithm is appropriate given the bounds in
  # lower and upper.  If not then change it, with a warning.
  if (d == 1 & a_algor == "optim" & any(is.infinite(c(lower,upper)))) {
    a_algor = "nlminb"
    warning("For d = 1 finite lower and upper bounds must be supplied when
            using a_algor = `optim'.  a_algor has been changed to `nlminb'")
  }
  if (d == 1 & b_algor == "optim" & any(is.infinite(c(lower,upper)))) {
    b_algor = "nlminb"
    warning("For d = 1 finite lower and upper bounds must be supplied when
            using b_algor = `optim'.  b_algor has been changed to `nlminb'")
  }
  if (b_algor == "optim") {
    if (b_method == "BFGS" | b_method == "CG") {
      warning("Using optim with b_method==`BFGS' or `CG' can produce the error
              message `non-finite finite-difference value'.  If you really want
              to use BFGS or CG try setting ep to be positive but small, e.g.
              ep=0.001.", immediate. = TRUE, noBreaks. = TRUE)
    }
  }
  # If d = 1 then set a_method and b_method to "Brent", just in case optim is
  # being used.
  if (d == 1) {
    a_method <- "Brent"
    b_method <- "Brent"
  }
  # Determine which values in lambda are not equal to 1.
  if (d == 1) {
    which_lam <- 1L
  } else {
    which_lam <- which(lambda != 1L)
  }
  # If no initial estimates have been supplied then use a vector of ones.
  if (is.null(init)) {
    init <- rep(1, d)
    warning("No initial estimate of the mode given: a vector of ones has
            been used", noBreaks. = TRUE)
  }
  len_init <- length(init)
  if (len_init == 1 & d > 1) {
    init <- rep(init, length.out = d)
    warning("d > 1 but init has length 1: a d-vector of inits has been used")
  }
  if (len_init != d & len_init != 1) {
      stop("the length of init is incompatible with d")
  }
  # The option rotate = TRUE is not relevant when d = 1
  if (d == 1 & rotate) {
    rotate <- FALSE
    warning("rotation is not relevant when d=1: no rotation is used")
  }
  ep <- abs(ep)
  # Objects to store the values underlying the ru box and convergence
  # indicators.
  #   vals: will contain the values of the variables at which
  #         the ru box dimensions occur.
  #   conv: will contain the corresponding covergence indicators returned by
  #         the optimisation algorithms.
  vals <- matrix(NA, ncol = d, nrow = 2 * d + 1)
  colnames(vals) <- paste("vals", 1:d, sep="")
  conv <- rep(NA, 2 * d + 1)
  # Large value to return when parameters are out of bounds
  big_val <- Inf
  # f will be scaled by exp(hscale).  ... but no scaling yet_
  hscale <- 0
  # Mode of (transformed) target density: set to zero initially
  psi_mode <- rep(0, d)
  if (trans == "none" & is.function(phi_to_theta)) {
    warning("phi_to_theta() not used when trans = ``none'': identity fn used")
  }
  if (!is.function(phi_to_theta) & !is.null(phi_to_theta)) {
    stop("phi_to_theta must be a function or NULL")
  }
  if (is.function(phi_to_theta) & is.null(log_j)) {
    log_j <- function(x) 0
    warning("No Jacobian for phi_to_theta(): constant Jacobian has been used")
  }
  if (is.null(phi_to_theta)) {
    phi_to_theta <- identity
    log_j <- function(x) 0
  }
  # variable rotation matrix: set to matrix of ones initially
  rot_mat <- diag(d)
#  if (trans == "BC") {
#    const <- lambda * gm ^ (lambda - 1)
#    if (any(lambda == 0L)) {
#      psi_to_phi <- function(psi) {
#        ifelse(lambda == 0, exp(psi / gm), (psi * const + 1) ^ (1 / lambda))
#      }
#    } else {
#      psi_to_phi <- function(psi) (psi * const + 1) ^ (1 / lambda)
#    }
#  }
  init_psi <- init
  #
#  if (trans == "none") {
#    rho_to_theta <- function(rho) rho_to_psi(rho)
#  }
#  if (trans == "BC") {
#    rho_to_theta <- function(rho) phi_to_theta(psi_to_phi(rho_to_psi(rho)))
#  }
#  if (trans == "user") {
#    rho_to_theta <- function(rho) phi_to_theta(rho_to_psi(rho))
#  }
  #
#  if (trans == "none") {
#    logf_rho <- function(rho,...) {
#      theta <- rho_to_psi(rho)
#      val <- logf(theta, ...) - hscale
#      structure(val, theta = theta)
#    }
#  }
#  if (trans == "BC") {
#    logf_rho <- function(rho, ...) {
#      psi <- rho_to_psi(rho)
#      # When lambda is not equal to one (psi * const + 1) must be positive
#      test <- (psi * const + 1)[which_lam]
#      if (any(test <= 0)) return(-Inf)
#      phi <- psi_to_phi(psi)
#      theta <- phi_to_theta(phi)
#      log_bc_jac <- sum((lambda - 1)[which_lam] * log(phi[which_lam]))
#      val <- logf(theta, ...) - log_bc_jac - log_j(theta) - hscale
#      structure(val, theta = theta)
#    }
#  }
  if (trans == "user") {
    logf_rho <- function(rho, ...) {
      phi <- rho_to_psi(rho)
      theta <- do.call(phi_to_theta, c(phi, user_args))
      if (!is.finite(theta)) return(-Inf)
      logj <- do.call(log_j, c(theta, user_args))
      val <- logf(theta, ...) - logj - hscale
      structure(val, theta = theta)
    }
  }
  #
  # Scale logf ---------------------------------
  #
  hscale <- cpp_logf_rho(init_psi, rep(0, d), diag(d), 0, logf,
                         pars = pars)
  if (is.infinite(hscale)) {
      stop("posterior density is zero at initial parameter values")
  }
  #
  # Calculate a(r) ----------------------------------
  # Create list of arguments for find_a()
  for_find_a <- list(logf = logf, init_psi = init_psi, d = d,
                     r = r, lower = lower, upper = upper, algor = a_algor,
                     method = a_method, control = a_control,
                     hscale = hscale, pars = pars)
  temp <- do.call("cpp_find_a", for_find_a)
  #
  # Check that logf is finite at 0
  #
  check_finite <- cpp_logf_rho(temp$par, rep(0, d), diag(d), hscale, logf,
                               pars = pars)
  if (!is.finite(check_finite)) {
    stop(paste("The target log-density is not finite at its mode: mode = ",
               temp$par, ", function value = ", check_finite, ".", sep=""))
  }
  #
  # Scale logf to have a maximum at 0, i.e. a=1 ------------
  #
  hscale <- check_finite + hscale
  a_box <- 1
  f_mode <- temp$par
  vals[1, ] <- rep(0, d)
  conv[1] <- temp$convergence
  pos_def <- TRUE
  if (class(temp$hessian) == "try-error") {
    pos_def <- FALSE
  } else {
    hess_mat <- temp$hessian
    e_vals <- eigen(hess_mat, symmetric = TRUE, only.values = TRUE)$values
    if (any(e_vals < 1e-6)) {
      pos_def <- FALSE
    }
  }
  # We check the eigenvalues of the estimated Hessian, hess_mat, of -logf at
  # its minimum.  This has two purposes:
  #
  #   1. If rotate = TRUE the transformation uses the Cholesky decomposition
  #      of hess_mat so hess_mat must be symmetric and positive definite, with
  #      all eigenvalues positive.
  #   2. Even if rotate = FALSE it is worth checking the positive-definiteness
  #      of hess_mat.  Lack of positive-definiteness of hess_mat could result
  #      from the mode, f_mode, of logf could being at or near a variable
  #      boundary.  This may be fine if logf is finite at the boundary, i.e.
  #      we have found the maximum of logf, but it is possible that logf is
  #      unbounded.
  #
  # In both cases, i.e. regardless of rotate, a warning is given.
  if (!pos_def) {
    warning("The Hessian of the target log-density at its mode is not positive
            definite. This may not be a problem, but it may be that a mode
            at/near a parameter boundary has been found and/or that the target
            function is unbounded.", immediate. = TRUE, noBreaks. = TRUE)
    if (trans != "BC") {
      cat("  It might be worth using the option trans = ``BC''.", "\n")
    }
    if (rotate) {
      rotate <- FALSE
      warning("rotate has been changed to FALSE.", immediate. = TRUE)
    }
  }
  if (rotate) {
    rot_mat <- solve(t(chol(hess_mat)))
    # We standardize so that the determinant of rot_mat is 1, i.e. the
    # transformation rotates but doesn't scale.  To do this we divide
    # by det(rot_mat)^(1/d).  The determinant is the product of the
    # eigenvalues of rot_mat.  These eigenvalues are the eigenvalues of
    # hess_mat in e_vals, raised to the power -1/2.
    rot_mat <- rot_mat / exp(-mean(log(e_vals)) / 2)
  }
  # In the C++ function cpp_rho_to_psi() the vectors are column vectors
  # so we need to transpose rot_mat.
  rot_mat <- t(rot_mat)
  psi_mode <- f_mode
  # Calculate biminus(r) and biplus(r), i = 1, ...d -----------
  # Create list of arguments for find_bs()
  for_find_bs <- list(logf = logf, d = d, r = r, lower = lower,
                      upper = upper, f_mode = f_mode, ep = ep, vals = vals,
                      conv = conv, algor = b_algor, method = b_method,
                      control = b_control, psi_mode = psi_mode,
                      rot_mat = rot_mat, hscale = hscale, pars = pars)
  temp <- do.call("cpp_find_bs", for_find_bs)
  vals <- temp$vals
  conv <- temp$conv
  l_box <- temp$l_box
  u_box <- temp$u_box
  #
  # Perform ratio-of-uniforms rejection samping ---------------
  # Call C++ function to do this.
  res <- ru_cpp(n = n, d = d, r = r, a_box = a_box, l_box = l_box,
                u_box = u_box, logf = logf, psi_mode= psi_mode,
                rot_mat = rot_mat, hscale = hscale, pars = pars)
  res$pa <- n / res$ntry
  box <- c(a_box, l_box, u_box)
  res$box <- cbind(box, vals, conv)
  bs <- paste(paste("b", 1:d, sep=""),rep(c("minus", "plus"), each=d), sep="")
  rownames(res$box) <- c("a", bs)
  if (any(conv != 0)) {
    warning("One or more convergence indicators are non-zero.",
            immediate.=TRUE)
    print(res$box)
  }
  res$d <- d
  res$logf <- logf
  res$logf <- cpp_logf
  res$logf_args <- list(pars = pars)
  res$logf_rho <- cpp_logf_rho
  res$cpp_logf_rho_args <- list(psi_mode = psi_mode, rot_mat = rot_mat,
                            hscale = hscale)
  res$f_mode <- f_mode
  class(res) <- "ru"
  return(res)
}

# =========================== cpp_find_a ===========================


cpp_find_a <-  function(logf, init_psi, d, r, lower, upper, algor,
                        method, control, hscale, pars) {
  #
  # Finds the value of a(r).
  #
  # Args:
  #   logf     : A pointer to the (original) target log-density function.
  #   init_psi     : A numeric scalar.  Initial value of psi.
  #   d            : A numeric scalar. Dimension of f.
  #   r            : A numeric scalar. Parameter of generalized
  #                  ratio-of-uniforms.
  #   lower        : A numeric vector.  Lower bounds on the arguments of logf.
  #   upper        : A numeric vector.  Upper bounds on the arguments of logf.
  #   algor        : A character scalar.  Algorithm ("optim" or "nlminb").
  #   method       : A character scalar.  Only relvant if algorithm = "optim".
  #   control      : A numeric list.  Control arguments to algor.
  #   hscale       : A numeric scalar.  Scales the target log-density.
  #
  # Returns: a list containing
  #   the standard returns from optim or nlminb
  #   hessian: the estimated hessian of -cpp_logf_rho/(d*r+1) at its minimum.
  #
  psi_mode <- rep(0, d)
  rot_mat <- diag(d)
  big_val <- 10 ^ 10
  #
  fn_args <- list(psi_mode = psi_mode, rot_mat = rot_mat, hscale = hscale,
                  logf = logf, d = d, r = r, pars = pars)
  #
  if (algor == "optim") {
    if (method == "L-BFGS-B" | method == "Brent") {
      add_args <- list(par = init_psi, fn = cpp_a_obj, method = method,
                       control = control, lower = lower, upper = upper,
                       big_val = big_val)
      temp <- do.call(stats::optim, c(fn_args, add_args))
    } else {
      add_args <- list(par = init_psi, fn = cpp_a_obj, method = method,
                       control = control, big_val = Inf)
      temp <- do.call(stats::optim, c(fn_args, add_args))
      # Sometimes Nelder-Mead fails if the initial estimate is too good.
      # ... so avoid non-zero convergence indicator by using BFGS instead.
      if (temp$convergence == 10) {
        add_args <- list(par = temp$par, fn = cpp_a_obj, method = "BFGS",
                         control = control, big_val = Inf)
        temp <- do.call(stats::optim, c(fn_args, add_args))
      }
    }
  } else {
    add_args <- list(start = init_psi, objective = cpp_a_obj, control = control,
                     lower = lower, upper = upper, big_val = Inf)
    temp <- do.call(stats::nlminb, c(fn_args, add_args))
    # Sometimes nlminb isn't sure that it has found the minimum when in fact
    # it has.  Try to check this, and avoid a non-zero convergence indicator
    # by using optim with method="BFGS", starting from nlminb's solution.
    if (temp$convergence > 0) {
      add_args <- list(par = temp$par, fn = cpp_a_obj, hessian = FALSE,
                       method = "BFGS", control = control, big_val = Inf)
      temp <- do.call(stats::optim, c(fn_args, add_args))
    }
  }
  # Try to calculate Hessian, but avoid problems if an error is produced.
  # An error may occur if the MAP estimate is very close to a parameter
  # boundary.
  add_args <- list(par = temp$par, fn = cpp_a_obj, big_val = Inf)
  temp$hessian <- try(do.call(stats::optimHess, c(fn_args, add_args)),
                      silent = TRUE)
  return(temp)
}

# =========================== cpp_find_bs ===========================

cpp_find_bs <-  function(logf, d, r, lower, upper, f_mode, ep, vals, conv,
                         algor, method, control, psi_mode, rot_mat, hscale,
                         pars) {
  # Finds the values of b-(r) and b+(r).
  #
  # Args:
  #   logf     : A pointer to the (original) target log-density function.
  #   init_psi     : A numeric scalar.  Initial value of psi.
  #   d            : A numeric scalar. Dimension of f.
  #   r            : A numeric scalar. Parameter of generalized
  #                  ratio-of-uniforms.
  #   lower        : A numeric vector.  Lower bounds on the arguments of logf.
  #   upper        : A numeric vector.  Upper bounds on the arguments of logf.
  #   f_mode       : A numeric scalar.  The estimated mode of the target
  #                  log-density logf.
  #   ep           : A numeric scalar.  Controls initial estimates for
  #                  optimizations to find the b-bounding box parameters.
  #                  The default (ep=0) corresponds to starting at the mode of
  #                  logf small positive values of ep move the constrained
  #                  variable slightly away from the mode in the correct
  #                  direction.  If ep is negative its absolute value is used,
  #                  with no warning given.
  #   vals         : A numeric matrix.  Will contain the values of the
  #                  variables at which the ru box dimensions occur.
  #                  Row 1 already contains the values for a(r).
  #   conv         : A numeric scalar.  Will contain the covergence
  #                  indicators returned by the optimisation algorithms.
  #                  Row 1 already contains the values for a(r).
  #   algor        : A character scalar. Algorithm ("optim" or "nlminb").
  #   method       : A character scalar.  Only relvant if algorithm = "optim".
  #   control      : A numeric list. Control arguments to algor.
  #   psi_mode     : A numeric vector.  Mode of the target log-density.
  #   rot_mat      : A numeric matrix.  Rotation matrix (equal to the identity
  #                  matrix if rotate = FALSE).
  #   hscale       : A numeric scalar.  Scales the target log-density.
  #
  # Returns: a list containing
  #   l_box : A numeric vector.  Values of biminus(r), i = 1, ...d.
  #   u_box : A numeric vector.  Values of biplus(r), i = 1, ...d.
  #   vals  : as described above in Args.
  #   conv  : as described above in Args.
  #
  big_val <- 10 ^ 10
  #
  fn_args <- list(psi_mode = psi_mode, rot_mat = rot_mat, hscale = hscale,
                  logf = logf, d = d, r = r, pars = pars)
  #
  l_box <- u_box <- NULL
  zeros <- rep(0, d)
  #
  # Find biminus(r) and biplus(s), i = 1, ...,d.
  #
  for (j in 1:d) {
    #
    # Find biminus(r) ----------
    #
    rho_init <- zeros
    rho_init[j] <- -ep
    t_upper <- upper - f_mode
    t_upper[j] <- 0
    if (algor == "nlminb") {
      add_args <- list(start = rho_init, objective = cpp_lower_box,
                       upper = t_upper, lower = lower - f_mode, j = j - 1,
                       control = control, big_val = Inf)
      temp <- do.call(stats::nlminb, c(fn_args, add_args))
      l_box[j] <- temp$objective
      # Sometimes nlminb isn't sure that it has found the minimum when in fact
      # it has.  Try to check this, and avoid a non-zero convergence indicator
      # by using optim with method="BFGS", starting from nlminb's solution.
      if (temp$convergence > 0) {
        add_args <- list(par = temp$par, fn = cpp_lower_box, j = j - 1,
                         method = "BFGS", big_val = Inf)
        temp <- do.call(stats::optim, c(fn_args, add_args))
        l_box[j] <- temp$value
      }
    }
    if (algor == "optim") {
      # L-BFGS-B and Brent don't like Inf or NA
      if (method == "L-BFGS-B" | method == "Brent") {
        add_args <- list(par = rho_init, fn = cpp_lower_box, upper = t_upper,
                         lower = lower - f_mode, j = j - 1, control = control,
                         method = method, big_val = big_val)
        temp <- do.call(stats::optim, c(fn_args, add_args))
      } else {
        add_args <- list(par = rho_init, fn = cpp_lower_box, j = j - 1,
                         control = control, method = method, big_val = Inf)
        temp <- do.call(stats::optim, c(fn_args, add_args))
        # Sometimes Nelder-Mead fails if the initial estimate is too good.
        # ... so avoid non-zero convergence indicator by using BFGS instead.
        if (temp$convergence == 10)
          add_args <- list(par = temp$par, fn = cpp_lower_box, j = j - 1,
                           control = control, method = "BFGS", big_val = Inf)
          temp <- do.call(stats::optim, c(fn_args, add_args))
      }
      l_box[j] <- temp$value
    }
    vals[j+1, ] <- temp$par
    conv[j+1] <- temp$convergence
    #
    # Find biplus(r) --------------
    #
    rho_init <- zeros
    rho_init[j] <- ep
    t_lower <- lower - f_mode
    t_lower[j] <- 0
    if (algor == "nlminb") {
      add_args <- list(start = rho_init, objective = cpp_upper_box,
                       lower = t_lower, upper = upper - f_mode, j = j - 1,
                       control = control, big_val = Inf)
      temp <- do.call(stats::nlminb, c(fn_args, add_args))
      u_box[j] <- -temp$objective
      # Sometimes nlminb isn't sure that it has found the minimum when in fact
      # it has.  Try to check this, and avoid a non-zero convergence indicator
      # by using optim with method="BFGS", starting from nlminb's solution.
      if (temp$convergence > 0) {
        add_args <- list(par = temp$par, fn = cpp_upper_box, j = j - 1,
                         method = "BFGS", big_val = Inf)
        temp <- do.call(stats::optim, c(fn_args, add_args))
        u_box[j] <- temp$value
      }
    }
    if (algor == "optim") {
      # L-BFGS-B and Brent don't like Inf or NA
      if (method == "L-BFGS-B" | method == "Brent") {
        add_args <- list(par = rho_init, fn = cpp_upper_box,
                         lower = t_lower, upper = upper - f_mode, j = j - 1,
                         control = control, method = method, big_val = big_val)
        temp <- do.call(stats::optim, c(fn_args, add_args))
      } else {
        add_args <- list(par = rho_init, fn = cpp_upper_box, j = j - 1,
                         control = control, method = method, big_val = Inf)
        temp <- do.call(stats::optim, c(fn_args, add_args))
        # Sometimes Nelder-Mead fails if the initial estimate is too good.
        # ... so avoid non-zero convergence indicator by using BFGS instead.
        if (temp$convergence == 10) {
          add_args <- list(par = temp$par, fn = cpp_upper_box, j = j - 1,
                           control = control, method = "BFGS", big_val = Inf)
          temp <- do.call(stats::optim, c(fn_args, add_args))
        }
      }
      u_box[j] <- -temp$value
    }
    vals[j+d+1, ] <- temp$par
    conv[j+d+1] <- temp$convergence
  }
  return(list(l_box = l_box, u_box = u_box, vals = vals, conv = conv))
}
