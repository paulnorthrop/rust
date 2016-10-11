# =========================== plot.ru ===========================

#' Plot diagnostics for an ru object
#'
#' \code{plot} method for class "ru".  Can only be used when d = 1 or d = 2.
#'   For d = 1 a histogram of the simulated values is plotted with a the
#'   density function superimposed.  The density is normalized crudely using
#'   the trapezium rule.  For d = 2 a scatter plot of the simulated values
#'   is produced with density contours superimposed.
#'
#' @param x an object of class "ru", a result of a call to \code{ru}.
#' @param y Not used.
#' @param ... Additional arguments passed on to \code{hist}, \code{lines},
#'   \code{contour} or \code{points}.
#' @param n A number.  The meaning depends on the value of x$d.
#' \itemize{
#'   \item {For d = 1 : n + 1 is the number of abscissae in the trapezium
#'      method used to normalize the density.}
#'   \item {For d = 2 : an n by n regular grid is used to contour the density.}
#' }
#' @param prob Numeric vector. Only relevant for d = 2.  The contour lines are
#'   drawn such that the respective probabilities that the variable lies
#'   within the contour are approximately prob.
#' @param ru_scale A logical scalar.  Should we plot data and density on the
#'   scale used in the ratio-of-uniforms algorithm (TRUE) or on the original
#'   scale (FALSE)?
#'   \code{hist} and \code{plot}.
#' @details
#' Note that \code{suppressWarnings} is used to avoid potential benign warnings
#'   caused by passing unused graphical parameters to \code{hist} and
#'   \code{lines} via \code{...}.
#' @examples
#' # Log-normal density ----------------
#' x <- ru(logf = dlnorm, log = TRUE, d = 1, n = 1000, lower = 0, init = 1)
#'
#' \dontrun{
#' plot(x)
#' }
#'
#' # Improve appearance using arguments to plot() and hist()
#' \dontrun{
#' plot(x, breaks = seq(0, ceiling(max(x$sim_vals)), by = 0.25),
#'   xlim = c(0, 10))
#' }
#'
#' # Two-dimensional normal with positive association ----------------
#' rho <- 0.9
#' covmat <- matrix(c(1, rho, rho, 1), 2, 2)
#' log_dmvnorm <- function(x, mean = rep(0, d), sigma = diag(d)) {
#'   x <- matrix(x, ncol = length(x))
#'   d <- ncol(x)
#'   - 0.5 * (x - mean) %*% solve(sigma) %*% t(x - mean)
#' }
#' x <- ru(logf = log_dmvnorm, sigma = covmat, d = 2, n = 1000, init = c(0, 0))
#'
#' \dontrun{
#' plot(x)
#' }
#' @seealso \code{\link{summary.ru}} for summaries of the simulated values
#'   and properties of the ratio-of-uniforms algorithm.
#' @export
plot.ru <- function(x, y, ..., n = ifelse(x$d == 1, 1001, 101),
                     prob = c(0.1, 0.25, 0.5, 0.75, 0.95, 0.99),
                     ru_scale = FALSE) {
  if (!inherits(x, "ru")) {
    stop("use only with \"ru\" objects")
  }
  if(x$d > 2) {
    stop("plot.ru is only available for d = 1 or d = 2")
  }
  if (n < 1) {
    stop("n must be no smaller than 1")
  }
  if (ru_scale) {
    plot_data <- x$sim_vals_rho
    plot_density <- x$logf_rho
  } else {
    plot_data <- x$sim_vals
    plot_density <- x$logf
  }
  if (x$d == 1) {
    temp <- suppressWarnings(graphics::hist(plot_data, prob = TRUE,
                                            plot = FALSE))
    a <- temp$breaks[1]
    b <- temp$breaks[length(temp$breaks)]
    h <- (b-a)/n
    xx <- seq(a, b, by = h)
    for_logf <- c(list(xx), x$logf_args)
    density_fun <- function(z) {
      density_list <- c(list(z), x$logf_args)
      exp(do.call(plot_density, density_list))
    }
    yy <- sapply(xx, density_fun)
    # Remove any infinite, missing, or undefined values
    cond <- is.finite(yy)
    yy <- yy[cond]
    xx <- xx[cond]
    n <- length(yy) - 1
    #
    area <- h * (yy[1] / 2 + sum(yy[2:n]) + yy[n + 1] / 2)
    yy <- yy / area
    max_y <- max(temp$density, yy)
    suppressWarnings(graphics::hist(plot_data, prob = TRUE, main="",
                                    ylim = c(0, max_y), ...))
    suppressWarnings(graphics::lines(xx, yy, ...))
  }
  if (x$d == 2) {
    r <- apply(plot_data, 2, range)
    xx <- seq(r[1,1], r[2,1], len = n)
    yy <- seq(r[1,2], r[2,2], len = n)
    xy <- cbind(xx, yy)
    zz <- matrix(NA, ncol = length(xx), nrow = length(yy))
    for (i in 1:length(xx)) {
      for (j in 1:length(yy)) {
        for_logf <- c(list(c(xx[i], yy[j])), x$logf_args)
        zz[i, j] <- exp(do.call(plot_density, for_logf))
      }
    }
    zz[zz == -Inf] <- NA
    dx <- diff(xx[1:2]); dy <- diff(yy[1:2]); sz <- sort(zz)
    c1 <- cumsum(sz) * dx * dy; c1 <- c1/max(c1)
    con.levs <- sapply(prob, function(x) stats::approx(c1, sz, xout = 1 - x)$y)
    #
    graphics::contour(xx, yy, zz, levels = con.levs, add = F, ann = F,
      labels = prob * 100, ...)
    graphics::points(plot_data, col = 8, ...)
    graphics::contour(xx, yy, zz, levels = con.levs, add = T, ann = T,
      labels = prob * 100, ...)
  }
}

# =========================== summary.ru ===========================

#' Summarizing ratio-of-uniforms samples
#'
#' \code{summary} method for class "ru"
#'
#' @param object an object of class "ru", a result of a call to
#'   \code{ru}.
#' @param ... Additional arguments passed on to \code{print} or \code{summary}.
#' @return Prints
#' \itemize{
#'   \item {a summary of the simulated values, via
#'     \code{summary(object$sim_vals)}}
#'   \item {an estimate of the probability of acceptance, i.e.
#'     \code{object$pa}}
#'   \item {information about the ratio-of-uniforms bounding box, i.e.
#'     \code{object$box}}
#' }
#' @examples
#' # one-dimensional standard normal ----------------
#' x <- ru(logf = function(x) -x ^ 2 / 2, d = 1, n = 1000, init = 0)
#' summary(x)
#'
#' # two-dimensional normal with positive association ----------------
#' rho <- 0.9
#' covmat <- matrix(c(1, rho, rho, 1), 2, 2)
#' log_dmvnorm <- function(x, mean = rep(0, d), sigma = diag(d)) {
#'   x <- matrix(x, ncol = length(x))
#'   d <- ncol(x)
#'   - 0.5 * (x - mean) %*% solve(sigma) %*% t(x - mean)
#' }
#' x <- ru(logf = log_dmvnorm, sigma = covmat, d = 2, n = 1000, init = c(0, 0))
#' summary(x)
#' @seealso \code{\link{ru}} for descriptions of \code{object$sim_vals} and
#'   \code{object$box}.
#' @seealso \code{\link{plot.ru}} for a diagnostic plot (for \code{d} = 1
#'   and \code{d} = 2 only).
#' @export
summary.ru <- function(object, ...) {
  cat("sample summary", "\n")
  print(summary(object$sim_vals, ...), ...)
  cat("\n")
  cat("estimated probability of acceptance: ", "\n")
  print(object$pa, ...)
  cat("\n")
  cat("ru bounding box: ", "\n")
  print(object$box, ...)
}

