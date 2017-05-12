#' rust: Ratio-of-Uniforms Simulation with Transformation
#'
#' Uses the multivariate generalized ratio-of-uniforms method to simulate from a
#' distribution with log-density \code{logf} (up to an additive constant).
#' \code{logf} must be bounded, perhaps after a transformation of variable.
#'
#' @details
#'
#' The main function in the rust package is \code{ru}, which implments the
#'   generalized ratio-of-uniforms algorithm.  Also provided are two functions,
#'   \code{find_lambda} and \code{find_lambda_one_d}, that may be used to
#'   set a suitable value for the parameter \code{lambda} if Box-Cox
#'   transformation is used prior to simulation.  Basic \code{plot} and
#'   \code{summary} methods are also provided.
#'
#' See \code{vignette("rust-vignette", package = "rust")} for an overview of
#' the package.
#'
#' @references Wakefield, J. C., Gelfand, A. E. and Smith, A. F. M. Efficient
#'  generation of random variates via the ratio-of-uniforms method. Statistics
#'  and Computing (1991) 1, 129-133.
#'  \url{http://dx.doi.org/10.1007/BF01889987}.
#' @references Box, G. and Cox, D. R. (1964) An Analysis of Transformations.
#'  Journal of the Royal Statistical Society. Series B (Methodological), 26(2),
#'  211-252, \url{http://www.jstor.org/stable/2984418}.
#'
#' @seealso \code{\link{ru}} to perform ratio-of-uniforms sampling.
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
#'
#' @docType package
#' @name rust
#' @import methods
#' @importFrom stats runif
#' @useDynLib rust
#' @importFrom Rcpp sourceCpp
NULL
