#context("Bounding box")

# Checks the values in the ratio-of-uniforms bounding box returned by ru()
# for some examples in which these values can be found explicitly.
# The examples are:
#   A: d-dimensional normal density for d = 1, 2, 3.
#   B: 1-dimensional log-normal density under log transformation back to the
#      normal density
#   C: 1-dimensional gamma density, with shape parameter not less than 1

# A. d-dimensional normal density

normal_box <- function(d, sigma = diag(d), rotate = TRUE, r = 1 / 2) {
  #
  # Calculates bounding box values in the case of a zero-mean unnormalized
  # d-dimensional normal density, with arbitrary covariance structure.
  # "Unnormalized" means that the function has a maximum of 1 attained
  # at the origin, which is what is required for comparability with the
  # output from the function ru(), which scales the input function to have
  # a maximum of 1 and relocates the mode to the origin.
  #
  # Args:
  #   d      : dimension of target density.
  #   sigma  : the covariance matrix of the normal density.
  #   rotate : should the rotation transformation used in the function ru()
  #            be applied?
  #   r      : ratio-of-uniforms tuning parameter.
  #
  # Returns:
  #   box : a  (2 * d + 1) by d + 2 matrix of ratio-of-uniforms bounding box
  #         information with the same structure as object$box returned by
  #         ru().
  #
  # Make sure that sigma is a matrix
  sigma <- as.matrix(sigma)
  # If d > 1 and if rotate = TRUE calculate the covariance matrix after
  # the rotation applied in the function ru.  The inverse of the
  # covariance matrix sigma plays the role of hess_mat in the function ru.
  if (rotate & d > 1) {
    l_mat <- t(chol(solve(sigma)))
    l_mat <- l_mat / det(l_mat) ^ (1/d)
    sigma <- t(l_mat) %*% sigma %*% l_mat
  }
  # Respective values at which bminus and bplus are obtained for each of
  # the d marginal variables in the N(0, 1), independent components case.
  val <- sqrt((r * d + 1) / r)
  # Calculate box
  ar <- 1
  bplus <- val * exp(-1 /2) * sqrt(diag(sigma))
  bminus <- -bplus
  box <- c(ar, bminus, bplus)
  # Adjustment of val for the marginal variances and for the correlation
  # between the components.
  adj_mat <- t(sqrt(diag(sigma)) * stats::cov2cor(sigma))
  vals <- rbind(rep(0, d), -val * adj_mat, val * adj_mat)
  colnames(vals) <- paste("vals", 1:d, sep="")
  # We hope to get perfect convegence indicators
  conv <- rep(0, 2 * d + 1)
  # Create matrix of bounding box information
  box <- cbind(box, vals, conv)
  bs <- paste(paste("b", 1:d, sep=""),rep(c("minus", "plus"), each=d), sep="")
  rownames(box) <- c("a", bs)
  #
  return(box)
}

my_tol <- 1e-5

# 1. 1-dimensional normal

# (a) N(0, 1)

x1a <- ru(logf = function(x) -x ^ 2 / 2, d = 1, n = 1, init = 0)
test_that("N(0,1)", {
  testthat::expect_equal(x1a$box, normal_box(d = 1), tolerance = my_tol)
})
x1a_mode <- ru(logf = function(x) -x ^ 2 / 2, d = 1, n = 1, mode = 0)
test_that("N(0,1), mode supplied", {
  testthat::expect_equal(x1a_mode$box, normal_box(d = 1), tolerance = my_tol)
})

# (b) N(0, 2)

sigma <- 2 # Note: sigma means variance here
x1b <- ru(logf = function(x) -x ^ 2 / (2 * sigma), d = 1, n = 1, init = 0)
test_that("N(0,2)", {
  testthat::expect_equal(x1b$box, normal_box(d = 1, sigma = 2),
                         tolerance = my_tol)
})
x1b_mode <- ru(logf = function(x) -x ^ 2 / (2 * sigma), d = 1, n = 1, mode = 0)
test_that("N(0,2), mode supplied", {
  testthat::expect_equal(x1b_mode$box, normal_box(d = 1, sigma = 2),
                         tolerance = my_tol)
})

# (c) N(1, 2)

sigma <- 2 # Note: sigma means variance here
x1c <- ru(logf = function(x) -(x - 1) ^ 2 / (2 * sigma), d = 1, n = 1,
          init = 0)
test_that("N(1,2)", {
  testthat::expect_equal(x1c$box, normal_box(d = 1, sigma = 2),
                       tolerance = my_tol)
})
x1c_mode <- ru(logf = function(x) -(x - 1) ^ 2 / (2 * sigma), d = 1, n = 1,
               mode = 1)
test_that("N(1,2), mode supplied", {
  testthat::expect_equal(x1c_mode$box, normal_box(d = 1, sigma = 2),
                       tolerance = my_tol)
})

# 2. 2-dimensional normal

# (a) Zero mean, unit variances and independent components
#     Note: rotating shouldn't make any difference in this example.

# (i) rotate = TRUE
x2ai <- ru(logf = function(x) -(x[1]^2 + x[2]^2) / 2, d = 2, n = 1000,
         init = c(0, 0), rotate = TRUE)
test_that("BVN, rotation", {
  testthat::expect_equal(x2ai$box, normal_box(d = 2), tolerance = my_tol)
})
x2ai_mode <- ru(logf = function(x) -(x[1]^2 + x[2]^2) / 2, d = 2, n = 1000,
                mode = c(0, 0), rotate = TRUE)
test_that("BVN, rotation, mode supplied", {
  testthat::expect_equal(x2ai_mode$box, normal_box(d = 2), tolerance = my_tol)
})

# (ii) rotate = FALSE
x2aii <- ru(logf = function(x) -(x[1]^2 + x[2]^2) / 2, d = 2, n = 1000,
           init = c(0, 0), rotate = FALSE)
test_that("BVN, no rotation", {
  testthat::expect_equal(x2aii$box, normal_box(d = 2), tolerance = my_tol)
})
x2aii_mode <- ru(logf = function(x) -(x[1]^2 + x[2]^2) / 2, d = 2, n = 1000,
                 mode = c(0, 0), rotate = FALSE)
test_that("BVN, no rotation, mode supplied", {
  testthat::expect_equal(x2aii_mode$box, normal_box(d = 2), tolerance = my_tol)
})

# Function to evaluate the unnormalized log-density of a d-dimensional
# normal distribution
log_dmvnorm <- function(x, mean = rep(0, d), sigma = diag(d)) {
  x <- matrix(x, ncol = length(x))
  d <- ncol(x)
  - 0.5 * (x - mean) %*% solve(sigma) %*% t(x - mean)
}

# (b) Zero mean, unit variances with positive association
rho <- 0.9
covmat <- matrix(c(1, rho, rho, 1), 2, 2)

# (i) rotate = FALSE
x2bi <- ru(logf = log_dmvnorm, sigma = covmat, d = 2, n = 1, init = c(0, 0),
           rotate = FALSE)
test_that("BVN, rho = 0.9, no rotation", {
  testthat::expect_equal(x2bi$box, normal_box(d = 2, sigma = covmat,
                                              rotate = FALSE),
                         tolerance = my_tol)
})
x2bi_mode <- ru(logf = log_dmvnorm, sigma = covmat, d = 2, n = 1,
                mode = c(0, 0), rotate = FALSE)
test_that("BVN, rho = 0.9, no rotation, mode supplied", {
  testthat::expect_equal(x2bi_mode$box, normal_box(d = 2, sigma = covmat,
                                                   rotate = FALSE),
                         tolerance = my_tol)
})

# (ii) rotate = TRUE
x2bii <- ru(logf = log_dmvnorm, sigma = covmat, d = 2, n = 1, init = c(0, 0),
           rotate = TRUE)
test_that("BVN, rho = 0.9, rotation", {
  testthat::expect_equal(x2bii$box, normal_box(d = 2, sigma = covmat,
                                               rotate = TRUE),
                         tolerance = my_tol)
})
x2bii_mode <- ru(logf = log_dmvnorm, sigma = covmat, d = 2, n = 1,
                 mode = c(0, 0), rotate = TRUE)
test_that("BVN, rho = 0.9, rotation, mode supplied", {
  testthat::expect_equal(x2bii_mode$box, normal_box(d = 2, sigma = covmat,
                                                    rotate = TRUE),
                         tolerance = my_tol)
})

# (c) Zero mean, different variances with positive association
covmat <- matrix(c(10, 3, 3, 2), 2, 2)

# (i) rotate = FALSE
x2ci <- ru(logf = log_dmvnorm, sigma = covmat, d = 2, n = 1,
           init = c(0, 0), rotate = FALSE)
test_that("BVN, general Sigma, no rotation", {
  testthat::expect_equal(x2ci$box, normal_box(d = 2, sigma = covmat,
                                              rotate = FALSE),
                         tolerance = my_tol)
})
x2ci_mode <- ru(logf = log_dmvnorm, sigma = covmat, d = 2, n = 1,
                mode = c(0, 0), rotate = FALSE)
test_that("BVN, general Sigma, no rotation, mode supplied", {
  testthat::expect_equal(x2ci_mode$box, normal_box(d = 2, sigma = covmat,
                                                   rotate = FALSE),
                         tolerance = my_tol)
})

# (ii) rotate = TRUE
x2cii <- ru(logf = log_dmvnorm, sigma = covmat, d = 2, n = 1,
            init = c(0, 0), rotate = TRUE)
test_that("BVN, general Sigma, rotation", {
  testthat::expect_equal(x2cii$box, normal_box(d = 2, sigma = covmat,
                                               rotate = TRUE),
                         tolerance = my_tol)
})
x2cii_mode <- ru(logf = log_dmvnorm, sigma = covmat, d = 2, n = 1,
                 mode = c(0, 0), rotate = TRUE)
test_that("BVN, general Sigma, rotation, mode supplied", {
  testthat::expect_equal(x2cii_mode$box, normal_box(d = 2, sigma = covmat,
                                                    rotate = TRUE),
                         tolerance = my_tol)
})

# (d) Mean (1,2), different variances with negative association
covmat <- matrix(c(10, -3, -3, 2), 2, 2)

# (i) rotate = FALSE
x2di <- ru(logf = log_dmvnorm, sigma = covmat, d = 2, n = 1,
           init = c(0, 0), rotate = FALSE, mean = c(1, 2))
test_that("BVN, non-zero mu, general Sigma, no rotation", {
  testthat::expect_equal(x2di$box, normal_box(d = 2, sigma = covmat,
                                              rotate = FALSE),
                         tolerance = my_tol)
})
x2di_mode <- ru(logf = log_dmvnorm, sigma = covmat, d = 2, n = 1,
                mode = c(1, 2), rotate = FALSE, mean = c(1, 2))
test_that("BVN, non-zero mu, general Sigma, no rotation, mode supplied", {
  testthat::expect_equal(x2di_mode$box, normal_box(d = 2, sigma = covmat,
                                                   rotate = FALSE),
                         tolerance = my_tol)
})

# (ii) rotate = TRUE
x2dii <- ru(logf = log_dmvnorm, sigma = covmat, d = 2, n = 1,
            init = c(0, 0), rotate = TRUE, mean = c(1, 2))
test_that("BVN, non-zero mu, general Sigma, rotation", {
  testthat::expect_equal(x2dii$box, normal_box(d = 2, sigma = covmat,
                                               rotate = TRUE),
                         tolerance = my_tol)
})
x2dii_mode <- ru(logf = log_dmvnorm, sigma = covmat, d = 2, n = 1,
                 mode = c(1, 2), rotate = TRUE, mean = c(1, 2))
test_that("BVN, non-zero mu, general Sigma, rotation, mode supplied", {
  testthat::expect_equal(x2dii_mode$box, normal_box(d = 2, sigma = covmat,
                                                    rotate = TRUE),
                         tolerance = my_tol)
})

# 3. 3-dimensional normal

# (a) Zero mean, unit variances with positive association
covmat <- matrix(rho, 3, 3) + diag(1 - rho, 3)

# (i) rotate = FALSE
x3ai <- ru(logf = log_dmvnorm, sigma = covmat, d = 3, n = 1,
           init = c(0, 0, 0), rotate = FALSE)
test_that("TVN, rho = 0.9, no rotation", {
  testthat::expect_equal(x3ai$box, normal_box(d = 3, sigma = covmat,
                                              rotate = FALSE),
                         tolerance = my_tol)
})
x3ai_mode <- ru(logf = log_dmvnorm, sigma = covmat, d = 3, n = 1,
                mode = c(0, 0, 0), rotate = FALSE)
test_that("TVN, rho = 0.9, no rotation, mode supplied", {
  testthat::expect_equal(x3ai_mode$box, normal_box(d = 3, sigma = covmat,
                                                   rotate = FALSE),
                         tolerance = my_tol)
})

# (ii) rotate = TRUE
x3aii <- ru(logf = log_dmvnorm, sigma = covmat, d = 3, n = 1,
            init = c(0, 0, 0), rotate = TRUE)
test_that("TVN, rho = 0.9, rotation", {
  testthat::expect_equal(x3aii$box, normal_box(d = 3, sigma = covmat,
                                               rotate = TRUE),
                         tolerance = my_tol)
})
x3aii_mode <- ru(logf = log_dmvnorm, sigma = covmat, d = 3, n = 1,
                 mode = c(0, 0, 0), rotate = TRUE)
test_that("TVN, rho = 0.9, rotation, mode supplied", {
  testthat::expect_equal(x3aii_mode$box, normal_box(d = 3, sigma = covmat,
                                                    rotate = TRUE),
                         tolerance = my_tol)
})

# B. 1-dimensional log-normal density

# Check that using a log transformation takes us back to the standard
# normal case
lambda <- 0
x <- ru(logf = stats::dlnorm, log = TRUE, d = 1, n = 1, init = 0.1,
        trans = "BC", lambda = lambda)
test_that("Log-normal", {
  testthat::expect_equal(x$box, normal_box(d = 1), tolerance = my_tol)
})
x_mode <- ru(logf = stats::dlnorm, log = TRUE, d = 1, n = 1, mode = 0,
             trans = "BC", lambda = lambda)
test_that("Log-normal, mode supplied", {
  testthat::expect_equal(x_mode$box, normal_box(d = 1), tolerance = my_tol)
})

# C: 1-dimensional gamma density, with shape parameter not less than 1

gamma_box <- function(shape = 1, rate = 1, r = 1 / 2) {
  #
  # Calculates bounding box values in the case of a zero-mean unnormalized
  # d-dimensional normal density, with arbitrary covariance structure.
  # "Unnormalized" means that the function has a maximum of 1 attained
  # at the origin, which is what is required for comparabilty with the
  # output from the function ru(), which scales the input function to have
  # a maximum of 1 and relocates the mode to the origin.
  #
  # Args:
  #   shape  : gamma shape parameter.
  #   rate   : gamma rate parameter.
  #   r      : ratio-of-uniforms tuning parameter.
  #
  # Returns:
  #   box : a  3 by 3 matrix of ratio-of-uniforms bounding box information
  #         with the same structure as object$box returned by ru().
  #
  if (shape < 1) {
    stop("A ratio-of-uniforms bounding box cannot be found when shape < 1")
  }
  # Mode of the gamma(shape, rate) density
  x_mode <- (shape - 1) / rate
  # Value of the (normalized) density at the mode
  x_max <- stats::dgamma(x_mode, shape = shape, rate = rate)
  #
  # Calculate box
  ar <- 1
  a_quad <- rate * r / (r + 1)
  b_quad <- -(1 + r * (shape - 1) / (r + 1) - a_quad * x_mode)
  c_quad <- -x_mode
  b_vals <- Re(polyroot(c(c_quad, b_quad, a_quad)))
  vals1 <- c(0, b_vals)
  b_fun <- function(x) {
    x * (stats::dgamma(x + x_mode, shape = shape) / x_max) ^ (r / (r + 1))
  }
  box <- c(ar, b_fun(b_vals))
  # Create matrix of bounding box information
  conv <- rep(0, 3)
  box <- cbind(box, vals1, conv)
  bs <- paste(paste("b", 1, sep=""),rep(c("minus", "plus"), each=1), sep="")
  rownames(box) <- c("a", bs)
  #
  return(box)
}

# Shape = 1
x1 <- suppressWarnings(ru(logf = dgamma, shape = 1, log = TRUE, d = 1, n = 1,
                          lower = 0, init = 1))
test_that("Gamma(1, 1)", {
  testthat::expect_equal(x1$box, gamma_box(shape = 1), tolerance = my_tol)
})
x1_mode <- suppressWarnings(ru(logf = dgamma, shape = 1, log = TRUE, d = 1,
                               n = 1, lower = 0, mode = 0))
test_that("Gamma(1, 1), mode supplied", {
  testthat::expect_equal(x1_mode$box, gamma_box(shape = 1), tolerance = my_tol)
})

# Shape = 10
x2 <- ru(logf = dgamma, shape = 10, log = TRUE, d = 1, n = 1, lower = 0,
         init = 10)
test_that("Gamma(10, 1)", {
  testthat::expect_equal(x2$box, gamma_box(shape = 10), tolerance = my_tol)
})
x2_mode <- ru(logf = dgamma, shape = 10, log = TRUE, d = 1, n = 1, lower = 0,
              mode = 9)
test_that("Gamma(10, 1), mode supplied", {
  testthat::expect_equal(x2_mode$box, gamma_box(shape = 10), tolerance = my_tol)
})
