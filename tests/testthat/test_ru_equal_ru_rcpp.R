context("rpost_rcpp vs rpost")

# We check that the values simulated using ru() and ru_rcpp() are
# (close enough to) identical when they are called using the same data,
# the same prior and starting from the same random number seed.

# Set a tolerance for the comparison of the simulated values

my_tol <- 1e-5
n <- 5
seed <- 27082017

# 1. One-dimensional standard normal ----------------

set.seed(seed)
x <- ru(logf = function(x) -x ^ 2 / 2, d = 1, n = n, init = 0.1)
ptr_N01 <- create_xptr("logdN01")
set.seed(seed)
x_rcpp <- ru_rcpp(logf = ptr_N01, d = 1, n = n, init = 0.1)
testthat::expect_equal(x$sim_vals, x_rcpp$sim_vals, tolerance = my_tol)

# 2. Posterior density of Generalized Pareto parameters

# Sample data from a GP(sigma, xi) distribution
seed <- 27082017

set.seed(seed)
gpd_data <- rgpd(m = 100, xi = 0, sigma = 1)
# Calculate summary statistics for use in the log-likelihood
ss <- gpd_sum_stats(gpd_data)
# Calculate an initial estimate
init <- c(mean(gpd_data), 0)

# 2a. Rotation of axes plus mode relocation ----------------
set.seed(seed)
x <- ru(logf = gpd_logpost, ss = ss, d = 2, n = n, init = init,
        lower = c(0, -Inf))

ptr_gp <- create_xptr("loggp")
for_ru_rcpp <- c(list(logf = ptr_gp, init = init, d = 2, n = n,
                      lower = c(0, -Inf)), ss)
set.seed(seed)
x_rcpp <- do.call(ru_rcpp, for_ru_rcpp)

testthat::expect_equal(x$sim_vals, x_rcpp$sim_vals, tolerance = my_tol)

# 2b. Box-Cox transformation, rotation of axes plus mode relocation -----------

# Calculate an initial estimate for phi = (phi1, phi2)
temp <- do.call(gpd_init, ss)
min_phi <- pmax(0, temp$init_phi - 2 * temp$se_phi)
max_phi <- pmax(0, temp$init_phi + 2 * temp$se_phi)

# find_lambda() -------------
phi_to_theta <- function(phi) c(phi[1], phi[2] - phi[1] / ss$xm)
log_j <- function(x) 0

lambda <- find_lambda(logf = gpd_logpost, ss = ss, d = 2, min_phi = min_phi,
                      max_phi = max_phi, phi_to_theta = phi_to_theta,
                      log_j = log_j)

# find_lambda_rcpp() -------------
ptr_gp <- create_xptr("loggp")
ptr_phi_to_theta_gp <- create_phi_to_theta_xptr("gp")
log_j <- create_log_jac_xptr("log_none_jac")
lambda_rcpp <- find_lambda_rcpp(logf = ptr_gp, ss = ss, d = 2,
                                min_phi = min_phi, max_phi = max_phi,
                                user_args = list(xm = ss$xm), log_j = log_j,
                                phi_to_theta = ptr_phi_to_theta_gp)

set.seed(seed)
x <- ru(logf = gpd_logpost, ss = ss, d = 2, n = n, trans = "BC",
        lambda = lambda)
set.seed(seed)
x_rcpp <- ru_rcpp(logf = ptr_gp, ss = ss, d = 2, n = n, trans = "BC",
                  lambda = lambda_rcpp)

testthat::expect_equal(x$sim_vals, x_rcpp$sim_vals, tolerance = my_tol)
