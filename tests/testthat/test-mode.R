#context("rpost_rcpp vs rpost")

# We check that setting the mode using the mode argument gives the same output
# as suppling this value as the initial vale init.
# Remember that if lambda = 1 then the effect is to subtract 1 from the variable

# Normal x log-normal: different Box-Cox parameters ==================
norm_lognorm <- function(x, ...) {
 dnorm(x[1], ...) + dlnorm(x[2], ...)
}
set.seed(12072022)
x1 <- ru(logf = norm_lognorm, log = TRUE, n = 100, d = 2, init = c(-1, 0),
         trans = "BC", lambda = c(1, 0))
set.seed(12072022)
x2 <- ru(logf = norm_lognorm, log = TRUE, n = 100, d = 2, mode = c(-1, 0),
         trans = "BC", lambda = c(1, 0))

my_tol <- 1e-5
test_that("init vs mode, sim_vals", {
  testthat::expect_equal(x1$sim_vals, x2$sim_vals, tolerance = my_tol)
})
test_that("init vs mode, box", {
  testthat::expect_equal(x1$box, x2$box, tolerance = my_tol)
})
