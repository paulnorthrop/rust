# rust 1.1.0.9000

## New features

* Package Rcpp (https://CRAN.R-project.org/package=Rcpp) can be used to speed up the computations.  New functions: `ru_rcpp`, `find_lambda_rcpp` and `find_lambda_one_d_rcpp`.  

* New vignette. "Rusting faster: Simulation using Rcpp".

## Bug fixes and minor improvements

* Bug fixed in `plot.ru()`: previously `plot.ru()` failed when `d > 2`  and no axis label names were provided.

* Bug fixed in `plot.ru` : previously, in the `d = 2` case, providing the graphical parameter `col` produced an error because `col = 8` was hard-coded in a call to `points`. Now the extra argument `points_par` enables the user to provide a list of arguments to `points`.

* "using `pairs()`" removed from the last sentence Description of `plot.ru()` because `pairs()` is not used when `d > 2`, rather a single plot is produced for each pair of variables.

# rust 1.1.0

## New features

* `plot.ru()` can now be used when `d > 2`: it produces pairwise plots of the simulated values.

* `find_lamba()`: argument `init_lambda` added to enable the user to supply an initial value for the Box-Cox transformation vector `lambda`.

## Bug fixes and minor improvements

* Unnecessary print statement `print(dim(phi))` removed from function `find_lambda()`.

* Unnecessary print statement `print(a_algor)` removed from function `ru()`.

* Correct `lambda$init` to `lambda$init_psi` in `ru()` when extracting Box-Cox information.
   
* Documentation of `ru()` updated to include a description of the returned function `logf_rho()` and simulated values `sim_vals_rho` and to clarify the meaning of the returned value of `f_mode`.

* `ru()`: the expression for the inverse Box-Cox transformation in the case where lambda is exactly 0 has been corrected. 

* `find_lambda()`: carry out calculation of the target on a shifted log scale to avoid underflow.

* Set up `plot.ru()` so that if the user supplies axis labels then they are used and otherwise the column name(s) of `ru_object$sim_vals` are used. Also enable plotmath symbols to be rendered in the axis labels.
