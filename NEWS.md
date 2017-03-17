# rust 1.1.0.9000

* Extra argument `bfgs_check` added `ru()` to allow the prevention of error messages resulting from calls to `stats::optim` with `method = "BFGS"`.  These calls occurred when other algorithms were not sure about the minumum they returned.  Often these minima are fine and starting very close to the minimum when calling `stats::optim` with `method = "BFGS"` can result in an error concerning finite differences.  The default is `bfgs_check = FALSE`, i.e. to avoid these extra calls to `stats:optim`.

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

# rust 1.1.0.9000

## New features

## Bug fixes and minor improvements

* "using `pairs()`" removed from the last sentence Description of `plot.ru()` because `pairs()` is not used when `d > 2`, rather a single plot is produced for each pair of variables.
