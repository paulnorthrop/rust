# rust 1.1.0

## New features

* `plot.ru()` can now be used when `d > 2`: it produces pairwise plots of the
  simulated values.

* `find_lamba()`: argument `init_lambda` added to enable the user to supply an
   initial value for the Box-Cox transformation vector `lambda`.

## Bug fixes and minor improvements

* Unnecessary print statement `print(dim(phi))` removed from function 
  `find_lambda()`.

* Unnecessary print statement `print(a_algor)` removed from function 
  `ru()`.

* Correct `lambda$init` to `lambda$init_psi` in `ru()` when extracting 
  Box-Cox information.
   
* Documentation of `ru()` updated to include a description of the returned 
  function `logf_rho()` and simulated values `sim_vals_rho` and to clarify 
  the meaning of the returned value of `f_mode`.

* `ru()`: the expression for the inverse Box-Cox transformation in the case
  where lambda is exactly 0 has been corrected. 

* `find_lambda()`: carry out calculation of the target on a shifted log scale
  to avoid underflow.

* Set up `plot.ru()` so that if the user supplies axis labels then they are 
  used and otherwise the column name(s) of `ru_object$sim_vals` are used.  
  Also enable plotmath symbols to be rendered in the axis labels.
