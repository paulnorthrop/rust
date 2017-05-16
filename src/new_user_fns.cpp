#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::interfaces(r, cpp)]]

// User-supplied C++ functions.
// Note that currently the only interface available in rust is
// double fun(const Rcpp::NumericVector& x, const Rcpp::List& pars)
// Each function must be prefaced by the line: // [[Rcpp::export]]

// [[Rcpp::export]]
double new_logdgamma(SEXP xS, SEXP env) {
  Rcpp::NumericVector x(xS) ;
  Rcpp::Environment e(env);
  double alpha = e["alpha"] ;
  if (x[0] > 0)
    return ((alpha - 1.0) * log(x[0]) - x[0]) ;
  if (alpha > 1)
    return R_NegInf;
  else if (alpha < 1)
    return R_PosInf ;
  else
    return 0.0 ;
}

// A function to create external pointers for any of the functions above.
// See http://gallery.rcpp.org/articles/passing-cpp-function-pointers/
// If you write a new function above called new_name then add the following
//
// else if (fstr == "new_name")
//   return(Rcpp::XPtr<funcPtr>(new funcPtr(&new_name))) ;

// [[Rcpp::export]]
SEXP new_create_xptr(std::string fstr) {
  typedef double (*funcPtr)(SEXP xS, SEXP env) ;
  if (fstr == "new_logdgamma")
    return(Rcpp::XPtr<funcPtr>(new funcPtr(&new_logdgamma))) ;
  else
    return(Rcpp::XPtr<funcPtr>(R_NilValue)) ;
}

// We could create the external pointers when this file is sourced using this
// embedded R code below and/or (re)create them using create_xptr() in an
// R session or R package..

/*** R
new_ptr_gam <- new_create_xptr("new_logdgamma")
*/
