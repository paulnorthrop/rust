## Resubmission

This is a resubmission, prompted by a request from Brian Ripley to correct the bugs that caused compilation errors in v1.2.0. I have

* Corrected ::pow to std::pow.

* Avoided using std::transform by not using it at all.

## R CMD check results

0 errors | 0 warnings | 0 notes

## Test environments

- Fedora Linux, clang, gfortran (on r-hub), R-devel 
- ubuntu 12.04 + GCC (on travis-ci), R-release, R-devel
- ubuntu 12.04 + clang (on travis-ci), R-release, R-devel
- osx (on travis-ci), R-oldrel, R-release
- win-builder (R-devel and R-release)

## Downstream dependencies

I have also run R CMD on downstream dependencies of rust.
All packages passed.
