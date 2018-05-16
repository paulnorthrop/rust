## Resubmission

This is a resubmission.  I have fixed a bug that caused an error in the 
building of a vignette of revdbayes package, a reverse dependency of rust.   

(I have kept the version number at 1.3.5 as the submission didn't make it to 
CRAN, but I'm happy to increment it if that's preferable.)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Test environments

- Fedora Linux, clang, gfortran (on r-hub), R-devel 
- ubuntu 12.04 + GCC (on travis-ci), R-release, R-devel
- ubuntu 12.04 + clang (on travis-ci), R-release, R-devel
- osx (on travis-ci), R-oldrel, R-release
- win-builder (R-devel and R-release)

## Downstream dependencies

All downstream dependencies of rust (bang, revdbayes and threshr) passed R CMD check.
