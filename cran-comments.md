## Resubmission

This is a resubmission, prompted by Brian identifying the remaining bug that prevented installation of v1.2.1 on Solaris.  There is only one change to the code in v1.2.1:

* An overloading ambiguity in function `loggp` has been corrected.

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
