Please may I submit a patch release to fix the installation ERROR of rust v1.3.3 on r-oldrel?

I have changed R (>= 3.4.1) back to R (>= 3.3.1) in DESCRIPTION to fix this.

Please accept my apologies for wasting your time by making such a stupid mistake.

## R CMD check results

0 errors | 0 warnings | 0 notes

## Test environments

- Fedora Linux, clang, gfortran (on r-hub), R-devel 
- ubuntu 12.04 + GCC (on travis-ci), R-release, R-devel
- ubuntu 12.04 + clang (on travis-ci), R-release, R-devel
- osx (on travis-ci), R-oldrel, R-release
- win-builder (R-devel and R-release)

## Downstream dependencies

All downstream dependencies of rust passed R CMD check.
