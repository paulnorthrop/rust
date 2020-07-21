This is a patch to avoid CRAN package check warnings on some platforms.

## R CMD check results

0 errors | 0 warnings | 0 notes

## Test environments

- ubuntu 12.04 + GCC (on travis-ci), R-release, R-devel
- ubuntu 12.04 + clang (on travis-ci), R-release, R-devel
- osx (on travis-ci), R-oldrel, R-release
- solaris-x86-patched using r-hub
- Debian Linux, R-devel, GCC, using r-hub
- win-builder (R-devel and R-release)

## Downstream dependencies

All downstream dependencies of rust (bang, revdbayes and threshr) passed R CMD check.
