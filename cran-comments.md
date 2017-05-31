This is a re-submission to fix a significant bug in the last release (v1.2.0), which caused compilation errors on some platforms. I apologize for the initial error and for having to trouble you again.

* Corrected ::pow to std::pow and avoided the incorrect use of std::transform,  as requested by Brian Ripley.

## R CMD check results

0 errors | 0 warnings | 0 notes

## Test environments

- ubuntu 12.04 + GCC (on travis-ci), R-oldrel, R-release, R-devel
- ubuntu 12.04 + clang (on travis-ci), R-oldrel, R-release, R-devel
- osx (on travis-ci), R-oldrel, R-release
- win-builder (devel and release)

## Downstream dependencies

I have also run R CMD on downstream dependencies of rust.
All packages passed.
