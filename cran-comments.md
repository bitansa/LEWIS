## Submission

This is the first submission of LEWIS.

## Test environments

* local macOS 27 (aarch64), R 4.6.1
* GitHub Actions: macOS (release), Windows (release), Ubuntu (devel, release, oldrel-1)
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new submission.
* "DAGs" in the Description is flagged as possibly misspelled. It is the
  plural of DAG (directed acyclic graph), which is defined in the same
  field.

## Notes for the CRAN team

* The method is described in a manuscript that is under review. There is no
  DOI to cite in the Description field yet.
* The package contains copies of code from three other sources. All
  copyright holders are listed in Authors@R.
  * The Fortran routines rqbr and rqfnb of the package quantreg (Roger
    Koenker, GPL (>= 2)). The routine rqbr is called through a small
    Fortran wrapper (src/rqbr_wrap.f) that passes its LOGICAL argument as
    an INTEGER.
  * The Shapiro-Wilk algorithm of R's file swilk.c (R Core Team,
    GPL (>= 2)).
  * The TL functions of the package TransGraph 1.1.0 (Mingyang Ren,
    Ruixuan Zhao, Xin He and Junhui Wang, GPL-2), in R/tl.R. TransGraph is
    archived on CRAN. LEWIS uses its TL method for the default start of
    the sampler. For this reason LEWIS is licensed GPL-2.
* The compiler warnings from src/rqbr.f and src/rqfnb.f (unused labels and
  similar) come from the unchanged quantreg code.
* The data set tcga_prostate is derived from open-access data of The Cancer
  Genome Atlas.
* The CRAN status badge in README.md points to the CRAN page of the
  package, which does not exist before acceptance.
* The examples run in a few seconds. The longer comparison with a pure R
  implementation in tests/testthat is skipped on CRAN.

## Downstream dependencies

There are currently no downstream dependencies.
