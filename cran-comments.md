## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.
* We are currently working on a manuscript describing the methods used in the package, so we don't have any references as of this moment.
* We wrap examples with `\donttest{}` due to long run times (i.e., >5 secs).
* We now use testthat with small toy examples to allow automatic testing.
* We made sure in R/generics.R; R/plot_data.R; R/trace_plot.R have an immediate call of on.exit().