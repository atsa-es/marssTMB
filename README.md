# marssTMB

## To use

Fit MARSS models as usual, but use `MARSS_tmb()` instead of `MARSS()`. Not R and Q can only be diagonal, unconstrained or fixed. `equalvarcov` is not available nor are other custom Q matrices available with the EM algorithm in `MARSS()`. But TMB is very fast. 

See the documentation at [marssTMB](https://atsa-es.github.io/marssTMB/)

## install

```
# Install marssTMB in R:
install.packages('marssTMB', repos = c('https://atsa-es.r-universe.dev', 'https://cloud.r-project.org'))
```

## Notes

I used the {[TMBtools](https://github.com/mlysy/TMBtools)} package to set-up the R package to work with TMB. See the documentation there.

* Add new models to `src` as hpp files. Then run `TMBtools::export_models()` to add them to the various TMB set-up files. See TMBtools documentation on how to make a cpp file into a hpp file.
* Add helper files (includes) to `inst/include/marssTMB` and then put a `#include marssTMB/newfile.hpp` in your hpp file.

Rebuilding:

* Build package and documentation
* `pkgdown::build_site()`
