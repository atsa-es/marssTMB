[![marssTMB status badge](https://atsa-es.r-universe.dev/badges/marssTMB)](https://atsa-es.r-universe.dev)

# marssTMB

## To use

Fit [MARSS](https://atsa-es.github.io/MARSS/) models as usual, but use [MARSS_tmb()] instead of [MARSS::MARSS()]. Note R and Q can only be diagonal, unconstrained or fixed. `equalvarcov` is not available nor are other custom Q matrices available with the EM algorithm in `MARSS()`. But TMB is very fast. 

See the [Quick Start Guide](https://atsa-es.github.io/marssTMB/articles/Quick_Start.html) to get started.

## install

```
# Install marssTMB in R:
install.packages('marssTMB', repos = c('https://atsa-es.r-universe.dev', 'https://cloud.r-project.org'))
```

## Notes

I used the {[TMBtools](https://github.com/mlysy/TMBtools)} package to set-up the R package to work with TMB. See the documentation there.

Rebuilding

* Add or edit models in `src` as hpp files. See TMBtools documentation on how to make a cpp file into a hpp file.
* Run `TMBtools::export_models()` whenever changes are made to any hpp files
* Rebuild the package
* Add helper files (includes) to `inst/include/marssTMB` and then put a `#include marssTMB/newfile.hpp` in your hpp file.
* Build package and documentation with `pkgdown::build_site()`


### NOAA Disclaimer

This repository is a scientific product and is not official communication of the National Oceanic and
Atmospheric Administration, or the United States Department of Commerce. All NOAA GitHub project code is
provided on an ‘as is’ basis and the user assumes responsibility for its use. Any claims against the Department of
Commerce or Department of Commerce bureaus stemming from the use of this GitHub project will be governed
by all applicable Federal law. Any reference to specific commercial products, processes, or services by service
mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recommendation or
favoring by the Department of Commerce. The Department of Commerce seal and logo, or the seal and logo of a
DOC bureau, shall not be used in any manner to imply endorsement of any commercial product or activity by
DOC or the United States Government.

