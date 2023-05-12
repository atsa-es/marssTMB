[![marssTMB status badge](https://atsa-es.r-universe.dev/badges/marssTMB)](https://atsa-es.r-universe.dev)

# marssTMB

## To use

Install MARSS 3.11.6 (development version)
```
install.packages('MARSS', repos = c('https://atsa-es.r-universe.dev', 'https://cloud.r-project.org'))
```

Install marssTMB
```
install.packages('marssTMB', repos = c('https://atsa-es.r-universe.dev', 'https://cloud.r-project.org'))
```

## Fit a MARSS model

See the [Quick Start Guide](https://atsa-es.github.io/marssTMB/articles/Quick_Start.html) to get started. Call `MARSS()` with `method="TMB"`
```
library(MARSS)
y <- t(harborSealWA); dat <- dat[2:4, ]
MARSS(y, method="TMB")
```

Note R and Q can only be diagonal, unconstrained or fixed. `equalvarcov` is not available nor are zeros allowed on the diagonal as is available with the Kalman filter + EM and BFGS algorithm in `MARSS()`. But TMB is very fast. Also be aware that {marssTMB} is still in early development and you should check your answers against `method="kem"` and `method="BFGS"`.

## Developer notes

Install {[TMBtools](https://github.com/mlysy/TMBtools)} package which is being used to help with set-up.

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

