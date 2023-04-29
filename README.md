# marssTMB

The beginning of a companion R package to MARSS that fits MARSS models with TMB.

See the documentation at [marssTMB]()

## install

```
remotes::install_github("atsa-es/marssTMB")
```

## Notes

I used the {[TMBtools](https://github.com/mlysy/TMBtools)} package to set-up the R package to work with TMB. See the documentation there.

* Add new models to `src` as hpp files. Then run `TMBtools::export_models()` to add them to the various TMB set-up files. See TMBtools documentation on how to make a cpp file into a hpp file.
* Add helper files (includes) to `inst/include/marssTMB` and then put a `#include marssTMB/newfile.hpp` in your hpp file.

Rebuilding:

* Build package and documentation
* `pkgdown::build_site()`
