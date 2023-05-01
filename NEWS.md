# marssTMB 0.0.4

* Eric added tests and covariates to DFA vignette.
* moved MARSS into Depends and adde to imports in `marssTMB-package.R`
* added some more info to `MARSStmb()` description

# marssTMB 0.0.3

* Working on `to_marssMLE()` which will convert the output to MARSS form.

**To do**

Fixing the parameter output from TMB to have the parameter names. Z and D are easy. For R, need to change how `marxss.hpp` deals with R. I think I don't need to separate into diag and non-diagonal elements. See how I dealt with this in `MARSSoptim()`. Also can add a `is.diag` flag to data if I need to id if R (and later Q) is diagonal to use faster code (diagonal matrices).

# marssTMB 0.0.2

* Added `MARSS_tmb()`. This sets up the model for `MARSStmb()` and for now ensures that the model will work with `marxss.hpp` which only allows DFA at the moment.
* Added `MARSStmb()`. This is symmetric to `MARSSoptim()`. Working to match the output to `MARSS()` so that all the MARSS functions work.
* Added `src/TMB/marxss.hpp` which is the MARSS model version with MARSS parameter names.

**To do**

* Currently `MARSStmb()` just fits a model. Next up is to convert this to MARSS output format.

# marssTMB 0.0.1

The first draft with the TMB package structure in place. Used the {[TMBtools](https://github.com/mlysy/TMBtools)} package to set-up the R package to work with TMB.

* Two functions `dfaTMB()` and `uniTMB()` which I will likely combine later.

<!--

## Breaking changes

* 

* 

## New features

* 

## Bug fixes

* 
-->