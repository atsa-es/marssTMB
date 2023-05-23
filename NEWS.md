# marssTMB 0.0.14

* Bug: if x0 had a column name (which it would if the data had column names), then the line to get the parvec would not have `x0.x0` as the name. It would have `x0` and this would cause a fail of the check that the parvec names match what `MARSSvectorizeparams()` expects.

# marssTMB 0.0.13

* Add check for data class in dfaTMB.
* Removed the kronecker products from the for loops to speed things up
* Don't create A x TT matrices when only need A x 1 because parameter is not time-varying
* Bug in `MARSS_TMB.R` that would not add one the obj when control$trace = 2.

# marssTMB 0.0.12

* Big clean. Move out uni.hpp (and R), estimate_marxss.hpp (and R), MARSS_tmb.R
* Keep dfaTMB for now while testing.

# marssTMB 0.0.11

* various bug fixes. Major one involving Z. Had (t(X) %x% I) instead of (I %x% t(X))
* bug with covariates since did not pass in form="marss" to coef()

# marssTMB 0.0.10

* added `marss2.hpp` and `estimate_marss2.R`. These use the chol() of the var-cov matrices and don't split into diag and corrmat. Seems more stable.

# marssTMB 0.0.9

Created a development version of {MARSS} with method TMB. `marssTMB()` will recognize if it was called with `MARSS()`. `MARSS_tmb()` also still works.

* Added   `marss.hpp`, a major update for MARSS models in vectorized "marss" form. This will allow time-varying parameters more easily and models with linear constraints. It's companion function is `estimate_marss.R`.
* Changed `MARSStmb.R` to `estimate_marxss.R`. It will become an internal function later.
* Added various helper functions in `inst/include/LOM.hpp`


# marssTMB 0.0.8

* Added a Quick Start vignette.
* Added check that no zeros on diagonal of Q or R (temporary)
* Added marssTMBCheckPackageVersions() to `zzz.R` and `.onLoad()` to deal with users getting warnings about TMB/Matrix version mismatch and not knowing what to do. https://glmmtmb.github.io/glmmTMB/#glmmtmbtmbmatrix-mismatches

# marssTMB 0.0.7

* updated `marxss.hpp` to be in MARSS format with X and Y as mxT and nxT.
* added Q, C, U, x0 and A estimation to `marxss.hpp`. Minimal testing so far.
* V0 = 0 is allowed.
* tinitx=1 or tinitx=0 allowed.

# marssTMB 0.0.6

* removed `to_marssTMB()` (not needed)
* completed `MARSStmb()` so that the marssMLE object is in proper form. All the {MARSS} helper functions should work.
* add more time comparisons to MARSS_tmb vignette and upped the maxit for MARSS()

# marssTMB 0.0.5

* Got `MARSStmb()` mostly working with marssMLE structure for output
* made default optimizer nlminb (faster)

# marssTMB 0.0.4

* Eric added tests and covariates to DFA vignette.
* moved MARSS into Depends and added to imports in `marssTMB-package.R`
* added some more info to `MARSStmb()` description
* started draft of vignette for `MARSS_tmb()`

# marssTMB 0.0.3

* Working on `to_marssMLE()` which will convert the output to MARSS form.

**To do**

* Fixing the parameter output from TMB to have the parameter names. *Done 0.0.6*
* Might also can add a `is.diag` flag to data if I need to id if R (and later Q) is diagonal to use faster code (diagonal matrices). *probably not needed?*

# marssTMB 0.0.2

* Added `MARSS_tmb()`. This sets up the model for `MARSStmb()` and for now ensures that the model will work with `marxss.hpp` which only allows DFA at the moment.
* Added `MARSStmb()`. This is symmetric to `MARSS::MARSSoptim()`. Working to match the output to `MARSS::MARSS()` so that all the MARSS functions work.
* Added `src/TMB/marxss.hpp` which is the MARSS model version with MARSS parameter names.

**To do**

* Currently `MARSStmb()` just fits a model. Next up is to convert this to MARSS output format.  *Done 0.0.5*

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