# HazReg: Parametric Hazard-Based Regression Models for Survival Data

[![R package](https://img.shields.io/badge/language-R-blue.svg)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

## Overview

`HazReg` is an R package for fitting parametric hazard-based regression models
for survival data in both the **overall survival** and **relative (net) survival**
frameworks. The package is built around the **General Hazard (GH)** structure,
which nests the most widely used hazard regression models as special cases:

| Model | Abbreviation | Special case of GH |
|---|---|---|
| General Hazard | GH | — |
| Proportional Hazards | PH | ✓ |
| Accelerated Failure Time | AFT | ✓ |
| Accelerated Hazards | AH | ✓ |

Models are fitted by maximum likelihood via `nlminb` and `optim`. Users should
specify initial values and verify convergence of the optimisation, as is standard
practice for these routines.

## Installation

```r
# install.packages("pak")
pak::pak("FJRubio67/HazReg")
library(HazReg)
```

## Main functions

| Function | Framework | Description |
|---|---|---|
| `GHMLE` | Overall survival | Fits GH, PH, AFT, and AH models |
| `GEHMLE` | Relative (excess) survival | Fits excess hazard versions of the above |
| `simGH` | — | Simulates survival times from a GH structure |

For full documentation: `?GHMLE`, `?GEHMLE`, `?simGH`

## Baseline hazard distributions

Both `GHMLE` and `GEHMLE` support the following parametric baseline hazards.
All positive parameters are log-transformed for unconstrained optimisation.

| Distribution | Key | GH | PH | AFT | AH |
|---|---|:---:|:---:|:---:|:---:|
| [Power Generalised Weibull](http://rpubs.com/FJRubio/PGW) | PGW | ✓ | ✓ | ✓ | ✓ |
| [Exponentiated Weibull](http://rpubs.com/FJRubio/EWD) | EW | ✓ | ✓ | ✓ | ✓ |
| [Generalised Gamma](http://rpubs.com/FJRubio/GG) | GenGamma | ✓ | ✓ | ✓ | ✓ |
| [Gamma](https://en.wikipedia.org/wiki/Gamma_distribution) | Gamma | ✓ | ✓ | ✓ | ✓ |
| [Log-normal](https://en.wikipedia.org/wiki/Log-normal_distribution) | LogNormal | ✓ | ✓ | ✓ | ✓ |
| [Log-logistic](https://en.wikipedia.org/wiki/Log-logistic_distribution) | LogLogistic | ✓ | ✓ | ✓ | ✓ |
| [Weibull](https://en.wikipedia.org/wiki/Weibull_distribution) | Weibull | — | ✓ | ✓ | ✓ |

Hazard and related functions (PDF, CDF, survival) are also exported directly,
e.g. `?hpgw`, `?hggama`.

## Tutorials and examples

- [Overall survival: HazReg models](https://fjrubio-hazreg.share.connect.posit.cloud/) — Cloud Connect
- [Relative survival: Excess hazard models](https://fjrubio-xhazreg.share.connect.posit.cloud/) — Cloud Connect
- [Simulating from a GH structure](https://fjrubio-simgh.share.connect.posit.cloud/) — Cloud Connect
- [Simulating survival times from a General Hazard structure with a flexible baseline hazard](https://fjrubio-ghsim.share.connect.posit.cloud/) - Cloud Connect

## Related resources

- [Short course on Parametric Survival Analysis](https://github.com/FJRubio67/ShortCourseParamSurvival)
- [HazReg.jl](https://github.com/FJRubio67/HazReg.jl) — Julia implementation
- [MEGH](https://github.com/FJRubio67/MEGH) — GH models for clustered survival data
- [SimLT](https://github.com/FJRubio67/SimLT) — Simulating survival times from life tables

## Citation

If you use `HazReg` in your work, please cite the package and the relevant
methodological papers linked in the tutorials above.
