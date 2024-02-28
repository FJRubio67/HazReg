# HazReg R package

## Overall Survival

The `HazReg` R package implements the following parametric hazard-based regression models for (overall) survival data.

- General Hazard (GH) model.

- Accelerated Failure Time (AFT) model.

- Proportional Hazards (PH) model.

- Accelerated Hazards (AH) model.


These models are fitted using the R commands `nlminb` and `optim`. Thus, the user needs to specify the initial points and to check the convergence of the optimisation step, as usual.


The current version of the `HazReg` R package implements the following parametric baseline hazards for the models discussed in the previous section, using the command `GHMLE`.

- [Power Generalised Weibull](http://rpubs.com/FJRubio/PGW) (PGW) distribution. 
 
- [Exponentiated Weibull](http://rpubs.com/FJRubio/EWD) (EW) distribution. 
 
- [Generalised Gamma](http://rpubs.com/FJRubio/GG) (GenGamma) distribuiton. 

- [Gamma](https://en.wikipedia.org/wiki/Gamma_distribution) (Gamma) distribution. 

- [Lognormal](https://en.wikipedia.org/wiki/Log-normal_distribution) (LogNormal) distribution. 

- [Log-logistic](https://en.wikipedia.org/wiki/Log-logistic_distribution) (LogLogistic) distribution. 

- [Weibull](https://en.wikipedia.org/wiki/Weibull_distribution) (Weibull) distribution. (only for AFT, PH, and AH models) 


All positive parameters are transformed into the real line using a `log` link (reparameterisation).

An Illustrative example and a description of the available models can be found at:

- [HazReg: Parametric Hazard-based regression models for survival data](https://rpubs.com/FJRubio/HazReg) [RPubs]
- [HazReg: Parametric Hazard-based regression models for survival data](https://fjrubio.quarto.pub/hazreg/) [quarto-pub]
- [simGH: simulating times to event from a general hazard structure](https://rpubs.com/FJRubio/simGH)

```
library(devtools)
install_github("FJRubio67/HazReg")

library(HazReg)

?GHMLE

?hpgw

?hggama

?simGH
```

## Relative (Net) Survival

Forthcoming

### See also: 
- [Simulating survival times from a General Hazard structure with a flexible baseline hazard](https://rpubs.com/FJRubio/GHSim)
- [Short course on Parametric Survival Analysis](https://github.com/FJRubio67/ShortCourseParamSurvival)
- [HazReg.jl](https://github.com/FJRubio67/HazReg.jl)
