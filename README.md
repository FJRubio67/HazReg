# HazReg

The `HazReg` R package implements the following parametric hazard-based regression models for (overall) survival data.

- General Hazard (GH) model.

- Accelerated Failure Time (AFT) model.

- Proportional Hazards (PH) model.

- Accelerated Hazards (AH) model.

The current version of the `HazReg` R package implements the following parametric baseline hazards for the models discussed in the previous section.

- [Power Generalised Weibull](http://rpubs.com/FJRubio/PGW) (PGW) distribution.

- [Exponentiated Weibull](http://rpubs.com/FJRubio/EWD) (EW) distribution.

- [Generalised Gamma](http://rpubs.com/FJRubio/GG) (GG) distribuiton.

- [Gamma](https://en.wikipedia.org/wiki/Gamma_distribution) (G) distribution.

- [Lognormal](https://en.wikipedia.org/wiki/Log-normal_distribution) (LN) distribution.

- [Log-logistic](https://en.wikipedia.org/wiki/Log-logistic_distribution) (LL) distribution.

- [Weibull](https://en.wikipedia.org/wiki/Weibull_distribution) (W) distribution. 

All positive parameters are transformed into the real line using a `log` reparameterisation.

An Illustrative example and a description of the available models can be found at:

[HazReg: Parametric Hazard-based regression models for survival data](https://rpubs.com/FJRubio/HazReg)

```
library(devtools)
install_github("FJRubio67/HazReg")

library(HazReg)
?PGWMLE
?GGMLE
?hpgw
?hggama
```

