# HazReg R package

The `HazReg` R package implements the following parametric hazard-based regression models for (overall) survival data.

- General Hazard (GH) model.

- Accelerated Failure Time (AFT) model.

- Proportional Hazards (PH) model.

- Accelerated Hazards (AH) model.


These models are fitted using the R commands `nlminb` and `optim`. Thus, the user needs to specify the initial points and to check the convergence of the optimisation step, as usual.


The current version of the `HazReg` R package implements the following parametric baseline hazards for the models discussed in the previous section.

- [Power Generalised Weibull](http://rpubs.com/FJRubio/PGW) (PGW) distribution. `PGWMLE`
 
- [Exponentiated Weibull](http://rpubs.com/FJRubio/EWD) (EW) distribution. `EWMLE`
 
- [Generalised Gamma](http://rpubs.com/FJRubio/GG) (GG) distribuiton. `GGMLE`

- [Gamma](https://en.wikipedia.org/wiki/Gamma_distribution) (G) distribution. `GMLE`

- [Lognormal](https://en.wikipedia.org/wiki/Log-normal_distribution) (LN) distribution. `LNMLE`

- [Log-logistic](https://en.wikipedia.org/wiki/Log-logistic_distribution) (LL) distribution. `LLMLE`

- [Weibull](https://en.wikipedia.org/wiki/Weibull_distribution) (W) distribution. (only for AFT, PH, and AH models) `WMLE`


All positive parameters are transformed into the real line using a `log` reparameterisation.

An Illustrative example and a description of the available models can be found at:

[HazReg: Parametric Hazard-based regression models for survival data](https://rpubs.com/FJRubio/HazReg) [RPubs]

[HazReg: Parametric Hazard-based regression models for survival data](https://fjrubio.quarto.pub/hazreg/) [quarto-pub]

```
library(devtools)
install_github("FJRubio67/HazReg")

library(HazReg)
?PGWMLE
?GGMLE
?hpgw
?hggama
```

See also: [GHSurv](https://github.com/FJRubio67/GHSurv), [LBANS](https://github.com/FJRubio67/LBANS)
