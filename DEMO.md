# SemibootFGLS.jl

A Julia package for a finite-sample heteroskedasticity-robust feasible generalized least squares (FGLS) test.

### Functions
---
```julia
CV(H, X, e; constant::Bool, nonparm::Function, KernelDensity=nothing)
```
**Description**

Apply cross validation to find the optimal smoothing parameter of the nonparametric estimator.

**Arguments**

```H```: a Vector, UnitRange, or StepRangeLen of grid points of bandwidth.

```X```: an n by k matrix of explanatory variables.

```e```: a vector of length n of the residuals.

```constant```: a boolean variable, ```=true``` if ```X``` contains a constant term, ```=false``` otherwise.

```nonparm```: a function of nonparametric method which must be chosen from ```NW()```, ```KNN()```, ```LL()```, ```SR()``` in this package.

```KernelDensity```: a kernel density which comes from [Distributions.jl](https://github.com/JuliaStats/Distributions.jl). ```nothing``` by default. If ```nonparm=NW()``` or ```LL()```, then this argument must be provided.

**Value**

A float or integer of the optimal bandwidth will be returned.

<br>
<br>

```julia
FGLS(X, y, Σ̂)
```
**Description**

Return FGLS estimates and variance-covariance matrix

**Arguments**

```X```: an n by k matrix of explanatory variables.

```y```: a vector of length n of the explained variable.

```Σ̂``` : an n by n diagonal matrix of the estimated skedastic function.

**Value**

A named tuple ```(estimates, vcov)``` will be returned. 

<br>
<br>

```julia
KNN(h, X, e; constant::Bool, CV::Bool=false)
```
**Description**

A k-nearest neighbors estimator for the skedastic function Σ.

**Arguments**

```h```: a scalar number of bandwidth which can either be a float or an integer.

```X```: an n by k matrix of explanatory variables.

```e```: a vector of length n of the residuals.

```constant```: a boolean variable, ```=true``` if ```X``` contains a constant term, ```=false``` otherwise.

```CV```: a boolean variable (```false``` by default), ```=true``` if this function is implemented in a cross-validation procedure, ```=false``` otherwise.

**Value**

A vector of length n of the diagonal elements in Σ will be returned. 

<br>
<br>

```julia
LL(h, X, e; constant::Bool, CV::Bool=false, kernel::Sampleable{Univariate,Continuous})
```
**Description**

A local linear estimator for the skedastic function Σ.

**Arguments**

```h```: a scalar number of bandwidth which can either be a float or an integer.

```X```: an n by k matrix of explanatory variables.

```e```: a vector of length n of the residuals.

```constant```: a boolean variable, ```=true``` if ```X``` contains a constant term, ```=false``` otherwise.

```CV```: a boolean variable (```false``` by default), ```=true``` if this function is implemented in a cross-validation procedure, ```=false``` otherwise.

```kernel```: a kernel density which comes from [Distributions.jl](https://github.com/JuliaStats/Distributions.jl).

**Value**

A vector of length n of the diagonal elements in Σ will be returned. 

<br>
<br>

```julia
NW(h, X, e; constant::Bool, CV::Bool=false, kernel::Sampleable{Univariate,Continuous})
```
**Description**

A Nadaraya-Watson kernel estimator for the skedastic function Σ.

**Arguments**

```h```: a scalar number of bandwidth which can either be a float or an integer.

```X```: an n by k matrix of explanatory variables.

```e```: a vector of length n of the residuals.

```constant```: a boolean variable, ```=true``` if ```X``` contains a constant term, ```=false``` otherwise.

```CV```: a boolean variable (```false``` by default), ```=true``` if this function is implemented in a cross-validation procedure, ```=false``` otherwise.

```kernel```: a kernel density which comes from [Distributions.jl](https://github.com/JuliaStats/Distributions.jl).

**Value**

A vector of length n of the diagonal elements in Σ will be returned. 

<br>
<br>

```julia
OLS(X, y)
```
**Description**

Return OLS estimates and residuals

**Arguments**

```X```: an n by k matrix of explanatory variables.

```y```: a vector of length n of the explained variable.

**Value**

A named tuple ```(estimates, resuduals)``` will be returned. 

<br>
<br>

```julia
RLS(X, y, R, q)
```
**Description**

Return RLS estimates and residuals

**Arguments**

```X```: an n by k matrix of explanatory variables.

```y```: a vector of length n of the explained variable.

```R```: a 1 by k matrix of restrictions on coefficients.

```q```: a one-element vector of restriction Rβ=q.

**Value**

A named tuple ```(estimates, resuduals)``` will be returned. 

<br>
<br>

```julia
semiparamFGLS(X, y, bandwidth, R, q; constant::Bool, nonparm::Function, KernelDensity=nothing, CrossValidation::Bool=true, restricted::Bool=true, bootstrap::Bool=true, twotails=nothing, α=nothing, numResample=nothing)
```
**Description**

Apply wild bootstrap to resample the test statistic.

**Arguments**

```X```: an n by k matrix of explanatory variables.

```y```: a vector of length n of the explained variable.

```bandwidth```: if ```CrossValidation=false```, then provide a scalar number of integer or float bandwidth. If ```CrossValidation=true```, then provide a Vector, UnitRange, or StepRangeLen of grid points of bandwidth.

```R```: a 1 by k matrix of restrictions on coefficients.

```q```: a one-element vector of restriction Rβ=q.

```constant```: a boolean variable, ```=true``` if ```X``` contains a constant term, ```=false``` otherwise.

```nonparm```: a function of nonparametric method which must be chosen from ```NW()```, ```KNN()```, ```LL()```, ```SR()``` in this package.

```KernelDensity```: a kernel density which comes from [Distributions.jl](https://github.com/JuliaStats/Distributions.jl). ```nothing``` by default. If ```nonparm=NW()``` or ```LL()```, then this argument must be provided.

```CrossValidation```: if ```false```, then a cross validation will not be implemented to find the optimal bandwidth. ```true``` by default.

```restricted```: if ```true```, then the RLS estimate is used. If ```false```, then the OLS estimate is used. ```true``` by default.

```bootstrap```: if ```true```, then a wild bootstrap is implemented. ```false``` otherwise. ```true``` by default.

```twotails```: a boolean variable, ```=true``` if the test is two-tail, ```=false``` if the test is one-tail. ```nothing``` by default. Only necessary if ```bootstrap=true```.

```α```: significance level of the test. ```nothing``` by default. Only necessary if ```bootstrap=true```.

```numresample```: an integer of the number of bootstrap resampling. ```nothing``` by default. Only necessary if ```bootstrap=true```.

**Value**

A named tuple ```(estimates, vcov, tstat, bootsample, uppercritval, lowercritval)``` will be returned. ```estimates``` is the FGLS estimates, ```vcov``` is the variance-covariance matrix of the FGLS estimates, ```tstat``` is the quasi-t statistic, ```bootsample``` is the bootstrapped sample, ```uppercritval``` is the upper critical value of the test, ```lowercritval``` is the lower critical value of the test.

<br>
<br>

```julia
SR(h, X, e; constant::Bool)
```
**Description**

A series estimator for the skedastic function Σ. A power series is applied to be the approximating function.

**Arguments**

```h```: a scalar number of the number of power terms which can either be a float or an integer.

```X```: an n by k matrix of explanatory variables.

```e```: a vector of length n of the residuals.

```constant```: a boolean variable, ```=true``` if ```X``` contains a constant term, ```=false``` otherwise.

**Value**

A vector of length n of the diagonal elements in Σ will be returned. 

<br>
<br>

```julia
tStatistic(β̂, varβ̂, R, q)
```
**Description**

Return a (quasi)-t statistic

**Arguments**

```β̂```: a vector of length k of the estimates.

```varβ̂```: a k by k variance-covariance matrix of the estimates.

```R```: a 1 by k matrix of restrictions on coefficients.

```q```: a one-element vector of restriction Rβ=q.

**Value**

A 1 by 1 matrix of the test statistic will be returned. 

<br>
<br>

```julia
WildBootstrap(numResample, X, e, β̂, Σ̂, R, q, α; twotails::Bool)
```
**Description**

Apply wild bootstrap to resample the test statistic.

**Arguments**

```numresample```: an integer of the number of bootstrap resampling.

```X```: an n by k matrix of explanatory variables.

```e```: a vector of length n of the residuals.

```β̂```: a vector of length k of the estimates.

```Σ̂```: an n by n diagonal matrix of the estimated skedastic function.

```R```: a 1 by k matrix of restrictions on coefficients.

```q```: a one-element vector of restriction Rβ=q.

```α```: significance level of the test.

```twotails```: a boolean variable, ```=true``` if the test is two-tail, ```=false``` if the test is one-tail.

**Value**

A named tuple ```(bootsample, uppercritval, lowercritval)``` will be returned. ```bootsample``` is the bootstrapped sample, ```uppercritval``` is the upper critical value of the test, ```lowercritval``` is the lower critical value of the test.