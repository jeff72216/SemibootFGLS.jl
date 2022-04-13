# SemibootFGLS.jl

A Julia package for a finite-sample heteroskedasticity-robust feasible generalized least squares (FGLS) test. Given a model:

![equation](https://latex.codecogs.com/svg.image?{\color{White}&space;y=X\beta&plus;\varepsilon})

where ![equation](https://latex.codecogs.com/svg.image?{\color{White}&space;\varepsilon}) has no serial correlation but potentially heteroskedastic, an FGLS estimator may be used:

![equation](https://latex.codecogs.com/svg.image?{\color{White}&space;\hat{\beta}^{FGLS}=(X'\hat{\Sigma}^{-1}X)^{-1}X'\hat{\Sigma}^{-1}y})

Here the diagonal elements in ![equation](https://latex.codecogs.com/svg.image?{\color{White}&space;\hat{\Sigma}}) may be estimated nonparametrically to avoid potential mis-specification. This ![equation](https://latex.codecogs.com/svg.image?{\color{White}&space;\hat{\beta}^{FGLS}}) is therefore a semi-parametric estimator.

Besides, when the sample size is small, inferences based on asymptotic theory may be unreliable. Thus, one may be willing to apply a bootstrap method to the test statistic of ![equation](https://latex.codecogs.com/svg.image?{\color{White}&space;\hat{\beta}^{FGLS}}) to construct a simulated distribution.

The testing procedure is as follows:

**Step 1.** Estimate the model under the null hypothesis ![equation](https://latex.codecogs.com/svg.image?{\color{White}&space;H_0:R\beta=q}), where R here is restricted to be a row vector and q is a scalar, using restricted least squares (RLS):

![equation](https://latex.codecogs.com/svg.image?{\color{White}&space;\tilde{\beta}^{RLS}=\hat{\beta}^{OLS}&plus;(X'X)^{-1}R'(R(X'X)^{-1}R')^{-1}(q-R\hat{\beta}^{OLS})})

Obtain residuals ![equation](https://latex.codecogs.com/svg.image?{\color{White}&space;\tilde{e}^{RLS}}).

**Step 2.** Implement a nonparametric method to estimate the (usually unknown, allowing heteroskedasticity) skedastic function ![equation](https://latex.codecogs.com/svg.image?{\color{White}&space;\Sigma}). This package contains four different nonparametric approaches to estimate the diagonal elements ![equation](https://latex.codecogs.com/svg.image?{\color{White}&space;\sigma^2_i}) in ![equation](https://latex.codecogs.com/svg.image?{\color{White}&space;\Sigma}):

1. Nadaraya-Watson kernel estimator:
![equation](https://latex.codecogs.com/svg.image?{\color{White}&space;\hat{\sigma}^2_i=\frac{\sum_{j=1}^{n}\tilde{e}^{2,RLS}_{j}\mathcal{K}_{h}(x_{-1,i}-x_{-1,j})}{\sum_{j=1}^{n}\mathcal{K}_{h}(x_{-1,i}-x_{-1,j})}})
where ![equation](https://latex.codecogs.com/svg.image?{\color{White}&space;\mathcal{K}_h}) is the kernel function with bandwidth h and ![equation](https://latex.codecogs.com/svg.image?{\color{White}&space;x_{-1,i}}) is ![equation](https://latex.codecogs.com/svg.image?{\color{White}&space;x_{i}}) without constant.

2. k-nearest neighbors estimator:
![equation](https://latex.codecogs.com/svg.image?{\color{White}&space;\hat{\sigma}^2_i=\frac{1}{k}\sum_{j=1}^{n}\mathbf{I}_{k,i}(x_{-1,j})\tilde{e}^{2,RLS}_{j}})
where ![equation](https://latex.codecogs.com/svg.image?{\color{White}&space;\mathbf{I}_{k,i}(x_{-1,j})}) is an indicator function equals one if ![equation](https://latex.codecogs.com/svg.image?{\color{White}&space;x_{-1,j}}) is within the k nearest observations of ![equation](https://latex.codecogs.com/svg.image?{\color{White}&space;x_{-1,i}}).

3. Local linear estimator:
![equation](https://latex.codecogs.com/svg.image?{\color{White}&space;\hat{\sigma}^2_i=\iota'_1(Z'_iW_iZ_i)^{-1}Z'_iW_i\tilde{e}^{2,RLS}})
where ![equation](https://latex.codecogs.com/svg.image?{\color{White}&space;\iota'_1=(1,0,...,0)}), ![equation](https://latex.codecogs.com/svg.image?{\color{White}&space;Z_{i}}) is an n by k matrix with j-th row equals ![equation](https://latex.codecogs.com/svg.image?{\color{White}&space;(1,x_{-1,j}-x_{-1,i})}), and ![equation](https://latex.codecogs.com/svg.image?{\color{White}&space;W_{i}}) is an n by n diagonal matrix of ![equation](https://latex.codecogs.com/svg.image?{\color{White}&space;\mathcal{K}_{h}(x_{-1,i}-x_{-1,j})}),.

4. Series estimator:
![equation](https://latex.codecogs.com/svg.image?{\color{White}&space;\big[\mathbf{\hat{\sigma}}_{1}^{2},...,\mathbf{\hat{\sigma}}_{n}^{2}\big]^{\top}=Q(Q'Q)&space;&space;^{-1}Q'\tilde{e}^{2,RLS}})
where ![equation](https://latex.codecogs.com/svg.image?{\color{White}&space;Q}) is a matrix of approximating function (in this package, a power series).

**Step 3.** Obtain the quasi-t statistic of the FGLS estimator:

![equation](https://latex.codecogs.com/svg.image?{\color{White}&space;t=\frac{R\hat{\beta}^{FGLS}-R\beta_{o}}{\sqrt{R\widehat{\text{var}}(\hat{\beta}^{FGLS})R'}}})

**Step 4** Perform a wild bootstrap to construct a bootstrap sample of the quasi-t statistic. Test the null hypothesis according to this sample. Specifically, the wild bootstrap DGP is

![equation](https://latex.codecogs.com/svg.image?{\color{White}&space;y^{*}_{i}=x'_i\tilde{\beta}^{RLS}&plus;\tilde{e}^{RLS}_i\cdot&space;u^{*}_i})

where ![equation](https://latex.codecogs.com/svg.image?{\color{White}&space;u_{i}^*}) equals 1 and 0 with equal probability. Use the artificially-constructed samples ![equation](https://latex.codecogs.com/svg.image?{\color{White}&space;(y^{*},X)}) to repeatedly calculate the quasi-t statistics. The inference is made by comparing t to the α/2 and 1-α/2 quantiles in the bootstrap sample (for a two-tail test).



[![Build Status](https://github.com/jeff72216/SemibootFGLS.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/jeff72216/SemibootFGLS.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/jeff72216/SemibootFGLS.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/jeff72216/SemibootFGLS.jl)
