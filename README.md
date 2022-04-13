# SemibootFGLS.jl

A Julia package for a finite-sample heteroskedasticity-robust feasible generalized least squares (FGLS) test. Given a model:
$$
y=X\beta+\varepsilon
$$
where $\varepsilon$ has no serial correlation but potentially heteroskedastic, an FGLS estimator may be used:
$$
\hat{\beta}^{FGLS}=(X'\hat{\Sigma}^{-1}X)^{-1}X'\hat{\Sigma}^{-1}y
$$
Here the diagonal elements in $\hat{\Sigma}$ may be estimated nonparametrically to avoid potential mis-specification. This $\hat{\beta}^{FGLS}$ is therefore a semi-parametric estimator.

Bestdes, when the sample size is small, inferences based on asymptotic theory may be unreliable. Thus, one may be willing to apply a bootstrap method to the test statistic of $\hat{\beta}^{FGLS}$ to construct a simulated distribution.

The testing procedure is as follows:

**Step 1.** Estimate the model under the null hypothesis $H_0:R\beta=q$, where $R$ here is restricted to be a row vector and $q$ is a scalar, using restricted least squares (RLS):
$$
\tilde{\beta}^{RLS}=\hat{\beta}^{OLS}+(X'X)^{-1}R'(R(X'X)^{-1}R')^{-1}(q-R\hat{\beta}^{OLS})
$$ 
Obtain residuals $\tilde{e}^{RLS}$.

**Step 2.** Implement a nonparametric method to estimate the (usually unknown, allowing heteroskedasticity) skedastic function $\Sigma$. This package contains four different nonparametric approaches to estimate the diagonal elements $\sigma^2_i$ in $\Sigma$:

1. Nadaraya-Watson kernel estimator:
$
\hat{\sigma}^2_i=\frac{\sum_{j=1}^{n}\tilde{e}^{2,RLS}_{j}\mathcal{K}_{h}(x_{-1,i}-x_{-1,j})}{\sum_{j=1}^{n}\mathcal{K}_{h}(x_{-1,i}-x_{-1,j})}
$
where $\mathcal{K}_h$ is the kernel function with bandwidth $h$ and $x_{-1,i}$ is $x_i$ without constant.

2. $k$-nearest neighbors estimator:
$
\hat{\sigma}^2_i=\frac{1}{k}\sum_{j=1}^{n}\mathbf{I}_{k,i}(x_{-1,j})\tilde{e}^{2,RLS}_{j}
$
where $\mathbf{I}_{k,i}(x_{-1,j})$ is an indicator function equals one if $x_{-1,j}$ is within the $k$ nearest observations of $x_{-1,i}$.

3. Local linear estimator:
$
\hat{\sigma}^2_i=\iota'_1(Z'_iW_iZ_i)^{-1}Z'_iW_i\tilde{e}^{2,RLS}
$
where $\iota'_1=(1,0,...,0)$, $Z_i$ is an $n\times k$ matrix with $j$-th row equals $(1,x_{-1,j}-x_{-1,i})$, and $W_i$ is an $n\times n$ diagonal matrix of $\mathcal{K}_{h}(x_{-1,i}-x_{-1,j})$,.

4. Series estimator:
$
\big[\mathbf{\hat{\sigma}}_{1}^{2},...,\mathbf{\hat{\sigma}}_{n}^{2}
\big]^{\top}=Q(Q'Q)  ^{-1}Q'\tilde{e}^{2,RLS}
$
where $Q$ is a matrix of approximating function (in this package, a power series).

**Step 3.** Obtain the quasi-$t$ statistic of the FGLS estimator:
$$
t=\frac{R\hat{\beta}^{FGLS}-R\beta_{o}}{\sqrt{R\widehat{\text{var}}(\hat{\beta}^{FGLS})R'}}
$$

**Step 4** Perform a wild bootstrap to construct a bootstrap sample of the quasi-$t$ statistic. Test the null hypothesis according to this sample. Specifically, the wild bootstrap DGP is
$$
y^{*}_{i}=x'_i\tilde{\beta}^{RLS}+\tilde{e}^{RLS}_i\cdot u^{*}_i
$$
where 
$$u^{*}_i=\left\{\begin{array}{ll}
1 & \text{with probability }0.5\\
0 & \text{with probability }0.5
\end{array}\right.$$
Use the artificially-constructed samples $(y^{*}, X)$ to repeatedly calculate the quasi-$t$ statistics. The inference is made by comparing $t$ to the $\alpha/2$ and $1-\alpha/2$ quantiles in the bootstrap sample (for a two-tail test).



[![Build Status](https://github.com/jeff72216/SemibootFGLS.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/jeff72216/SemibootFGLS.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/jeff72216/SemibootFGLS.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/jeff72216/SemibootFGLS.jl)
