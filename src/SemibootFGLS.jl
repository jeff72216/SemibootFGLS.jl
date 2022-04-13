module SemibootFGLS

using LinearAlgebra, Distributions

#-------------------------------------------------------------------
#  OLS/RLS estimates & residuals
#-------------------------------------------------------------------
function OLS(X, y)
    β̂_OLS = inv(transpose(X) * X) * transpose(X) * y
    e_OLS = y - X * β̂_OLS
    return (estimates=β̂_OLS, resuduals=e_OLS)
end

function RLS(X, y, R, q)
    β̂_OLS = inv(transpose(X) * X) * transpose(X) * y
    β̂_RLS = β̂_OLS + inv(transpose(X) * X) * transpose(R) * inv(R * inv(transpose(X) * X) * transpose(R)) * (q - R * β̂_OLS)
    e_RLS = y - X * β̂_RLS
    return (estimates=β̂_RLS, resuduals=e_RLS)
end

#-------------------------------------------------------------------
#  FGLS estimates & variance-covariance matrix
#-------------------------------------------------------------------
function FGLS(X, y, Σ̂)
    β̂_FGLS = inv(transpose(X) * inv(Σ̂) * X) * transpose(X) * inv(Σ̂) * y
    varβ̂_FGLS = inv(transpose(X) * inv(Σ̂) * X)
    return (estimates=β̂_FGLS, vcov=varβ̂_FGLS)
end

#-------------------------------------------------------------------
#  test statistic
#-------------------------------------------------------------------
function tStatistic(β̂, varβ̂, R, q)
    t = (R * β̂ - q) / sqrt(R * varβ̂ * transpose(R))
    return t
end

#-------------------------------------------------------------------
#  nonparatric estimates of Σ̂
#-------------------------------------------------------------------
function NW(h, X, e; constant::Bool, CV::Bool=false, kernel::Sampleable{Univariate,Continuous})
    if constant==true X=X[:, 2:size(X)[2]] end
    N = size(X)[1]
    w = zeros(N, N)
    for i in 1:N
        Xi = transpose(X[i, :])
        density = pdf.(kernel, (X .- Xi) ./ h)
        k = prod(density, dims=2)
        if CV==true k[i] = 0 end
        w[:, i] = k / sum(k)
    end
    σ̂² = transpose(w) * e.^2
    small = mean(e.^2) / 10
    σ̂²[σ̂².<small] .= small
    return σ̂²
end

function KNN(h, X, e; constant::Bool, CV::Bool=false)
    if constant==true X=X[:, 2:size(X)[2]] end
    N = size(X)[1]
    w = zeros(N, N)
    for i in 1:N
        Xi = transpose(X[i, :])
        dist = sqrt.(sum((X .- Xi).^2, dims=2))
        indx = partialsortperm(vec(dist), 1:h)[1:h]
        k = zeros(N)
        k[indx] .= 1
        if CV==true k[i] = 0 end
        w[:, i] = k ./ h
    end
    σ̂² = transpose(w) * e.^2
    small = mean(e.^2) / 10
    σ̂²[σ̂².<small] .= small
    return σ̂²
end

function LL(h, X, e; constant::Bool, CV::Bool=false, kernel::Sampleable{Univariate,Continuous})
    if constant==true X=X[:, 2:size(X)[2]] end
    N = size(X)[1]
    σ̂² = zeros(N)
    for i in 1:N
        Xi = transpose(X[i, :])
        density = pdf.(kernel, (X .- Xi) ./ h)
        k = prod(density, dims=2)
        k[k.<0.0001] .= 0.0001
        if CV==true k[i] = 0 end
        w = Diagonal(vec(k))
        constant==true ? z=hcat(ones(N,1),X.-Xi) : z=X.-Xi
        σ̂²[i] = (inv(transpose(z) * w * z) * transpose(z) * w * e.^2)[1, 1]
    end
    small = mean(e.^2) / 10
    σ̂²[σ̂².<small] .= small
    return σ̂²
end

function SR(h, X, e; constant::Bool)
    if constant==true X=X[:, 2:size(X)[2]] end
    N = size(X)[1]
    constant==true ? Q=hcat(ones(N,1),X) : Q=X
    for g in 2:h
        Q = hcat(Q, X.^g)
    end
    σ̂² = vec(Q * inv(transpose(Q) * Q) * transpose(Q) * e.^2)
    small = mean(e.^2) / 10
    σ̂²[σ̂².<small] .= small
    return σ̂²
end

#-------------------------------------------------------------------
#  cross validation
#-------------------------------------------------------------------
function CV(H, X, e; constant::Bool, nonparm::Function, KernelDensity=nothing)
    N = size(X)[1]
    MSE = zeros(length(H))
    for r in 1:length(H)
        h = H[r]
        if nonparm==NW || nonparm==LL
            cv_σ̂² = nonparm(h, X, e, constant=constant, CV=true, kernel=KernelDensity)
        elseif nonparm==SR
            cv_σ̂² = nonparm(h, X, e, constant=constant)
        else
            cv_σ̂² = nonparm(h, X, e, constant=constant, CV=true)
        end
        errorsq = (e.^2 - cv_σ̂²).^2
        nonparm==SR ? MSE[r]=sum(errorsq)/(N-r-1) : MSE[r]=mean(errorsq)
    end
    h_optimal = H[argmin(MSE)]
    return h_optimal
end

#-------------------------------------------------------------------
#  wild Bootstrap on the test statistic
#-------------------------------------------------------------------
function WildBootstrap(numResample, X, e, β̂, Σ̂, R, q, α; twotails::Bool)
    N = size(X)[1]
    t_memory = zeros(numResample)
    for b in 1:numResample
        u = float(rand(Bernoulli(0.5), N))
        u[u.==0] .= -1
        ỹ = X * β̂ + e.*u
        β̃ = FGLS(X, ỹ, Σ̂)[1]
        varβ̃ = FGLS(X, ỹ, Σ̂)[2]
        t_memory[b] = tStatistic(β̃, varβ̃, R, q)[1]
    end
    if twotails == true
        uppercrit = quantile(t_memory, 1 - α/2)
        lowercrit = quantile(t_memory, α/2)
    else
        uppercrit = quantile(t_memory, 1 - α)
        lowercrit = quantile(t_memory, α)
    end
    return (bootsample=t_memory, uppercritval=uppercrit, lowercritval=lowercrit)
end

#-------------------------------------------------------------------
#  the heteroskedasticity-robust semiparametric FGLS test 
#-------------------------------------------------------------------
function semiparamFGLS(X, y, bandwidth, R, q; constant::Bool, nonparm::Function, KernelDensity=nothing, CrossValidation::Bool=true, restricted::Bool=true, bootstrap::Bool=true, twotails=nothing, α=nothing, numResample=nothing)
    if restricted==false
        OLSresult = OLS(X, y)
        β̂ = OLSresult[1]
        e = OLSresult[2]
    else
        RLSresult = RLS(X, y, R, q)
        β̂ = RLSresult[1]
        e = RLSresult[2]
    end

    if CrossValidation==false
        if length(bandwidth)==1
            h = bandwidth
            if !isnothing(KernelDensity)
                Σ̂ = Diagonal(nonparm(h, X, e, constant=constant, kernel=KernelDensity))
            else
                Σ̂ = Diagonal(nonparm(h, X, e, constant=constant))
            end
        else
            error("bandwidth must be a scalar value")
        end
    else
        if length(bandwidth)==1
            error("bandwidth candidates must be a vector")
        else
            if !isnothing(KernelDensity)
                h = CV(bandwidth, X, e, nonparm=nonparm, constant=constant, KernelDensity=KernelDensity)
                Σ̂ = Diagonal(nonparm(h, X, e, constant=constant, kernel=KernelDensity))
            else
                h = CV(bandwidth, X, e, nonparm=nonparm, constant=constant)
                Σ̂ = Diagonal(nonparm(h, X, e, constant=constant))
            end
        end
    end

    FGLSresult = FGLS(X, y, Σ̂)
    β̂_FGLS = FGLSresult[1]
    varβ̂_FGLS = FGLSresult[2]
    t_FGLS = tStatistic(β̂_FGLS, varβ̂_FGLS, R, q)

    if bootstrap==false
        return (estimates=β̂_FGLS, vcov=varβ̂_FGLS, tstat=t_FGLS)
    else
        if isnothing(α) && isnothing(twotails) && isnothing(numResample)
            error("must provide details in the bootstrap function (see menu)")
        else
            boot = WildBootstrap(numResample, X, e, β̂, Σ̂, R, q, α, twotails=twotails)
            t_memory = boot[1]
            uppercrit = boot[2]
            lowercrit = boot[3]
            return (estimates=β̂_FGLS, vcov=varβ̂_FGLS, tstat=t_FGLS, bootsample=t_memory, uppercritval=uppercrit, lowercritval=lowercrit)
        end
    end
end

export OLS, RLS, FGLS, tStatistic, NW, KNN, LL, SR, CV, WildBootstrap, semiparamFGLS

end