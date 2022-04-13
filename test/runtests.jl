using SemibootFGLS
using LinearAlgebra, Distributions
using Test

@testset "SemibootFGLS.jl" begin
    # test dataset1
    x1 = ones(100)
    x25 = rand(Normal(3, 1), 100, 4)
    X = hcat(x1, x25)
    β = [1, 1, 1, 1, 0]
    ε = rand(100)
    y = X * β + ε
    R = [0 0 0 0 1] # linear restriction: Rβ=q
    q = [0]

    # test dataset2
    X2 = x25
    β2 = [1, 1, 1, 0]
    y2 = X2 * β2 + ε
    R2 = [0 0 0 1]
    
    # OLS/RLS test
    @test length(OLS(X, y)[1]) == length(β)
    @test length(OLS(X, y)[2]) == size(X)[1]
    @test length(RLS(X, y, R, q)[1]) == length(β)
    @test length(RLS(X, y, R, q)[2]) == size(X)[1]
    
    # nonparametric test
    e = RLS(X, y, R, q)[2]
    e2 = RLS(X2, y2, R2, q)[2]
    @test length(NW(5, X, e, constant=true, kernel=Normal(0, 1))) == length(e)
    @test length(NW(5, X2, e2, constant=false, kernel=Epanechnikov(0, 1))) == length(e2)
    @test length(KNN(5, X, e, constant=true)) == length(e)
    @test length(KNN(5, X2, e2, constant=false)) == length(e2)
    @test length(LL(5, X, e, constant=true, kernel=Normal(0, 1))) == length(e)
    @test length(LL(5, X2, e2, constant=false, kernel=Epanechnikov(0, 1))) == length(e2)
    @test length(SR(5, X, e, constant=true)) == length(e)
    @test length(SR(5, X2, e2, constant=false)) == length(e2)
    bw_NW = 1.0:0.5:50.0
    @test typeof(CV(bw_NW, X, e, constant=true, nonparm=NW, KernelDensity=Normal(0, 1))) == Float64
    @test typeof(CV(bw_NW, X2, e2, constant=false, nonparm=NW, KernelDensity=Epanechnikov(0, 1))) == Float64
    bw_KNN = 1:50
    @test typeof(CV(bw_KNN, X, e, constant=true, nonparm=KNN)) == Int64
    @test typeof(CV(bw_KNN, X2, e2, constant=false, nonparm=KNN)) == Int64
    bw_LL = 1.0:0.5:50.0
    @test typeof(CV(bw_LL, X, e, constant=true, nonparm=LL, KernelDensity=Normal(0, 1))) == Float64
    @test typeof(CV(bw_LL, X2, e2, constant=false, nonparm=LL, KernelDensity=Epanechnikov(0, 1))) == Float64
    bw_SR = 1:10
    @test typeof(CV(bw_SR, X, e, constant=true, nonparm=SR)) == Int64
    @test typeof(CV(bw_SR, X2, e2, constant=false, nonparm=SR)) == Int64

    # FGLS test
    Σ̂ = Diagonal(NW(5, X, e, constant=true, kernel=Normal(0, 1)))
    @test length(FGLS(X, y, Σ̂)[1]) == length(β)
    @test size(FGLS(X, y, Σ̂)[2]) == (length(β), length(β))
    
    # t-statistic test
    β̂ = FGLS(X, y, Σ̂)[1]
    varβ̂ = FGLS(X, y, Σ̂)[2]
    @test size(tStatistic(β̂, varβ̂, R, q)) == (length(q), length(q))

    # wild bootstrap test
    bootsample = 999
    @test length(WildBootstrap(bootsample, X, e, β̂, Σ̂, R, q, 0.05, twotails=true)[1]) == bootsample
    @test length(WildBootstrap(bootsample, X, e, β̂, Σ̂, R, q, 0.025, twotails=false)[1]) == bootsample
    @test typeof(WildBootstrap(bootsample, X, e, β̂, Σ̂, R, q, 0.05, twotails=true)[2]) == Float64
    @test typeof(WildBootstrap(bootsample, X, e, β̂, Σ̂, R, q, 0.025, twotails=false)[2]) == Float64
    @test typeof(WildBootstrap(bootsample, X, e, β̂, Σ̂, R, q, 0.05, twotails=true)[3]) == Float64
    @test typeof(WildBootstrap(bootsample, X, e, β̂, Σ̂, R, q, 0.025, twotails=false)[3]) == Float64
    @test WildBootstrap(bootsample, X, e, β̂, Σ̂, R, q, 0.05, twotails=true)[2] > WildBootstrap(bootsample, X, e, β̂, Σ̂, R, q, 0.05, twotails=true)[3]

    # testing procedure test
    @test length(semiparamFGLS(X, y, bw_NW, R, q, constant=true, nonparm=NW, KernelDensity=Normal(0,1), twotails=true, α=0.05, numResample=bootsample)[1]) == length(β)
    @test size(semiparamFGLS(X, y, bw_KNN, R, q, constant=true, nonparm=KNN, bootstrap=false)[2]) == (length(β), length(β))
    @test size(semiparamFGLS(X, y, bw_LL, R, q, constant=true, nonparm=LL, KernelDensity=Epanechnikov(0,1), twotails=false, α=0.05, numResample=bootsample)[3]) == (length(q), length(q))
    @test length(semiparamFGLS(X, y, 3, R, q, constant=true, nonparm=SR, CrossValidation=false, twotails=false, α=0.05, numResample=bootsample)[4]) == bootsample
    @test typeof(semiparamFGLS(X2, y2, bw_NW, R2, q, constant=false, nonparm=NW, KernelDensity=Normal(0,1), restricted=false, twotails=true, α=0.05, numResample=bootsample)[5]) == Float64
    @test typeof(semiparamFGLS(X2, y2, bw_KNN, R2, q, constant=false, nonparm=KNN, restricted=false, twotails=true, α=0.05, numResample=bootsample)[6]) == Float64

end
