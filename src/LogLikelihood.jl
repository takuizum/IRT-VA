using Optim, Distributions, StatsFuns, StatsBase, CategoricalArrays, Random, OrdinalMultinomialModels, GLM, Statistics, LinearAlgebra

# function Distributions.logpdf(d::VAGraded, i::Int64, j::Int64, k::Int64)
#     log(Distributions.cdf(Normal(0, 1), d.ζ[j, k] - d.η[i, j]) - Distributions.cdf(Normal(0, 1), d.ζ[j, k+1] - d.η[i, j]))
# end

# function Distributions.pdf(d::VAGraded, i::Int64, j::Int64, k::Int64)
#     Distributions.cdf(Normal(0, 1), d.ζ[j, k] - d.η[i, j]) - Distributions.cdf(Normal(0, 1), d.ζ[j, k+1] - d.η[i, j])
# end

# Samejima Graded pmf with GLLVM parameterization.
# How to transform?
pmf(η, ζ, k) = cdf(Normal(0,1), ζ[k+1] - η) - cdf(Normal(0,1), ζ[k] - η)
∂pmf(η, ζ, k) = (pdf(Normal(0, 1), ζ[k+1] - η) - pdf(Normal(0, 1), ζ[k] - η)) / pmf(η, ζ, k)

function loglikelihood(Model, Item, Person)
    N, J = size(Model.y)
    η = CalcEta(Item, Person, Model.X)
    # Calc Loglikelihood
    lnp = zero(Float64)
    for j in 1:J, i in 1:N
        lnp += log(pmf(η[i, j], Item.ζ[j,:], Model.y[i, j]))
    end
    # Calc quad term
    lnp += _QuadTerm(Item.λ, Person.Σ, Person.μ)
    return lnp
end

function loglikelihood_β(τ, β₀, β, λ, ζ, μ, y, X)
    K = size(ζ)
    N = size(μ, 1)
    # Calc eta
    η = _CalcEta(τ, β₀, β, λ, ζ, μ[:,:], X)
    # Calc Loglikelihood
    lnp = 0.0
    for i in 1:N
        lnp += log(pmf(η[i], ζ, y[i]))
    end
    return -lnp # minimize
end

function loglikelihood_λ(τ, β₀, β, λ, ζ, μ, Σ, y, X)
    K = size(ζ)
    N = size(y, 1)
    # Calc eta
    η = _CalcEta(τ, β₀, β, λ, ζ, μ[:,:], X)    # Calc Loglikelihood
    lnp = 0.0
    for i in 1:N
        lnp += log(pmf(η[i], ζ, y[i]))
    end
    # Calc quad term
    lnp += _QuadTerm(λ, Σ, μ[:,:])
    return -lnp # minimize
end

function loglikelihood_ζ(ζ, η, y)
    N = size(y, 1)
    ζ = [-Inf; ζ; Inf]
    lnp = zero(eltype(η))
    for i in 1:N
        lnp += log(pmf(η[i], ζ, y[i]))
    end
    return -lnp
end

function ∂loglikelihood_ζ(ζ, η, y)
    N = size(y, 1)
    ζ = [-Inf; ζ; Inf]
    ∂lnp = zeros(eltype(η), length(ζ)+2)
    for i in 1:N
        ∂lnp[y[i]+1] += ∂pmf(η[i], ζ, y[i])
        ∂lnp[y[i]] -= ∂pmf(η[i], ζ, y[i]-1)
    end
    return -∂lnp[2:end-1]
end

function loglikelihood_μ(τ, β₀, β, λ, ζ, μ, Σ, y, X)
    J = size(β₀, 1)
    η = _CalcEta(τ, β₀, β, λ, ζ, μ, X)
    # Calc Loglikelihood
    lnp = 0.0
    for j in 1:J
        lnp += log(pmf(η[j], ζ[j,:], y[j]))
    end
    # Calc quad term
    lnp += 0.5*(logdet(Σ) - tr(Σ) - μ'μ) - (0.5μ'Σ*μ)[1]
    return -lnp # minimize
end
