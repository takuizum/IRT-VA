using Optim, Distributions, StatsFuns, StatsBase, CategoricalArrays, Random, OrdinalMultinomialModels, GLM, Statistics, LinearAlgebra

function Distributions.logpdf(d::VAGraded, i::Int64, j::Int64, k::Int64)
    log(Distributions.cdf(Normal(0, 1), d.ζ[j, k] - d.η[i, j]) - Distributions.cdf(Normal(0, 1), d.ζ[j, k+1] - d.η[i, j]))
end

function Distributions.pdf(d::VAGraded, i::Int64, j::Int64, k::Int64)
    Distributions.cdf(Normal(0, 1), d.ζ[j, k] - d.η[i, j]) - Distributions.cdf(Normal(0, 1), d.ζ[j, k+1] - d.η[i, j])
end

# Samejima Graded pmf with GLLVM parameterization.
# How to transform?
pmf(η, ζ, k) = cdf(Normal(0,1), ζ[k] - η) - cdf(Normal(0,1), ζ[k+1] - η)

function loglikelihood(τ, β₀, β, λ, ζ, μ, Σ, y, X)
    J, K = size(ζ)
    N = size(μ, 1)
    # Calc eta
    η = isnothing(X) ? μ * λ : X * β + μ * λ
    for i in 1:N, j in 1:J
        η[i, j] += τ[i] + β₀[j]
    end
    # Calc Loglikelihood
    lnp = 0.0
    for i in 1:N, j in 1:J
        lnp += log(pmf(η[i, j], ζ[j,:], y[i, j]))
    end
    # Calc quad term
    lnp += QuadTerm(λ, Σ, μ)
    return lnp
end


function loglikelihoodⱼ(τ, β₀, β, λ, ζ, μ, y, X, qt)
    J, K = size(ζ)
    N = size(μ, 1)
    # Calc eta
    η = isnothing(X) ? μ * λ : X * β + μ * λ
    for i in 1:N
        η[i] += τ[i] + β₀
    end
    # Calc Loglikelihood
    lnp = 0.0
    for i in 1:N
        lnp += log(pmf(η[i], ζ, y[i]))
    end
    # Calc quad term
    # lnp += QuadTerm(λ, Σ, μ)
    lnp += qt
    return lnp
end

# Calculation 
QuadTermSub(μ::AbstractArray{Float64, 1}, Σ::AbstractArray{Float64, 2}) = logdet(Σ) - tr(Σ) - μ'μ

function QuadTerm(λ::AbstractArray{Float64, 2}, Σ::AbstractArray{Float64, 3}, μ::AbstractArray{Float64, 2})
    N = size(Σ, 1) # N of subjects
    J = size(λ, 1)
    # unstructure
    term1 = zeros(Float64)
    term2 = zeros(Float64)
    res = zeros(Float64)
    for i in 1:N
        for j in 1:J
            term1 += λ[j, :]' * Σ[i, :, :] * λ[j, :] * 0.5
        end
        term2 += QuadTermSub(μ[i, :], Σ[i, :, :])
    end
    return -0.5term1 + 0.5term2
end

function QuadTerm(λ::AbstractArray{Float64, 2}, Σ::AbstractArray{Float64, 2}, μ::AbstractArray{Float64, 2})
    term1 = λ' * Σ * λ
    term2 = zero(Float64)
    for i in 1:size(μ, 2)
        term2 += sum(log(Σ[i,:])) - sum(diag(Σ[i,:])) -sum(μ[i,:] .^2)
    end
    return -0.5term1 + 0.5term2
end


