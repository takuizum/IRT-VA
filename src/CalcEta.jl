# Calcuration η matrix
function _CalcEta(τ, β₀, β, λ, ζ::AbstractArray{Float64, 2}, μ::AbstractArray{Float64, 2}, X)
    J = size(ζ, 1)
    N = size(μ, 1)
    # Calc eta
    η = isnothing(X) ? λ*μ' : β*X' + λ*μ'
    for j in 1:J, i in 1:N
        η[j, i] += τ[i] + β₀[j]
    end
    return η'
end

# In item j
function _CalcEta(τ, β₀, β, λ, ζ::AbstractArray{Float64, 1}, μ::AbstractArray{Float64, 2}, X)
    N = size(μ, 1)
    # Calc eta
    η = isnothing(X) ? λ*μ' : β*X' + λ*μ'
    for i in 1:N
        η[i] += τ[i] + β₀
    end
    return η'
end

function _CalcEta(τ, β₀, β, λ, ζ::AbstractArray{Float64, 2}, μ, X)
    J = size(ζ, 1)
    # Calc eta
    η = isnothing(X) ? λ*μ' : β*X' + λ*μ'
    for j in 1:J
        η[j] += τ + β₀[j]
    end
    return η
end

function CalcEta(i::ItemParameters, p::PersonParameters, X)
    _CalcEta(p.τ, i.β₀, i.β, i.λ, i.ζ, p.μ[:,:], X)
end