# struct StartingValues
#     τ
#     β₀
#     β
#     ζ
#     λ
#     μ
#     Σ
# end

using DataFrames, Distributions, CategoricalArrays, Random, OrdinalMultinomialModels, LinearAlgebra, Statistics

function CalcStartingValues(y, Numθ; X = nothing)
    N, p = size(y)
    if isnothing(X)
        NumCovariates = 0
    else
        NumCovariates = size(X)
    end
    # Starting values of a(VA mean = mean of latent trauts)
    unique_ind = findall(x -> !x, nonunique(y)) # y is a dataframe that dosen't contain ID col
    y = convert(Matrix{Int64}, y)
    rs = sum(y; dims = 2)[:]
    len_uni = length(unique(rs))
    rs = levelcode.(CategoricalArray(rs, ordered = true)) # ranking of raw score
    sortperm(rs)
    μ = range(-3, stop = 3, length = len_uni)[rs]
    if Numθ > 1
        μ = hcat(μ, rand(MvNormal(ones(Numθ-1)), length(μ))')
    end
    # unique_μ = μ[unique_ind, :]
    # fit GLM
    max_v = maximum(y)
    params = zeros(Float64, p, Numθ + NumCovariates + 1) # τ, β, λ
    β₀ = zeros(Float64, p)
    β = zeros(Float64, p, NumCovariates)
    λ = zeros(Float64, p, Numθ)
    ζ = zeros(Float64, p, max_v+1)
    ζ[:, 1] .= -Inf
    D = 1/1.702
    for j in 1:p
        println("Item", j)
        # Logit Link, which is used in `polr` to fit orderd regression analysis, is more stable than ProbitLink in my experience.
        Kⱼ = length(unique(y[:, j]))
        if isnothing(X)
            plfit = polr(μ, y[:, j], LogitLink())
            λ[j, :] .= plfit.β .* D
        else
            plfit = polr([μ X], y[:, j], LogitLink())
            λ[j, :] .= plfit.β[1:Numθ] .* D
            β[j, :] .= plfit.β[(Numθ+1):end] .* D
        end
        β₀[j] = plfit.β[1] .* D
        λ[j, :] .= plfit.β .* D
        ζ[j, 2:Kⱼ] .= plfit.θ .* D #  .- mean(plfit.θ)
        ζ[j, Kⱼ+1:max_v+1] .= Inf
    end
    σ = diagm(Numθ, Numθ, fill(1.0, Numθ))
    Σ = zeros(Float64, N, Numθ, Numθ)
    map(i -> Σ[i, :, :] = σ, 1:N)
    τ = fill(0.0, N)
    return τ, β₀, β, ζ, λ, μ, Σ
end

# res_l = polr(μ, y[:, 25], LogitLink()); # NLoptsolverだと計算されない。
# res_p = polr(μ, y[:, 25], ProbitLink()); # NLoptsolverだと計算されない。

# res_p |> coeftable
# res_l |> coeftable
# res_p.β  # regression coefficient
# 0.588res_l.β  # regression coefficient
# res_p.θ # intercept, satisfying `θ[1]≤...≤θ[J-1]`
# 0.588res_l.θ # intercept, satisfying `θ[1]≤...≤θ[J-1]`
# res_p.η # 
# 0.588res_l.η # 
# res_p.α # unconstrained parameterization of θ
# 0.588res_l.α # unconstrained parameterization of θ


# length(unique(y[:, 25]))