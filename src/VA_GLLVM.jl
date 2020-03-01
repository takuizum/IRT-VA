using Optim, Distributions, StatsFuns, StatsBase, CategoricalArrays, Random, OrdinalMultinomialModels, GLM, Statistics

# ]add https://github.com/OpenMendel/OrdinalMultinomialModels.jl

struct VAGradedItem{T1} <: DiscreteUnivariateDistribution
   β₀::T1
   ζ::Vector{T1}
   λ::AbstractArray{T1}
   β::AbstractArray{T1}
end

struct VAGradedPerson{T1}
    a::AbstractArray{T1}
    A::AbstractArray{T1}
end

# function Distributions.rand(s::VAGraded)
#     K = length(s.ξ)+1
# 	ζ = [-Inf; s.ζ; Inf]
# 	p = zeros(K)
# 	p[1] = Distributions.cdf(Normal(0, 1), (ζ[1] - s.η)) - Distributions.cdf(Normal(0, 1), ζ[2] - s.η)
# 	for k in 2:K
# 		p[k] = Distributions.cdf(Normal(0, 1), (ζ[k] - s.η)) - Distributions.cdf(Normal(0, 1), ζ[k+1] - s.η)
# 	end
# 	sum(p .< rand()) + 1
# end

function Distributions.logpdf(d::VAGraded, k::Int)
	ζ = [-Inf; d.ζ; Inf]
    log(Distributions.cdf(Normal(0, 1), ζ[k] - d.η) - Distributions.cdf(Normal(0, 1), ζ[k+1] - d.η))
end

function Distributions.pdf(d::VAGraded, k::Int)
    ζ = [-Inf; d.ζ; Inf]
    Distributions.cdf(Normal(0, 1), ζ[k] - d.η) - Distributions.cdf(Normal(0, 1), ζ[k+1] - d.η)
end

StatsBase.loglikelihood(d::Vector{VAGraded{Float64}}, X) = sum(logpdf.(d, X))

testPara = VAGraded(1.0, [-1.0, 0.0, 1.0])
pdf(testPara, 3)

testVec = [
    VAGraded(1.0, [-1.0, 0.0, 1.0], 0.0), 
    VAGraded(2.0, [1.0, 1.5, 2.0], 0.0), 
    VAGraded(0.5, [-2.0, 0.0, 2.0], 0.0)
]

logpdf.(testVec, [1, 2, 1])
loglikelihood(testVec, [1, 2, 1])

function StartingValues(y, Numθ; X = nothing)
    # y in integer matrix
    N, p = size(y)
    if isnothing(X)
        NumCovariates = 0
    else
        NumCovariates = size(X)
    end
    # Starting values of a(VA mean = mean of latent trauts)
    unique_ind = findall(x -> !x, nonunique(y)) # y is a dataframe that dosen't contain ID col
    y = convert(Matrix{Int64}, y[!, Not(:ID)])
    rs = sum(y; dims = 2)[:]
    len_uni = length(unique(rs))
    rs = levelcode.(CategoricalArray(rs, ordered = true)) # ranking of raw score
    sortperm(rs)
    a = range(-3, stop = 3, length = len_uni)[rs]
    if Numθ > 1
        a = hcat(a, rand(MvNormal(ones(Numθ-1)), length(a))')
    end
    unique_a = a[unique_ind, :]
    # fit GLM
    max_v = maximum(y)
    params = zeros(Float64, p, Numθ + NumCovariates + 1) # τ, β, λ
    β₀ = zeros(Float64, p)
    β = zeros(Float64, p, NumCovariates)
    λ = zeros(Float64, p, Numθ)
    ζ = zeors(Float64, p, max_v+1)
    ζ[:, 1] = -Inf
    # df_y = DataFrame(y)
    # rename!(df_y, Symbol.("Y" .* string.(1:p)))
    # if !isnothing(X)
    #     df_x = DataFrame(X)
    #     rename!(df_x, Symbol.("Y" .* string.(1:NumCovariates)))
    # end
    for j in 1:p
        Kⱼ = length(unique(y[:, j]))
        if isnothing(X)
            plfit = polr(a, y[:, j], ProbitLink())
            λ[j, :] = plfit.β
        else
            plfit = polr([a X], y[:, j], ProbitLink())
            λ[j, :] = plfit.β[1:Numθ]
            β[j, :] = plfit.β[(Numθ+1):end]
        end
        β₀[j] = plfit.β[1]
        λ[j, :] = plfit.β
        ζ[j, 2:Kⱼ] = plfit.θ .- mean(plfit.θ)
        ζ[j, Kⱼ:max_v+1] = Inf
    end
    return VAGradedItem(β₀, ζ, λ, β)
end

res = polr(a, y[:, 2], ProbitLink(), IpoptSolver()) # NLoptsolverだと計算されない。
tab = res |> coeftable
res.β # regression coefficient
res.θ # intercept, satisfying `θ[1]≤...≤θ[J-1]`
res.η # 
res.α # unconstrained parameterization of θ