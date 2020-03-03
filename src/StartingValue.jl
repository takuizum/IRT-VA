struct StartingValues
    τ
    β₀
    β
    ζ
    λ
    μ
    Σ
end

function CalcStartingValues(y, Numθ; X = nothing)
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
    μ = range(-3, stop = 3, length = len_uni)[rs]
    if Numθ > 1
        μ = hcat(μ, rand(MvNormal(ones(Numθ-1)), length(a))')
    end
    # unique_μ = μ[unique_ind, :]
    # fit GLM
    max_v = maximum(y)
    params = zeros(Float64, p, Numθ + NumCovariates + 1) # τ, β, λ
    β₀ = zeros(Float64, p)
    β = zeros(Float64, p, NumCovariates)
    λ = zeros(Float64, p, Numθ)
    ζ = zeors(Float64, p, max_v+1)
    ζ[:, 1] = -Inf
    for j in 1:p
        Kⱼ = length(unique(y[:, j]))
        if isnothing(X)
            plfit = polr(μ, y[:, j], ProbitLink())
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
    σ = diagm(Numθ, Numθ, fill(1.0, Numθ))
    Σ = fill(σ, N)
    τ = fill(0.0, N)
    return VAGradedItem(τ, β₀, β, ζ, λ, μ, Σ)
end

res = polr(a, y[:, 2], ProbitLink(), IpoptSolver()) # NLoptsolverだと計算されない。
tab = res |> coeftable
res.β # regression coefficient
res.θ # intercept, satisfying `θ[1]≤...≤θ[J-1]`
res.η # 
res.α # unconstrained parameterization of θ

# 