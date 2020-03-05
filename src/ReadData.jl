using TableReader, DataFrames
include.("src/" .* ["LogLikelihood.jl", "StartingValue.jl", "unstractured.jl"])

grat = readcsv("data/grat.csv")
rename!(grat, :UNNAMED_1 => :ID)
y = grat[!, Not(:ID)]
X = nothing
# CreateModel
model = GradedModel(y)

# Initial values
τ, β₀, β, ζ, λ, μ, Σ = CalcStartingValues(y, model.d);

# QuadTerm
qt = QuadTerm(λ, Σ, μ)
qt = QuadTermⱼ(λ[j,:], Σ, μ)

# η matrix
η = CalcEta(τ, β₀, β, λ, ζ, μ[:,:], X)

# optimize step
j = 1
using Plots
# β₀
CalcEta(τ, β₀[j], β[j,:], λ[j,:], ζ[j,:], μ[:,:], X)
loglikelihood_β(τ, β₀[j], β[j,:], λ[j,:], ζ[j,:], μ, Σ, y[:,j], X)
plot([-4:0.1:4], map(x -> loglikelihood_β(τ, x, β[j,:], λ[j,:], ζ[j,:], μ, Σ, y[:,j], X), [-4:0.1:4;]))
res_β₀ = Optim.optimize(x -> loglikelihood_β(τ, x, β[j,:], λ[j,:], ζ[j,:], μ, y[:,j], X, Σ), [β₀[j]], Newton()) # More efficient if the cost of evalueation of optimize function is high and derivatives are easy to calculate.
res_β₀ = Optim.optimize(x -> loglikelihood_β(τ, x, β[j,:], λ[j,:], ζ[j,:], μ, y[:,j], X, Σ), -10, 10, Brent())
res_β₀.minimizer

# λ
plot([0:0.1:4], map(x -> loglikelihood_λ(τ, β₀[j], β[j,:], [x], ζ[j,:], μ, Σ, y[:,j], X), [0:0.1:4;]))
res_λ = Optim.optimize(x -> loglikelihood_λ(τ, β₀[j], β[j,:], x, ζ[j,:], μ, Σ, y[:,j], X), λ[j, :], LBFGS())
res_λ.minimizer

# ζ
η = CalcEta(τ, β₀, β, λ, ζ, μ[:,:], X);
_ζ = ζ[j,:]
loglikelihood_ζ(_ζ[-Inf .< _ζ .< Inf], η[:,j], y[:,j])
res_ζ = Optim.optimize(x -> loglikelihood_ζ(x, η[:,j], y[:,j]), _ζ[-Inf .< _ζ .< Inf], LBFGS())
# It is necessary for some proper constraints.