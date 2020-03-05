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


res_λ = optimize(x -> loglikelihoodⱼ(τ, β₀[j], β[j,:], x, ζ[j,:], μ, y[:,j], X, qt), λ[j, :], BFGS(); autodiff=:forward)
# optimize step
j = 6
using Plots
# β₀
loglikelihoodⱼ(τ, β₀[j], β[j,:], λ[j,:], ζ[j,:], μ, y[:,j], X, Σ)
plot([-4:0.1:4], map(x -> loglikelihoodⱼ(τ, x, β[j,:], λ[j,:], ζ[j,:], μ, y[:,j], X, Σ), [-4:0.1:4;]))
res_β₀ = Optim.optimize(x -> loglikelihoodⱼ(τ, x, β[j,:], λ[j,:], ζ[j,:], μ, y[:,j], X, Σ), [β₀[j]], Newton()) # More efficient if the cost of evalueation of optimize function is high and derivatives are easy to calculate.
res_β₀ = Optim.optimize(x -> loglikelihoodⱼ(τ, x, β[j,:], λ[j,:], ζ[j,:], μ, y[:,j], X, Σ), -10, 10, Brent())
res_β₀.minimizer

# λ
res_λ = Optim.optimize(x -> loglikelihoodⱼ(τ, β₀[j], β[j,:], x, ζ[j,:], μ, y[:,j], X, Σ), λ[j, :], BFGS())
res_λ.minimizer

# ζ
_ζ = ζ[j,:]
up = fill(Inf, 3)
lp = -Inf
res_ζ = Optim.optimize(x -> loglikelihoodⱼ(τ, β₀[j], β[j,:], λ[j, :], [lp; x; up], μ, y[:,j], X, Σ), _ζ[-Inf .< _ζ .< Inf], LBFGS())
# It is necessary for some proper constraints.