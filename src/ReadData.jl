using TableReader, DataFrames
using Plots

include.("src/" .* ["LogLikelihood.jl", "StartingValue.jl", "unstractured.jl", "CalcEta.jl", "QuadTerm.jl"])

grat = readcsv("data/grat.csv")
rename!(grat, :UNNAMED_1 => :ID)
y = grat[!, Not(:ID)]
X = nothing
# CreateModel
model = GradedModel(y)

# Initial values
oldItem, oldPerson = CalcStartingValues(model.y, model.d);
newItem, newPerson = copy.([oldItem, oldPerson])

# QuadTerm
qt = QuadTerm(oldItem, oldPerson)
qt = QuadTerm(λ[j,:], Σ, μ)

# η matrix
η = CalcEta(τ, β₀, β, λ, ζ, μ[:,:], X);

# UpdateModel
@time tmp = UpdateModelParameters(oldItem, oldPerson, η, model; debug = true)
tmp.β₀
newItem.β₀
tmp.ζ
oldItem.ζ

# Update variational parameters
i = 1
Optim.optimize(x -> loglikelihood_μ(oldPerson.τ[i], oldItem.β₀, oldItem.β, oldItem.λ, oldItem.ζ, [x], oldPerson.Σ[i,:,:], y[i,:], model.X), -10, 10, Brent())
y = Matrix(model.y)
loglikelihood_μ(oldPerson.τ[i], oldItem.β₀, oldItem.β, oldItem.λ, oldItem.ζ, [1.0], oldPerson.Σ[i,:,:], Matrix, y[i,:], model.X)