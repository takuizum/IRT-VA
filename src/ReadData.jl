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

# UpdateModel
@time tmp = UpdateModelParameters(newItem, newPerson, η, model; debug = true)
tmp.β₀
newItem.β₀
tmp.ζ
oldItem.ζ

# Update variational parameters
UpdateVariationalParameters(newItem, newPerson, model)
oldPerson.Σ[1,:,:]
newPerson.Σ[1,:,:]
oldPerson.μ[1]
newPerson.μ[1]
