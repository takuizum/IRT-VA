using TableReader, DataFrames
using Plots

include.("src/" .* ["LogLikelihood.jl", "StartingValue.jl", "unstractured.jl", "CalcEta.jl", "QuadTerm.jl"]);
include.(["LogLikelihood.jl", "StartingValue.jl", "unstractured.jl", "CalcEta.jl", "QuadTerm.jl"]);

grat = readcsv("data/grat.csv")
rename!(grat, :UNNAMED_1 => :ID)
y = grat[!, Not(:ID)]
X = nothing
# CreateModel
model = GradedModel(y; d = 3);

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

# Update the both of parameters
mod = GradedModel(grat[!, Not(:ID)]; MaxIter = 50, d = 3);
item, person = gradedVA(mod);

item.β₀
item.λ .* 1.702
item.ζ

CSV.write("data/vs_zeta.csv", DataFrame(item.ζ))
CSV.write("data/vs_lambda.csv", DataFrame(item.λ))
CSV.write("data/vs_beta0.csv", DataFrame(β = item.β₀))

histogram(person.μ; nbins = 50)
first(person.μ)
using CSV
CSV.write("data/va_theta.csv", DataFrame(person.μ))
CSV.write("data/va_theta_init.csv", DataFrame(oldPerson.μ))
person.μ
person.Σ[1, :, :]

map(i -> zscore(person.μ[:,i]), 1:size(person.μ, 2))


# Code optimization
grat2 = grat[!, Not(:ID)];
@code_warntype GradedModel(grat2; MaxIter = 2, d = 3);
mod = GradedModel(grat2; MaxIter = 10, d = 3);
@code_warntype gradedVA(mod);