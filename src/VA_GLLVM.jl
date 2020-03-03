using Optim, Distributions, StatsFuns, StatsBase, CategoricalArrays, Random, OrdinalMultinomialModels, GLM, Statistics, LinearAlgebra

# ]add https://github.com/OpenMendel/OrdinalMultinomialModels.jl

# Read test dataset
include("src/ReadData.jl")

struct VAGradedItem{T1}
   β₀::AbstractArray{T1}
   ζ::AbstractArray{T1}
   λ::AbstractArray{T1}
   β::AbstractArray{T1}
end

struct VAGradedPerson{T1}
    a::AbstractArray{T1}
    A::AbstractArray{T1}
end
struct VAGraded{T1} <: DiscreteUnivariateDistribution
    ζ::AbstractArray{T1}
    η::AbstractArray{T1}
end


