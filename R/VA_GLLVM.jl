using Optim, Distributions, StatsFuns, StatsBase

struct GradedLogistic{T1} <: DiscreteUnivariateDistribution
   a::T1
   b::Vector{T1}
   θ::T1
end

function Distributions.rand(s::GradedLogistic)
    D = 1.702 # Scaling constants to be appriximated probit function
    K = length(s.b)+1
	b = [-Inf; s.b; Inf]
	p = zeros(K)
	p[1] = logistic(D * s.a*(s.θ - b[1])) - logistic(D * s.a*(s.θ - b[2]))
	for k in 2:K
		p[k] = p[k-1] + logistic(D * s.a*(s.θ - b[k])) - logistic(D * s.a*(s.θ - b[k+1]))
	end
	sum(p .< rand()) + 1
end

function Distributions.logpdf(d::GradedLogistic, k::Int)
    D = 1.702 # Scaling constants to be appriximated probit function
	b = [-Inf; d.b; Inf]
    log(logistic(D * d.a*(d.θ - b[k])) - logistic(D * d.a*(d.θ - b[k+1])))
end

function Distributions.pdf(d::GradedLogistic, k::Int)
    D = 1.702 # Scaling constants to be appriximated probit function
	b = [-Inf; d.b; Inf]
    logistic(D * d.a*(d.θ - b[k])) - logistic(D * d.a*(d.θ - b[k+1]))
end

StatsBase.loglikelihood(d::Vector{GradedLogistic{Float64}}, X) = sum(logpdf.(d, X))

testPara = GradedLogistic(1.0, [-1.0, 0.0, 1.0], 0.0)
pdf(testPara, 2)

testVec = [
    GradedLogistic(1.0, [-1.0, 0.0, 1.0], 0.0), 
    GradedLogistic(2.0, [1.0, 1.5, 2.0], 0.0), 
    GradedLogistic(0.5, [-2.0, 0.0, 2.0], 0.0)
]

logpdf.(testVec, [1, 2, 1])
loglikelihood(testVec, [1, 2, 1])

