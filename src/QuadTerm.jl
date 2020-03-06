QuadTermSub(μ::AbstractArray{Float64, 1}, Σ::AbstractArray{Float64, 2}) = logdet(Σ) - tr(Σ) - μ'μ

function _QuadTerm(λ::AbstractArray{Float64, 2}, 
                   Σ::AbstractArray{Float64, 3}, 
                   μ::AbstractArray{Float64, 2}
                   )
    N = size(Σ, 1) # N of subjects
    J = size(λ, 1)
    # unstructure
    term1 = zero(Float64)
    term2 = zero(Float64)
    for i in 1:N
        for j in 1:J
            term1 += λ[j, :]' * Σ[i, :, :] * λ[j, :]
        end
        term2 += QuadTermSub(μ[i, :], Σ[i, :, :])
    end
    return -0.5term1 + 0.5term2
end

function _QuadTerm(λ::AbstractArray{Float64, 2}, Σ::AbstractArray{Float64, 2}, μ::AbstractArray{Float64, 2})
    term1 = λ' * Σ * λ
    term2 = zero(Float64)
    for i in 1:size(μ, 2)
        term2 += sum(log(Σ[i,:])) - sum(diag(Σ[i,:])) -sum(μ[i,:] .^2)
    end
    return -0.5term1 + 0.5term2
end

# For λ optimization
function _QuadTerm(λ::AbstractArray{Float64, 1}, 
                   Σ::AbstractArray{Float64, 3}, 
                   μ::AbstractArray{Float64, 2}
                   )
    N = size(Σ, 1) # N of subjects
    # unstructure
    term1 = zero(Float64)
    term2 = zero(Float64)
    for i in 1:N
        term1 += λ' * Σ[i, :, :] * λ
        term2 += QuadTermSub(μ[i, :], Σ[i, :, :])
    end
    return -0.5term1 + 0.5term2
end

function _QuadTerm(λ, 
                   Σ::AbstractArray{Float64, 3}, 
                   μ::AbstractArray{Float64, 2}
                   )
    N = size(Σ, 1) # N of subjects
    # unstructure
    term1 = zero(Float64)
    term2 = zero(Float64)
    for i in 1:N
        term1 += λ' * Σ[i, :, :] * λ
        term2 += QuadTermSub(μ[i, :], Σ[i, :, :])
    end
    return -0.5term1 + 0.5term2
end

function QuadTerm(i::ItemParameters, p::PersonParameters)
    _QuadTerm(i.λ, p.Σ, p.μ[:,:])
end