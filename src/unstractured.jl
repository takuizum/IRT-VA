abstract type unstractured end

struct GradedModel <: unstractured
    y
    X
    d
    MaxIter::Int64
    RowEffect::Bool
    ϵ
end

function gradedVA(model::unstractured)
    #Starting values
    τ, β₀, β, ζ, λ, μ, Σ = StartingValues(model.y, model.d; X = model.X)
    # Initialize
    N, J = size(model.y)
    y = model.y
    X = model.X
    while(true)
        qt = QuadTerm(λ, Σ, μ)
        for j in 1:J
            # update β₀
            res_β₀ = optimize(x -> loglikelihoodⱼ(τ, x, β[j,:], λ[j,:], ζ[j,:], μ, y[:,j], X, qt), β₀[j], Brent(); autodiff=:forward)
            # update β
            if !isnothing(X)
                res_β = optimize(x -> loglikelihoodⱼ(τ, β₀[j], x, λ[j,:], ζ[j,:], μ, y[:,j], X, qt), β[j, :], BFGS(); autodiff=:forward)
            end
            # update ζ (impose proper constraints not to reverse the adjacent category boundary parameters.)
            # In Baker & Kim (2004), all of the ζ were estimated concurrently 
            # update λ
            res_λ = optimize(x -> loglikelihoodⱼ(τ, β₀[j], β[j,:], x, ζ[j,:], μ, y[:,j], X, qt), λ[j, :], BFGS(); autodiff=:forward)
        end
        # update τ and μ
        res_person = optimize(x -> loglikelihood(x, β₀, β, λ, ζ, μ, y, X), τ, LBFGS())

        # update Σ(closed form)
        
        if(model.iter == model.MaxIter || diff < model.ϵ) break
    end

end

