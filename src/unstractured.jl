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
    sv = StartingValues(model.y, model.d; X = model.X)
    # Initialize
    
    while(true)
        # update β
        
        # update λ

        # update τ and μ

        # update Σ(closed form)
        
        if(iter == model.MaxIter || diff < ϵ) break
    end

end

