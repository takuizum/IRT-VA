abstract type diagonal end

struct GradedModelDiag <: diagonal
    y
    d
    MaxIter::Int64
    RowEffect::Bool
end

function gradedVA(model::diagonal)
    #
end