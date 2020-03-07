function CalcDiffItem(old::ItemParameters, new::ItemParameters)
    β = old.β .- new.β
    ζ = old.ζ .- new.ζ
    λ = old.λ .- new.λ
end

function CalcDiffPerson(old::PersonParameters, new::PersonParameters)
    #
end
