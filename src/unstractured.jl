<<<<<<< HEAD
struct GradedModel
=======
struct GradedModel 
>>>>>>> b48efc19327ea337f04fad793bbb7609fd80ae05
    y
    X
    d
    MaxIter::Int64
    RowEffect::Bool
    ϵ
end
function GradedModel(y; d = 1, X = nothing, MaxIter = 100, RowEffect = true, ϵ = 1e-3)
    GradedModel(y, X, d, MaxIter, RowEffect, ϵ)
end

function gradedVA(Model::GradedModel)
    #Starting values
    oldItem, oldPerson = CalcStartingValues(Model.y, Model.d);
    newItem, newPerson = copy.([oldItem, oldPerson])
    # Initialize
    N, J = size(model.y)
<<<<<<< HEAD
    for iter in 1:Model.MaxIter
        η = CalcEta(newItem, newPerson, Model.X)
        # Update the item parameters
        UpdateModelParameters(newItem, newPerson, η, Model)
        # update τ and μ
        UpdateVariationalParameters(newItem, newPerson, Model)
        if(diff < Model.ϵ)
=======
    while(true)
        η = CalcEta(newItem, newPerson, Model.X)
        # Update the item parameters
        UpdateModelParameters(newItem, newPerson, η, Model)

        # update τ and μ
        res_person = optimize(x -> loglikelihood(x, β₀, β, λ, ζ, μ, y, X), τ, LBFGS())

        # update Σ(closed form)
        
        if(Model.iter == Model.MaxIter || diff < Model.ϵ) 
>>>>>>> b48efc19327ea337f04fad793bbb7609fd80ae05
            break
        end
    end

end

function UpdateModelParameters(Item, Person, η, Model; debug = false)
    if debug
        Item = copy(Item) # only use in debugging
    end
    J = size(Model.y, 2)
    for j in 1:J
        # println("Item ", j)
        # println("Update β")
        res_β₀ = Optim.optimize(x -> loglikelihood_β(Person.τ, x, Item.β[j,:], Item.λ[j,:], Item.ζ[j,:], Person.μ, Model.y[:,j], Model.X), Item.β₀[j]-1.0, Item.β₀[j]+1.0, Brent())
        Item.β₀[j] = res_β₀.minimizer[1]
        if !isnothing(Model.X)
            res_β = Optim.optimize(x -> loglikelihood_β(Person.τ, Item.β₀[j], x, Item.λ[j,:], Item.ζ[j,:], Person.μ, Model.y[:,j], Model.X), Item.β[j,:], BFGS(); autodiff = :forward)
            Item.β[j,:] = res_β.minimizer[:]
        end
        # println("Update ζ")
        # (impose proper constraints not to reverse the adjacent category boundary parameters.)
        update_location = -Inf .< Item.ζ[j,:] .< Inf
        _ζ = Item.ζ[j, update_location]
        # temporary constraints
        lower = [-Inf            ; _ζ[1:end-1] .* 0.9]
        upper = [_ζ[2:end] .* 0.9; Inf]
        try
            res_ζ = Optim.optimize(x -> loglikelihood_ζ(x, η[:,j], Model.y[:,j]), lower, upper, _ζ, Fminbox(LBFGS()); autodiff = :forward)
            Item.ζ[j, update_location] = res_ζ.minimizer[:]
        catch
            @warn "Failed to Optimize in Item $j"
        finally
            Item
        end
        # println("Update λ")
        res_λ = Optim.optimize(x -> loglikelihood_λ(Person.τ, Item.β₀[j], Item.β[j,:], x, Item.ζ[j,:], Person.μ, Person.Σ, Model.y[:,j], Model.X), Item.λ[j, :], LBFGS(); autodiff = :forward)
        Item.λ[j,:] = res_λ.minimizer[:]
    end
    return Item
end

<<<<<<< HEAD
function UpdateVariationalParameters(Item::ItemParameters, Person::PersonParameters, Model::GradedModel)
    N, J = size(Model.y)
    # λ square term
    λλ = zeros(eltype(Item.λ), Model.d, Model.d)
    for j in 1:J
        λλ += Item.λ[j,:] * Item.λ[j,:]'
    end
    # Update
    for i in 1:N
        Person.μ[i,:] = Optim.optimize(x -> loglikelihood_μ(Person.τ[i], Item.β₀, Item.β, Item.λ, Item.ζ, x, Person.Σ[i,:,:], Model.y[i,:], Model.X), Person.μ[i,:], BFGS()).minimizer
        Person.Σ[i,:,:] = inv(diagm(Model.d, Model.d, ones(Float64, Model.d)) + λλ)
    end
end
=======
function UpdateVariatinalParameters
    #
end
>>>>>>> b48efc19327ea337f04fad793bbb7609fd80ae05
