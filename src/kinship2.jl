function _previous_generation(pedigree::Pedigree{T}, next_generation::Vector{Int64}) where T <: AbstractIndividual
    previous_generation = Int64[]
    for ID ∈ next_generation
        individual = pedigree[ID]
        father = individual.father
        mother = individual.mother
        if isnothing(father) && isnothing(mother)
            push!(previous_generation, ID)
        else
            if !isnothing(father)
                push!(previous_generation, father.ID)
            end
            if !isnothing(mother)
                push!(previous_generation, mother.ID)
            end
        end
    end
    sort(unique(previous_generation))
end

function phi(pedigree::Pedigree, Ψ::Matrix{Float64}, topIDs::Vector{Int64}, bottomIDs::Vector{Int64})
    Φ = ones(length(pedigree), length(pedigree)) .* -1
    indices = [pedigree[ID].rank for ID ∈ topIDs]
    Φ[indices, indices] = copy(Ψ)
    for individualᵢ in values(pedigree)
        i = individualᵢ.rank
        for individualⱼ in values(pedigree)
            j = individualⱼ.rank
            if Φ[i, j] > -1
                continue
            elseif i == j
                father = individualᵢ.father
                mother = individualᵢ.mother
                if !isnothing(father) && !isnothing(mother)
                    coefficient = (1 + Φ[mother.rank, father.rank]) / 2
                    Φ[i, i] = coefficient
                else
                    coefficient = 0.5
                    Φ[i, i] = coefficient
                end
            elseif i < j
                fatherᵢ = individualᵢ.father
                motherᵢ = individualᵢ.mother
                coefficientᵢ = 0.
                if !isnothing(fatherᵢ)
                    coefficientᵢ += Φ[fatherᵢ.rank, individualⱼ.rank] / 2
                end
                if !isnothing(motherᵢ)
                    coefficientᵢ += Φ[motherᵢ.rank, individualⱼ.rank] / 2
                end
                fatherⱼ = individualⱼ.father
                motherⱼ = individualⱼ.mother
                coefficientⱼ = 0.
                if !isnothing(fatherⱼ)
                    coefficientⱼ += Φ[fatherⱼ.rank, individualᵢ.rank] / 2
                end
                if !isnothing(motherⱼ)
                    coefficientⱼ += Φ[motherⱼ.rank, individualᵢ.rank] / 2
                end
                coefficient = max(coefficientᵢ, coefficientⱼ)
                Φ[i, j] = coefficient
                Φ[j, i] = coefficient
            end
        end
    end
    indices = [pedigree[ID].rank for ID ∈ bottomIDs]
    Φ[indices, indices]
end

function recursive_cut(pedigree::Pedigree{T}, probandIDs::Vector{Int64} = pro(pedigree)) where T <: AbstractIndividual
    global Ψ, next_generation
    isolated_pedigree = branching(pedigree, pro = probandIDs)
    cut_vertices = [probandIDs]
    founderIDs = founder(isolated_pedigree)
    previous_generation = probandIDs
    while true
        previous_generation = _previous_generation(isolated_pedigree, previous_generation)
        if previous_generation != founderIDs
            pushfirst!(cut_vertices, previous_generation)
        else
            break
        end
    end
    pushfirst!(cut_vertices, founderIDs)
    for i ∈ 1:length(cut_vertices)-1
        previous_generation = cut_vertices[i]
        next_generation = cut_vertices[i+1]
        if i == 1
            Ψ = zeros(length(previous_generation), length(previous_generation))
            for j in eachindex(previous_generation)
                Ψ[j, j] = 0.5
            end
        end
        isolated_pedigree = branching(pedigree, pro = next_generation, ancestors = previous_generation)
        Ψ = phi(isolated_pedigree, Ψ, previous_generation, next_generation)
    end
    indices = [findfirst(ID .== next_generation) for ID ∈ probandIDs]
    Φ = Ψ[indices, indices]
end