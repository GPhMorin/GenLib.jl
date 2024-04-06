"""
    struct PossibleFounder <: AbstractIndividual
        ID::Int64
        father::Union{Nothing, PossibleFounder}
        mother::Union{Nothing, PossibleFounder}
        children::Vector{PossibleFounder}
        sex::Int64
        rank::Int64
        index::Int64
    end

An individual with an index to access the founder's kinships ∈ the Ψ matrix.
"""
struct PossibleFounder <: AbstractIndividual
    ID::Int64
    father::Union{Nothing, PossibleFounder}
    mother::Union{Nothing, PossibleFounder}
    children::Vector{PossibleFounder}
    sex::Int64
    rank::Int64
    index::Int64
end

"""
    phi(individualᵢ::T, individualⱼ::T) where T <: AbstractIndividual

Return the kinship coefficient between two individuals.

Adapted from [Karigl, 1981](@ref).

# Example

```julia
import GenLib as gen
genea140 = gen.genea140
ped = gen.genealogy(genea140)
pro1 = ped[10033]
pro2 = ped[113470]
gen.phi(pro1, pro2)
```
"""
function phi(individualᵢ::T, individualⱼ::T) where T <: AbstractIndividual
    value = 0.
    if individualᵢ.rank > individualⱼ.rank # From the genealogical order, i cannot be an ancestor of j
        # Φᵢⱼ = (Φₚⱼ + Φₘⱼ) / 2, if i is not an ancestor of j (Karigl, 1981)
        if !isnothing(individualᵢ.father)
            value += phi(individualᵢ.father, individualⱼ) / 2
        end
        if !isnothing(individualᵢ.mother)
            value += phi(individualᵢ.mother, individualⱼ) / 2
        end
    elseif individualⱼ.rank > individualᵢ.rank # Reverse the order since a > b
        # Φⱼᵢ = (Φₚⱼ + Φₘⱼ) / 2, if j is not an ancestor of i (Karigl, 1981)
        if !isnothing(individualⱼ.father)
            value += phi(individualⱼ.father, individualᵢ) / 2
        end
        if !isnothing(individualⱼ.mother)
            value += phi(individualⱼ.mother, individualᵢ) / 2
        end
    elseif individualᵢ.rank == individualⱼ.rank # Same individual
        # Φₐₐ = (1 + Φₚₘ) / 2 (Karigl, 1981)
        value += 1/2
        if !isnothing(individualᵢ.father) & !isnothing(individualᵢ.mother)
            value += phi(individualᵢ.father, individualᵢ.mother) / 2
        end
    end
    return value
end

"""
    phi(individualᵢ::PossibleFounder, individualⱼ::PossibleFounder, Ψ::Matrix{Float64})

Return the kinship coefficient between two individuals given a matrix of the founders' kinships.

Adapted from [Karigl, 1981](@ref), and [Kirkpatrick et al., 2019](@ref).
"""
function phi(individualᵢ::PossibleFounder, individualⱼ::PossibleFounder, Ψ::Matrix{Float64})
    if individualᵢ.index != 0 && individualⱼ.index != 0 # They are both founders
        return Ψ[individualᵢ.index, individualⱼ.index]
    end
    value = 0.
    if individualᵢ.rank == individualⱼ.rank # Same individual
        # Φₐₐ = (1 + Φₚₘ) / 2 (Karigl, 1981)
        value += 1/2
        if !isnothing(individualᵢ.father) & !isnothing(individualᵢ.mother)
            value += phi(individualᵢ.father, individualᵢ.mother, Ψ) / 2
        end
    else
        valueᵢ = 0.
        if !isnothing(individualᵢ.father)
            valueᵢ += phi(individualᵢ.father, individualⱼ, Ψ) / 2
        end
        if !isnothing(individualᵢ.mother)
            valueᵢ += phi(individualᵢ.mother, individualⱼ, Ψ) / 2
        end
        valueⱼ = 0.
        if !isnothing(individualⱼ.father)
            valueⱼ += phi(individualⱼ.father, individualᵢ, Ψ) / 2
        end
        if !isnothing(individualⱼ.mother)
            valueⱼ += phi(individualⱼ.mother, individualᵢ, Ψ) / 2
        end
        value = max(valueᵢ, valueⱼ)
    end
    return value
end

"""
    function _lowest_founders(pedigree::Pedigree{T}) where T <: AbstractIndividual

Return the lowest founders of a given pedigree.

The lowest founder is defined as an only child who is either a
founder or as the only child of a lineage of only children.
"""
function _lowest_founders(pedigree::Pedigree{T}) where T <: AbstractIndividual
    founderIDs = founder(pedigree)
    deepestIDs = Int64[]
    for founderID ∈ founderIDs
        deepest = pedigree[founderID]
        while length(deepest.children) == 1
            deepest = deepest.children[1]
        end
        push!(deepestIDs, deepest.ID)
    end
    sort(unique(deepestIDs))
end

"""
    _previous_generation(pedigree::Pedigree{T}, next_generation::Vector{Int64}) where T <: AbstractIndividual

Return the previous generation of a given set of individuals.
"""
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

"""
    phi(pedigree::Pedigree, Ψ::Matrix{Float64}, topIDs::Vector{Int64}, bottomIDs::Vector{Int64})

Return a square matrix of pairwise kinship coefficients
between all probands given the founders' kinships.

An implementation of the recursive-cut algorithm presented ∈ [Kirkpatrick et al., 2019](@ref).
"""
function phi(pedigree::Pedigree, Ψ::Matrix{Float64}, topIDs::Vector{Int64}, bottomIDs::Vector{Int64})
    Φ = ones(length(pedigree), length(pedigree)) .* -1
    indices = [pedigree[ID].rank for ID ∈ topIDs]
    Φ[indices, indices] = copy(Ψ)
    for individualᵢ ∈ values(pedigree)
        i = individualᵢ.rank
        for individualⱼ ∈ values(pedigree)
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

"""
    phi(pedigree::Pedigree, probandIDs::Vector{Int64} = pro(pedigree); MT::Bool = false, verbose::Bool = false)

Return a square matrix of pairwise kinship coefficients between probands.

If no probands are given, return the square matrix
for all probands ∈ the pedigree.

If MT = false: an implementation of the recursive-cut algorithm
presented ∈ [Kirkpatrick et al., 2019](@ref).

If MT = true: pairwise kinships ∈ parallel, like ∈ GENLIB.

# Example

```julia
import GenLib as gen
geneaJi = gen.geneaJi
ped = gen.genealogy(geneaJi)
gen.phi(ped)
```
"""
function phi(pedigree::Pedigree, probandIDs::Vector{Int64} = pro(pedigree); MT::Bool = false, verbose::Bool = false)
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
        if verbose
            println("Step $i / $(length(cut_vertices)-1): $(length(previous_generation)) founders, $(length(next_generation)) probands.")
        end
    end
    for i ∈ 1:length(cut_vertices)-1
        previous_generation = cut_vertices[i]
        next_generation = cut_vertices[i+1]
        if verbose
            println("Running step $i / $(length(cut_vertices)-1) ($(length(previous_generation)) founders, $(length(next_generation)) probands)")
        end
        if i == 1
            Ψ = zeros(length(previous_generation), length(previous_generation))
            for j ∈ eachindex(previous_generation)
                Ψ[j, j] = 0.5
            end
        end
        if MT
            temporary_pedigree = branching(pedigree, pro = next_generation, ancestors = previous_generation)
            isolated_pedigree = Pedigree{PossibleFounder}()
            for individual ∈ values(temporary_pedigree)
                isolated_pedigree[individual.ID] = PossibleFounder(
                    individual.ID,
                    !isnothing(individual.father) ? isolated_pedigree[individual.father.ID] : nothing,
                    !isnothing(individual.mother) ? isolated_pedigree[individual.mother.ID] : nothing,
                    [],
                    individual.sex,
                    individual.rank,
                    individual.ID ∈ previous_generation ? findfirst(previous_generation .== individual.ID) : 0
                )
            end
            for individual ∈ values(isolated_pedigree)
                if !isnothing(individual.father)
                    push!(individual.father.children, individual)
                end
                if !isnothing(individual.mother)
                    push!(individual.mother.children, individual)
                end
            end
            Φ = Matrix{Float64}(undef, length(next_generation), length(next_generation))
            probands = [isolated_pedigree[ID] for ID ∈ next_generation]
            Threads.@threads for i ∈ eachindex(probands)
                Threads.@threads for j ∈ eachindex(probands)
                    if i ≤ j
                        Φ[i, j] = Φ[j, i] = phi(probands[i], probands[j], Ψ)
                    end
                end
            end
            Ψ = Φ
        else
            isolated_pedigree = branching(pedigree, pro = next_generation, ancestors = previous_generation)
            Ψ = phi(isolated_pedigree, Ψ, previous_generation, next_generation)
        end
    end
    Ψ
end

"""
    phi(pedigree::Pedigree, rowIDs::Vector{Int64}, columnIDs::Vector{Int64})

Return a rectangle matrix of kinship coefficients,
as defined by a list or row IDs and column IDs.
"""
function phi(pedigree::Pedigree, rowIDs::Vector{Int64}, columnIDs::Vector{Int64})
    Φ = Matrix{Float64}(undef, length(rowIDs), length(columnIDs))
    Threads.@threads for i ∈ eachindex(rowIDs)
        individualᵢ = pedigree[rowIDs[i]]
        Threads.@threads for j ∈ eachindex(columnIDs)
            individualⱼ = pedigree[columnIDs[j]]
            Φ[i, j] = phi(individualᵢ, individualⱼ)
        end
    end
    Φ
end

"""
    _cleanup(kinship_matrix::Matrix{Float64, Float64}, threshold::Float64 = 0.0625)

Return a vector of booleans of probands whose kinships
between them never exceed a given `threshold`.
    
For instance, a threshold of 0.0625 (the default) removes individuals
who are first-degree cousins or closer.
"""
function _cleanup(kinship_matrix::Matrix{Float64}, threshold::Float64 = 0.0625)
    visited = [false for i ∈ 1:size(kinship_matrix, 1)]
    to_keep = [true for i ∈ 1:size(kinship_matrix, 1)]
    for i ∈ eachindex(to_keep)
        row = kinship_matrix[i, :]
        row[i] = 0
        row[visited] .= 0
        row[.!to_keep] .= 0
        to_reject = findall(row .≥ threshold)
        to_keep[to_reject] .= false
        visited[i] = true
    end
    to_keep
end

"""
    f(pedigree::Pedigree, ID::Int64)

Return the coefficient of inbreeding of an individual.
"""

function f(pedigree::Pedigree, ID::Int64)
    individual = pedigree[ID]
    if isnothing(individual.father) || isnothing(individual.mother)
        0.
    else
        phi(individual.father, individual.mother)
    end
end
