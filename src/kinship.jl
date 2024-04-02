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

An individual with an index to access the founder's kinships in the Ψ matrix.
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
    if individualᵢ.rank > individualⱼ.rank # From the genealogical order, i cannot be an ancestor of j
        # Φᵢⱼ = (Φₚⱼ + Φₘⱼ) / 2, if i is not an ancestor of j (Karigl, 1981)
        if !isnothing(individualᵢ.father)
            value += phi(individualᵢ.father, individualⱼ, Ψ) / 2
        end
        if !isnothing(individualᵢ.mother)
            value += phi(individualᵢ.mother, individualⱼ, Ψ) / 2
        end
    elseif individualⱼ.rank > individualᵢ.rank # Reverse the order since a > b
        # Φⱼᵢ = (Φₚⱼ + Φₘⱼ) / 2, if j is not an ancestor of i (Karigl, 1981)
        if !isnothing(individualⱼ.father)
            value += phi(individualⱼ.father, individualᵢ, Ψ) / 2
        end
        if !isnothing(individualⱼ.mother)
            value += phi(individualⱼ.mother, individualᵢ, Ψ) / 2
        end
    elseif individualᵢ.rank == individualⱼ.rank # Same individual
        # Φₐₐ = (1 + Φₚₘ) / 2 (Karigl, 1981)
        value += 1/2
        if !isnothing(individualᵢ.father) & !isnothing(individualᵢ.mother)
            value += phi(individualᵢ.father, individualᵢ.mother, Ψ) / 2
        end
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
    _cut_vertex(individual::T, candidateID::Int64) where T <: AbstractIndividual

Return whether an individual can be used as a cut vertex.

A cut vertex is an individual that "when removed,
disrupt every path from any source [founder]
to any sink [proband]" ([Kirkpatrick et al., 2019](@ref)).
"""
function _cut_vertex(individual::T, candidateID::Int64) where T <: AbstractIndividual
    value = true
    for child in individual.children
        # Check if going down the pedigree
        # while avoiding the candidate ID
        # reaches a proband (sink) anyway
        if isempty(child.children) # The child is a proband
            return false
        elseif child.ID != candidateID
            value = value && _cut_vertex(child, candidateID)
        end
    end
    value
end

"""
    _bottlenecks(pedigree::Pedigree{T}) where T <: AbstractIndividual
Return the bottleneck ancestors of a given pedigree.

A bottleneck ancestor is an individual who is the only transmitter of 
their ancestors' genetic contributions. In graph terms, it may be defined as
a cut vertex, i. e. an individual that "when removed, disrupt every path from
any source [founder] to any sink [proband]" ([Kirkpatrick et al., 2019](@ref)).
"""
function _bottlenecks(pedigree::Pedigree{T}) where T <: AbstractIndividual
    bottlenecks = Int64[]
    probandIDs = pro(pedigree)
    founderIDs = _lowest_founders(pedigree)
    isolated_pedigree = branching(pedigree, ancestors=founderIDs)
    candidateIDs = [ID for ID ∈ keys(isolated_pedigree)
                    if !isnothing(pedigree[ID].father)
                    && !isnothing(pedigree[ID].mother)]
    candidateIDs = [ID for ID ∈ candidateIDs
                     if length(pedigree[ID].father.children) == 1
                     && length(pedigree[ID].mother.children) == 1]
    for candidateID ∈ setdiff(candidateIDs, ∪(probandIDs, founderIDs))
        ancestors = ancestor(pedigree, candidateID)
        sourceIDs = ∩(ancestors, founderIDs)
        if all(_cut_vertex(pedigree[sourceID], candidateID) for sourceID ∈ sourceIDs)
            push!(bottlenecks, candidateID)
        end
    end
    sort(unique(bottlenecks))
end


"""
    phi(pedigree::Pedigree, Ψ::Matrix{Float64})

Return a square matrix of pairwise kinship coefficients
between all probands given the founders' kinships.

An implementation of the recursive-cut algorithm presented in [Kirkpatrick et al., 2019](@ref).
"""
function phi(pedigree::Pedigree, Ψ::Matrix{Float64})
    probandIDs = filter(x -> isempty(pedigree[x].children), collect(keys(pedigree)))
    probands = [pedigree[ID] for ID in probandIDs]
    founderIDs = filter(x -> isnothing(pedigree[x].father) && isnothing(pedigree[x].mother), collect(keys(pedigree)))
    founders = [pedigree[ID] for ID in founderIDs]
    Φ = zeros(length(pedigree), length(pedigree))
    for (index₁, founder₁) in enumerate(founders)
        for (index₂, founder₂) in enumerate(founders)
            if index₁ ≤ index₂
                coefficient = Ψ[index₁, index₂]
                Φ[founder₁.rank, founder₂.rank] = coefficient
                Φ[founder₂.rank, founder₁.rank] = coefficient
            end
        end
    end
    for individualᵢ in values(pedigree)
        i = individualᵢ.rank
        for individualⱼ in values(pedigree)
            j = individualⱼ.rank
            if Φ[i, j] > 0
                continue
            elseif i > j # i cannot be an ancestor of j
                father = individualᵢ.father
                mother = individualᵢ.mother
                if !isnothing(father)
                    coefficient = Φ[father.rank, individualⱼ.rank] / 2
                    Φ[i, j] += coefficient
                    Φ[j, i] += coefficient
                end
                if !isnothing(mother)
                    coefficient = Φ[mother.rank, individualⱼ.rank] / 2
                    Φ[i, j] += coefficient
                    Φ[j, i] += coefficient
                end
            elseif j > i # j cannot be an ancestor of i
                father = individualⱼ.father
                mother = individualⱼ.mother
                if !isnothing(father)
                    coefficient = Φ[father.rank, individualᵢ.rank] / 2
                    Φ[i, j] += coefficient
                    Φ[j, i] += coefficient
                end
                if !isnothing(mother)
                    coefficient = Φ[mother.rank, individualᵢ.rank] / 2
                    Φ[i, j] += coefficient
                    Φ[j, i] += coefficient
                end
            else # i == j, same individual
                father = individualᵢ.father
                mother = individualᵢ.mother
                if !isnothing(father) && !isnothing(mother)
                    Φ[i, i] += (1 + Φ[mother.rank, father.rank]) / 2
                else
                    Φ[i, i] += 0.5
                end
            end
        end
    end
    indices = [proband.rank for proband in probands]
    Φ[indices, indices]
end

"""
    phi(pedigree::Pedigree, probandIDs::Vector{Int64} = pro(pedigree); MT::Bool = false)

Return a square matrix of pairwise kinship coefficients between probands.

If no probands are given, return the square matrix for all probands in the pedigree.

An implementation of the recursive-cut algorithm presented in [Kirkpatrick et al., 2019](@ref).

# Example

```julia
import GenLib as gen
geneaJi = gen.geneaJi
ped = gen.genealogy(geneaJi)
gen.phi(ped)
```
"""
function phi(pedigree::Pedigree, probandIDs::Vector{Int64} = pro(pedigree); MT::Bool = false)
    isolated_pedigree = branching(pedigree, pro = probandIDs)
    founderIDs = _lowest_founders(isolated_pedigree)
    isolated_pedigree = branching(pedigree, ancestors = founderIDs)
    Ψ = zeros(length(founderIDs), length(founderIDs))
    for f in eachindex(founderIDs)
        Ψ[f, f] = 0.5
    end
    if MT
        index_pedigree = Pedigree{PossibleFounder}()
        index = 0
        for individual in collect(values(isolated_pedigree))
            father = individual.father
            mother = individual.mother
            if isnothing(father) && isnothing(mother)
                index += 1
                individual_index = copy(index)
            else
                individual_index = 0
            end
            index_pedigree[individual.ID] = PossibleFounder(
                individual.ID,
                isnothing(father) ? nothing : index_pedigree[father.ID],
                isnothing(mother) ? nothing : index_pedigree[mother.ID],
                PossibleFounder[],
                individual.sex,
                individual.rank,
                individual_index
            )
            if !isnothing(father)
                push!(index_pedigree[father.ID].children, index_pedigree[individual.ID])
            end
            if !isnothing(mother)
                push!(index_pedigree[mother.ID].children, index_pedigree[individual.ID])
            end
        end
        probands = filter(x -> isempty(x.children), collect(values(index_pedigree)))
        Φ = Matrix{Float64}(undef, length(probands), length(probands))
        Threads.@threads for i in eachindex(probands)
            Threads.@threads for j in eachindex(probands)
                if i ≤ j
                    Φ[i, j] = Φ[j, i] = phi(probands[i], probands[j], Ψ)
                end
            end
        end
    else
        Φ = phi(isolated_pedigree, Ψ)
    end
    probandIDs = filter(x -> x ∈ probandIDs, collect(keys(isolated_pedigree)))
    order = sortperm(probandIDs)
    Φ[order, order]
end

"""
    phi(pedigree::Pedigree, rowIDs::Vector{Int64}, columnIDs::Vector{Int64})

Return a rectangle matrix of kinship coefficients,
as defined by a list or row IDs and column IDs.
"""
function phi(pedigree::Pedigree, rowIDs::Vector{Int64}, columnIDs::Vector{Int64})
    Φ = Matrix{Float64}(undef, length(rowIDs), length(columnIDs))
    probandIDs = union(rowIDs, columnIDs)
    isolated_pedigree = branching(pedigree, pro = probandIDs)
    founderIDs = _lowest_founders(isolated_pedigree)
    isolated_pedigree = branching(pedigree, ancestors = founderIDs)
    Ψ = zeros(length(founderIDs), length(founderIDs))
    for f in eachindex(founderIDs)
        Ψ[f, f] = 0.5
    end
    phi_pedigree = Pedigree{PossibleFounder}()
    founder_index = 1
    for (rank, individual) in enumerate(collect(values(isolated_pedigree)))
        father = individual.father
        mother = individual.mother
        if isnothing(father) && isnothing(mother)
            index = copy(founder_index)
            founder_index += 1
        else
            index = 0
        end
        phi_pedigree[individual.ID] = PossibleFounder(
            individual.ID,
            isnothing(father) ? nothing : phi_pedigree[father.ID],
            isnothing(mother) ? nothing : phi_pedigree[mother.ID],
            PossibleFounder[],
            individual.sex,
            rank,
            index
        )
        if !isnothing(father)
            push!(phi_pedigree[father.ID].children, phi_pedigree[individual.ID])
        end
        if !isnothing(mother)
            push!(phi_pedigree[mother.ID].children, phi_pedigree[individual.ID])
        end
    end
    Threads.@threads for i in eachindex(rowIDs)
        individualᵢ = phi_pedigree[rowIDs[i]]
        Threads.@threads for j in eachindex(columnIDs)
            individualⱼ = phi_pedigree[columnIDs[j]]
            Φ[i, j] = phi(individualᵢ, individualⱼ, Ψ)
        end
    end
    Φ
end

"""
    _cleanup(kinship_matrix::Matrix{Float64, Float64}, threshold::Float64 = 0.0625)

Return a vector of booleans of probands whose kinships
between them never exceed a minimum `threshold`.
    
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