"""
    mutable struct IndexedIndividual <: AbstractIndividual
        ID::Int64
        father::Union{Nothing, IndexedIndividual}
        mother::Union{Nothing, IndexedIndividual}
        sex::Int64
        children::Vector{IndexedIndividual}
        max_height::Int64
        rank::Int64
        founder_index::Int64
    end

An individual with an index to access the founder's kinships in the Ψ matrix.

`max_height` is also used to do a topological sort of the pedigree.
"""
mutable struct IndexedIndividual <: AbstractIndividual
    ID::Int64
    father::Union{Nothing, IndexedIndividual}
    mother::Union{Nothing, IndexedIndividual}
    sex::Int64
    children::Vector{IndexedIndividual}
    max_height::Int64
    rank::Int64
    founder_index::Int64
end

"""
    phi(individualᵢ::Individual, individualⱼ::Individual)

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
function phi(individualᵢ::Individual, individualⱼ::Individual)
    value = 0.
    if individualᵢ.rank > individualⱼ.rank
        # From the genealogical order, i cannot be an ancestor of j
        # Φᵢⱼ = (Φₚⱼ + Φₘⱼ) / 2, if i is not an ancestor of j (Karigl, 1981)
        if !isnothing(individualᵢ.father)
            value += phi(individualᵢ.father, individualⱼ) / 2
        end
        if !isnothing(individualᵢ.mother)
            value += phi(individualᵢ.mother, individualⱼ) / 2
        end
    elseif individualⱼ.rank > individualᵢ.rank
        # Reverse the order since i can be an ancestor of j
        # Φⱼᵢ = (Φₚᵢ + Φₘᵢ) / 2, if j is not an ancestor of i (Karigl, 1981)
        if !isnothing(individualⱼ.father)
            value += phi(individualⱼ.father, individualᵢ) / 2
        end
        if !isnothing(individualⱼ.mother)
            value += phi(individualⱼ.mother, individualᵢ) / 2
        end
    elseif individualᵢ.rank == individualⱼ.rank
        # Same individual
        # Φₐₐ = (1 + Φₚₘ) / 2 (Karigl, 1981)
        value += 0.5
        if !isnothing(individualᵢ.father) & !isnothing(individualᵢ.mother)
            value += phi(individualᵢ.father, individualᵢ.mother) / 2
        end
    end
    return value
end

"""
    phi(pedigree::Pedigree, rowIDs::Vector{Int64}, columnIDs::Vector{Int64})

Return a rectangle matrix of kinship coefficients,
as defined by a list or row IDs and column IDs.
"""
function phi(pedigree::Pedigree, rowIDs::Vector{Int64}, columnIDs::Vector{Int64})
    ϕ = Matrix{Float64}(undef, length(rowIDs), length(columnIDs))
    if rowIDs == columnIDs
        Threads.@threads for i in eachindex(rowIDs)
            Threads.@threads for j in eachindex(columnIDs)
                if i ≤ j
                    individualᵢ = pedigree[rowIDs[i]]
                    individualⱼ = pedigree[columnIDs[j]]
                    coefficient = phi(individualᵢ, individualⱼ)
                    ϕ[i, j] = coefficient
                    ϕ[j, i] = coefficient
                end
            end
        end
    else
        Threads.@threads for i in eachindex(rowIDs)
            individualᵢ = pedigree[rowIDs[i]]
            Threads.@threads for j in eachindex(columnIDs)
                individualⱼ = pedigree[columnIDs[j]]
                ϕ[i, j] = phi(individualᵢ, individualⱼ)
            end
        end
    end
    ϕ
end

"""
    phi(individualᵢ::IndexedIndividual, individualⱼ::IndexedIndividual, Ψ::Matrix{Float64})

Return the kinship coefficient between two individuals given a matrix of the founders' kinships.

Adapted from [Karigl, 1981](@ref), and [Kirkpatrick et al., 2019](@ref).
"""
function phi(individualᵢ::IndexedIndividual, individualⱼ::IndexedIndividual, Ψ::Matrix{Float64})
    value = 0.
    if individualᵢ.founder_index != 0 && individualⱼ.founder_index != 0
        # Both individuals are founders, so we already know their kinship coefficient
        value += Ψ[individualᵢ.founder_index, individualⱼ.founder_index]
    elseif individualᵢ.founder_index != 0
        # Individual i is a founder, so we climb the pedigree on individual j's side
        if !isnothing(individualⱼ.father)
            value += phi(individualᵢ, individualⱼ.father, Ψ) / 2
        end
        if !isnothing(individualⱼ.mother)
            value += phi(individualᵢ, individualⱼ.mother, Ψ) / 2
        end
    elseif individualⱼ.founder_index != 0
        # Individual j is a founder, so we climb the pedigree on individual i's side
        if !isnothing(individualᵢ.father)
            value += phi(individualⱼ, individualᵢ.father, Ψ) / 2
        end
        if !isnothing(individualᵢ.mother)
            value += phi(individualⱼ, individualᵢ.mother, Ψ) / 2
        end
    else
        # None of the individuals are founders, so we climb on the side
        # of the individual that appears lowest in the pedigree
        if individualᵢ.rank > individualⱼ.rank
            # From the genealogical order, i cannot be an ancestor of j
            # Φᵢⱼ = (Φₚⱼ + Φₘⱼ) / 2, if i is not an ancestor of j (Karigl, 1981)
            if !isnothing(individualᵢ.father)
                value += phi(individualᵢ.father, individualⱼ, Ψ) / 2
            end
            if !isnothing(individualᵢ.mother)
                value += phi(individualᵢ.mother, individualⱼ, Ψ) / 2
            end
        elseif individualⱼ.rank > individualᵢ.rank
            # Reverse the order since i can be an ancestor of j
            # Φⱼᵢ = (Φₚᵢ + Φₘᵢ) / 2, if j is not an ancestor of i (Karigl, 1981)
            if !isnothing(individualⱼ.father)
                value += phi(individualⱼ.father, individualᵢ, Ψ) / 2
            end
            if !isnothing(individualⱼ.mother)
                value += phi(individualⱼ.mother, individualᵢ, Ψ) / 2
            end
        elseif individualᵢ.rank == individualⱼ.rank
            # Same individual
            # Φₐₐ = (1 + Φₚₘ) / 2 (Karigl, 1981)
            value += 0.5
            if !isnothing(individualᵢ.father) && !isnothing(individualᵢ.mother)
                value += phi(individualᵢ.father, individualᵢ.mother, Ψ) / 2
            end
        end
    end
    value
end

"""
    function _lowest_founders(pedigree::Pedigree)

Return the lowest founders of a given pedigree.

The lowest founder is defined as an only child who is either a
founder or as the only child of a lineage of only children.
"""
function _lowest_founders(pedigree::Pedigree)
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
    _max_height!(individual::IndexedIndividual)

Set and return the maximum height of an individual in the pedigree.
"""
function _max_height!(individual::IndexedIndividual)
    if individual.max_height == -1
        if isempty(individual.children)
            individual.max_height = 0
        else
            individual.max_height = maximum([_max_height!(child) for child ∈ individual.children]) + 1
        end
    end
    individual.max_height
end

"""
    _topological_sort(pedigree::Pedigree)

Return a reordered pedigree where the individuals are in chronological order,
i.e. any individual's parents appear before them. In this topological sort,
the parents of the probands appear as far in the pedigree as possible, etc.
"""
function _topological_sort(pedigree::Pedigree)
    indexed_pedigree = Pedigree{IndexedIndividual}()
    for individual ∈ values(pedigree)
        indexed_pedigree[individual.ID] = IndexedIndividual(
            individual.ID,
            !isnothing(individual.father) ? indexed_pedigree[individual.father.ID] : nothing,
            !isnothing(individual.mother) ? indexed_pedigree[individual.mother.ID] : nothing,
            individual.sex, [], -1, 0, 0
        )
    end
    for individual ∈ values(indexed_pedigree)
        if !isnothing(individual.father)
            push!(individual.father.children, individual)
        end
        if !isnothing(individual.mother)
            push!(individual.mother.children, individual)
        end
    end
    IDs = collect(keys(pedigree))
    heights = [_max_height!(individual) for individual ∈ values(indexed_pedigree)]
    order = sortperm(heights, rev=true)
    sortedIDs = IDs[order]
    for (rank, ID) ∈ enumerate(sortedIDs)
        indexed_pedigree[ID].rank = rank
    end
    indexed_pedigree
end

"""
    _previous_generation(pedigree::Pedigree, next_generation::Vector{Int64})

Return the previous generation of a given set of individuals.
"""
function _previous_generation(pedigree::Pedigree, next_generationIDs::Vector{Int64})
    previous_generationIDs = Int64[]
    # We extract all the parents of the next generation,
    # and keep the founders of the next generation
    for ID ∈ next_generationIDs
        individual = pedigree[ID]
        father = individual.father
        mother = individual.mother
        if isnothing(father) && isnothing(mother)
            # The individual is a founder, so we keep them
            push!(previous_generationIDs, ID)
        else
            if !isnothing(father)
                push!(previous_generationIDs, father.ID)
            end
            if !isnothing(mother)
                push!(previous_generationIDs, mother.ID)
            end
        end
    end
    minimum_rank = minimum([pedigree[ID].rank for ID ∈ previous_generationIDs])
    candidateIDs = [individual.ID for individual ∈ values(pedigree) if individual.rank ≥ minimum_rank]
    parentIDs = Set{Int64}()
    for ID ∈ candidateIDs
        individual = pedigree[ID]
        if length(individual.children) == 0
            push!(parentIDs, ID)
        else
            for child ∈ individual.children
                if !isnothing(child.father)
                    push!(parentIDs, child.father.ID)
                end
                if !isnothing(child.mother)
                    push!(parentIDs, child.mother.ID)
                end
            end
        end
    end
    descendantIDs = Set{Int64}()
    for ID ∈ candidateIDs
        descendants = descendant(pedigree, ID)
        for descendant ∈ descendants
            push!(descendantIDs, descendant)
        end
    end
    sort(collect(setdiff(parentIDs, descendantIDs)))
end

"""
    phi(pedigree::Pedigree, Ψ::Matrix{Float64}, topIDs::Vector{Int64}, bottomIDs::Vector{Int64})

Return a square matrix of pairwise kinship coefficients
between all probands given the founders' kinships.

An implementation of the recursive-cut algorithm presented in [Kirkpatrick et al., 2019](@ref).
"""
function phi(pedigree::Pedigree, Ψ::Matrix{Float64}, topIDs::Vector{Int64}, bottomIDs::Vector{Int64})
    A = Dict{Int64, Set{Int64}}()
    for individual ∈ values(pedigree)
        A[individual.ID] = Set(individual.ID)
        if !isnothing(individual.father)
            union!(A[individual.ID], A[individual.father.ID])
        end
        if !isnothing(individual.mother)
            union!(A[individual.ID], A[individual.mother.ID])
        end
    end
    ϕ = ones(length(pedigree), length(pedigree)) .* -1
    indices = [pedigree[ID].rank for ID ∈ topIDs]
    ϕ[indices, indices] = copy(Ψ)
    for individualᵢ ∈ values(pedigree)
        i = individualᵢ.rank
        for individualⱼ ∈ values(pedigree)
            j = individualⱼ.rank
            if ϕ[i, j] > -1
                continue
            elseif i == j
                coefficient = 0.5
                father = individualᵢ.father
                mother = individualᵢ.mother
                if !isnothing(father) && !isnothing(mother)
                    coefficient += ϕ[mother.rank, father.rank] / 2
                end
                ϕ[i, i] = coefficient
            else
                if ϕ[individualⱼ.rank, individualⱼ.rank] > -1 && individualᵢ.ID ∉ A[individualⱼ.ID]
                    father = individualᵢ.father
                    mother = individualᵢ.mother
                    individual = individualⱼ
                else
                    father = individualⱼ.father
                    mother = individualⱼ.mother
                    individual = individualᵢ
                end
                coefficient = 0.
                if !isnothing(father)
                    coefficient += ϕ[father.rank, individual.rank] / 2
                end
                if !isnothing(mother)
                    coefficient += ϕ[mother.rank, individual.rank] / 2
                end
                ϕ[i, j] = coefficient
                ϕ[j, i] = coefficient
            end
        end
    end
    indices = [pedigree[ID].rank for ID ∈ bottomIDs]
    ϕ[indices, indices]
end

"""
    phi(pedigree::Pedigree, probandIDs::Vector{Int64} = pro(pedigree); MT::Bool = false, verbose::Bool = false)

Return a square matrix of pairwise kinship coefficients between probands.

If no probands are given, return the square matrix
for all probands in the pedigree.

If MT = false: an implementation of the recursive-cut algorithm
presented in [Kirkpatrick et al., 2019](@ref).

If MT = true: pairwise kinships in parallel, a hybrid between
the algorithms of [Karigl, 1981](@ref), and [Kirkpatrick et al., 2019](@ref).

# Example

```julia
import GenLib as gen
geneaJi = gen.geneaJi
ped = gen.genealogy(geneaJi)
gen.phi(ped)
```
"""
function phi(pedigree::Pedigree, probandIDs::Vector{Int64} = pro(pedigree); MT::Bool = false, verbose::Bool = false)
    global Ψ
    isolated_pedigree = branching(pedigree, pro = probandIDs)
    indexed_pedigree = _topological_sort(isolated_pedigree)
    cut_vertices = [probandIDs]
    founderIDs = founder(indexed_pedigree)
    previous_generation = _previous_generation(indexed_pedigree, probandIDs)
    while previous_generation != founderIDs
        pushfirst!(cut_vertices, previous_generation)
        previous_generation = _previous_generation(indexed_pedigree, previous_generation)
    end
    pushfirst!(cut_vertices, founderIDs)
    if verbose
        for i ∈ 1:length(cut_vertices)-1
            previous_generation = cut_vertices[i]
            next_generation = cut_vertices[i+1]
            verbose_pedigree = branching(isolated_pedigree, pro = next_generation, ancestors = previous_generation)
            println("Step $i / $(length(cut_vertices)-1): $(length(previous_generation)) founders, $(length(next_generation)) probands (n = $(length(verbose_pedigree))).")
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
            for (index, ID) ∈ enumerate(previous_generation)
                indexed_pedigree[ID].founder_index = index
            end
            ϕ = Matrix{Float64}(undef, length(next_generation), length(next_generation))
            probands = [indexed_pedigree[ID] for ID ∈ next_generation]
            Threads.@threads for i ∈ eachindex(probands)
                Threads.@threads for j ∈ eachindex(probands)
                    if i ≤ j
                        ϕ[i, j] = ϕ[j, i] = phi(probands[i], probands[j], Ψ)
                    end
                end
            end
            Ψ = ϕ
        else
            isolated_pedigree = branching(pedigree, pro = next_generation, ancestors = previous_generation)
            Ψ = phi(isolated_pedigree, Ψ, previous_generation, next_generation)
        end
    end
    Ψ
end

"""
    _trim_kinships(kinship_matrix::Matrix{Float64, Float64}, threshold::Float64 = 0.0625)

Return a vector of booleans of probands whose kinships
between them never exceed a given `threshold`.
    
For instance, a threshold of 0.0625 (the default) removes individuals
who are first-degree cousins or closer.
"""
function _trim_kinships(kinship_matrix::Matrix{Float64}, threshold::Float64 = 0.0625)
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

"""
    mutable struct PossibleDescendant <: AbstractIndividual

An individual with its ancestor's contribution.
"""
mutable struct PossibleDescendant <: AbstractIndividual
    ID::Int64
    father::Union{Nothing, PossibleDescendant}
    mother::Union{Nothing, PossibleDescendant}
    children::Vector{PossibleDescendant}
    contribution::Float64
end

"""
    _contribute!(individual::Individual, depth::Int64 = 0)

Recursively compute the genetic contributions of an individual.
"""
function _contribute!(individual::PossibleDescendant, depth::Int64 = 0)
    # Ported from GENLIB's ExploreConGenProposant
    if isempty(individual.children)
        individual.contribution += 0.5 ^ depth
    else
        for child ∈ individual.children
            _contribute!(child, depth + 1)
        end
    end
end

"""
    gc(pedigree::Pedigree; pro::Vector{Int64} = pro(pedigree), ancestors::Vector{Int64} = founder(pedigree))

Return a matrix of the genetic contribution
of each ancestor (columns) to each proband (rows).

# Example

```julia
import GenLib as gen
genea140 = gen.genea140
ped = gen.genealogy(genea140)
contributions = gen.gc(ped)
```
"""
function gc(
    pedigree::Pedigree;
    pro::Vector{Int64} = pro(pedigree),
    ancestors::Vector{Int64} = founder(pedigree))
    
    # Ported from GENLIB's Congen
    matrix = zeros(length(pro), length(ancestors))
    contribution_pedigree = Pedigree{PossibleDescendant}()
    for individual ∈ collect(values(pedigree))
        father = individual.father
        mother = individual.mother
        contribution_pedigree[individual.ID] = PossibleDescendant(
            individual.ID,
            isnothing(father) ? nothing : contribution_pedigree[father.ID],
            isnothing(mother) ? nothing : contribution_pedigree[mother.ID],
            PossibleDescendant[],
            0.
        )
        if !isnothing(father)
            push!(contribution_pedigree[father.ID].children, contribution_pedigree[individual.ID])
        end
        if !isnothing(mother)
            push!(contribution_pedigree[mother.ID].children, contribution_pedigree[individual.ID])
        end
    end
    for (index₁, ancestorID) ∈ enumerate(ancestors)
        ancestor = contribution_pedigree[ancestorID]
        _contribute!(ancestor)
        for (index₂, probandID) ∈ enumerate(pro)
            proband = contribution_pedigree[probandID]
            matrix[index₂, index₁] = proband.contribution
            proband.contribution = 0.
        end
    end
    matrix
end