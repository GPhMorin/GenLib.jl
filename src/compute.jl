"""
    mutable struct IndexedIndividual <: AbstractIndividual
        ID::Int64
        father::Union{Nothing, IndexedIndividual}
        mother::Union{Nothing, IndexedIndividual}
        rank::Int64
        founder_index::Int64
    end

An individual with an index to access the founder's kinships in the Ψ matrix.
"""
mutable struct IndexedIndividual <: AbstractIndividual
    ID::Int64
    father::Union{Nothing, IndexedIndividual}
    mother::Union{Nothing, IndexedIndividual}
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
    _phi(individualᵢ::Individual, individualⱼ::Individual)

Return the kinship coefficient between two individuals.

Adapted from the recursive function of [Karigl, 1981](@ref).
Actually slower than the recursive version, so not in use.

function _phi(individualᵢ::Individual, individualⱼ::Individual)
    coefficient = 0.
    if individualᵢ.rank < individualⱼ.rank
        # First one in the tuple is always ranked lower.
        # This is done so that once the second individual has no parents,
        # we know both of them are founders.
        stack = [(individualᵢ, individualⱼ, 1)]
    else
        stack = [(individualⱼ, individualᵢ, 1)]
    end
    while !isempty(stack)
        individualᵢ, individualⱼ, depth = pop!(stack)
        if individualᵢ.rank == individualⱼ.rank # Same individual
            coefficient += 0.5 ^ depth
            if !isnothing(individualᵢ.father) && !isnothing(individualⱼ.mother)
                if individualᵢ.father.rank < individualᵢ.mother.rank
                    pushfirst!(stack, (individualᵢ.father, individualⱼ.mother, depth + 1))
                else
                    pushfirst!(stack, (individualᵢ.mother, individualᵢ.father, depth + 1))
                end
            end
        else
            if !isnothing(individualⱼ.father)
                if individualⱼ.father.rank < individualᵢ.rank
                    pushfirst!(stack, (individualⱼ.father, individualᵢ, depth + 1))
                else
                    pushfirst!(stack, (individualᵢ, individualⱼ.father, depth + 1))
                end
            end
            if !isnothing(individualⱼ.mother)
                if individualⱼ.mother.rank < individualᵢ.rank
                    pushfirst!(stack, (individualⱼ.mother, individualᵢ, depth + 1))
                else
                    pushfirst!(stack, (individualᵢ, individualⱼ.mother, depth + 1))
                end
            end
        end
    end
    coefficient
end
"""

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
    _max_height(individual::Individual)

Return the maximum height of an individual's pedigree.
"""
function _max_height(individual::Individual)
    if isempty(individual.children)
        max_children_height = 0
    else
        max_children_height = maximum([_max_height(child) for child ∈ individual.children])
    end
    max_children_height + 1
end

"""
    _topological_sort(pedigree::Pedigree)

Return a reordered pedigree where the individuals are in chronological order,
i.e. any individual's parents appear before them. In this topological sort,
the parents of the probands appear as far in the pedigree as possible, etc.
"""
function _topological_sort(pedigree::Pedigree)
    IDs = collect(keys(pedigree))
    heights = [_max_height(individual) for individual ∈ values(pedigree)]
    order = sortperm(heights, rev=true)
    sortedIDs = IDs[order]
    ordered_pedigree = Pedigree{Individual}()
    rank = 1
    for ID ∈ sortedIDs
        individual = pedigree[ID]
        ordered_pedigree[ID] = Individual(
            ID,
            !isnothing(individual.father) ? ordered_pedigree[individual.father.ID] : nothing,
            !isnothing(individual.mother) ? ordered_pedigree[individual.mother.ID] : nothing,
            [],
            individual.sex,
            rank
        )
        rank += 1
    end
    for individual ∈ values(ordered_pedigree)
        if !isnothing(individual.father)
            push!(individual.father.children, individual)
        end
        if !isnothing(individual.mother)
            push!(individual.mother.children, individual)
        end
    end
    ordered_pedigree
end

"""
    _previous_generation(pedigree::Pedigree, next_generation::Vector{Int64})

Return the previous generation of a given set of individuals.
"""
function _previous_generation(pedigree::Pedigree, next_generationIDs::Vector{Int64})
    intermediateIDs = Int64[]
    # We extract all the parents of the next generation,
    # and keep the founders of the next generation
    for ID ∈ next_generationIDs
        individual = pedigree[ID]
        father = individual.father
        mother = individual.mother
        if isnothing(father) && isnothing(mother)
            # The individual is a founder, so we keep them
            push!(intermediateIDs, ID)
        else
            if !isnothing(father)
                push!(intermediateIDs, father.ID)
            end
            if !isnothing(mother)
                push!(intermediateIDs, mother.ID)
            end
        end
    end
    # The minimum rank indicates where we find the highest parent (reference
    # individual) so that all the other individuals will be below them
    minimum_rank = minimum([pedigree[ID].rank for ID ∈ intermediateIDs])
    previous_generationIDs = Int64[]
    for ID ∈ intermediateIDs
        # We replace the intermediate candidates with the highest ancestors
        # that are ranked equal or below the reference individual
        individual = pedigree[ID]
        stack = [individual]
        while !isempty(stack)
            candidate = pop!(stack)
            if candidate.rank == minimum_rank
                # This is the reference individual with the highest position
                # among the parents, so we keep that individual
                push!(previous_generationIDs, candidate.ID)
            elseif isnothing(candidate.father) && isnothing(candidate.mother)
                # The individual is a founder, so we keep them
                push!(previous_generationIDs, candidate.ID)
            elseif !isnothing(candidate.father) && !isnothing(candidate.mother)
                if candidate.father.rank ≥ minimum_rank && candidate.mother.rank ≥ minimum_rank
                    # We need to check both parents' ranks because we don't want to keep
                    # both a parent and their child, as that would cause a conflict
                    push!(stack, candidate.father)
                    push!(stack, candidate.mother)
                else
                    # At least one parent is ranked higher
                    # than the reference individual,
                    # so we keep the child
                    push!(previous_generationIDs, candidate.ID)
                end
            elseif !isnothing(candidate.father)
                # The individual is a half founder,
                # so we just check the father
                if candidate.father.rank ≥ minimum_rank
                    push!(stack, candidate.father)
                else
                    push!(previous_generationIDs, candidate.ID)
                end
            elseif !isnothing(candidate.mother)
                # The individual is a half founder,
                # so we just check the mother
                if candidate.mother.rank ≥ minimum_rank
                    push!(stack, candidate.mother)
                else
                    push!(previous_generationIDs, candidate.ID)
                end
            end
        end
    end
    sort(unique(previous_generationIDs))
end

"""
    phi(pedigree::Pedigree, Ψ::Matrix{Float64}, topIDs::Vector{Int64}, bottomIDs::Vector{Int64})

Return a square matrix of pairwise kinship coefficients
between all probands given the founders' kinships.

An implementation of the recursive-cut algorithm presented in [Kirkpatrick et al., 2019](@ref).
"""
function phi(pedigree::Pedigree, Ψ::Matrix{Float64}, topIDs::Vector{Int64}, bottomIDs::Vector{Int64})
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
                father = individualᵢ.father
                mother = individualᵢ.mother
                if !isnothing(father) && !isnothing(mother)
                    coefficient = (1 + ϕ[mother.rank, father.rank]) / 2
                    ϕ[i, i] = coefficient
                else
                    coefficient = 0.5
                    ϕ[i, i] = coefficient
                end
            elseif i < j
                fatherᵢ = individualᵢ.father
                motherᵢ = individualᵢ.mother
                coefficient₁ = 0.
                if !isnothing(fatherᵢ)
                    coefficient₁ += ϕ[fatherᵢ.rank, individualⱼ.rank] / 2
                end
                if !isnothing(motherᵢ)
                    coefficient₁ += ϕ[motherᵢ.rank, individualⱼ.rank] / 2
                end
                fatherⱼ = individualⱼ.father
                motherⱼ = individualⱼ.mother
                coefficient₂ = 0.
                if !isnothing(fatherⱼ)
                    coefficient₂ += ϕ[fatherⱼ.rank, individualᵢ.rank] / 2
                end
                if !isnothing(motherⱼ)
                    coefficient₂ += ϕ[motherⱼ.rank, individualᵢ.rank] / 2
                end
                coefficient = max(coefficient₁, coefficient₂)
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
    isolated_pedigree = _topological_sort(isolated_pedigree)
    if MT
        indexed_pedigree = Pedigree{IndexedIndividual}()
        for individual ∈ values(isolated_pedigree)
            indexed_pedigree[individual.ID] = IndexedIndividual(
                individual.ID,
                !isnothing(individual.father) ? indexed_pedigree[individual.father.ID] : nothing,
                !isnothing(individual.mother) ? indexed_pedigree[individual.mother.ID] : nothing,
                individual.rank,
                0
            )
        end
    end
    cut_vertices = [probandIDs]
    founderIDs = founder(isolated_pedigree)
    previous_generation = _previous_generation(isolated_pedigree, probandIDs)
    while previous_generation != founderIDs
        pushfirst!(cut_vertices, previous_generation)
        previous_generation = _previous_generation(isolated_pedigree, previous_generation)
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