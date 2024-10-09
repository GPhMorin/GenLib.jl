"""
    mutable struct IndexedIndividual <: AbstractIndividual
        ID::Int
        father::Union{Nothing, IndexedIndividual}
        mother::Union{Nothing, IndexedIndividual}
        sex::Int
        children::Vector{IndexedIndividual}
        rank::Int
        founder_index::Int
    end

An individual with an index to access the founder's kinships in the Ψ matrix.
"""
mutable struct IndexedIndividual <: AbstractIndividual
    ID::Int
    father::Union{Nothing, IndexedIndividual}
    mother::Union{Nothing, IndexedIndividual}
    sex::Int
    children::Vector{IndexedIndividual}
    rank::Int
    founder_index::Int
    is_proband::Bool
    children_to_process::Int
end

"""
    struct KinshipMatrix

A minimal structure wrapping an `Dict` with kinships of individuals accessed by IDs.
"""
struct KinshipMatrix
    dict::Dict{Int32, Dict{Int32, Float32}}
    ID_to_rank::Dict{Int32, Int}
end

function Base.getindex(ϕ::KinshipMatrix, ID₁::Int, ID₂::Int)
    (rank₁, rank₂) = (ϕ.ID_to_rank[ID₁], ϕ.ID_to_rank[ID₂])
    (rank₁, rank₂) = rank₁ < rank₂ ? (rank₁, rank₂) : (rank₂, rank₁)
    haskey(ϕ.dict[rank₁], rank₂) ? ϕ.dict[rank₁][rank₂] : 0.
end

function Base.show(io::IO, ::MIME"text/plain", ϕ::KinshipMatrix)
    nz = sum(length(kinships) for kinships ∈ values(ϕ.dict))
    print(io, "$(length(ϕ.dict))×$(length(ϕ.dict)) KinshipMatrix " *
        "with $nz stored entries.")
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
    phi(individualᵢ::IndexedIndividual, individualⱼ::IndexedIndividual, Ψ::Matrix{Float32})

Return the kinship coefficient between two individuals given a matrix of the founders'
kinships.

Adapted from [Karigl, 1981](@ref), and [Kirkpatrick et al., 2019](@ref).
"""
function phi(individualᵢ::IndexedIndividual, individualⱼ::IndexedIndividual,
    Ψ::Matrix{Float32})
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
        # None of the individuals are founders, so we climb on the side of the individual
        # that appears lowest in the pedigree
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
    _index_pedigree(pedigree::Pedigree)

Return a pedigree with an index that indicates the individual's position among the founders.
"""
function _index_pedigree(pedigree::Pedigree)
    indexed_pedigree = Pedigree{IndexedIndividual}()
    for individual ∈ values(pedigree)
        indexed_pedigree[individual.ID] = IndexedIndividual(
            individual.ID,
            !isnothing(individual.father) ?
                indexed_pedigree[individual.father.ID] : nothing,
            !isnothing(individual.mother) ?
                indexed_pedigree[individual.mother.ID] : nothing,
            individual.sex, [], individual.rank, 0, false, 0
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
    indexed_pedigree
end

"""
    _next_generation(pedigree::Pedigree, previous_generationIDs::Vector{Int})

Return the next generation of a given set of individuals.
"""
function _previous_generation(pedigree::Pedigree, next_generationIDs::Vector{Int})
    previous_generationIDs = Int[]
    for ID ∈ next_generationIDs
        individual = pedigree[ID]
        father = individual.father
        if !isnothing(father)
            push!(previous_generationIDs, father.ID)
        end
        mother = individual.mother
        if !isnothing(mother)
            push!(previous_generationIDs, mother.ID)
        end
    end
    unique!(previous_generationIDs)
end

"""
    phi(pedigree::Pedigree, probandIDs::Vector{Int} = pro(pedigree);
        verbose::Bool = false)

Return a square matrix of pairwise kinship coefficients between probands.

If no probands are given, return the square matrix for all probands in the pedigree.

The algorithm computes pairwise kinships in parallel, a hybrid between the algorithms of
[Karigl, 1981](@ref), and [Kirkpatrick et al., 2019](@ref).

If `verbose` is `true`, print the information about the cut vertices. If `compute` is
`true` (the default), compute the kinship matrix. If it is `false`, only print the
information about the cut vertices.

# Example

```julia
import GenLib as gen
geneaJi = gen.geneaJi
ped = gen.genealogy(geneaJi)
gen.phi(ped)
```
"""
function phi(pedigree::Pedigree, probandIDs::Vector{Int} = pro(pedigree);
    verbose::Bool = false, compute::Bool = true)
    # Start from the probands and go up until the highest founder(s)
    cut_vertices = [probandIDs]
    previous_generationIDs = _previous_generation(pedigree, probandIDs)
    while !isempty(previous_generationIDs)
        pushfirst!(cut_vertices, previous_generationIDs)
        previous_generationIDs = _previous_generation(pedigree, previous_generationIDs)
    end
    # Drag the individuals down the generations for as long as they are required
    top_down = copy(cut_vertices)
    bottom_up = copy(cut_vertices)
    reverse!(bottom_up)
    for i ∈ 1:length(cut_vertices)-1
        top_down[i + 1] = union(top_down[i + 1], top_down[i])
        bottom_up[i + 1] = union(bottom_up[i + 1], bottom_up[i])
    end
    reverse!(bottom_up)
    cut_vertices = [∩(i, j) for (i, j) ∈ zip(top_down, bottom_up)]
    # Describe each pair of generations, if desired
    if verbose || !compute
        for i ∈ 1:length(cut_vertices)-1
            previous_generationIDs = cut_vertices[i]
            next_generationIDs = cut_vertices[i+1]
            println("Step $i of $(length(cut_vertices)-1): " *
                "$(length(previous_generationIDs)) founders, " *
                "$(length(next_generationIDs)) probands, " *
                "$(length(∩(previous_generationIDs, next_generationIDs))) both.")
        end
    end
    # Stop here if the user only wants to print information about each pair of generations
    if !compute
        return
    end
    # Add the `founder_index` attribute to the individuals and make them mutable so we can
    # quickly track the location of the founders in their kinship matrix.
    indexed_pedigree = _index_pedigree(pedigree)
    # Initialize the kinship matrix of the top founders
    Ψ = zeros(Float32, length(first(cut_vertices)), length(first(cut_vertices)))
    for i ∈ axes(Ψ, 1)
        Ψ[i, i] = 0.5
    end
    # For each pair of generations…
    for k ∈ 1:length(cut_vertices)-1
        previous_generationIDs = cut_vertices[k]
        next_generationIDs = cut_vertices[k+1]
        # Describe each pair of generations, if desired
        if verbose
            println("Running step $k of $(length(cut_vertices)-1) " *
                "($(length(previous_generationIDs)) founders, " *
                "$(length(next_generationIDs)) probands, " *
                "$(length(∩(previous_generationIDs, next_generationIDs))) both).")
        end
        # Assign the index to each individual from the previous generation
        for (index, ID) ∈ enumerate(previous_generationIDs)
            indexed_pedigree[ID].founder_index = index
        end
        # Fill the matrix in parallel, using the adapted algorithm from Karigl, 1981
        ϕ = Matrix{Float32}(undef, length(next_generationIDs), length(next_generationIDs))
        probands = [indexed_pedigree[ID] for ID ∈ next_generationIDs]
        Threads.@threads for i ∈ eachindex(probands)
            Threads.@threads for j ∈ eachindex(probands)
                if i ≤ j
                    ϕ[i, j] = ϕ[j, i] = phi(probands[i], probands[j], Ψ)
                end
            end
        end
        # At each iteration, the next generation becomes the previous generation
        Ψ = ϕ
    end
    Ψ
end

"""
    function sparse_phi(pedigree::Pedigree, probandIDs::Vector{Int} = pro(pedigree))

An implementation of Kirkpatrick et al.'s algorithm to compute the kinship matrix.

Based on this interpretation: https://lineagekit.github.io/lineagekit/use_cases/kinship.html

# Example

```julia
import GenLib as gen
geneaJi = gen.geneaJi
ped = gen.genealogy(geneaJi)
gen.sparse_phi(ped)
"""
function sparse_phi(pedigree::Pedigree, probandIDs::Vector{Int} = pro(pedigree))
    # Remove unrelated individuals
    isolated_pedigree = branching(pedigree, pro = probandIDs)
    # Index the pedigree for faster access
    indexed_pedigree = _index_pedigree(isolated_pedigree)
    # Mark the probands for faster access
    for ID ∈ probandIDs
        indexed_pedigree[ID].is_proband = true
    end
    # Initialize the kinship matrix
    ϕ = Dict{Int32, Dict{Int32, Float32}}()
    # Initialize the queue with the probands
    ranks_to_visit = Set{Int}()
    sizehint!(ranks_to_visit, length(isolated_pedigree))
    queue = Deque{IndexedIndividual}()
    for ID ∈ founder(isolated_pedigree)
        individual = indexed_pedigree[ID]
        push!(queue, individual)
    end
    rank = 1
    while !isempty(queue)
        individualᵢ = popfirst!(queue)
        rankᵢ = individualᵢ.rank
        father_rank = isnothing(individualᵢ.father) ? 0 : individualᵢ.father.rank
        mother_rank = isnothing(individualᵢ.mother) ? 0 : individualᵢ.mother.rank
        # Initialize the kinship dictionary for the individual
        ϕ[rankᵢ] = Dict{Int32, Float32}()
        # Kinship with self
        coefficient = 0.5
        if father_rank != 0 && mother_rank != 0
            if father_rank < mother_rank
                if haskey(ϕ[father_rank], mother_rank)
                    coefficient += ϕ[father_rank][mother_rank] / 2
                end
            else
                if haskey(ϕ[mother_rank], father_rank)
                    coefficient += ϕ[mother_rank][father_rank] / 2
                end
            end
        end
        ϕ[rankᵢ][rankᵢ] = coefficient
        # Kinship with previous individuals
        for rankⱼ ∈ ranks_to_visit
            coefficient = 0.
            if father_rank != 0
                # In order to make the kinships as sparse as possible,
                # we only store the kinship with the lowest ranked founder
                # and the kinship is inserted in rank order for faster lookup
                if rankⱼ < father_rank
                    if haskey(ϕ[rankⱼ], father_rank)
                        coefficient += ϕ[rankⱼ][father_rank] / 2
                    end
                else
                    if haskey(ϕ[father_rank], rankⱼ)
                        coefficient += ϕ[father_rank][rankⱼ] / 2
                    end
                end
            end
            if mother_rank != 0
                # Same thing but on the mother's side
                if rankⱼ < mother_rank
                    if haskey(ϕ[rankⱼ], mother_rank)
                        coefficient += ϕ[rankⱼ][mother_rank] / 2
                    end
                else
                    if haskey(ϕ[mother_rank], rankⱼ)
                        coefficient += ϕ[mother_rank][rankⱼ] / 2
                    end
                end
            end
            if coefficient > 0.
                # Store the non-zero kinship with the lowest ID
                ϕ[rankⱼ][rankᵢ] = coefficient
            end
        end
        # Mark the individual as processed
        push!(ranks_to_visit, rankᵢ)
        individualᵢ.founder_index = 1
        individualᵢ.children_to_process = length(individualᵢ.children)
        # If all of a parent's children are processed, we can remove the parent
        if father_rank != 0
            if !individualᵢ.father.is_proband
                individualᵢ.father.children_to_process -= 1
                if individualᵢ.father.children_to_process == 0
                    delete!(ranks_to_visit, father_rank)
                    empty!(ϕ[father_rank])
                    delete!(ϕ, father_rank)
                    for rankⱼ ∈ ranks_to_visit
                        if rankⱼ < father_rank
                            delete!(ϕ[rankⱼ], father_rank)
                        end
                    end
                end
            end
        end
        if mother_rank != 0
            if !individualᵢ.mother.is_proband
                individualᵢ.mother.children_to_process -= 1
                if individualᵢ.mother.children_to_process == 0
                    delete!(ranks_to_visit, mother_rank)
                    empty!(ϕ[mother_rank])
                    delete!(ϕ, mother_rank)
                    for rankⱼ ∈ ranks_to_visit
                        if rankⱼ < mother_rank
                            delete!(ϕ[rankⱼ], mother_rank)
                        end
                    end
                end
            end
        end
        for child ∈ individualᵢ.children
            if !isnothing(child.father) && !isnothing(child.mother)
                if child.father.founder_index != 0 && child.mother.founder_index != 0
                    push!(queue, child)
                end
            else
                push!(queue, child)
            end
        end
        rank += 1
    end
    ID_to_rank = Dict{Int32, Int}()
    for ID ∈ probandIDs
        ID_to_rank[ID] = indexed_pedigree[ID].rank
    end
    KinshipMatrix(ϕ, ID_to_rank)
end

"""
    function phiMean(ϕ::Matrix{Float32})

Return the mean kinship from a given kinship matrix.
"""
function phiMean(ϕ::Matrix{Float32})::Float32
    total = sum(ϕ)
    diagonal = sum([ϕ[i, i] for i ∈ axes(ϕ, 1)])
    total -= diagonal
    total / (length(ϕ) - size(ϕ, 1))
end

"""
    function phiMean(ϕ::KinshipMatrix)

Return the mean kinship from a given sparse kinship matrix.
"""
function phiMean(ϕ::KinshipMatrix)::Float32
    total = sum(sum(values(kinships)) for kinships ∈ values(ϕ.dict))
    diagonal = sum([ϕ.dict[ID][ID] for ID ∈ keys(ϕ.dict)])
    total -= diagonal
    count = length(ϕ.dict) * (length(ϕ.dict) - 1) / 2
    total / count
end

"""
    function _lowest_founders(pedigree::Pedigree)

Return the lowest founders of a given pedigree.

The lowest founder is defined as an only child who is either a founder or as the only child
of a lineage of only children.
"""
function _lowest_founders(pedigree::Pedigree)
    founderIDs = founder(pedigree)
    deepestIDs = Int[]
    for founderID ∈ founderIDs
        deepest = pedigree[founderID]
        while length(deepest.children) == 1
            deepest = first(deepest.children)
        end
        push!(deepestIDs, deepest.ID)
    end
    sort(unique(deepestIDs))
end

"""
    f(pedigree::Pedigree, IDs::Vector{Int})

Return the coefficients of inbreeding of a vector of individuals.
"""
function f(pedigree::Pedigree, IDs::Vector{Int})
    coefficients = Float32[]
    for ID ∈ IDs
        individual = pedigree[ID]
        if isnothing(individual.father) || isnothing(individual.mother)
            push!(coefficients, 0.)
        else
            push!(coefficients, phi(individual.father, individual.mother))
        end
    end
    coefficients
end

"""
    mutable struct PossibleDescendant <: AbstractIndividual

An individual with its ancestor's contribution.
"""
mutable struct PossibleDescendant <: AbstractIndividual
    ID::Int
    father::Union{Nothing, PossibleDescendant}
    mother::Union{Nothing, PossibleDescendant}
    children::Vector{PossibleDescendant}
    contribution::Float32
end

"""
    _contribute!(individual::Individual, depth::Int = 0)

Recursively compute the genetic contributions of an individual.
"""
function _contribute!(individual::PossibleDescendant, depth::Int = 0)
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
    gc(pedigree::Pedigree; pro::Vector{Int} = pro(pedigree),
        ancestors::Vector{Int} = founder(pedigree))

Return a matrix of the genetic contribution of each ancestor (columns) to each proband
(rows).

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
    pro::Vector{Int} = pro(pedigree),
    ancestors::Vector{Int} = founder(pedigree))
    
    # Ported from GENLIB's Congen
    matrix = zeros(Float32, length(pro), length(ancestors))
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
            push!(contribution_pedigree[father.ID].children,
                contribution_pedigree[individual.ID])
        end
        if !isnothing(mother)
            push!(contribution_pedigree[mother.ID].children,
                contribution_pedigree[individual.ID])
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
