"""
    mutable struct IndexedIndividual <: AbstractIndividual
        ID::Int64
        father::Union{Nothing, IndexedIndividual}
        mother::Union{Nothing, IndexedIndividual}
        sex::Int64
        children::Vector{IndexedIndividual}
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

Return the kinship coefficient between two individuals given a matrix of the founders'
kinships.

Adapted from [Karigl, 1981](@ref), and [Kirkpatrick et al., 2019](@ref).
"""
function phi(individualᵢ::IndexedIndividual, individualⱼ::IndexedIndividual,
    Ψ::Matrix{Float64})
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
            individual.sex, [], individual.rank, 0
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
    _next_generation(pedigree::Pedigree, previous_generationIDs::Vector{Int64})

Return the next generation of a given set of individuals.
"""
function _previous_generation(pedigree::Pedigree, next_generationIDs::Vector{Int64})
    previous_generationIDs = Int64[]
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
    phi(pedigree::Pedigree, probandIDs::Vector{Int64} = pro(pedigree);
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
function phi(pedigree::Pedigree, probandIDs::Vector{Int64} = pro(pedigree);
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
    if compute
        # Add the `founder_index` attribute to the individuals and make them mutable so we can
        # quickly track the location of the founders in their kinship matrix.
        indexed_pedigree = _index_pedigree(pedigree)
        # Initialize the kinship matrix of the top founders
        Ψ = zeros(length(cut_vertices[1]), length(cut_vertices[1]))
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
            ϕ = Matrix{Float64}(undef, length(next_generationIDs), length(next_generationIDs))
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
end

"""
    struct KinshipMatrix

A minimal structure wrapping an `Dict` with kinships of individuals accessed by IDs.
"""
struct KinshipMatrix
    dict::Dict{Int64, Dict{Int64, Float64}}
end

KinshipMatrix() = KinshipMatrix(Dict{Int64, Dict{Int64, Float64}}())

Base.length(ϕ::KinshipMatrix) = length(ϕ.dict)
Base.keys(ϕ::KinshipMatrix) = keys(ϕ.dict)
Base.values(ϕ::KinshipMatrix) = values(ϕ.dict)
Base.iterate(ϕ::KinshipMatrix) = iterate(ϕ.dict)
Base.iterate(ϕ::KinshipMatrix, i) = iterate(ϕ.dict, i)
Base.pop!(ϕ::KinshipMatrix) = pop!(ϕ.dict)

function Base.setindex!(ϕ::KinshipMatrix, value::Float64, ID₁::Int64, ID₂::Int64)
    if !haskey(ϕ.dict, ID₁)
        ϕ.dict[ID₁] = Dict{Int64, Float64}()
    end
    ϕ.dict[ID₁][ID₂] = value
end

function Base.getindex(ϕ::KinshipMatrix, ID₁::Int64, ID₂::Int64)
    (ID₁, ID₂) = ID₁ < ID₂ ? (ID₁, ID₂) : (ID₂, ID₁)
    if haskey(ϕ.dict, ID₁)
        return haskey(ϕ.dict[ID₁], ID₂) ? ϕ.dict[ID₁][ID₂] : 0.
    else
        return 0.
    end
end

function Base.delete!(ϕ::KinshipMatrix, ID::Int64)
    empty!(ϕ.dict[ID])
    delete!(ϕ.dict, ID)
    for (ID₂, kinships) ∈ ϕ.dict
        if ID₂ < ID
            delete!(kinships, ID)
        end
    end
end

function Base.show(io::IO, ::MIME"text/plain", ϕ::KinshipMatrix)
    nz = 0
    for kinships ∈ values(ϕ)
        nz += length(kinships)
    end
    print(io, "$(length(ϕ))×$(length(ϕ)) KinshipMatrix with $nz stored entries.")
end

"""
    phi(individualᵢ::IndexedIndividual, individualⱼ::IndexedIndividual,
    ϕ::Dict{Int64, Vector{Pair{Int64, Float64}}})

Return the kinship coefficient between two individuals given a dictionary of the
individuals' kinships.

Adapted from [Karigl, 1981](@ref), and [Kirkpatrick et al., 2019](@ref).
"""
function phi(individualᵢ::IndexedIndividual, individualⱼ::IndexedIndividual,
    ϕ::KinshipMatrix)
    value = 0.
    if individualᵢ.founder_index != 0 && individualⱼ.founder_index != 0
        # Both individuals are founders, so we already know their kinship coefficient
        value += ϕ[individualᵢ.ID, individualⱼ.ID]
    elseif individualᵢ.founder_index != 0
        # Individual i is a founder, so we climb the pedigree on individual j's side
        if !isnothing(individualⱼ.father)
            value += phi(individualᵢ, individualⱼ.father, ϕ) / 2
        end
        if !isnothing(individualⱼ.mother)
            value += phi(individualᵢ, individualⱼ.mother, ϕ) / 2
        end
    elseif individualⱼ.founder_index != 0
        # Individual j is a founder, so we climb the pedigree on individual i's side
        if !isnothing(individualᵢ.father)
            value += phi(individualⱼ, individualᵢ.father, ϕ) / 2
        end
        if !isnothing(individualᵢ.mother)
            value += phi(individualⱼ, individualᵢ.mother, ϕ) / 2
        end
    else
        # None of the individuals are founders, so we climb on the side of the individual
        # that appears lowest in the pedigree
        if individualᵢ.rank > individualⱼ.rank
            # From the genealogical order, i cannot be an ancestor of j
            # Φᵢⱼ = (Φₚⱼ + Φₘⱼ) / 2, if i is not an ancestor of j (Karigl, 1981)
            if !isnothing(individualᵢ.father)
                value += phi(individualᵢ.father, individualⱼ, ϕ) / 2
            end
            if !isnothing(individualᵢ.mother)
                value += phi(individualᵢ.mother, individualⱼ, ϕ) / 2
            end
        elseif individualⱼ.rank > individualᵢ.rank
            # Reverse the order since i can be an ancestor of j
            # Φⱼᵢ = (Φₚᵢ + Φₘᵢ) / 2, if j is not an ancestor of i (Karigl, 1981)
            if !isnothing(individualⱼ.father)
                value += phi(individualⱼ.father, individualᵢ, ϕ) / 2
            end
            if !isnothing(individualⱼ.mother)
                value += phi(individualⱼ.mother, individualᵢ, ϕ) / 2
            end
        elseif individualᵢ.rank == individualⱼ.rank
            # Same individual
            # Φₐₐ = (1 + Φₚₘ) / 2 (Karigl, 1981)
            value += 0.5
            if !isnothing(individualᵢ.father) && !isnothing(individualᵢ.mother)
                value += phi(individualᵢ.father, individualᵢ.mother, ϕ) / 2
            end
        end
    end
    value
end

"""
    sparse_phi(pedigree::Pedigree, probandIDs::Vector{Int64} = pro(pedigree);
        verbose::Bool = false)

Return a sparse matrix of pairwise kinship coefficients between probands.

If no probands are given, return the sparse matrix for all probands in the pedigree.

The algorithm is a hybrid between the algorithms of
[Karigl, 1981](@ref), and [Kirkpatrick et al., 2019](@ref).

If `verbose` is `true`, print the information about the cut vertices. If `compute` is
`true` (the default), compute the kinship matrix. If it is `false`, only print the
information about the cut vertices.

# Example

```julia
import GenLib as gen
geneaJi = gen.geneaJi
ped = gen.genealogy(geneaJi)
gen.sparse_phi(ped)
```
"""
function sparse_phi(pedigree::Pedigree, probandIDs::Vector{Int64} = pro(pedigree);
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
    if compute
        # Add the `founder_index` attribute to the individuals and make them mutable so we can
        # quickly track the location of the founders in their kinship vector.
        indexed_pedigree = _index_pedigree(pedigree)
        # Initialize the kinship vector of the top founders
        ϕ = KinshipMatrix()
        for ID ∈ cut_vertices[1]
            ϕ[ID, ID] = 0.5
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
            for ID ∈ previous_generationIDs
                indexed_pedigree[ID].founder_index = 1
            end
            # Make a parallel copy of the kinships
            ϕs = [KinshipMatrix() for _ ∈ 1:Threads.nthreads()]
            # Fill the dictionary in parallel, using the adapted algorithm from Karigl, 1981
            individuals = [indexed_pedigree[ID] for ID ∈ next_generationIDs]
            Threads.@threads :static for i ∈ eachindex(individuals)
                j = Threads.threadid()
                individualᵢ = individuals[i]
                for individualⱼ ∈ individuals
                    if individualᵢ.ID ≤ individualⱼ.ID &&
                        (individualᵢ.founder_index == 0 || individualⱼ.founder_index == 0)
                        coefficient = phi(individualᵢ, individualⱼ, ϕ)
                        if coefficient > 0
                            ϕs[j][individualᵢ.ID, individualⱼ.ID] = coefficient
                        end
                    end
                end
            end
            # Remove the individuals that are not probands
            non_recurringIDs = setdiff(previous_generationIDs, next_generationIDs)
            for ID ∈ non_recurringIDs
                delete!(ϕ, ID)
            end
            # Merge the kinships
            while !(isempty(ϕs))
                ϕᵢ = pop!(ϕs)
                while !(isempty(ϕᵢ))
                    (IDᵢ, kinships) = pop!(ϕᵢ)
                    if !haskey(ϕ.dict, IDᵢ)
                        ϕ.dict[IDᵢ] = Dict{Int64, Float64}()
                    end
                    while !(isempty(kinships))
                        (IDⱼ, coefficient) = pop!(kinships)
                        ϕ.dict[IDᵢ][IDⱼ] = coefficient
                    end
                end
            end
        end
        ϕ
    end
end

"""
    function phiMean(ϕ::Matrix{Float64})

Return the mean kinship from a given kinship matrix.
"""
function phiMean(ϕ::Matrix{Float64})
    total = sum(ϕ)
    diagonal = sum([ϕ[i, i] for i ∈ axes(ϕ, 1)])
    total -= diagonal
    total / (length(ϕ) - size(ϕ, 1))
end

"""
    function phiMean(ϕ::KinshipMatrix)

Return the mean kinship from a given sparse kinship matrix.
"""
function phiMean(ϕ::KinshipMatrix)
    total = 0.
    for (ID, kinships) ∈ ϕ
        for coefficient ∈ values(kinships)
            total += coefficient
        end
        total -= kinships[ID]
    end
    count = length(ϕ) * (length(ϕ) - 1) / 2
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
    f(pedigree::Pedigree, IDs::Vector{Int64})

Return the coefficients of inbreeding of a vector of individuals.
"""
function f(pedigree::Pedigree, IDs::Vector{Int64})
    coefficients = Float64[]
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
    gc(pedigree::Pedigree; pro::Vector{Int64} = pro(pedigree),
        ancestors::Vector{Int64} = founder(pedigree))

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
