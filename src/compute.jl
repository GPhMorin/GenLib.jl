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
    proband_index::Int
end

"""
    struct KinshipMatrix

A minimal structure wrapping an `Dict` with kinships of individuals accessed by IDs.
"""
struct KinshipMatrix
    values::Vector{Vector{Pair{Int, Float64}}}
    ID_to_index::Dict{Int, Int}
    ID_to_rank::Dict{Int, Int}
end

function Base.getindex(ϕ::KinshipMatrix, ID₁::Int, ID₂::Int)
    (ID₁, ID₂) = ϕ.ID_to_rank[ID₁] < ϕ.ID_to_rank[ID₂] ? (ID₁, ID₂) : (ID₂, ID₁)
    kinship = searchsorted(
        ϕ.values[ϕ.ID_to_index[ID₁]],
        ϕ.ID_to_rank[ID₂] => 0, by = first
    )
    isempty(kinship) ? 0 : ϕ.values[ϕ.ID_to_index[ID₁]][first(kinship)].second
end

function Base.show(io::IO, ::MIME"text/plain", ϕ::KinshipMatrix)
    nz = sum(length(kinships) for kinships ∈ ϕ.values)
    print(io, "$(length(ϕ.values))×$(length(ϕ.values)) KinshipMatrix " *
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
            individual.sex, [], individual.rank, 0, 0
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
    for (index, ID) ∈ enumerate(probandIDs)
        indexed_pedigree[ID].proband_index = index
    end
    # Initialize the kinship matrix
    ϕ = Vector{Vector{Pair{Int, Float32}}}()
    index_to_ID = Dict{Int, Int}()
    # Initialize the queue with the probands
    IDs = Set{Int}()
    queue = founder(pedigree)
    # Print information about the steps of the algorithm
    depths = [_max_depth(individual) for individual ∈ values(isolated_pedigree)]
    n_per_depth = [count(i -> i == depth, depths) for depth ∈ depths]
    max_depth = maximum(depths)
    for depth ∈ unique(depths)
        println("Depth $depth / $max_depth: $(count(i -> i == depth, depths)) individuals")
    end
    current_depth = 0
    # Reassign a rank to the individuals
    rank = 1
    while !isempty(queue)
        IDᵢ = popfirst!(queue)
        # Print information about the steps of the algorithm
        depth = popfirst!(depths)
        n = popfirst!(n_per_depth)
        if depth != current_depth
            current_depth = depth
            println("Running: $current_depth / $max_depth ($n individuals)")
        end
        individualᵢ = indexed_pedigree[IDᵢ]
        # Initialize the kinship array for the individual
        push!(ϕ, Vector{Pair{Int, Float32}}())
        # The founder index is used to retrieve the kinship of the individual
        individualᵢ.founder_index = lastindex(ϕ)
        index_to_ID[lastindex(ϕ)] = IDᵢ
        # Assign a new rank to the individual
        individualᵢ.rank = rank
        rank += 1
        father = individualᵢ.father
        mother = individualᵢ.mother
        # Kinship with previous individuals
        for IDⱼ ∈ IDs
            individualⱼ = indexed_pedigree[IDⱼ]
            coefficient = 0.
            if !isnothing(father)
                # In order to make the kinships as sparse as possible,
                # we only store the kinship with the lowest ranked founder
                # and the kinship is inserted in rank order for faster lookup
                (individual₁, individual₂) = father.rank < individualⱼ.rank ?
                    (father, individualⱼ) : (individualⱼ, father)
                kinship = searchsorted(
                    ϕ[individual₁.founder_index],
                    individual₂.rank => 0,
                    by = first
                )
                if !isempty(kinship)
                    coefficient += ϕ[individual₁.founder_index][first(kinship)].second / 2
                end
            end
            if !isnothing(mother)
                # Same thing but on the mother's side
                (individual₁, individual₂) = mother.rank < individualⱼ.rank ?
                    (mother, individualⱼ) : (individualⱼ, mother)
                kinship = searchsorted(
                    ϕ[individual₁.founder_index],
                    individual₂.rank => 0,
                    by = first
                )
                if !isempty(kinship)
                    coefficient += ϕ[individual₁.founder_index][first(kinship)].second / 2
                end
            end
            if coefficient > 0.
                # Store the non-zero kinship with the lowest ranked individual
                push!(ϕ[individualⱼ.founder_index], individualᵢ.rank => coefficient)
            end
        end
        # Kinship with self
        coefficient = 0.5
        if !isnothing(father) && !isnothing(mother)
            (individual₁, individual₂) = father.rank < mother.rank ?
                (father, mother) : (mother, father)
            kinship = searchsorted(
                ϕ[individual₁.founder_index],
                individual₂.rank => 0,
                by = first
            )
            if !isempty(kinship)
                coefficient += ϕ[individual₁.founder_index][first(kinship)].second / 2
            end
        end
        push!(ϕ[individualᵢ.founder_index], individualᵢ.rank => coefficient)
        # If all of a parent's children are processed, we can remove the parent
        if !isnothing(father) && father.proband_index == 0
            if all(child.founder_index != 0 for child ∈ father.children)
                delete!(IDs, father.ID)
                last_ID = index_to_ID[lastindex(ϕ)]
                individual = indexed_pedigree[last_ID]
                new_index = father.founder_index
                father.founder_index = -1
                kinships = pop!(ϕ)
                if individual.ID != father.ID
                    ϕ[new_index] = kinships
                    individual.founder_index = new_index
                    index_to_ID[new_index] = last_ID
                end
                for ID ∈ IDs
                    individual = indexed_pedigree[ID]
                    if individual.rank < father.rank
                        kinship = searchsorted(
                            ϕ[individual.founder_index],
                            father.rank => 0,
                            by = first
                        )
                        if !isempty(kinship)
                            deleteat!(ϕ[individual.founder_index], first(kinship))
                        end
                    end
                end
            end
        end
        if !isnothing(mother) && mother.proband_index == 0
            if all(child.founder_index != 0 for child ∈ mother.children)
                delete!(IDs, mother.ID)
                last_ID = index_to_ID[lastindex(ϕ)]
                individual = indexed_pedigree[last_ID]
                new_index = mother.founder_index
                mother.founder_index = -1
                kinships = pop!(ϕ)
                if individual.ID != mother.ID
                    ϕ[new_index] = kinships
                    individual.founder_index = new_index
                    index_to_ID[new_index] = last_ID
                end
                for ID ∈ IDs
                    individual = indexed_pedigree[ID]
                    if individual.rank < mother.rank
                        kinship = searchsorted(
                            ϕ[individual.founder_index],
                            mother.rank => 0,
                            by = first
                        )
                        if !isempty(kinship)
                            deleteat!(ϕ[individual.founder_index], first(kinship))
                        end
                    end
                end
            end
        end
        # Add the individual to the list of IDs that are already processed
        push!(IDs, individualᵢ.ID)
        for child ∈ individualᵢ.children
            if !isnothing(child.father) && !isnothing(child.mother)
                if child.father.founder_index != 0 && child.mother.founder_index != 0
                    push!(queue, child.ID)
                end
            else
                push!(queue, child.ID)
            end
        end
    end
    ID_to_index = Dict{Int, Int}(
        individual.ID => individual.founder_index
        for individual ∈ values(indexed_pedigree)
    )
    ID_to_rank = Dict{Int, Int}(
        individual.ID => individual.rank
        for individual ∈ values(indexed_pedigree)
    )
    KinshipMatrix(ϕ, ID_to_index, ID_to_rank)
end

"""
    function phiMean(ϕ::Matrix{Float32})

Return the mean kinship from a given kinship matrix.
"""
function phiMean(ϕ::Matrix{Float32})
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
    total = sum(sum(kinship.second for kinship ∈ kinships) for kinships ∈ ϕ.values)
    diagonal = sum([first(kinships).second for kinships ∈ ϕ.values])
    total -= diagonal
    count = length(ϕ.values) * (length(ϕ.values) - 1) / 2
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
