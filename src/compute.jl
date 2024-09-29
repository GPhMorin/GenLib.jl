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
    kinship::Vector{Pair{Int32, Float64}}
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
            individual.sex, [], individual.rank, 0, Pair{Int32, Float64}[]
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
    phi(individualᵢ::IndexedIndividual, individualⱼ::IndexedIndividual,
    ϕ::Dict{Int64, Vector{Pair{Int64, Float64}}})

Return the kinship coefficient between two individuals given a dictionary of the
individuals' kinships.

Adapted from [Karigl, 1981](@ref), and [Kirkpatrick et al., 2019](@ref).
"""
function phi(individualᵢ::IndexedIndividual, individualⱼ::IndexedIndividual,
    ϕ::Dict{Int64, Vector{Pair{Int64, Float64}}})
    value = 0.
    (individualᵢ, individualⱼ) = individualᵢ.ID ≤ individualⱼ.ID ?
        (individualᵢ, individualⱼ) : (individualⱼ, individualᵢ)
    if individualᵢ.founder_index != 0 && individualⱼ.founder_index != 0
        # Both individuals are founders, so we already know their kinship coefficient
        if haskey(ϕ, individualᵢ.ID)
            slice = searchsorted(ϕ[individualᵢ.ID], individualⱼ.ID => 0, by = first)
            if !isempty(slice)
                value += ϕ[individualᵢ.ID][slice[1]].second
            end
        end
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
    probands_sparse_phi(pedigree::Pedigree, probandIDs::Vector{Int64} = pro(pedigree);
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
gen.probands_sparse_phi(ped)
```
"""
function probands_sparse_phi(pedigree::Pedigree, probandIDs::Vector{Int64} = pro(pedigree);
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
        ϕ = Dict{Int64, Vector{Pair{Int64, Float64}}}()
        for ID ∈ cut_vertices[1]
            ϕ[ID] = [ID => 0.5]
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
            ϕs = [Dict{Int64, Vector{Pair{Int64, Float64}}}() for _ ∈ 1:Threads.nthreads()]
            # Fill the dictionary in parallel, using the adapted algorithm from Karigl, 1981
            individuals = [indexed_pedigree[ID] for ID ∈ sort(next_generationIDs)]
            Threads.@threads :static for i ∈ eachindex(individuals)
                individualᵢ = individuals[i]
                has_key = false
                for individualⱼ ∈ individuals
                    if individualᵢ.ID ≤ individualⱼ.ID
                        coefficient = phi(individualᵢ, individualⱼ, ϕ)
                        if coefficient > 0
                            if has_key
                                push!(ϕs[Threads.threadid()][individualᵢ.ID],
                                    individualⱼ.ID => coefficient)
                            else
                                ϕs[Threads.threadid()][individualᵢ.ID] =
                                    [individualⱼ.ID => coefficient]
                                has_key = true
                            end
                        end
                    end
                end
            end
            # Merge the dictionaries
            empty!(ϕ)
            for ϕᵢ ∈ ϕs
                while !isempty(ϕᵢ)
                    (ID, kinships) = pop!(ϕᵢ)
                    ϕ[ID] = kinships
                end
            end
        end
        if verbose
            println("Transforming the kinships into a sparse CSC matrix.")
        end
        # Convert ID to index
        ID_to_index = Dict{Int64, Int64}()
        for (index, ID) ∈ enumerate(probandIDs)
            ID_to_index[ID] = index
        end
        # Convert the dictionary to a sparse CSC matrix
        matrix = SparseMatrixCSC{Float64, Int64}(
            length(probandIDs), 
            length(probandIDs),
            [1 for _ ∈ 1:length(probandIDs) + 1],
            Int64[],
            Float64[]
        )
        while(!isempty(ϕ))
            (IDᵢ, kinships) = pop!(ϕ)
            while(!isempty(kinships))
                (IDⱼ, value) = pop!(kinships)
                matrix[ID_to_index[IDᵢ], ID_to_index[IDⱼ]] = value
            end
        end
        matrix
    end
end

"""
    complete_sparse_phi(pedigree::Pedigree; verbose::Bool = false)

Return a sparse matrix of pairwise kinship coefficients between all individuals of the given
pedigree.

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
gen.complete_sparse_phi(ped)
```
"""
function complete_sparse_phi(pedigree::Pedigree)    
    # Index the individuals with integrated kinships
    indexed_pedigree = _index_pedigree(pedigree)
    # Fill the vectors, using the adapted algorithm from Kirkpatrick, 2019
    total_count = length(pedigree) * (length(pedigree) - 1) / 2 + length(pedigree)
    decile = ceil(total_count / 10)
    count = 0
    percentage = 0
    for individualᵢ ∈ values(indexed_pedigree)
        for individualⱼ ∈ values(indexed_pedigree)
            if individualᵢ.rank ≥ individualⱼ.rank
                count += 1
                if count % decile == 0
                    percentage += 10
                    println("$percentage% processed.")
                end
            end
            if individualᵢ.rank > individualⱼ.rank
                coefficient = 0.
                if !isnothing(individualᵢ.father)
                    i = individualⱼ.rank < individualᵢ.father.rank ?
                        individualⱼ : individualᵢ.father
                    j = individualⱼ.rank > individualᵢ.father.rank ?
                        individualⱼ : individualᵢ.father
                    slice = searchsorted(i.kinship, j.rank => 0, by = first)
                    if !isempty(slice)
                        coefficient += i.kinship[slice[1]].second / 2
                    end
                end
                if !isnothing(individualᵢ.mother)
                    i = individualⱼ.rank < individualᵢ.mother.rank ?
                        individualⱼ : individualᵢ.mother
                    j = individualⱼ.rank > individualᵢ.mother.rank ?
                        individualⱼ : individualᵢ.mother
                    slice = searchsorted(i.kinship, j.rank => 0, by = first)
                    if !isempty(slice)
                        coefficient += i.kinship[slice[1]].second / 2
                    end
                end
                if coefficient > 0
                    push!(individualⱼ.kinship, individualᵢ.rank => coefficient)
                end
            elseif individualᵢ.rank == individualⱼ.rank
                coefficient = 0.5
                if !isnothing(individualᵢ.father) && !isnothing(individualᵢ.mother)
                    parentᵢ = individualᵢ.father.rank < individualᵢ.mother.rank ?
                        individualᵢ.father : individualᵢ.mother
                    parentⱼ = individualᵢ.father.rank > individualᵢ.mother.rank ?
                        individualᵢ.father : individualᵢ.mother
                    slice = searchsorted(parentᵢ.kinship, parentⱼ.rank => 0, by = first)
                    if !isempty(slice)
                        coefficient += parentᵢ.kinship[slice[1]].second / 2
                    end
                end
                push!(individualᵢ.kinship, individualᵢ.rank => coefficient)
            else
                break
            end
        end
    end
    println("Transforming the kinships into a sparse CSC matrix.")
    # Convert the vector to COO matrix
    rows = Int64[]
    cols = Int64[]
    vals = Float64[]
    for individual ∈ values(indexed_pedigree)
        while(!isempty(individual.kinship))
            kinship = popfirst!(individual.kinship)
            push!(rows, individual.rank)
            push!(cols, kinship.first)
            push!(vals, kinship.second)
        end
    end
    # Convert COO format to CSC sparse matrix
    sparse(rows, cols, vals)
end

"""
    function phiMean(phiMatrix::Matrix{Float64})

Return the mean kinship from a given kinship matrix.
"""
function phiMean(phiMatrix::Matrix{Float64})
    total = sum(phiMatrix)
    diagonal = sum([phiMatrix[i, i] for i ∈ axes(phiMatrix, 1)])
    total -= diagonal
    total / (length(phiMatrix) - size(phiMatrix, 1))
end

"""
    function phiMean(phiMatrix::SparseMatrixCSC{Float64, Int64})

Return the mean kinship from a given sparse kinship matrix.
"""
function phiMean(phiMatrix::SparseMatrixCSC{Float64, Int64})
    total = sum(phiMatrix)
    diagonal = sum([phiMatrix[i, i] for i ∈ axes(phiMatrix, 1)])
    total -= diagonal
    total * 2 / (length(phiMatrix) - size(phiMatrix, 1))
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
