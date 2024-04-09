"""
    nomen(pedigree::Pedigree)

Return the number of men in the pedigree.
"""
function nomen(pedigree::Pedigree)
    number = 0
    for individual ∈ values(pedigree)
        if individual.sex == 1
            number += 1
        end
    end
    number
end

"""
    nowomen(pedigree::Pedigree)

Return the number of women in the pedigree.
"""
function nowomen(pedigree::Pedigree)
    number = 0
    for individual ∈ values(pedigree)
        if individual.sex == 2
            number += 1
        end
    end
    number
end

"""
    noind(pedigree::Pedigree)

Return the number of individuals in the pedigree.
"""
noind(pedigree::Pedigree) = length(pedigree)

"""
    _max_depth(individual::Individual)

Return the maximum depth of an individual's pedigree.
"""
function _max_depth(individual::Individual)
    father_depth = 1
    mother_depth = 1
    if !isnothing(individual.father)
        father_depth += _max_depth(individual.father)
    end
    if !isnothing(individual.mother)
        mother_depth += _max_depth(individual.mother)
    end
    max(father_depth, mother_depth)
end

"""
    depth(pedigree::Pedigree)

Return the maximum depth of a pedigree.
"""
function depth(pedigree::Pedigree)
    max_depth = 0
    for individual ∈ values(pedigree)
        max_depth = max(max_depth, _max_depth(individual))
    end
    max_depth
end

"""
    rec(pedigree::Pedigree, probandIDs::Vector{Int64} = pro(genealogy), ancestorIDs::Vector{Int64} = founder(genealogy))

Return the number of descendants of each ancestor.
"""
function rec(
    pedigree::Pedigree,
    probandIDs::Vector{Int64} = pro(pedigree),
    ancestorIDs::Vector{Int64} = founder(pedigree))
    
    coverage = Vector{Int64}()
    for ancestorID ∈ ancestorIDs
        descendantIDs = descendant(pedigree, ancestorID)
        descendantIDs = filter!(x -> x ∈ probandIDs, descendantIDs)
        append!(coverage, length(descendantIDs))
    end
    coverage
end

"""
    mutable struct Occurrent <: AbstractIndividual
        ID::Int64
        father::Union{Nothing, Occurrent}
        mother::Union{Nothing, Occurrent}
        children::Vector{Occurrent}
        sex::Int64
        rank::Int64
        is_ancestor::Bool
        occurrence::Int64
    end

An individual with a number of occurrences and whether they are an ancestor.
"""
mutable struct Occurrent <: AbstractIndividual
    ID::Int64
    father::Union{Nothing, Occurrent}
    mother::Union{Nothing, Occurrent}
    children::Vector{Occurrent}
    sex::Int64
    rank::Int64
    is_ancestor::Bool
    occurrence::Int64
end

"""
    occ(pedigree::Pedigree; pro::Vector{Int64} = pro(genealogy), ancestors::Vector{Int64} = founder(genealogy), typeOcc::String = "IND")

Return a matrix of ancestors' occurrences.

If `typeOcc` is "IND" (default), then the matrix corresponds to the occurrence per individual.
If `typeOcc` is "TOTAL", then the matrix corresponds to the total occurrence.

# Example

```julia
import GenLib as gen
geneaJi = gen.geneaJi
ped = gen.genealogy(geneaJi)
occ = gen.occ(ped, typeOcc = "TOTAL")
```
"""
function occ(
    pedigree::Pedigree;
    pro::Vector{Int64} = pro(pedigree),
    ancestors::Vector{Int64} = founder(pedigree),
    typeOcc::String = "IND")
    
    occurrence_matrix = Matrix{Int64}(undef, length(ancestors), length(pro))
    occurrence_pedigree = Pedigree{Occurrent}()
    for individual ∈ collect(values(pedigree))
        father = individual.father
        mother = individual.mother
        occurrence_pedigree[individual.ID] = Occurrent(
            individual.ID,
            isnothing(father) ? nothing : occurrence_pedigree[father.ID],
            isnothing(mother) ? nothing : occurrence_pedigree[mother.ID],
            Occurrent[],
            individual.sex,
            individual.rank,
            individual.ID ∈ ancestors ? true : false,
            0
        )
    end
    for (j, probandID) ∈ enumerate(pro)
        proband = occurrence_pedigree[probandID]
        occur!(proband)
        for (i, ancestorID) ∈ enumerate(ancestors)
            ancestor = occurrence_pedigree[ancestorID]
            occurrence_matrix[i, j] = ancestor.occurrence
            ancestor.occurrence = 0
        end
    end
    if typeOcc == "IND"
        return occurrence_matrix
    elseif typeOcc == "TOTAL"
        return sum(occurrence_matrix, dims=2)
    end
end

"""
    occur!(individual::Occurrent)

Recursively increment the occurrence of an `individual` if they are an ancestor.
"""
function occur!(individual::Occurrent)
    if individual.is_ancestor
        individual.occurrence += 1
    end
    if !isnothing(individual.father)
        occur!(individual.father)
    end
    if !isnothing(individual.mother)
        occur!(individual.mother)
    end
end

"""
    _get_paths(pedigree::Pedigree, individual::Individual)

Return the paths from an individual to their ancestors.
"""
function _get_paths(pedigree::Pedigree, individual::Individual)
    paths = Vector{Vector{Int64}}([[individual.ID]])
    if !isnothing(individual.father)
        fathers_paths = _get_paths(pedigree, individual.father)
        for path ∈ fathers_paths
            push!(path, individual.ID)
        end
        append!(paths, fathers_paths)
    end
    if !isnothing(individual.mother)
        mothers_paths = _get_paths(pedigree, individual.mother)
        for path ∈ mothers_paths
            push!(path, individual.ID)
        end
        append!(paths, mothers_paths)
    end
    paths
end

"""
    _findDistance(pedigree::Pedigree, descendantID::Int64, ancestorID::Int64)

Return a vector of distances between an individual and their ancestor.
"""
function _findDistance(
    pedigree::Pedigree,
    descendantID::Int64,
    ancestorID::Int64)
    
    descendant = pedigree[descendantID]
    paths = _get_paths(pedigree, descendant)
    lengths = Vector{Int64}()
    for path ∈ paths
        if path[1] ≡ ancestorID
            push!(lengths, length(path) - 1)
        end
    end
    lengths
end

"""
    _findMinDistance(pedigree::Pedigree, descendantID::Int64, ancestorID::Int64)

Return the minimum distance between an individual and their ancestor.
"""
function _findMinDistance(pedigree::Pedigree, descendantID::Int64, ancestorID::Int64)
    lengths = _findDistance(pedigree, descendantID, ancestorID)
    minimum(lengths)
end

"""
    findDistance(pedigree::Pedigree, IDs::Vector{Int64}, ancestorID::Int64)

Return the distance between two individuals and their ancestor.
"""
function findDistance(pedigree::Pedigree, IDs::Vector{Int64}, ancestorID::Int64)
    distance₁ = _findMinDistance(pedigree, IDs[1], ancestorID)
    distance₂ = _findMinDistance(pedigree, IDs[2], ancestorID)
    distance₁ + distance₂
end