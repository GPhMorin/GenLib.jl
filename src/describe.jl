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
    _completeness!(completeness::Vector{Int}, individual::Individual, depth::Int)

Return the completeness of an individual filled recursively.
"""
function _completeness!(completeness::Vector{Int}, individual::Individual, depth::Int)
    if !isnothing(individual.father) || !isnothing(individual.mother)
        if length(completeness) < depth + 1
            push!(completeness, 0)
        end
    end
    if !isnothing(individual.father)
        completeness[depth + 1] += 1
        _completeness!(completeness, individual.father, depth + 1)
    end
    if !isnothing(individual.mother)
        completeness[depth + 1] += 1
        _completeness!(completeness, individual.mother, depth + 1)
    end
    return completeness
end

"""
    completeness(pedigree::Pedigree; pro::Vector{Int} = pro(pedigree),
    genNo::Vector{Int} = Int[], type::String = "MEAN")

Return a dataframe with the completeness at each generation (one row per generation).

genNo: A vector of the generations to output. The probands are at generation 0.

type: If ```"MEAN"```, the mean completeness for each generation. If ```"IND"```, the
completeness for each generation for each proband.
"""
function completeness(pedigree::Pedigree, pro::Vector{Int} = pro(pedigree);
    genNo::Vector{Int} = Int[], type::String = "MEAN")
    completenesses = Vector{Vector{Int}}()
    for ID ∈ pro
        proband = pedigree[ID]
        push!(completenesses, _completeness!(Int[], proband, 0))
    end
    max_depth = maximum([length(completeness) for completeness ∈ completenesses])
    matrix = zeros(max_depth, length(pro))
    for column ∈ eachindex(pro)
        completeness = completenesses[column]
        for row ∈ eachindex(completeness)
            matrix[row, column] = completeness[row] ./ (2^row) * 100
        end
    end
    matrix = vcat(ones(1, length(pro)) .* 100, matrix)
    if !isempty(genNo)
        matrix = matrix[genNo .+ 1, :]
    end
    if type == "IND"
        return matrix
    elseif type == "MEAN"
        return sum(matrix, dims=2) ./ length(pro)
    end
end

"""
    rec(pedigree::Pedigree, probandIDs::Vector{Int} = pro(genealogy),
    ancestorIDs::Vector{Int} = founder(genealogy))

Return the number of descendants of each ancestor.
"""
function rec(
    pedigree::Pedigree,
    probandIDs::Vector{Int} = pro(pedigree),
    ancestorIDs::Vector{Int} = founder(pedigree))
    
    coverage = Vector{Int}()
    for ancestorID ∈ ancestorIDs
        descendantIDs = descendant(pedigree, ancestorID)
        descendantIDs = filter!(x -> x ∈ probandIDs, descendantIDs)
        append!(coverage, length(descendantIDs))
    end
    coverage
end

"""
    mutable struct Occurrent <: AbstractIndividual
        ID::Int
        father::Union{Nothing, Occurrent}
        mother::Union{Nothing, Occurrent}
        is_ancestor::Bool
        occurrence::Int
    end

An individual with a number of occurrences and whether they are an ancestor.
"""
mutable struct Occurrent <: AbstractIndividual
    ID::Int
    father::Union{Nothing, Occurrent}
    mother::Union{Nothing, Occurrent}
    is_ancestor::Bool
    occurrence::Int
end

"""
    occ(pedigree::Pedigree; pro::Vector{Int} = pro(genealogy),
    ancestors::Vector{Int} = founder(genealogy), typeOcc::String = "IND")

Return a matrix of ancestors' occurrences.

If `typeOcc` is "IND" (default), then the matrix corresponds to the occurrence per
individual. If `typeOcc` is "TOTAL", then the matrix corresponds to the total occurrence.

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
    pro::Vector{Int} = pro(pedigree),
    ancestors::Vector{Int} = founder(pedigree),
    typeOcc::String = "IND")
    
    occurrence_matrix = Matrix{Int}(undef, length(ancestors), length(pro))
    occurrence_pedigree = Pedigree{Occurrent}()
    for individual ∈ collect(values(pedigree))
        father = individual.father
        mother = individual.mother
        occurrence_pedigree[individual.ID] = Occurrent(
            individual.ID,
            isnothing(father) ? nothing : occurrence_pedigree[father.ID],
            isnothing(mother) ? nothing : occurrence_pedigree[mother.ID],
            individual.ID ∈ ancestors ? true : false,
            0
        )
    end
    for (j, probandID) ∈ enumerate(pro)
        proband = occurrence_pedigree[probandID]
        _occur!(proband)
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
function _occur!(individual::Occurrent)
    if individual.is_ancestor
        individual.occurrence += 1
    end
    if !isnothing(individual.father)
        _occur!(individual.father)
    end
    if !isnothing(individual.mother)
        _occur!(individual.mother)
    end
end

"""
    _get_paths(pedigree::Pedigree, individual::Individual)

Return the paths from an individual to their ancestors.
"""
function _get_paths(pedigree::Pedigree, individual::Individual)
    paths = Vector{Vector{Int}}([[individual.ID]])
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
    _findDistance(pedigree::Pedigree, descendantID::Int, ancestorID::Int)

Return a vector of distances between an individual and their ancestor.
"""
function _findDistance(
    pedigree::Pedigree,
    descendantID::Int,
    ancestorID::Int)
    
    descendant = pedigree[descendantID]
    paths = _get_paths(pedigree, descendant)
    lengths = Vector{Int}()
    for path ∈ paths
        if path[1] ≡ ancestorID
            push!(lengths, length(path) - 1)
        end
    end
    lengths
end

"""
    _findMinDistance(pedigree::Pedigree, descendantID::Int, ancestorID::Int)

Return the minimum distance between an individual and their ancestor.
"""
function _findMinDistance(pedigree::Pedigree, descendantID::Int, ancestorID::Int)
    lengths = _findDistance(pedigree, descendantID, ancestorID)
    minimum(lengths)
end

"""
    findDistance(pedigree::Pedigree, IDs::Vector{Int}, ancestorID::Int)

Return the distance between two individuals and their ancestor.
"""
function findDistance(pedigree::Pedigree, IDs::Vector{Int}, ancestorID::Int)
    distance₁ = _findMinDistance(pedigree, IDs[1], ancestorID)
    distance₂ = _findMinDistance(pedigree, IDs[2], ancestorID)
    distance₁ + distance₂
end