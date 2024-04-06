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