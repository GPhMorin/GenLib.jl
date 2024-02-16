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
    is_ancestors = fill(false, length(pedigree))
    occurrences = fill(0, length(pedigree))
    for ID in ancestors
        ancestor = pedigree[ID]
        is_ancestors[ancestor.index] = true
    end
    for (j, probandID) in enumerate(pro)
        proband = pedigree[probandID]
        occur!(proband, is_ancestors, occurrences)
        for (i, ancestorID) in enumerate(ancestors)
            ancestor = pedigree[ancestorID]
            occurrence_matrix[i, j] = occurrences[ancestor.index]
            occurrences[ancestor.index] = 0
        end
    end
    if typeOcc == "IND"
        return occurrence_matrix
    elseif typeOcc == "TOTAL"
        return sum(occurrence_matrix, dims=2)
    end
end

"""
    occur!(individual::Individual)

Recursively increment the occurrence of an `individual` if they are an ancestor.
"""
function occur!(individual::Individual, is_ancestors::Vector{Bool}, occurrences::Vector{Int64})
    if is_ancestors[individual.index]
        occurrences[individual.index] += 1
    end
    if !isnothing(individual.father)
        occur!(individual.father, is_ancestors, occurrences)
    end
    if !isnothing(individual.mother)
        occur!(individual.mother, is_ancestors, occurrences)
    end
end