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
    for (ID, individual) in pedigree
        if ID in ancestors
            individual.stats["is_ancestor"] = true
            individual.stats["occurrence"] = 0
        else
            individual.stats["is_ancestor"] = false
        end
    end
    for (j, probandID) in enumerate(pro)
        proband = pedigree[probandID]
        occur!(proband)
        for (i, ancestorID) in enumerate(ancestors)
            ancestor = pedigree[ancestorID]
            occurrence_matrix[i, j] = ancestor.stats["occurrence"]
            ancestor.stats["occurrence"] = 0
        end
    end
    for (_, individual) in pedigree
        empty!(individual.stats)
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
function occur!(individual::Individual)
    if individual.stats["is_ancestor"] == true # is an ancestor
        individual.stats["occurrence"] += 1
    end
    if !isnothing(individual.father)
        occur!(individual.father)
    end
    if !isnothing(individual.mother)
        occur!(individual.mother)
    end
end