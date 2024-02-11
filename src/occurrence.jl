"""
    occ(pedigree::OrderedDict{Int64, Individual}; pro::Vector{Int64} = pro(genealogy), ancestors::Vector{Int64} = founder(genealogy), typeOcc::String = "IND")

Return a matrix of ancestors' occurrences.

If `typeOcc` is `:ind` (default), then the matrix corresponds to the occurrence per individual.
If `typeOcc` is `:total`, then the matrix corresponds to the total occurrence.
"""
function occ(
    pedigree::OrderedDict{Int64, Individual};
    pro::Vector{Int64} = pro(pedigree),
    ancestors::Vector{Int64} = founder(pedigree),
    typeOcc::String = "IND")
    
    occurrence_matrix = Matrix{Int64}(undef, length(pro), length(ancestors))
    for ID in ancestors
        pedigree[ID].ancestor = true
    end
    for (i, probandID) in enumerate(pro)
        proband = pedigree[probandID]
        occur!(proband)
        for (j, ancestorID) in enumerate(ancestors)
            ancestor = pedigree[ancestorID]
            occurrence_matrix[i, j] = ancestor.occurrence
            ancestor.occurrence = 0
        end
    end
    for (_, individual) in pedigree
        individual.ancestor = false
        individual.occurrence = 0
    end
    if typeOcc == "IND"
        return occurrence_matrix
    elseif typeOcc == "TOTAL"
        return sum(occurrence_matrix, dims=1)'
    end
end

"""
    occur!(individual::Individual)

Recursively increment the occurrence of an `individual` if they are an ancestor.
"""
function occur!(individual::Individual)
    if individual.ancestor
        individual.occurrence += 1
    end
    if !isnothing(individual.father)
        occur!(individual.father)
    end
    if !isnothing(individual.mother)
        occur!(individual.mother)
    end
end