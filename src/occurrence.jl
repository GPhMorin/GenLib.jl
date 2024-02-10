"""
    occ(genealogy::OrderedDict{Int64, Individual}; pro::Vector{Int64} = pro(genealogy), ancestors::Vector{Int64} = founder(genealogy), typeOcc::String = "IND")

Return a matrix of ancestors' occurrences.

If `typeOcc` is `:ind` (default), then the matrix corresponds to the occurrence per individual.
If `typeOcc` is `:total`, then the matrix corresponds to the total occurrence.
"""
function occ(
    genealogy::OrderedDict{Int64, Individual};
    pro::Vector{Int64} = pro(genealogy),
    ancestors::Vector{Int64} = founder(genealogy),
    typeOcc::String = "IND")
    
    occurrence_matrix = Matrix{Int64}(undef, length(pro), length(ancestors))
    reference = refer(genealogy)
    for ID in ancestors
        reference[ID].ancestor = true
    end
    for (i, probandID) in enumerate(pro)
        proband = reference[probandID]
        occur!(proband)
        for (j, ancestorID) in enumerate(ancestors)
            ancestor = reference[ancestorID]
            occurrence_matrix[i, j] = ancestor.occurrence
            ancestor.occurrence = 0
        end
    end
    if typeOcc == "IND"
        return occurrence_matrix
    elseif typeOcc == "TOTAL"
        return sum(occurrence_matrix, dims=1)'
    end
end

"""
    occur!(individual::ReferenceIndividual)

Recursively increment the occurrence of an `individual` if they are an ancestor.
"""
function occur!(individual::ReferenceIndividual)
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