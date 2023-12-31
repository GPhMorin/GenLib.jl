"""
occ(genealogy::Dict{Int64, Individual}, probandIDs::Vector{Int64}, ancestorIDs::Vector{Int64}; type::Symbol = :ind)

Takes a `genealogy`, a vector of `probandIDs`, a vector of `ancestorIDs`
and optionally a `type` (`:ind` or `:total`) and returns a matrix of occurrences.
If the type is `:ind` (default), then the matrix corresponds to *probandIDs* × *ancestorIDs*.
If the type is `:total`, then the matrix corresponds to *ancestorIDs* × 1.
"""
function occ(
    genealogy::Dict{Int64, Individual},
    probandIDs::Vector{Int64},
    ancestorIDs::Vector{Int64};
    type::Symbol = :ind)
    
    occurrence_matrix = Matrix{Int64}(undef, length(probandIDs), length(ancestorIDs))
    reference = refer(genealogy)
    for ancestorID in ancestorIDs
        reference[ancestorID].state = ANCESTOR
    end
    for (i, probandID) in enumerate(probandIDs)
        proband = reference[probandID]
        occur!(proband)
        for (j, ancestorID) in enumerate(ancestorIDs)
            ancestor = reference[ancestorID]
            occurrence_matrix[i, j] = ancestor.occurrence
            ancestor.occurrence = 0
        end
    end
    if type == :ind
        return occurrence_matrix
    elseif type == :total
        return sum(occurrence_matrix, dims=1)'
    end
end

"""
occur!(individual::ReferenceIndividual)

A recursive function that increments the occurrence of a reference `individual`
if they are an ancestor.
"""
function occur!(individual::ReferenceIndividual)
    if individual.state == ANCESTOR
        individual.occurrence += 1
    end
    if !isnothing(individual.father)
        occur!(individual.father)
    end
    if !isnothing(individual.mother)
        occur!(individual.mother)
    end
end