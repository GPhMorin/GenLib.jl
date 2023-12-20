function occ(genealogy::Dict{Int64, Individual}, probandIDs::Vector{Int64}, ancestorIDs::Vector{Int64}, type::Symbol = :ind)::Matrix{Int64}
    occurence_matrix = Matrix{Int64}(undef, length(probandIDs), length(ancestorIDs))
    pointer = point(genealogy)
    for ancestorID in ancestorIDs
        pointer[ancestorID].state = ANCESTOR
    end
    for (i, probandID) in enumerate(probandIDs)
        proband = pointer[probandID]
        occur!(proband)
        for (j, ancestorID) in enumerate(ancestorIDs)
            ancestor = pointer[ancestorID]
            occurence_matrix[i, j] = ancestor.occurrence
            ancestor.occurrence = 0
        end
    end
    if type == :ind
        return occurence_matrix
    elseif type == :total
        return sum(occurrence_matrix, dims=2)
    end
end

function occur!(individual::PointerIndividual)::Nothing
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