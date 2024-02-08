"""
remove_relatives!(probandIDs::Vector{Int64}, genealogy::OrderedDict{Int64, Individual})

Takes a list of `probandIDs` and, according to a given `genealogy`,
removes IDs of individuals who are first cousins or closer in the genealogy.
"""
function remove_relatives!(probandIDs::Vector{Int64}, genealogy::OrderedDict{Int64, Individual})
    reference = refer(genealogy)
    candidateIDs = copy(probandIDs)
    empty!(probandIDs)
    for probandID in candidateIDs
        proband = reference[probandID]
        relativeIDs = Vector{Int64}()
        add_relatives!(relativeIDs, proband, Int8(0))
        if isempty(relativeIDs ∩ probandIDs)
            push!(probandIDs, proband.ID)
        end
    end
    probandIDs
end

"""
add_relatives!(relativeIDs::Vector{Int64}, individual::ReferenceIndividual, depth::Int64)

Recursively add IDs of individuals who are first cousins or closer.
"""
function add_relatives!(relativeIDs::Vector{Int64}, individual::ReferenceIndividual, depth::Int64)
    push!(relativeIDs, individual.ID)
    if depth < 4
        if !isnothing(individual.father)
            add_relatives!(relativeIDs, individual.father, depth += 1)
        end
        if !isnothing(individual.mother)
            add_relatives!(relativeIDs, individual.mother, depth += 1)
        end
        for child in individual.children
            add_relatives!(relativeIDs, child, depth += 1)
        end
    end
end