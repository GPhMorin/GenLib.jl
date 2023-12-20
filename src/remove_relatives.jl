"""
remove_relatives!(probandIDs::Vector{Int64}, genealogy::Dict{Int64, Individual})::Vector{Int64}

Takes a list of `probandIDs` and, according to a given `genealogy`,
removes IDs of individuals who are first cousins or closer in the genealogy.
"""
function remove_relatives!(probandIDs::Vector{Int64}, genealogy::Dict{Int64, Individual})::Vector{Int64}
    pointer = point(genealogy)
    candidateIDs = copy(probandIDs)
    empty!(probandIDs)
    for probandID in candidateIDs
        proband = pointer[probandID]
        relativeIDs = Vector{Int64}()
        add_relatives!(relativeIDs, proband, Int8(0))
        if isempty(relativeIDs ∩ probandIDs)
            push!(probandIDs, proband.ID)
        end
    end
    probandIDs
end

"""
add_relatives!(relativeIDs::Vector{Int64}, individual::PointerIndividual, depth::Int8)::Nothing

Recursively add IDs of individuals who are first cousins or closer.
"""
function add_relatives!(relativeIDs::Vector{Int64}, individual::PointerIndividual, depth::Int8)::Nothing
    push!(relativeIDs, individual.ID)
    if depth < 4
        if !isnothing(individual.father)
            add_relatives!(relativeIDs, individual.father, depth += Int8(1))
        end
        if !isnothing(individual.mother)
            add_relatives!(relativeIDs, individual.mother, depth += Int8(1))
        end
        for child in individual.children
            add_relatives!(relativeIDs, child, depth += Int8(1))
        end
    end
end