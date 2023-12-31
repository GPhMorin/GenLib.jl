"""
branching(genealogy::Dict{Int64, Individual}, probands::Vector{Int64}, ancestors::Vector{Int64})

Takes a `genealogy` and removes individuals who are not in the paths between select `probands` and `ancestors`.
"""
function branching(genealogy::Dict{Int64, Individual}, probandIDs::Vector{Int64}, ancestorIDs::Vector{Int64})::Dict{Int64, Individual}
    pointer = point(genealogy)
    for ID in probandIDs
    end
end

function mark_ancestors!(individual::PointerIndividual)
    for (ID, individual) in genealogy
        if ID == ancestor
            individual.ancestor = true
        elseif ID ∈ descendant(genealogy, ancestor)
            individual.ancestor = true
        else
            individual.ancestor = false
        end
    end
end

function mark_descendants!(gen::Dict{Int64, Individual}, ancestor::Int64)
    for (ID, individual) in genealogy
        if ID == ancestor
            individual.descendant = true
        elseif ID ∈ descendant(genealogy, ancestor)
            individual.descendant = true
        else
            individual.descendant = false
        end
    end
end