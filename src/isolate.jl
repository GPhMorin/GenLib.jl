"""
mark_ancestors!(individual::ReferenceIndividual)

A recursive function that marks the ancestors of a `individual`.
"""
function mark_ancestors!(individual::ReferenceIndividual)
    individual.ancestor = true
    if !isnothing(individual.father)
        mark_ancestors!(individual.father)
    end
    if !isnothing(individual.mother)
        mark_ancestors!(individual.mother)
    end
end

"""
mark_descendants!(individual::ReferenceIndividual)

A recursive function that marks the descendants of a `individual`.
"""
function mark_descendants!(individual::ReferenceIndividual)
    individual.descendant = true
    for child in individual.children
        mark_descendants!(child)
    end
end

"""
branching(genealogy::Dict{Int64, Individual}, probandIDs::Vector{Int64}, ancestorIDs::Vector{Int64})

Takes a `genealogy` and removes individuals who are not in the paths between select `probands` and `ancestors`.
"""
function branching(genealogy::Dict{Int64, Individual}; probandIDs::Vector{Int64} = pro(genealogy), ancestorIDs::Vector{Int64} = founder(genealogy))
    isolated_genealogy = Dict{Int64, Individual}()
    reference = point(genealogy)
    for ID in probandIDs
        proband = reference[ID]
        mark_ancestors!(proband)
    end
    for ID in ancestorIDs
        ancestor = reference[ID]
        mark_descendants!(ancestor)
    end
    for (ID, individual) in genealogy
        if reference[ID].ancestor | reference[ID].descendant
            new_genealogy[ID] = individual
        end
    end
    isolated_genealogy
end