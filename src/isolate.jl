"""
mark_ancestors!(individual::ReferenceIndividual)

A recursive function that marks the ancestors of an `individual`.
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

A recursive function that marks the descendants of an `individual`.
"""
function mark_descendants!(individual::ReferenceIndividual)
    individual.descendant = true
    for child in individual.children
        mark_descendants!(child)
    end
end

"""
branching(genealogy::Dict{Int64, Individual}, probandIDs::Vector{Int64}, ancestorIDs::Vector{Int64})

Takes a `genealogy` and removes individuals who are not in the paths between select probands and ancestors.
"""
function branching(genealogy::Dict{Int64, Individual}; probandIDs::Vector{Int64} = pro(genealogy), ancestorIDs::Vector{Int64} = founder(genealogy))
    isolated_genealogy = Dict{Int64, Individual}()
    ref = refer(genealogy)
    for ID in probandIDs
        proband = ref[ID]
        mark_ancestors!(proband)
    end
    for ID in ancestorIDs
        ancestor = ref[ID]
        mark_descendants!(ancestor)
    end
    for (ID, individual) in genealogy
        if ref[ID].ancestor & ref[ID].descendant
            father = 0
            if individual.father != 0
                father = ref[individual.father].ancestor & ref[individual.father].descendant ? individual.father : 0
            end
            mother = 0
            if individual.mother != 0
                mother = ref[individual.mother].ancestor & ref[individual.mother].descendant ? individual.mother : 0
            end
            children = filter(x -> ref[x].ancestor & ref[x].descendant, individual.children)
            isolated_genealogy[ID] = Individual(father, mother, individual.index, children, individual.sex)
        end
    end
    isolated_genealogy
end

"""
branching(genealogy, IDs)

Takes a `genealogy` and isolates individuals with one of the given `IDs`.
"""
function branching(genealogy, IDs)
    isolated_genealogy = Dict()
    for (ID, individual) in genealogy
        if ID ∈ IDs
            father = 0
            if individual.father != 0
                father = individual.father ∈ IDs ? individual.father : 0
            end
            mother = 0
            if individual.mother != 0
                mother = individual.mother ∈ IDs ? individual.mother : 0
            end
            children = filter(x -> x ∈ IDs, individual.children)
            isolated_genealogy[ID] = Individual(father, mother, individual.index, children, individual.sex)
        end
    end
    isolated_genealogy
end

function ablate(genealogy, IDs)
    ablated_genealogy = Dict{Int64, Individual}()
    for (ID, individual) in genealogy
        if ID ∉ IDs
            father = 0
            if individual.father != 0
                father = individual.father ∈ IDs ? individual.father : 0
            end
            mother = 0
            if individual.mother != 0
                mother = individual.mother ∈ IDs ? individual.mother : 0
            end
            children = filter(x -> x ∈ IDs, individual.children)
            ablated_genealogy[ID] = Individual(father, mother, individual.index, children, individual.sex)
        end
    end
    ablated_genealogy
end