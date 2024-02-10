"""
    mark_ancestors!(individual::ReferenceIndividual)

Recursively mark the ancestors of an `individual`.
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

Recursively mark the descendants of an `individual`.
"""
function mark_descendants!(individual::ReferenceIndividual)
    individual.descendant = true
    for child in individual.children
        mark_descendants!(child)
    end
end

"""
    branching(genealogy::OrderedDict{Int64, Individual}; pro::Vector{Int64} = pro(genealogy), ancestors::Vector{Int64} = founder(genealogy))

Return a pedigree that filters individuals who are in the paths
between select probands and ancestors.
"""
function branching(genealogy::OrderedDict{Int64, Individual}; pro::Vector{Int64} = pro(genealogy), ancestors::Vector{Int64} = founder(genealogy))
    isolated_genealogy = OrderedDict{Int64, Individual}()
    ref = refer(genealogy)
    for ID in pro
        proband = ref[ID]
        mark_ancestors!(proband)
    end
    for ID in ancestors
        ancestor = ref[ID]
        mark_descendants!(ancestor)
    end
    index = 0
    for (ID, individual) in genealogy
        if (ref[ID].ancestor & ref[ID].descendant)
            index += 1
            father = 0
            if individual.father != 0
                father = ref[individual.father].ancestor && ref[individual.father].descendant ? individual.father : 0 
            end
            mother = 0
            if individual.mother != 0
                mother = ref[individual.mother].ancestor && ref[individual.mother].descendant ? individual.mother : 0 
            end
            children = filter(x -> ref[x].ancestor && ref[x].descendant, individual.children)
            isolated_genealogy[ID] = Individual(father, mother, index, children, individual.sex)
        end
    end
    isolated_genealogy
end