"""
    mark_ancestors!(individual::Individual)

Recursively mark the ancestors of an `individual`.
"""
function mark_ancestors!(individual::Individual)
    individual.ancestor = true
    if !isnothing(individual.father)
        mark_ancestors!(individual.father)
    end
    if !isnothing(individual.mother)
        mark_ancestors!(individual.mother)
    end
end

"""
    mark_descendants!(individual::Individual)

Recursively mark the descendants of an `individual`.
"""
function mark_descendants!(individual::Individual)
    individual.descendant = true
    for child in individual.children
        mark_descendants!(child)
    end
end

"""
    branching(pedigree::OrderedDict{Int64, Individual}; pro::Vector{Int64} = pro(pedigree), ancestors::Vector{Int64} = founder(pedigree))

Return a pedigree that filters individuals who are in the paths
between select probands and ancestors.
"""
function branching(pedigree::OrderedDict{Int64, Individual}; pro::Vector{Int64} = pro(pedigree), ancestors::Vector{Int64} = founder(pedigree))
    isolated_pedigree = OrderedDict{Int64, Individual}()
    for ID in pro
        proband = pedigree[ID]
        mark_ancestors!(proband)
    end
    for ID in ancestors
        ancestor = pedigree[ID]
        mark_descendants!(ancestor)
    end
    for (index, (ID, individual)) in enumerate(pedigree)
        if (individual.ancestor && individual.descendant)
            isolated_pedigree[ID] = Individual(ID, nothing, nothing, index, [],
                                                individual.sex, UNEXPLORED, 0., 0,
                                                false, false, 0)
        end
    end
    for (ID, individual) in isolated_pedigree
        father = pedigree[ID].father
        if !isnothing(father)
            if father.ID ∈ keys(isolated_pedigree)
                individual.father = isolated_pedigree[father.ID]
                push!(isolated_pedigree[father.ID].children, individual)
            end
        end
        mother = pedigree[ID].mother
        if !isnothing(mother)
            if mother.ID ∈ keys(isolated_pedigree)
                individual.mother = isolated_pedigree[mother.ID]
                push!(isolated_pedigree[mother.ID].children, individual)
            end
        end
    end
    for (_, individual) in pedigree
        individual.ancestor = false
        individual.descendant = false
    end
    isolated_pedigree
end