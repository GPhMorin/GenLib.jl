"""
    _mark_ancestors!(individual::Individual)

Recursively mark the ancestors of an `individual`.
"""
function _mark_ancestors!(individual::Individual)
    individual.stats[1] = true
    if !isnothing(individual.father)
        _mark_ancestors!(individual.father)
    end
    if !isnothing(individual.mother)
        _mark_ancestors!(individual.mother)
    end
end

"""
    _mark_descendants!(individual::Individual)

Recursively mark the descendants of an `individual`.
"""
function _mark_descendants!(individual::Individual)
    individual.stats[2] = true
    for child in individual.children
        _mark_descendants!(child)
    end
end

"""
    branching(pedigree::Pedigree; pro::Vector{Int64} = pro(pedigree), ancestors::Vector{Int64} = founder(pedigree))

Return a pedigree that filters individuals who are in the paths
between select probands and ancestors.
"""
function branching(pedigree::Pedigree; pro::Vector{Int64} = pro(pedigree), ancestors::Vector{Int64} = founder(pedigree))
    isolated_pedigree = Pedigree()
    for (_, individual) in pedigree
        individual.stats = [false, false]
    end
    for ID in pro
        proband = pedigree[ID]
        _mark_ancestors!(proband)
    end
    for ID in ancestors
        ancestor = pedigree[ID]
        _mark_descendants!(ancestor)
    end
    index = 0
    for (ID, individual) in pedigree
        if (individual.stats[1] && individual.stats[2])
            index += 1
            isolated_pedigree[ID] = Individual(ID, nothing, nothing, Int64[], individual.sex, index, [])
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
        empty!(individual.stats)
    end
    isolated_pedigree
end