"""
    _mark_ancestors!(individual::Individual)

Recursively mark the ancestors of an `individual`.
"""
function _mark_ancestors!(individual::Individual, is_ancestors::Vector{Bool})
    is_ancestors[individual.index] = true
    if !isnothing(individual.father)
        _mark_ancestors!(individual.father, is_ancestors)
    end
    if !isnothing(individual.mother)
        _mark_ancestors!(individual.mother, is_ancestors)
    end
end

"""
    _mark_descendants!(individual::Individual)

Recursively mark the descendants of an `individual`.
"""
function _mark_descendants!(individual::Individual, is_descendants::Vector{Bool})
    is_descendants[individual.index] = true
    for child in individual.children
        _mark_descendants!(child, is_descendants)
    end
end

"""
    branching(pedigree::Pedigree; pro::Vector{Int64} = pro(pedigree), ancestors::Vector{Int64} = founder(pedigree))

Return a pedigree that filters individuals who are in the paths
between select probands and ancestors.
"""
function branching(pedigree::Pedigree; pro::Vector{Int64} = pro(pedigree), ancestors::Vector{Int64} = founder(pedigree))
    isolated_pedigree = Pedigree()
    is_ancestors = fill(false, length(pedigree))
    is_descendants = fill(false, length(pedigree))
    for ID in pro
        proband = pedigree[ID]
        _mark_ancestors!(proband, is_ancestors)
    end
    for ID in ancestors
        ancestor = pedigree[ID]
        _mark_descendants!(ancestor, is_descendants)
    end
    index = 0
    for ((ID, individual), is_ancestor, is_descendant) in zip(pedigree, is_ancestors, is_descendants)
        if is_ancestor && is_descendant
            index += 1
            isolated_pedigree[ID] = Individual(ID, nothing, nothing, Int64[], individual.sex, index, Dict{String, Any}())
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
    isolated_pedigree
end