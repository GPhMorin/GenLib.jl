"""
    gc(genealogy::OrderedDict{Int64, Individual}; pro::Vector{Int64} = pro(genealogy), ancestors::Vector{Int64} = founder(genealogy))

Takes a `genealogy` dictionary, computes the genetic contribution of each ancestor to each proband using a vector of `probandIDs` and a vector of `ancestorIDs` and returns a matrix.
"""
function gc(
    genealogy::OrderedDict{Int64, Individual};
    probands::Vector{Int64} = pro(genealogy),
    ancestors::Vector{Int64} = founder(genealogy))
    
    # Ported from GENLIB's Congen
    matrix = zeros(length(probands), length(ancestors))
    reference = refer(genealogy)
    ref_probands = [reference[ID] for ID in probands]
    ref_ancestors = [reference[ID] for ID in ancestors]
    for ref_proband in ref_probands
        ref_proband.probability = 0.
        ref_proband.state = PROBAND
    end
    for (index₁, ref_ancestor) in enumerate(ref_ancestors)
        contribute!(ref_ancestor)
        for (index₂, ref_proband) in enumerate(ref_probands)
            matrix[index₂, index₁] = ref_proband.probability
            ref_proband.probability = 0.
        end
    end
    for ref_proband in ref_probands
        ref_proband.state = UNEXPLORED
    end
    matrix
end

"""
    contribute!(individual::ReferenceIndividual, depth::Int64 = 0)

Recursively computes the genetic contribution of an `individual` using a reference at a certain `depth`.
"""
function contribute!(individual::ReferenceIndividual, depth::Int64 = 0)
    # Ported from GENLIB's ExploreConGenProposant
    if individual.state == PROBAND
        individual.probability += 0.5 ^ depth
    end
    for child in individual.children
        contribute!(child, depth + 1)
    end
end