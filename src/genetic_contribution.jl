"""
gc(genealogy::Dict{Int64, Individual}, probandIDs::Vector{Int64} = pro(genealogy), ancestorIDs::Vector{Int64} = founder(genealogy))

Takes a `genealogy` dictionary, computes the genetic contribution of each ancestor to each proband using a vector of `probandIDs` and a vector of `ancestorIDs` and returns a matrix.
"""
function gc(
    genealogy::Dict{Int64, Individual};
    probandIDs::Vector{Int64} = pro(genealogy),
    ancestorIDs::Vector{Int64} = founder(genealogy))
    
    # Ported from GENLIB's Congen
    matrix = zeros(length(probandIDs), length(ancestorIDs))
    reference = refer(genealogy)
    probands = [reference[ID] for ID in probandIDs]
    ancestors = [reference[ID] for ID in ancestorIDs]
    for proband in probands
        proband.probability = 0.
        proband.state = PROBAND
    end
    for (index₁, ancestor) in enumerate(ancestors)
        contribute!(ancestor)
        for (index₂, proband) in enumerate(probands)
            matrix[index₂, index₁] = proband.probability
            proband.probability = 0.
        end
    end
    for proband in probands
        proband.state = UNEXPLORED
    end
    matrix
end

"""
contribute!(individual::ReferenceIndividual, depth::Int8 = Int8(0))

Recursively computes the genetic contribution of an `individual` using a reference at a certain `depth`.
"""
function contribute!(individual::ReferenceIndividual, depth::Int8 = Int8(0))
    # Ported from GENLIB's ExploreConGenProposant
    if individual.state == PROBAND
        individual.probability += 0.5 ^ depth
    end
    for child in individual.children
        contribute!(child, depth + Int8(1))
    end
end