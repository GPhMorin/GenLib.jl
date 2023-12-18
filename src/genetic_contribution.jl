"""
gc(genealogy::Dict{Int64, Individual}, probandIDs::Vector{Int64} = pro(genealogy), ancestorIDs::Vector{Int64} = founder(genealogy))::Matrix{Float64}

Takes a `genealogy` dictionary, computes the genetic contribution of each ancestor to each proband using a vector of `probandIDs` and a vector of `ancestorIDs` and returns a matrix.
"""
function gc(genealogy::Dict{Int64, Individual}, probandIDs::Vector{Int64} = pro(genealogy), ancestorIDs::Vector{Int64} = founder(genealogy))::Matrix{Float64}
    # Ported from GENLIB's Congen
    matrix = zeros(length(probandIDs), length(ancestorIDs))
    pointer = point(genealogy)
    probands = [pointer[ID] for ID in probandIDs]
    ancestors = [pointer[ID] for ID in ancestorIDs]
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
contribute!(individual::PointerIndividual, depth::Int8 = Int8(0))::Nothing

Recursively computes the genetic contribution of an `individual` using a pointer at a certain `depth`.
"""
function contribute!(individual::PointerIndividual, depth::Int8 = Int8(0))::Nothing
    # Ported from GENLIB's ExploreConGenProposant
    if individual.state == PROBAND
        individual.probability += 0.5 ^ depth
    end
    for child in individual.children
        contribute!(child, depth + Int8(1))
    end
end