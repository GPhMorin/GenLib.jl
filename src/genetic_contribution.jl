"""
    gc(pedigree::Pedigree; pro::Vector{Int64} = pro(pedigree), ancestors::Vector{Int64} = founder(pedigree))

Return a matrix of the genetic contribution
of each ancestor (columns) to each proband (rows).

# Example

```julia
import GenLib as gen
genea140 = gen.genea140
ped = gen.genealogy(genea140)
contributions = gen.gc(ped)
```
"""
function gc(
    pedigree::Pedigree;
    pro::Vector{Int64} = pro(pedigree),
    ancestors::Vector{Int64} = founder(pedigree))
    
    # Ported from GENLIB's Congen
    matrix = zeros(length(pro), length(ancestors))
    for (ID, individual) in pedigree
        individual.stats = [0., UNEXPLORED]
        if ID in pro
            individual.stats[2] = PROBAND
        end
    end
    for (index₁, ancestorID) in enumerate(ancestors)
        ancestor = pedigree[ancestorID]
        _contribute!(ancestor)
        for (index₂, probandID) in enumerate(pro)
            proband = pedigree[probandID]
            matrix[index₂, index₁] = proband.stats[1]
            proband.stats[1] = 0.
        end
    end
    for (_, individual) in pedigree
        individual.stats = nothing
    end
    matrix
end

"""
    _contribute!(individual::Individual, depth::Int64 = 0)

Recursively compute the genetic contributions of an individual.
"""
function _contribute!(individual::Individual, depth::Int64 = 0)
    # Ported from GENLIB's ExploreConGenProposant
    if individual.stats[2] == PROBAND
        individual.stats[1] += 0.5 ^ depth
    end
    for child in individual.children
        _contribute!(child, depth + 1)
    end
end