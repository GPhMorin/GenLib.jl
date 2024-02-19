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
    contributions = fill(0., length(pedigree))
    states = fill(UNEXPLORED, length(pedigree))
    for ID in pro
        proband = pedigree[ID]
        states[proband.rank] = PROBAND
    end
    for (index₁, ancestorID) in enumerate(ancestors)
        ancestor = pedigree[ancestorID]
        _contribute!(ancestor, contributions, states)
        for (index₂, probandID) in enumerate(pro)
            proband = pedigree[probandID]
            matrix[index₂, index₁] = contributions[proband.rank]
            contributions[proband.rank] = 0.
        end
    end
    matrix
end

"""
    _contribute!(individual::Individual, depth::Int64 = 0)

Recursively compute the genetic contributions of an individual.
"""
function _contribute!(individual::Individual, contributions::Vector{Float64}, states::Vector{State}, depth::Int64 = 0)
    # Ported from GENLIB's ExploreConGenProposant
    if states[individual.rank] == PROBAND
        contributions[individual.rank] += 0.5 ^ depth
    end
    for child in individual.children
        _contribute!(child, contributions, states, depth + 1)
    end
end