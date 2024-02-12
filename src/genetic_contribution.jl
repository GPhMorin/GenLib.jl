"""
    gc(pedigree::OrderedDict{Int64, Individual}; pro::Vector{Int64} = pro(pedigree), ancestors::Vector{Int64} = founder(pedigree))

Return a matrix of the genetic contribution
of each ancestor (columns) to each proband (rows).

# Example

```@repl
import GenLib as gen
genea140 = gen.genea140;
ped = gen.genealogy(genea140);
contributions = gen.gc(ped)
sum(contributions, dims=2)
```
"""
function gc(
    pedigree::OrderedDict{Int64, Individual};
    probandIDs::Vector{Int64} = pro(pedigree),
    ancestorIDs::Vector{Int64} = founder(pedigree))
    
    # Ported from GENLIB's Congen
    matrix = zeros(length(probandIDs), length(ancestorIDs))
    probands = [pedigree[ID] for ID in probandIDs]
    ancestors = [pedigree[ID] for ID in ancestorIDs]
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
    for (_, individual) in pedigree
        individual.probability = 0.
        individual.state = UNEXPLORED
    end
    matrix
end

"""
    contribute!(individual::Individual, depth::Int64 = 0)

Recursively compute the genetic contributions of an individual.
"""
function contribute!(individual::Individual, depth::Int64 = 0)
    # Ported from GENLIB's ExploreConGenProposant
    if individual.state == PROBAND
        individual.probability += 0.5 ^ depth
    end
    for child in individual.children
        contribute!(child, depth + 1)
    end
end