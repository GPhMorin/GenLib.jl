"""
    pro(pedigree::OrderedDict{Int64, Individual})

Return a vector of proband IDs in alphabetical order.

# Example

```julia
import GenLib as gen
genea140 = gen.genea140
ped = gen.genealogy(genea140)
probands = gen.pro(ped)
```
"""
function pro(pedigree::OrderedDict{Int64, Individual})
    probands = [individual.ID for (_, individual) in pedigree if isempty(individual.children)]
    sort(probands)
end

"""
    founder(pedigree::OrderedDict{Int64, Individual})

Return a vector of founder IDs in alphabetical order.

# Example

```julia
import GenLib as gen
genea140 = gen.genea140
ped = gen.genealogy(genea140)
founders = gen.founder(ped)
```
"""
function founder(pedigree::OrderedDict{Int64, Individual})
    founders = [individual.ID for (_, individual) in pedigree if isnothing(individual.father) && isnothing(individual.mother)]
    sort(founders)
end

"""
    findFounders(pedigree::OrderedDict{Int64}, IDs::Vector{Int64})

Return a vector of founders from whom the `IDs` descend.
"""
function findFounders(pedigree::OrderedDict{Int64}, IDs::Vector{Int64})
    ancestorIDs = [ancestor(pedigree, ID) for ID in IDs]
    common_ancestorIDs = ∩(ancestorIDs...)
    founderIDs = [ancestorID for ancestorID in common_ancestorIDs
                  if isnothing(pedigree[ancestorID].father) && isnothing([ancestorID].mother)]
    sort(founderIDs)
end

"""
    get_paths(pedigree::OrderedDict{Int64, Individual}, individual::Individual)

Return the paths from an individual to their ancestors.
"""
function get_paths(pedigree::OrderedDict{Int64, Individual}, individual::Individual)
    paths = Vector{Vector{Int64}}([[individual.ID]])
    if !isnothing(individual.father)
        fathers_paths = get_paths(pedigree, individual.father)
        for path in fathers_paths
            push!(path, individual.ID)
        end
        append!(paths, fathers_paths)
    end
    if !isnothing(individual.mother)
        mothers_paths = get_paths(pedigree, individual.mother)
        for path in mothers_paths
            push!(path, individual.ID)
        end
        append!(paths, mothers_paths)
    end
    paths
end

"""
    children(pedigree::OrderedDict{Int64, Individual}, ID::Int64)

Return the children of an individual.
"""
function children(pedigree::OrderedDict{Int64, Individual}, ID::Int64)
    individual = pedigree[ID]
    sort([child.ID for child in individual.children])
end

"""
    descendant(pedigree::OrderedDict{Int64, Individual}, ID::Int64)

Return the descendants of an individual.
"""
function descendant(pedigree::OrderedDict{Int64, Individual}, ID::Int64)
    descendantIDs = Set{Int64}()
    stack = Stack{Int64}()
    push!(stack, ID)
    while !isempty(stack)
        ID = pop!(stack)
        individual = pedigree[ID]
        for child in individual.children
            push!(stack, child.ID)
            push!(descendantIDs, child.ID)
        end
    end
    sort(collect(descendantIDs))
end

"""
    rec(pedigree::OrderedDict{Int64, Individual}, probandIDs::Vector{Int64} = pro(genealogy), ancestorIDs::Vector{Int64} = founder(genealogy))

Return the number of descendants of each ancestor.
"""
function rec(
    pedigree::OrderedDict{Int64, Individual},
    probandIDs::Vector{Int64} = pro(pedigree),
    ancestorIDs::Vector{Int64} = founder(pedigree))
    
    coverage = Vector{Int64}()
    for ancestorID in ancestorIDs
        descendantIDs = descendant(pedigree, ancestorID)
        descendantIDs = filter!(x -> x in probandIDs, descendantIDs)
        append!(coverage, length(descendantIDs))
    end
    coverage
end

"""
    genout(pedigree::OrderedDict{Int64, Individual}, sorted::Bool = false)

Return a pedigree as a `DataFrame`.

If `sorted` is `false` (the default), then the individuals
will appear in the same order as in the pedigree.

If `sorted` is `true`, then the individuals
will appear in alphabetical ID order.

```julia
import GenLib as gen
geneaJi = gen.geneaJi
ped = gen.genealogy(geneaJi)
gen.genout(ped)
```
"""
function genout(pedigree::OrderedDict{Int64, Individual}; sorted::Bool = false)
    inds = Int64[]
    fathers = Int64[]
    mothers = Int64[]
    sexes = Int64[]
    if !sorted
        for (ID, individual) in pedigree
            push!(inds, ID)
            push!(fathers, !isnothing(individual.father) ? individual.father.ID : 0)
            push!(mothers, !isnothing(individual.mother) ? individual.mother.ID : 0)
            push!(sexes, individual.sex)
        end
    else # if sorted
        inds = sort(collect(keys(pedigree)))
        for ID in inds
            individual = pedigree[ID]
            push!(fathers, !isnothing(individual.father) ? individual.father.ID : 0)
            push!(mothers, !isnothing(individual.mother) ? individual.mother.ID : 0)
            push!(sexes, individual.sex)
        end
    end
    df = DataFrame([inds, fathers, mothers, sexes], ["ind", "father", "mother", "sex"])
    df
end