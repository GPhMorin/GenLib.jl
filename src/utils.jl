"""
    pro(pedigree::Pedigree)

Return a vector of proband IDs ∈ alphabetical order.

# Example

```julia
import GenLib as gen
genea140 = gen.genea140
ped = gen.genealogy(genea140)
probands = gen.pro(ped)
```
"""
function pro(pedigree::Pedigree)
    probands = [individual.ID for individual ∈ values(pedigree) if isempty(individual.children)]
    sort(probands)
end

"""
    founder(pedigree::Pedigree)

Return a vector of founder IDs ∈ alphabetical order.

# Example

```julia
import GenLib as gen
genea140 = gen.genea140
ped = gen.genealogy(genea140)
founders = gen.founder(ped)
```
"""
function founder(pedigree::Pedigree)
    founders = [individual.ID for individual ∈ values(pedigree) if isnothing(individual.father) && isnothing(individual.mother)]
    sort(founders)
end

"""
    findFounders(pedigree::Pedigree, IDs::Vector{Int64})

Return a vector of founders from whom the `IDs` descend.
"""
function findFounders(pedigree::Pedigree, IDs::Vector{Int64})
    ancestorIDs = [ancestor(pedigree, ID) for ID ∈ IDs]
    common_ancestorIDs = ∩(ancestorIDs...)
    founderIDs = [ancestorID for ancestorID ∈ common_ancestorIDs
                  if isnothing(pedigree[ancestorID].father) && isnothing(pedigree[ancestorID].mother)]
    sort(founderIDs)
end

"""
    get_paths(pedigree::Pedigree{T}, individual::T) where T <: AbstractIndividual

Return the paths from an individual to their ancestors.
"""
function get_paths(pedigree::Pedigree{T}, individual::T) where T <: AbstractIndividual
    paths = Vector{Vector{Int64}}([[individual.ID]])
    if !isnothing(individual.father)
        fathers_paths = get_paths(pedigree, individual.father)
        for path ∈ fathers_paths
            push!(path, individual.ID)
        end
        append!(paths, fathers_paths)
    end
    if !isnothing(individual.mother)
        mothers_paths = get_paths(pedigree, individual.mother)
        for path ∈ mothers_paths
            push!(path, individual.ID)
        end
        append!(paths, mothers_paths)
    end
    paths
end

"""
    children(pedigree::Pedigree, ID::Int64)

Return the children of an individual.
"""
function children(pedigree::Pedigree, ID::Int64)
    individual = pedigree[ID]
    sort([child.ID for child ∈ individual.children])
end

"""
    descendant(pedigree::Pedigree, ID::Int64)

Return the descendants of an individual.
"""
function descendant(pedigree::Pedigree, ID::Int64)
    descendantIDs = Set{Int64}()
    stack = Int64[]
    push!(stack, ID)
    while !isempty(stack)
        ID = pop!(stack)
        individual = pedigree[ID]
        for child ∈ individual.children
            push!(stack, child.ID)
            push!(descendantIDs, child.ID)
        end
    end
    sort(collect(descendantIDs))
end

"""
    rec(pedigree::Pedigree, probandIDs::Vector{Int64} = pro(genealogy), ancestorIDs::Vector{Int64} = founder(genealogy))

Return the number of descendants of each ancestor.
"""
function rec(
    pedigree::Pedigree,
    probandIDs::Vector{Int64} = pro(pedigree),
    ancestorIDs::Vector{Int64} = founder(pedigree))
    
    coverage = Vector{Int64}()
    for ancestorID ∈ ancestorIDs
        descendantIDs = descendant(pedigree, ancestorID)
        descendantIDs = filter!(x -> x ∈ probandIDs, descendantIDs)
        append!(coverage, length(descendantIDs))
    end
    coverage
end

"""
    genout(pedigree::Pedigree, sorted::Bool = false)

Return a pedigree as a `DataFrame`.

If `sorted` is `false` (the default), then the individuals
will appear ∈ the same order as ∈ the pedigree.

If `sorted` is `true`, then the individuals
will appear ∈ alphabetical ID order.

```julia
import GenLib as gen
geneaJi = gen.geneaJi
ped = gen.genealogy(geneaJi)
gen.genout(ped)
```
"""
function genout(pedigree::Pedigree; sorted::Bool = false)
    inds = Int64[]
    fathers = Int64[]
    mothers = Int64[]
    sexes = Int64[]
    if !sorted
        for (ID, individual) ∈ pedigree
            push!(inds, ID)
            push!(fathers, !isnothing(individual.father) ? individual.father.ID : 0)
            push!(mothers, !isnothing(individual.mother) ? individual.mother.ID : 0)
            push!(sexes, individual.sex)
        end
    else # if sorted
        inds = sort(collect(keys(pedigree)))
        for ID ∈ inds
            individual = pedigree[ID]
            push!(fathers, !isnothing(individual.father) ? individual.father.ID : 0)
            push!(mothers, !isnothing(individual.mother) ? individual.mother.ID : 0)
            push!(sexes, individual.sex)
        end
    end
    df = DataFrame([inds, fathers, mothers, sexes], ["ind", "father", "mother", "sex"])
    df
end

"""
    _max_depth(individual::T) where T <: AbstractIndividual

Return the maximum depth of an individual's pedigree.
"""
function _max_depth(individual::T) where T <: AbstractIndividual
    father_depth = 1
    mother_depth = 1
    if !isnothing(individual.father)
        father_depth += _max_depth(individual.father)
    end
    if !isnothing(individual.mother)
        mother_depth += _max_depth(individual.mother)
    end
    max(father_depth, mother_depth)
end