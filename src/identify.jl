"""
    founder(pedigree::Pedigree)

Return a vector of founder IDs in alphabetical order.

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
    pro(pedigree::Pedigree)

Return a vector of proband IDs in alphabetical order.

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
    children(pedigree::Pedigree, ID::Int64)

Return the children of an individual.
"""
function children(pedigree::Pedigree, ID::Int64)
    individual = pedigree[ID]
    sort([child.ID for child ∈ individual.children])
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
    _findMRCA(pedigree::Pedigree, probandIDs::Vector{Int64})

Return the most recent common ancestors of a list of probands.
"""
function _findMRCA(pedigree::Pedigree, individuals::Vector{Int64})
    ancestorIDs = [ancestor(pedigree, ID) for ID ∈ individuals]
    common_ancestorIDs = ∩(ancestorIDs...)
    older_common_ancestorIDs = ancestor(pedigree, common_ancestorIDs)
    mrcaIDs = sort(collect(setdiff(common_ancestorIDs, older_common_ancestorIDs)))
    mrcaIDs
end

"""
    _findMinDistanceMRCA(pedigree::Pedigree, individuals::Vector{Int64})

Return the minimum distance (meioses) between two individuals.
"""
function _findMinDistanceMRCA(pedigree::Pedigree, individuals::Vector{Int64})
    mrcas = _findMRCA(pedigree, individuals)
    distances = [findDistance(pedigree, individuals, mrca) for mrca ∈ mrcas]
    minimum(distances)
end

"""
    struct GenMatrix
        individuals::Vector{Int64}
        ancestors::Vector{Int64}
        meioses::Matrix{Int64}
    end

A matrix that goes with individuals as rows and ancestors as columns.
"""
struct GenMatrix
    individuals::Vector{Int64}
    ancestors::Vector{Int64}
    meioses::Matrix{Int64}
end

"""
    findMRCA(pedigree::Pedigree, IDs::Vector{Int64})

Return a [`GenMatrix`](@ref) of meioses between individuals and their
most recent common ancestors (MRCAs).

# Example

```julia
import GenLib as gen
genea140 = gen.genea140
ped = gen.genealogy(genea140)
pro = gen.pro(ped)
pro1 = pro[1]
pro2 = pro[2]
genMatrix = gen.findMRCA(ped, [pro1, pro2])
```
"""
function findMRCA(pedigree::Pedigree, individuals::Vector{Int64})
    mrcaIDs = _findMRCA(pedigree, individuals)
    meioses_matrix = Matrix{Int64}(undef, length(individuals), length(mrcaIDs))
    for (i, ID) ∈ enumerate(individuals), (j, mrcaID) ∈ enumerate(mrcaIDs)
        meioses_matrix[i, j] = _findMinDistance(pedigree, ID, mrcaID)
    end
    genMatrix = GenMatrix(individuals, mrcaIDs, meioses_matrix)
    genMatrix
end

"""
    ancestor(pedigree::Pedigree, ID::Int64)

Return a vector of an individual's ancestors.
"""
function ancestor(pedigree::Pedigree, ID::Int64)
    ancestorIDs = Set{Int64}()
    stack = Int64[]
    push!(stack, ID)
    while !isempty(stack)
        ID = pop!(stack)
        individual = pedigree[ID]
        if !isnothing(individual.father)
            push!(stack, individual.father.ID)
            push!(ancestorIDs, individual.father.ID)
        end
        if !isnothing(individual.mother)
            push!(stack, individual.mother.ID)
            push!(ancestorIDs, individual.mother.ID)
        end
    end
    ancestorIDs = collect(ancestorIDs)
    sort!(ancestorIDs)
end

"""
    ancestor(pedigree::Pedigree, IDs::Vector{Int64})

Return a vector of several individual's ancestors.
"""
function ancestor(pedigree::Pedigree, IDs::Vector{Int64})
    ancestorIDs = union([ancestor(pedigree, ID) for ID ∈ IDs]...)
    sort!(ancestorIDs)
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