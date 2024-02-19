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
    ancestorIDs = union([ancestor(pedigree, ID) for ID in IDs]...)
    sort!(ancestorIDs)
end

"""
    _findMRCA(pedigree::Pedigree, probandIDs::Vector{Int64})

Return the most recent common ancestors of a list of probands.
"""
function _findMRCA(pedigree::Pedigree, individuals::Vector{Int64})
    ancestorIDs = [ancestor(pedigree, ID) for ID in individuals]
    common_ancestorIDs = ∩(ancestorIDs...)
    older_common_ancestorIDs = ancestor(pedigree, common_ancestorIDs)
    mrcaIDs = sort(collect(setdiff(common_ancestorIDs, older_common_ancestorIDs)))
    mrcaIDs
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
    for (i, ID) in enumerate(individuals), (j, mrcaID) in enumerate(mrcaIDs)
        meioses_matrix[i, j] = _findMinDistance(pedigree, ID, mrcaID)
    end
    genMatrix = GenMatrix(individuals, mrcaIDs, meioses_matrix)
    genMatrix
end

"""
    _findMinDistanceMRCA(pedigree::Pedigree, individuals::Vector{Int64})

Return the minimum distance (meioses) between two individuals.
"""
function _findMinDistanceMRCA(pedigree::Pedigree, individuals::Vector{Int64})
    mrcas = _findMRCA(pedigree, individuals)
    distances = [findDistance(pedigree, individuals, mrca) for mrca in mrcas]
    minimum(distances)
end