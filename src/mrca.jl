"""
    ancestor(pedigree::OrderedDict{Int64, Individual}, ID::Int64)

Return a vector of an individual's ancestors.
"""
function ancestor(pedigree::OrderedDict{Int64, Individual}, ID::Int64)
    ancestorIDs = Set{Int64}()
    stack = Stack{Int64}()
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
    ancestor(pedigree::OrderedDict{Int64, Individual}, IDs::Vector{Int64})

Return a vector of several individual's ancestors.
"""
function ancestor(pedigree::OrderedDict{Int64, Individual}, IDs::Vector{Int64})
    ancestorIDs = union([ancestor(pedigree, ID) for ID in IDs]...)
    sort!(ancestorIDs)
end

"""
    findMRCA(pedigree::OrderedDict{Int64, Individual}, IDs::Vector{Int64})

Return a vector of individuals' most recent common ancestors (MRCAs).
"""
function findMRCA(pedigree::OrderedDict{Int64, Individual}, IDs::Vector{Int64})
    ancestorIDs = [ancestor(pedigree, ID) for ID in IDs]
    common_ancestorIDs = ∩(ancestorIDs...)
    older_common_ancestorIDs = ancestor(pedigree, common_ancestorIDs)
    most_recent_commmon_ancestorIDs = setdiff(common_ancestorIDs, older_common_ancestorIDs)
    most_recent_commmon_ancestorIDs
end