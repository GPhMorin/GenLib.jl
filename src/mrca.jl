"""
    ancestor(genealogy::OrderedDict{Int64, Individual}, ID::Int64)

Return a vector of an individual's ancestors.
"""
function ancestor(genealogy::OrderedDict{Int64, Individual}, ID::Int64)
    ancestors = Set{Int64}()
    stack = Stack{Int64}()
    push!(stack, ID)
    while !isempty(stack)
        ID = pop!(stack)
        individual = genealogy[ID]
        if individual.father != 0
            push!(stack, individual.father)
            push!(ancestors, individual.father)
        end
        if individual.mother != 0
            push!(stack, individual.mother)
            push!(ancestors, individual.mother)
        end
    end
    ancestors = collect(ancestors)
    sort!(ancestors)
end

"""
    ancestor(genealogy::OrderedDict{Int64, Individual}, IDs::Vector{Int64})

Return a vector of several individual's ancestors.
"""
function ancestor(genealogy::OrderedDict{Int64, Individual}, IDs::Vector{Int64})
    ancestors = union([ancestor(genealogy, ID) for ID in IDs]...)
    sort!(ancestors)
end

"""
    findMRCA(genealogy::OrderedDict{Int64, Individual}, IDs::Vector{Int64})

Return a vector of individuals' most recent common ancestors (MRCAs).
"""
function findMRCA(genealogy::OrderedDict{Int64, Individual}, IDs::Vector{Int64})
    ancestors = [ancestor(genealogy, ID) for ID in IDs]
    common_ancestors = ∩(ancestors...)
    older_common_ancestors = ancestor(genealogy, common_ancestors)
    most_recent_commmon_ancestors = setdiff(common_ancestors, older_common_ancestors)
    most_recent_commmon_ancestors
end