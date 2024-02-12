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

Return a tuple of type `::Tuple{::Vector{Int64}, ::Matrix{Int64}}` consisting in
the vector of individuals' most recent common ancestors (MRCAs)
the matrix of meioses between each proband and each MRCA.
"""
function findMRCA(pedigree::OrderedDict{Int64, Individual}, probandIDs::Vector{Int64})
    ancestorIDs = [ancestor(pedigree, ID) for ID in probandIDs]
    common_ancestorIDs = ∩(ancestorIDs...)
    older_common_ancestorIDs = ancestor(pedigree, common_ancestorIDs)
    mrcaIDs = sort(collect(setdiff(common_ancestorIDs, older_common_ancestorIDs)))
    meioses_matrix = Matrix{Int64}(undef, length(probandIDs), length(mrcaIDs))
    for (i, probandID) in enumerate(probandIDs), (j, mrcaID) in enumerate(mrcaIDs)
        meioses_matrix[i, j] = findDistance(pedigree, probandID, mrcaID)
    end
    mrcaIDs, meioses_matrix
end