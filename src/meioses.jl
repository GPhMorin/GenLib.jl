"""
    _findDistance(pedigree::Pedigree, descendantID::Int64, ancestorID::Int64)

Return a vector of distances between an individual and their ancestor.
"""
function _findDistance(
    pedigree::Pedigree,
    descendantID::Int64,
    ancestorID::Int64)
    
    descendant = pedigree[descendantID]
    paths = get_paths(pedigree, descendant)
    lengths = Vector{Int8}()
    for path in paths
        if path[1] ≡ ancestorID
            push!(lengths, length(path) - 1)
        end
    end
    lengths
end

"""
    _findMinDistance(pedigree::Pedigree, descendantID::Int64, ancestorID::Int64)

Return the minimum distance between an individual and their ancestor.
"""
function _findMinDistance(pedigree::Pedigree, descendantID::Int64, ancestorID::Int64)
    lengths = _findDistance(pedigree, descendantID, ancestorID)
    minimum(lengths)
end

"""
    findDistance(pedigree::Pedigree, IDs::Vector{Int64}, ancestorID::Int64)

Return the distance between two individuals and their ancestor.
"""
function findDistance(pedigree::Pedigree, IDs::Vector{Int64}, ancestorID::Int64)
    distance₁ = _findMinDistance(pedigree, IDs[1], ancestorID)
    distance₂ = _findMinDistance(pedigree, IDs[2], ancestorID)
    distance₁ + distance₂
end