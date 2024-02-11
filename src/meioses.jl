"""
    findDistances(pedigree::OrderedDict{Int64, Individual}, descendantID::Int64, ancestorID::Int64)

Return a vector of distances between an individual and their ancestor.
"""
function findDistances(
    pedigree::OrderedDict{Int64, Individual},
    descendantID::Int64,
    ancestorID::Int64)
    
    paths = get_paths(pedigree, descendantID)
    lengths = Vector{Int8}()
    for path in paths
        if path[1] ≡ ancestorID
            push!(lengths, length(path) - 1)
        end
    end
    lengths
end

"""
    findDistance(pedigree::OrderedDict{Int64, Individual}, descendantID::Int64, ancestorID::Int64)

Return the minimum distance between an individual and their ancestor.
"""
function findDistance(pedigree::OrderedDict{Int64, Individual}, descendantID::Int64, ancestorID::Int64)
    lengths = findDistances(pedigree, descendantID, ancestorID)
    minimum(lengths)
end

"""
    findDistance(genealogy::OrderedDict{Int64, Individual}, ID₁::Int64, ID₂::Int64, ancestorID::Int64)

Return the distance between two individuals and their ancestor.
"""
function findDistance(pedigree::OrderedDict{Int64, Individual}, ID₁::Int64, ID₂::Int64, ancestorID::Int64)
    distance₁ = findDistance(pedigree, ID₁, ancestorID)
    distance₂ = findDistance(pedigree, ID₂, ancestorID)
    distance₁ + distance₂
end