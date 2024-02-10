"""
    findDistances(genealogy::OrderedDict{Int64, Individual}, descendantID::Int64, ancestorID::Int64)

Return a vector of distances between an individual and their ancestor.
"""
function findDistances(
    genealogy::OrderedDict{Int64, Individual},
    descendantID::Int64,
    ancestorID::Int64)
    
    paths = get_paths(genealogy, descendantID)
    lengths = Vector{Int8}()
    for path in paths
        if path[1] ≡ ancestorID
            push!(lengths, length(path) - 1)
        end
    end
    lengths
end

"""
    findDistance(genealogy::OrderedDict{Int64, Individual}, descendantID::Int64, ancestorID::Int64)

Return the minimum distance between an individual and their ancestor.
"""
function findDistance(genealogy::OrderedDict{Int64, Individual}, descendantID::Int64, ancestorID::Int64)
    lengths = findDistances(genealogy, descendantID, ancestorID)
    minimum(lengths)
end

"""
    findDistance(genealogy::OrderedDict{Int64, Individual}, ID₁::Int64, ID₂::Int64, ancestorID::Int64)

Return the distance between two individuals and their ancestor.
"""
function findDistance(genealogy::OrderedDict{Int64, Individual}, ID₁::Int64, ID₂::Int64, ancestorID::Int64)
    distance₁ = findDistance(genealogy, ID₁, ancestorID)
    distance₂ = findDistance(genealogy, ID₂, ancestorID)
    distance₁ + distance₂
end