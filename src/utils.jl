"""
pro(genealogy::Dict{Int64, Individual})::Vector{Int64}

Takes a `genealogy` dictionary and returns a vector of proband IDs.
"""
function pro(genealogy::Dict{Int64, Individual})::Vector{Int64}
    probands = [ID for (ID, individual) in genealogy if isempty(individual.children)]
    sort!(probands)
end

"""
founder(genealogy::Dict{Int64, Individual})::Vector{Int64}

Takes a `genealogy` dictionary and returns a vector of founder IDs.
"""
function founder(genealogy::Dict{Int64, Individual})::Vector{Int64}
    founders = [ID for (ID, individual) in genealogy if individual.father == 0]
    sort!(founders)
end

"""
findFounders(genealogy::Dict{Int64}, IDs::Vector{Int64})::Vector{Int64}

Takes a `genealogy` and returns a vector of founders from whom the `IDs` descend.
"""
function findFounders(genealogy::Dict{Int64}, IDs::Vector{Int64})::Vector{Int64}
    ancestorIDs = [ancestor(genealogy, ID) for ID in IDs]
    common_ancestorIDs = ∩(ancestorIDs...)
    founderIDs = [ancestorID for ancestorID in common_ancestorIDs
                  if (genealogy[ancestorID].father == 0) & (genealogy[ancestorID].mother == 0)]
    founderIDs
end

"""
get_paths(genealogy::Dict{Int64, Individual}, ID::Int64)::Vector{Vector{Int64}}

Takes a `genealogy` dictionary and an `ID` and returns the paths from an individual to their ancestors.
"""
function get_paths(genealogy::Dict{Int64, Individual}, ID::Int64)::Vector{Vector{Int64}}
    paths = Vector{Vector{Int64}}([[ID]])
    individual = genealogy[ID]
    if individual.father != 0
        fathers_paths = get_paths(genealogy, individual.father)
        for path in fathers_paths
            push!(path, ID)
        end
        append!(paths, fathers_paths)
    end
    if individual.mother != 0
        mothers_paths = get_paths(genealogy, individual.mother)
        for path in mothers_paths
            push!(path, ID)
        end
        append!(paths, mothers_paths)
    end
    paths
end

"""
children(genealogy::Dict{Int64, Individual}, ID::Int64)::Vector{Int64}

Takes a `genealogy` dictionary and an `ID` and returns the children of an individual.
"""
function children(genealogy::Dict{Int64, Individual}, ID::Int64)::Vector{Int64}
    individual = genealogy[ID]
    individual.children
end

"""
descendant(genealogy::Dict{Int64, Individual}, ID::Int64)::Vector{Int64}

Takes a `genealogy` dictionary and an `ID` and returns the descendants of an individual.
"""
function descendant(genealogy::Dict{Int64, Individual}, ID::Int64)::Vector{Int64}
    descendants = Set{Int64}()
    stack = Stack{Int64}()
    push!(stack, ID)
    while !isempty(stack)
        ID = pop!(stack)
        individual = genealogy[ID]
        for child in individual.children
            push!(stack, child)
            push!(descendants, child)
        end
    end
    sort(collect(descendants))
end

"""
rec(genealogy::Dict{Int64, Individual}, probandIDs::Vector{Int64} = pro(genealogy), ancestorIDs::Vector{Int64} = founder(genealogy))::Vector{Int64}

Takes a `genealogy` dictionary, a vector of `probandIDs` and a vector of `ancestorIDs` and returns the number of descendants of each ancestor.
"""
function rec(
    genealogy::Dict{Int64, Individual},
    probandIDs::Vector{Int64} = pro(genealogy),
    ancestorIDs::Vector{Int64} = founder(genealogy)
    )::Vector{Int64}
    
    coverage = Vector{Int64}()
    for ancestorID in ancestorIDs
        descendants = descendant(genealogy, ancestorID)
        descendants = filter!(x -> x in probandIDs, descendants)
        append!(coverage, length(descendants))
    end
    coverage
end

"""
distance_matrix(similarity_matrix::Matrix{Float64})::Matrix{Float64}

Converts a similarity matrix into a distance matrix normalized with values [0, 1].
"""
function distance_matrix(matrix::Matrix{Float64})::Matrix{Float64}
    if maximum(matrix) - minimum(matrix) == 0
        distance_matrix = zeros(size(matrix))
    else
        distance_matrix = (maximum(matrix .- minimum(matrix)) .- (matrix .- minimum(matrix)))/maximum(matrix .- minimum(matrix))
    end
end