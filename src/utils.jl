"""
pro(genealogy::OrderedDict{Int64, Individual})

Takes a `genealogy` dictionary and returns a vector of proband IDs.
"""
function pro(genealogy::OrderedDict{Int64, Individual})
    probands = [ID for (ID, individual) in genealogy if isempty(individual.children)]
    sort(probands)
end

"""
founder(genealogy::OrderedDict{Int64, Individual})

Takes a `genealogy` dictionary and returns a vector of founder IDs.
"""
function founder(genealogy::OrderedDict{Int64, Individual})
    founders = [ID for (ID, individual) in genealogy if (individual.father == 0) && (individual.mother == 0)]
    sort(founders)
end

"""
findFounders(genealogy::OrderedDict{Int64}, IDs::Vector{Int64})

Takes a `genealogy` and returns a vector of founders from whom the `IDs` descend.
"""
function findFounders(genealogy::OrderedDict{Int64}, IDs::Vector{Int64})
    ancestorIDs = [ancestor(genealogy, ID) for ID in IDs]
    common_ancestorIDs = ∩(ancestorIDs...)
    founderIDs = [ancestorID for ancestorID in common_ancestorIDs
                  if (genealogy[ancestorID].father == 0) & (genealogy[ancestorID].mother == 0)]
    founderIDs
end

"""
get_paths(genealogy::OrderedDict{Int64, Individual}, ID::Int64)

Takes a `genealogy` dictionary and an `ID` and returns the paths from an individual to their ancestors.
"""
function get_paths(genealogy::OrderedDict{Int64, Individual}, ID::Int64)
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
children(genealogy::OrderedDict{Int64, Individual}, ID::Int64)

Takes a `genealogy` dictionary and an `ID` and returns the children of an individual.
"""
function children(genealogy::OrderedDict{Int64, Individual}, ID::Int64)
    individual = genealogy[ID]
    individual.children
end

"""
descendant(genealogy::OrderedDict{Int64, Individual}, ID::Int64)

Takes a `genealogy` dictionary and an `ID` and returns the descendants of an individual.
"""
function descendant(genealogy::OrderedDict{Int64, Individual}, ID::Int64)
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
rec(genealogy::OrderedDict{Int64, Individual}, probandIDs::Vector{Int64} = pro(genealogy), ancestorIDs::Vector{Int64} = founder(genealogy))

Takes a `genealogy` dictionary, a vector of `probandIDs` and a vector of `ancestorIDs` and returns the number of descendants of each ancestor.
"""
function rec(
    genealogy::OrderedDict{Int64, Individual},
    probandIDs::Vector{Int64} = pro(genealogy),
    ancestorIDs::Vector{Int64} = founder(genealogy))
    
    coverage = Vector{Int64}()
    for ancestorID in ancestorIDs
        descendants = descendant(genealogy, ancestorID)
        descendants = filter!(x -> x in probandIDs, descendants)
        append!(coverage, length(descendants))
    end
    coverage
end

"""
distance_matrix(matrix::Matrix{Float64})

Converts a similarity matrix into a distance matrix normalized with values [0, 1].
"""
function distance_matrix(matrix::Matrix{Float64})
    if maximum(matrix) - minimum(matrix) == 0
        distance_matrix = zeros(size(matrix))
    else
        distance_matrix = (maximum(matrix .- minimum(matrix)) .- (matrix .- minimum(matrix)))/maximum(matrix .- minimum(matrix))
    end
end

"""
explore_tree(individual::ReferenceIndividual)

A recursive function that labels `individual`s on whether
they lead from relevant ancestors to relevant probands.
"""
function explore_tree(individual::ReferenceIndividual)
    # Ported from GENLIB's ExploreArbre
    state = individual.state
    if state == EXPLORED
        0
    elseif state == NODE
        1
    elseif state == PROBAND
        1
    elseif state == START
        for child in individual.children
            explore_tree(child)
        end
        1
    elseif state == EXPLOREDPROBAND
        state = PROBAND
        for child in individual.children
            explore_tree(child)
        end
        1
    elseif state == UNEXPLORED
        useful = 0
        for child in individual.children
            useful += explore_tree(child)
        end
        if useful > 0
            state = NODE
            1
        else
            state = EXPLORED
            0
        end
    end
    99
end

"""
prepare_priority_sort(nodes::Vector{ReferenceIndividual})

Initializes the value of the individuals' `sort` value.
"""
function prepare_priority_sort(nodes::Vector{ReferenceIndividual})
    # Ported from GENLIB's PrepareSortPrioriteArbre
    for node in nodes
        if isnothing(node.father)
            node.sort = -1
        elseif node.father.state == UNEXPLORED | node.father.state == EXPLORED
            node.sort = -1
        elseif isnothing(node.mother)
            node.sort = -1
        elseif node.mother.state == UNEXPLORED | node.mother.state == EXPLORED
            node.sort = -1
        elseif node.father.state == EXPLOREDPROBAND | node.mother.state == EXPLOREDPROBAND
            node.sort = -1
        else
            node.sort = 0
        end
    end
end

"""
start_priority_sort(node::ReferenceIndividual,
                    order::Dict{Int64, Int64},
                    index::Int64,
                    jumps::Dict{Int64, Int64})
"""
function start_priority_sort(node::ReferenceIndividual,
                             order::Dict{Int64, Int64},
                             index::Int64,
                             jumps::Dict{Int64, Int64})


    # Ported from GENLIB's StartSortPrioriteArbre

    node.sort = 5
    nodelist = []

    for child in node.children
        if child.sort == -1
            priority_sort!(child, order, index, jumps, nodelist)
        end
    end
    empty!(nodelist)

    for child in node.children
        if child.sort == -1
        elseif child.sort == 0
            child.sort = 1
        elseif child.sort == 1
            priority_sort!(child, order, index, jumps, nodelist)
        end
    end
    empty!(nodelist)
end

"""
priority_sort!(node::ReferenceIndividual, 
               order::Dict{Int64, ReferenceIndividual},
               index::Int64,
               jumps::Dict{Int64, Int64},
               nodelist::Vector{ReferenceIndividual})
"""
function priority_sort!(node::ReferenceIndividual, 
                        order::Dict{Int64, ReferenceIndividual},
                        index::Int64,
                        jumps::Dict{Int64, Int64},
                        nodelist::Vector{ReferenceIndividual})

    # Ported from GENLIB's SortPrioriteArbre

    jump = 0

    if !isempty(nodelist)
        return 0
    end

    if (node.state != PROBAND & node.state != NODE) | node.sort == 5
        return 0
    end

    order[index] = node
    index += 1

    former_flag = copy(node.sort)
    node.sort = 5

    for child in node.children
        if child.sort == -1
            jump += priority_sort!(child, order, index, jumps, nodelist)
        end
    end

    jumps[index] = jump

    if former_flag == -1
        jump += 1
    end

    for child in node.children
        if child.sort == -1
        elseif child.sort == 0
            child.sort = 1
        elseif child.sort == 1
            push!(nodelist, child)
        end
    end

    jump
end