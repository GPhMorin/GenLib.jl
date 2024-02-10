"""
    pro(genealogy::OrderedDict{Int64, Individual})

Return a vector of proband IDs in alphabetical order.
"""
function pro(genealogy::OrderedDict{Int64, Individual})
    probands = [ID for (ID, individual) in genealogy if isempty(individual.children)]
    sort(probands)
end

"""
    founder(genealogy::OrderedDict{Int64, Individual})

Return a vector of founder IDs in alphabetical order.
"""
function founder(genealogy::OrderedDict{Int64, Individual})
    founders = [ID for (ID, individual) in genealogy if (individual.father == 0) && (individual.mother == 0)]
    sort(founders)
end

"""
    findFounders(genealogy::OrderedDict{Int64}, IDs::Vector{Int64})

Return a vector of founders from whom the `IDs` descend.
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

Return the paths from an individual to their ancestors.
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

Return the children of an individual.
"""
function children(genealogy::OrderedDict{Int64, Individual}, ID::Int64)
    individual = genealogy[ID]
    individual.children
end

"""
    descendant(genealogy::OrderedDict{Int64, Individual}, ID::Int64)

Return the descendants of an individual.
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

Return the number of descendants of each ancestor.
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
    genout(genealogy::OrderedDict, sorted::Bool = false)

Return a pedigree as a DataFrame.

If `sorted` is `false` (the default), then the individuals
will appear in the same order as in the pedigree.

If `sorted` is `true`, then the individuals
will appear in alphabetical ID order.
"""
function genout(genealogy::OrderedDict; sorted::Bool = false)
    inds = Int64[]
    fathers = Int64[]
    mothers = Int64[]
    sexes = Int64[]
    if !sorted
        for (ID, individual) in genealogy
            push!(inds, ID)
            push!(fathers, individual.father)
            push!(mothers, individual.mother)
            push!(sexes, individual.sex)
        end
    else # if sorted
        inds = sort(collect(keys(genealogy)))
        for ID in inds
            individual = genealogy[ID]
            push!(fathers, individual.father)
            push!(mothers, individual.mother)
            push!(sexes, individual.sex)
        end
    end
    df = DataFrame([inds, fathers, mothers, sexes], ["ind", "father", "mother", "sex"])
    df
end