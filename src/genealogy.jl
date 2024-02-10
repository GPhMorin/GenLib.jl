struct Individual
    father::Int64
    mother::Int64
    index::Int64
    children::Vector{Int64}
    sex::Int64
end

"""
    genealogy(dataframe::DataFrame)

Return an ordered dictionary of individuals from a `DataFrame`.
"""
function genealogy(dataframe::DataFrame)
    genealogy::OrderedDict{Int64, Individual} = Dict()
    for (index, row) in enumerate(eachrow(dataframe))
        genealogy[row.ind] = Individual(row.father, row.mother, index, [], row.sex)
    end
    for row in eachrow(dataframe)
        individual = genealogy[row.ind]
        if individual.father != 0
            push!(genealogy[individual.father].children, row.ind)
        end
        if individual.mother != 0
            push!(genealogy[individual.mother].children, row.ind)
        end
    end
    if !check_order(genealogy)
        println("Reordering the genealogy...")
        return order_genealogy(genealogy)
    else
        return genealogy
    end
end

"""
    genealogy(filename::String)

Return an ordered dictionary of individuals from a CSV file.
"""
function genealogy(filename::String)
    dataset = CSV.read(filename, DataFrame, delim='\t', types=Dict(:ind => Int64, :father => Int64, :mother => Int64, :sex => Int64))
    genealogy::OrderedDict{Int64, Individual} = Dict()
    for (index, row) in enumerate(eachrow(dataset))
        genealogy[row.ind] = Individual(row.father, row.mother, index, [], row.sex)
    end
    for row in eachrow(dataset)
        individual = genealogy[row.ind]
        if individual.father != 0
            push!(genealogy[individual.father].children, row.ind)
        end
        if individual.mother != 0
            push!(genealogy[individual.mother].children, row.ind)
        end
    end
    if !check_order(genealogy)
        println("Reordering the genealogy...")
        return order_genealogy(genealogy)
    else
        return genealogy
    end
end

"""
    check_order(genealogy::OrderedDict{Int64, Individual})

Return whether the individuals appear in chronological order in the dictionary,
i.e. any individual's parents appear before them.
"""
function check_order(genealogy::OrderedDict{Int64, Individual})
    value = true
    for (_, individual) in genealogy
        father = individual.father
        mother = individual.mother
        if (father != 0)
            if individual.index < genealogy[father].index
                value = false
                break
            end
        end
        if (mother != 0)
            if individual.index < genealogy[mother].index
                value = false
                break
            end
        end
    end
    value
end

"""
    order_genealogy(genealogy::OrderedDict{Int64, Individual})

Return a reordered pedigree where the individuals are in chronological order,
i.e. any individual's parents appear before them.
"""
function order_genealogy(genealogy::OrderedDict{Int64, Individual})
    founderIDs = founder(genealogy)
    sortedIDs = []
    queue = founderIDs
    while !isempty(queue)
        ID = popfirst!(queue)
        individual = genealogy[ID]
        if (individual.father ∈ sortedIDs) && (individual.mother ∈ sortedIDs)
            push!(sortedIDs, ID)
            filter!(x -> x ≠ ID, queue)
            push!(queue, individual.children...)
        elseif (individual.father == 0) && (individual.mother == 0)
            push!(sortedIDs, ID)
            filter!(x -> x ≠ ID, queue)
            push!(queue, individual.children...)
        elseif (individual.father == 0) && (individual.mother ∈ sortedIDs)
            push!(sortedIDs, ID)
            filter!(x -> x ≠ ID, queue)
            push!(queue, individual.children...)
        elseif (individual.mother == 0) && (individual.father ∈ sortedIDs)
            push!(sortedIDs, ID)
            filter!(x -> x ≠ ID, queue)
            push!(queue, individual.children...)
        else
            push!(queue, ID)
        end
    end
    ordered_genealogy = OrderedDict{Int64, Individual}()
    for (index, ID) in enumerate(sortedIDs)
        individual = genealogy[ID]
        ordered_genealogy[ID] = Individual(individual.father, individual.mother, index, individual.children, individual.sex)
    end
    ordered_genealogy
end

"""
    save_genealogy(genealogy::OrderedDict, path::String, sorted::Bool = false)

Export the pedigree as a CSV file at a given `path`.

If `sorted` is `false` (the default), then the individuals
will appear in the same order as in the genealogy.

If `sorted` is `true`, then the individuals
will appear in alphabetical ID order.
"""
function save_genealogy(genealogy::OrderedDict, path::String, sorted::Bool = false)
    df = genout(genealogy, sorted = sorted)
    CSV.write(path, df, delim="\t")
end