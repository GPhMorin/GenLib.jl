@enum SEX begin
    MALE = 1
    FEMALE = 2
end

"""
An `Individual` is an immutable object containing
the ID to their `father` and `mother` (0 if unknown),
the `index` in which it appeared in the ancestry file,
the IDs of their `children`,
and their `sex` (1 for a male, 2 for a female).
"""
struct Individual
    father::Int64
    mother::Int64
    index::Int64
    children::Vector{Int64}
    sex::SEX
end

"""
genealogy(filename::String)

Reads a DataFrame and returns a dictionary of individuals.
"""
function genealogy(dataframe::DataFrame)
    genealogy::OrderedDict{Int64, Individual} = Dict()
    for (index, row) in enumerate(eachrow(dataframe))
        genealogy[row.ind] = Individual(row.father, row.mother, index, [], SEX(row.sex))
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

Reads a file with `filename` and returns a dictionary of individuals.
"""
function genealogy(filename::String)
    dataset = CSV.read(filename, DataFrame, delim='\t', types=Dict(:ind => Int64, :father => Int64, :mother => Int64, :sex => Int64))
    genealogy::OrderedDict{Int64, Individual} = Dict()
    for (index, row) in enumerate(eachrow(dataset))
        genealogy[row.ind] = Individual(row.father, row.mother, index, [], SEX(row.sex))
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

Takes a given `genealogy` dictionary and checks whether
the individuals appear in chronological order,
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

Given a `genealogy` dictionary, returns a reordered genealogy
where the individuals are in chronological order,
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