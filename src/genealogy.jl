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
genealogy(filename::String)::Dict{Int64, Individual}

Reads a DataFrame and returns a dictionary of individuals.
"""
function genealogy(dataframe::DataFrame)::Dict{Int64, Individual}
    genealogy::Dict{Int64, Individual} = Dict()
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
    genealogy
end

"""
genealogy(filename::String)::Dict{Int64, Individual}

Reads a file with `filename` and returns a dictionary of individuals.
"""
function genealogy(filename::String)::Dict{Int64, Individual}
    dataset = CSV.read(filename, DataFrame, delim='\t', types=Dict(:ind => Int64, :father => Int64, :mother => Int64, :sex => Int64))
    genealogy::Dict{Int64, Individual} = Dict()
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
    genealogy
end

"""
isolate_genealogy(genealogy::Dict{Int64, Individual}, ancestor::Int64)::Dict{Int64, Individual}

Takes a `genealogy` and ablates parents from individuals who do not descend from the `ancestor`.
"""
function isolate_genealogy(genealogy::Dict{Int64, Individual}, ancestor::Int64)::Dict{Int64, Individual}
    descendants = descendant(genealogy, ancestor)
    isolated_genealogy = Dict{Int64, Individual}()
    for (ID, individual) in genealogy
        if ID ∈ descendants
            isolated_genealogy[ID] = individual
        else
            isolated_genealogy[ID] = Individual(0, 0, individual.index, individual.children, individual.sex)
        end
    end
    isolated_genealogy
end