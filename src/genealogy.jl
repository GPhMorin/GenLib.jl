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
    sex::Int8
end

"""
genealogy(filename::String)::Dict{Int64, Individual}

Reads a file with `filename` and returns a dictionary of individuals.
"""
function genealogy(filename::String)::Dict{Int64, Individual}
    dataset = CSV.read(filename, DataFrame, delim='\t', types=Dict(:ind => Int64, :father => Int64, :mother => Int64, :sex => Int8))
    genealogy::Dict{Int64, Individual} = Dict()
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
    genealogy
end