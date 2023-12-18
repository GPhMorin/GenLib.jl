"""
population(filename::String)::Dict{Int64, String}

Reads a file with `filename` and returns a dictionary of populations.
"""
function population(filename::String)::Dict{Int64, String}
    dataset = CSV.read(filename, DataFrame, delim='\t', types=Dict(:ind => Int64, :father => Int64, :mother => Int64, :sex => Int8))
    population::Dict{Int64, String} = Dict()
    for row in eachrow(dataset)
        population[row.ind] = row.pop
    end
    population
end