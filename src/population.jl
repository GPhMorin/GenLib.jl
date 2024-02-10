"""
    population(filename::String)

Return a dictionary of populations.
"""
function population(filename::String)
    dataset = CSV.read(filename, DataFrame, delim='\t', types=Dict(:ind => Int64, :pop => String))
    population::Dict{Int64, String} = Dict()
    for row in eachrow(dataset)
        population[row.ind] = row.pop
    end
    population
end