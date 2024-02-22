"""
    population(filename::String)

Return a dictionary of populations.

# Example

```@example
import GenLib as gen
pop140 = gen.pop140
pop = gen.population(pop140)
```
"""
function population(filename::String)
    population::Dict{Int64, String} = Dict()
    open(filename) do file
        rank = 0
        while !eof(file)
            line = readline(file)
            rank += 1
            if rank == 1
                continue
            end
            (ind, pop) = split(line)
            ind = parse(Int64, ind)
            population[ind] = pop
        end
    end
    population
end