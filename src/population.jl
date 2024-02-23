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
        is_firstline = true
        while !eof(file)
            line = readline(file)
            if is_firstline
                is_firstline = false
                continue
            end
            (ind, pop) = split(line)
            ind = parse(Int64, ind)
            population[ind] = pop
        end
    end
    population
end