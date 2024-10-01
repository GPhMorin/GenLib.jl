"""
GenLib.jl: An unofficial, pure Julia port of R's GENLIB genetics and genealogical library.
"""
module GenLib

using DataFrames: DataFrame
using DataStructures: OrderedDict, SwissDict

# GENLIB datasets
"""
Genealogical information for 140 individuals from the Quebec Reference Sample.

According to the R GENLIB documentation, `genea140` corresponds to "a genealogical corpus
made of 41523 individuals from the province of Quebec, Canada. A total of 140 individuals
have been sampled in seven sub-populations, listed in pop140, and their genealogies were
reconstructed as far back as possible using the BALSAC population register and the Early
Quebec Population Register."

# Loading the Dataset

```julia
import GenLib as gen
genea140 = gen.genea140
ped = gen.genealogy(genea140)
```
"""
const genea140 = "$(chop(pathof(GenLib), tail=13))data/genea140.csv"

"""
A highly inbred pedigree.

According to the R GENLIB documentation, `geneaJi` corresponds to "a modified version of a
pedigree of two Jicaque Indians studied by Chapman & Jacquard (1971)."

# Loading the Dataset

```julia
import GenLib as gen
geneaJi = gen.geneaJi
ped = gen.genealogy(geneaJi)
```
"""
const geneaJi = "$(chop(pathof(GenLib), tail=13))data/geneaJi.csv"

"""
Population of origin of the 140 Quebec samples.

According to the R GENLIB documentation, `pop140` corresponds to "140 individuals from the
genealogical corpus from Quebec (…) sampled from 7 different populations from 5 regions:
Quebec City, Montreal, Saguenay, North Shore, Gaspesia. In Gaspesia we find 3 different
populations: French-Canadians, Acadians and Loyalists."

# Loading the Dataset

```julia
import GenLib as gen
pop140 = gen.pop140
pop = gen.population(pop140)
```
"""
const pop140 = "$(chop(pathof(GenLib), tail=13))data/pop140.csv"

"""
    _pop(filename::String)

Return a dictionary of populations.

# Example

```@example
import GenLib as gen
pop140 = gen.pop140
pop = gen.population(pop140)
```
"""
function _pop(filename::String)
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

include("create.jl")
include("extract.jl")
include("output.jl")
include("identify.jl")
include("describe.jl")
include("compute.jl")

end