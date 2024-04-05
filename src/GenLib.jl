"""
GenLib.jl: An unofficial, pure Julia port of R's GENLIB genetics and genealogical library.
"""
module GenLib

using DataFrames: DataFrame
using DataStructures: OrderedDict

# The unit structure for pedigrees
export Individual,
       Pedigree

# GENLIB types and functions
export GenMatrix,
       ancestor,
       branching,
       children,
       f,
       findDistance,
       findFounders,
       findMRCA,
       founder,
       gc,
       genealogy,
       genout,
       occ,
       phi,
       pro,
       rec

# GENLIB datasets
"""
Genealogical information for 140 individuals from the Quebec Reference Sample.

According to the R GENLIB documentation, `genea140` corresponds to
"a genealogical corpus made of 41523 individuals from the province
of Quebec, Canada. A total of 140 individuals have been sampled in
seven sub-populations, listed in pop140, and their genealogies
were reconstructed as far back as possible using the BALSAC population
register and the Early Quebec Population Register."

# Loading the Dataset

```julia
import GenLib as gen
genea140 = gen.genea140
ped = gen.genealogy(genea140)
```
"""
const genea140 = "$(chop(pathof(GenLib), tail=13))data/genea140.asc"

"""
A highly inbred pedigree.

According to the R GENLIB documentation, `geneaJi` corresponds to
"a modified version of a pedigree of two Jicaque Indians studied by
Chapman & Jacquard (1971)."

# Loading the Dataset

```julia
import GenLib as gen
geneaJi = gen.geneaJi
ped = gen.genealogy(geneaJi)
```
"""
const geneaJi = "$(chop(pathof(GenLib), tail=13))data/geneaJi.asc"

"""
Population of origin of the 140 Quebec samples.

According to the R GENLIB documentation, `pop140` corresponds to
"140 individuals from the genealogical corpus from Quebec (…) sampled
from 7 different opulations from 5 regions: Quebec City, Montreal, Saguenay,
North Shore, Gaspesia. In Gaspesia we find 3 different populations:
French-Canadians, Acadians and Loyalists."

# Loading the Dataset

```julia
import GenLib as gen
pop140 = gen.pop140
pop = gen.population(pop140)
```
"""
const pop140 = "$(chop(pathof(GenLib), tail=13))data/pop140.csv"

export genea140,
       geneaJi,
       pop140

include("genealogy.jl")
include("isolate.jl")
include("utils.jl")
include("genetic_contribution.jl")
include("kinship.jl")
include("kinship2.jl")
include("mrca.jl")
include("meioses.jl")
include("population.jl")
include("occurrence.jl")

end