"""
GenLib.jl: An unofficial, pure Julia port of R's GENLIB genetics and genealogical library.
"""
module GenLib

using CSV
using DataFrames: DataFrame
using DataStructures: Stack, OrderedDict

# The unit structure for pedigrees
export Individual,
       STATE

# GENLIB functions
export ancestor,
       branching,
       children,
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

# Custom functions
export check_order,
       descendant,
       findDistances,
       get_paths,
       population,
       refer,
       remove_relatives!,
       save_genealogy,
       ϕ

# GENLIB datasets
"""
Genealogical information for 140 individuals from the Quebec Reference Sample.

According to the R GENLIB documentation, `genea140` corresponds to
"a genealogical corpus made of 41523 individuals from the province
of Quebec, Canada. A total of 140 individuals have been sampled in
seven sub-populations, listed in pop140, and their genealogies
were reconstructed as far back as possible using the BALSAC population
register and the Early Quebec Population Register."
"""
const genea140 = "$(chop(pathof(GenLib), tail=13))data/genea140.asc"

"""
A highly inbred pedigree.

According to the R GENLIB documentation, `geneaJi` corresponds to
"a modified version of a pedigree of two Jicaque Indians studied by
Chapman & Jacquard (1971)."
"""
const geneaJi = "$(chop(pathof(GenLib), tail=13))data/geneaJi.asc"

"""
Population of origin of the 140 Quebec samples.

According to the R GENLIB documentation, `pop140` corresponds to
"140 individuals from the genealogical corpus from Quebec (…) sampled
from 7 different opulations from 5 regions: Quebec City, Montreal, Saguenay,
North Shore, Gaspesia. In Gaspesia we find 3 different populations:
French-Canadians, Acadians and Loyalists."
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
include("mrca.jl")
include("meioses.jl")
include("population.jl")
include("remove_relatives.jl")
include("occurrence.jl")

end
