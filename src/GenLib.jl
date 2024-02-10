"""
GenLib.jl: An unofficial, pure Julia port of R's GENLIB genetics and genealogical library.
"""
module GenLib

using CSV
using DataFrames: DataFrame
using DataStructures: Stack, OrderedDict

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
       phi,
       pro,
       rec

# Custom functions
export check_order,
       descendant,
       distance_matrix,
       findDistances,
       get_paths,
       occ,
       population,
       refer,
       remove_relatives!,
       save_genealogy,
       ϕ

# GENLIB datasets
const genea140 = "$(chop(pathof(GenLib), tail=13))data/genea140.asc"
const geneaJi = "$(chop(pathof(GenLib), tail=13))data/geneaJi.asc"
const pop140 = "$(chop(pathof(GenLib), tail=13))data/pop140.csv"

export genea140,
       geneaJi,
       pop140

include("genealogy.jl")
include("reference.jl")
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
