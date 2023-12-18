@doc read(joinpath(dirname(@__DIR__), "README.md"), String)

module GenLib

using CSV
using DataFrames: DataFrame
using DataStructures: Stack

# For creating new genealogies manually
export Individual

# GENLIB functions
export genealogy,
       gc,
       phi,
       findDistances,
       findDistance,
       ancestor,
       findMRCA,
       point,
       pro,
       founder,
       get_paths,
       children,
       descendant,
       rec

# GENLIB datasets
const genea140 = "$(chop(pathof(GenLib), tail=13))/data/genea140.asc"
const genJi = "$(chop(pathof(GenLib), tail=13))/data/geneaJi.asc"
const pop140 = "$(chop(pathof(GenLib), tail=13))/data/pop140.asc"

export genea140,
       genJi,
       pop140

include("genealogy.jl")
include("pointers.jl")
include("genetic_contribution.jl")
include("utils.jl")
include("kinship.jl")
include("mrca.jl")
include("meioses.jl")

end
