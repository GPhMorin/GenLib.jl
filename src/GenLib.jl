module GenLib

using CSV
using DataFrames: DataFrame
using DataStructures: Stack
using Distributions

# For creating new genealogies manually
export Individual

# GENLIB functions
export genealogy,
       isolate_genealogy,
       gc,
       phi,
       ϕ,
       findDistances,
       findDistance,
       ancestor,
       findMRCA,
       point,
       pro,
       founder,
       findFounders,
       get_paths,
       children,
       descendant,
       rec,
       distance_matrix,
       population,
       remove_relatives!,
       occ,
       simulate

# GENLIB datasets
const genea140 = "$(chop(pathof(GenLib), tail=13))data/genea140.asc"
const geneaJi = "$(chop(pathof(GenLib), tail=13))data/geneaJi.asc"
const pop140 = "$(chop(pathof(GenLib), tail=13))data/pop140.csv"

export genea140,
       geneaJi,
       pop140

include("genealogy.jl")
include("pointers.jl")
include("genetic_contribution.jl")
include("utils.jl")
include("kinship.jl")
include("mrca.jl")
include("meioses.jl")
include("population.jl")
include("remove_relatives.jl")
include("occurrence.jl")
include("simulations.jl")

end
