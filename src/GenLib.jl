module GenLib

using CSV
using DataFrames: DataFrame
using DataStructures: Stack, OrderedDict

# For evaluating the sex manually
export SEX,
       MALE,
       FEMALE

# For creating new genealogies manually
export Individual

# GENLIB functions
export genealogy,
       gc,
       phi,
       Φ,
       findDistances,
       findDistance,
       ancestor,
       findMRCA,
       refer,
       branching,
       pro,
       founder,
       findFounders,
       get_paths,
       children,
       descendant,
       rec

# Custom functions
export distance_matrix,
       population,
       remove_relatives!,
       occ,
       check_order,
       save_genealogy

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
