@doc read(joinpath(dirname(@__DIR__), "README.md"), String)

using CSV
using DataFrames: DataFrame
using DataStructures: Stack

module GenLib

include("genealogy.jl")
include("genetic_contribution.jl")
include("pointers.jl")
include("utils.jl")
include("kinship.jl")
include("mrca.jl")
include("meioses.jl")

export Individual # For creating new genealogies manually

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

end