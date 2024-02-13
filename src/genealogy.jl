"""
    @enum STATE begin
        PROBAND
        FOUNDER
        UNEXPLORED
    end

An enumeration of the [`Individual`](@ref)'s state.
"""
@enum STATE begin
    PROBAND
    FOUNDER
    UNEXPLORED
end

"""
    mutable struct Individual
        ID::Int64
        father::Union{Nothing, Individual}
        mother::Union{Nothing, Individual}
        index::Int64
        children::Vector{Individual}
        sex::Int64
        state::STATE
        probability::Float64
        sort::Int64
        ancestor::Bool
        descendant::Bool
        occurrence::Int64
    end

The unit structure of a pedigree.
"""
mutable struct Individual
    ID::Int64
    father::Union{Nothing, Individual}
    mother::Union{Nothing, Individual}
    index::Int64
    children::Vector{Individual}
    sex::Int64
    state::STATE
    probability::Float64
    sort::Int64
    ancestor::Bool
    descendant::Bool
    occurrence::Int64
end

function Base.show(io::IO, individual::Individual)
    println(io, "ind: ", individual.ID)
    println(io, "father: ", !isnothing(individual.father) ? individual.father.ID : 0)
    println(io, "mother: ", !isnothing(individual.mother) ? individual.mother.ID : 0)
    print(io, "sex: ", individual.sex)
end

"""
    Pedigree = OrderedDict{Int64, Individual}

A particular case of an `OrderedDict` containing individuals accessed by ID.
"""
Pedigree = OrderedDict{Int64, Individual}

function Base.show(io::IO, ::MIME"text/plain", pedigree::Pedigree)
    n = 0
    parent_child = 0
    men = 0
    women = 0
    probands = 0
    depth = 0
    for (_, individual) in pedigree
        n += 1
        children = length(individual.children)
        if children > 0
            parent_child += children
        else
            probands += 1
        end
        if individual.sex == 1
            men += 1
        else
            women += 1
        end
        depth = max(depth, max_depth(individual))
    end
    println(io, "A pedigree with:")
    s = n == 1 ? "" : "s"
    println(io, n, " individual", s, ";")
    s = parent_child == 1 ? "" : "s"
    println(io, parent_child, " parent-child relation", s, ";")
    man = men == 1 ? "man" : "men"
    println(io, men, " ", man, ";")
    woman = women == 1 ? "woman" : "women"
    println(io, women, " ", woman, ";")
    s = probands == 1 ? "" : "s"
    println(io, probands, " subject", s, ";")
    s = depth == 1 ? "" : "s"
    print(io, depth, " generation", s, ".")
end

"""
    genealogy(dataframe::DataFrame)

Return an ordered pedigree of individuals from a `DataFrame`.

# Example

```julia
import GenLib as gen
using DataFrames
inds = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
fathers = [0, 0, 0, 1, 1, 0, 3, 3, 6, 6]
mothers = [0, 0, 0, 2, 2, 0, 4, 4, 5, 5]
sexes = [1, 2, 1, 2, 2, 1, 2, 1, 1, 2]
df = DataFrame([inds, fathers, mothers, sexes], [:ind, :father, :mother, :sex])
ped = gen.genealogy(df)
```
"""
function genealogy(dataframe::DataFrame)
    pedigree = Pedigree()
    for (index, row) in enumerate(eachrow(dataframe))
        pedigree[row.ind] = Individual(
            row.ind, # ID
            nothing, # father
            nothing, # mother
            index, # index in the genealogy
            [], # children
            row.sex, # sex (1 for male, 2 for female)
            UNEXPLORED, # status
            0., # probability
            0, # sort
            false, # whether the individual is an ancestor
            false, # whether the individual is a descendant
            0) # occurrence
    end
    for row in eachrow(dataframe)
        father = row.father
        if father != 0
            pedigree[row.ind].father = pedigree[father]
            push!(pedigree[father].children, pedigree[row.ind])
        end
        mother = row.mother
        if mother != 0
            pedigree[row.ind].mother = pedigree[mother]
            push!(pedigree[mother].children, pedigree[row.ind])
        end
    end
    order_pedigree!(pedigree)
end

"""
    genealogy(filename::String)

Return an ordered pedigree of individuals from a CSV file.

# Example

```julia
import GenLib as gen
genea140 = gen.genea140
ped = gen.genealogy(genea140)
```
"""
function genealogy(filename::String)
    dataset = CSV.read(filename, DataFrame, delim='\t', types=Dict(:ind => Int64, :father => Int64, :mother => Int64, :sex => Int64))
    pedigree = Pedigree()
    for (index, row) in enumerate(eachrow(dataset))
        pedigree[row.ind] = Individual(
            row.ind, # ID
            nothing, # father
            nothing, # mother
            index, # index in the genealogy
            [], # children
            row.sex, # sex (1 for male, 2 for female)
            UNEXPLORED, # status
            0., # probability
            0, # sort
            false, # whether the individual is an ancestor
            false, # whether the individual is a descendant
            0) # occurrence
    end
    for row in eachrow(dataset)
        father = row.father
        if father != 0
            pedigree[row.ind].father = pedigree[father]
            push!(pedigree[father].children, pedigree[row.ind])
        end
        mother = row.mother
        if mother != 0
            pedigree[row.ind].mother = pedigree[mother]
            push!(pedigree[mother].children, pedigree[row.ind])
        end
    end
    order_pedigree!(pedigree)
end

"""
    order_pedigree(pedigree::Pedigree)

Return a reordered pedigree where the individuals are in chronological order,
i.e. any individual's parents appear before them.
"""
function order_pedigree!(pedigree::Pedigree)
    individuals = [individual for (ID, individual) in pedigree]
    depths = [max_depth(individual) for individual in individuals]
    order = sortperm(depths)
    sorted_individuals = individuals[order]
    empty!(pedigree)
    for (index, individual) in enumerate(sorted_individuals)
        individual.index = index
        pedigree[individual.ID] = individual
    end
    pedigree
end

"""
    save_genealogy(pedigree::Pedigree, path::String, sorted::Bool = false)

Export the pedigree as a CSV file at a given `path`.

If `sorted` is `false` (the default), then the individuals
will appear in the same order as in the genealogy.

If `sorted` is `true`, then the individuals
will appear in alphabetical ID order.
"""
function save_genealogy(pedigree::Pedigree, path::String, sorted::Bool = false)
    df = genout(pedigree, sorted = sorted)
    CSV.write(path, df, delim="\t")
end