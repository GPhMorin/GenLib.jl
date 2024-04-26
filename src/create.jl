"""
    abstract type AbstractIndividual end

The abstract type of `Individual` and its subtypes.
"""
abstract type AbstractIndividual end

"""
    struct IntIndividual <: AbstractIndividual
        ID::Int64
        father::Int64
        mother::Int64
        sex::Int64
        max_depth::Int64
    end

The temporary unit structure of a [`GenLib.Pedigree`](@ref).
"""
mutable struct IntIndividual <: AbstractIndividual
    ID::Int64
    father::Int64
    mother::Int64
    sex::Int64
    max_depth::Int64
end

"""
    struct Individual <: AbstractIndividual
        ID::Int64
        father::Union{Nothing, Individual}
        mother::Union{Nothing, Individual}
        children::Vector{Individual}
        sex::Int64
        rank::Int64
    end

The unit structure of a [`GenLib.Pedigree`](@ref).
"""
struct Individual <: AbstractIndividual
    ID::Int64
    father::Union{Nothing, Individual}
    mother::Union{Nothing, Individual}
    children::Vector{Individual}
    sex::Int64
    rank::Int64
end

function Base.show(io::IO, individual::T) where T <: AbstractIndividual
    println(io, "ind: ", individual.ID)
    println(io, "father: ", !isnothing(individual.father) ? individual.father.ID : 0)
    println(io, "mother: ", !isnothing(individual.mother) ? individual.mother.ID : 0)
    print(io, "sex: ", individual.sex)
end

"""
    struct Pedigree{T}

A minimal structure wrapping an `OrderedDict` with individuals accessed by ID.
"""
struct Pedigree{T}
    dict::OrderedDict{Int64, T}
end

Pedigree{T}() where T <: AbstractIndividual = Pedigree(OrderedDict{Int64, T}())

Base.copy(p::Pedigree{T}) where T <: AbstractIndividual = Pedigree{T}(copy(p.dict))
Base.length(p::Pedigree{T}) where T <: AbstractIndividual = length(p.dict)
Base.setindex!(p::Pedigree{T}, value::T, key::Int64) where T <: AbstractIndividual = setindex!(p.dict, value, key)
Base.getindex(p::Pedigree{T}, key::Int64) where T <: AbstractIndividual = p.dict[key]
Base.keys(p::Pedigree{T}) where T <: AbstractIndividual = keys(p.dict)
Base.values(p::Pedigree{T}) where T <: AbstractIndividual = values(p.dict)
Base.iterate(p::Pedigree{T}) where T <: AbstractIndividual = iterate(p.dict)
Base.iterate(p::Pedigree{T}, i) where T <: AbstractIndividual = iterate(p.dict, i)

function Base.show(io::IO, ::MIME"text/plain", pedigree::Pedigree)
    n = 0
    parent_child = 0
    men = 0
    women = 0
    probands = 0
    depth = 0
    for (_, individual) ∈ pedigree
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
        depth = max(depth, _max_depth(individual))
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
    genealogy(dataframe::DataFrame; sort::Bool = true)

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
function genealogy(dataframe::DataFrame; sort::Bool = true)
    pedigree = Pedigree{IntIndividual}()
    for row ∈ eachrow(dataframe)
        pedigree[row.ind] = IntIndividual(
            row.ind,
            row.father,
            row.mother,
            row.sex,
            -1
        )
    end
    if sort
        pedigree = _ordered_pedigree(pedigree)
    end
    _finalize_pedigree(pedigree)
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
function genealogy(filename::String; sort = true)
    pedigree = Pedigree{IntIndividual}()
    open(filename) do file
        is_firstline = true
        while !eof(file)
            line = readline(file)
            if is_firstline
                is_firstline = false
                continue
            end
            (ind, father, mother, sex) = split(line)
            ind = parse(Int64, ind)
            father = parse(Int64, father)
            mother = parse(Int64, mother)
            sex = parse(Int64, sex)
            pedigree[ind] = IntIndividual(
                ind,
                father,
                mother,
                sex,
                -1
            )
        end
    end
    if sort
        pedigree = _ordered_pedigree(pedigree)
    end
    _finalize_pedigree(pedigree)
end

"""
    _max_depth!(individual::IntIndividual, pedigree::Pedigree{IntIndividual})

Set and return the maximum depth of an individual's pedigree.
"""
function _max_depth!(individual::IntIndividual, pedigree::Pedigree{IntIndividual})
    if individual.max_depth == -1
        father_depth = 0
        mother_depth = 0
        if individual.father != 0
            father_depth += _max_depth!(pedigree[individual.father], pedigree)
        end
        if individual.mother != 0
            mother_depth += _max_depth!(pedigree[individual.mother], pedigree)
        end
        individual.max_depth = max(father_depth, mother_depth) + 1
    end
    individual.max_depth
end

"""
    _ordered_pedigree(pedigree::Pedigree)

Return a reordered pedigree where the individuals are in chronological order,
i.e. any individual's parents appear before them.
"""
function _ordered_pedigree(pedigree::Pedigree{IntIndividual})
    IDs = collect(keys(pedigree))
    depths = [_max_depth!(pedigree[ID], pedigree) for ID ∈ IDs]
    order = sortperm(depths)
    sortedIDs = IDs[order]
    ordered_pedigree = Pedigree{IntIndividual}()
    for ID ∈ sortedIDs
        ordered_pedigree[ID] = pedigree[ID]
    end
    ordered_pedigree
end

"""
    _finalize_pedigree(pedigree::Pedigree)

Return a pedigree of immutable individuals.
"""
function _finalize_pedigree(pedigree::Pedigree)
    temporary_pedigree = copy(pedigree)
    pedigree = Pedigree{Individual}()
    for (rank, individual) ∈ enumerate(values(temporary_pedigree))
        pedigree[individual.ID] = Individual(
            individual.ID,
            individual.father == 0 ? nothing : pedigree[individual.father],
            individual.mother == 0 ? nothing : pedigree[individual.mother],
            Individual[],
            individual.sex,
            rank
        )
        if individual.father != 0
            push!(pedigree[individual.father].children, pedigree[individual.ID])
        end
        if individual.mother != 0
            push!(pedigree[individual.mother].children, pedigree[individual.ID])
        end
    end
    pedigree
end

"""
    _save(path::String, pedigree::Pedigree, sorted::Bool = false)

Export the pedigree as a CSV file at a given `path`.

If `sorted` is `false` (the default), then the individuals
will appear in the same order as in the genealogy.

If `sorted` is `true`, then the individuals
will appear in alphabetical ID order.
"""
function _save(path::String, pedigree::Pedigree; sorted::Bool = false)
    df = genout(pedigree, sorted = sorted)
    open(path, "w") do file
        firstline = "ind\tfather\tmother\tsex"
        println(file, firstline)
        for row ∈ eachrow(df)
            line = join(row, "\t")
            println(file, line)
        end
    end
end