"""
    abstract type AbstractIndividual end

The abstract type of `Individual` and its subtypes.
"""
abstract type AbstractIndividual end

"""
    mutable struct MutableIndividual <: AbstractIndividual
        ID::Int64
        father::Union{Nothing, Individual}
        mother::Union{Nothing, Individual}
        children::Vector{Individual}
        sex::Int64
        rank::Int64
    end

The temporary unit structure of a [`Pedigree`](@ref).
"""
mutable struct MutableIndividual <: AbstractIndividual
    ID::Int64
    father::Union{Nothing, MutableIndividual}
    mother::Union{Nothing, MutableIndividual}
    children::Vector{MutableIndividual}
    sex::Int64
    rank::Int64
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

The unit structure of a [`Pedigree`](@ref).
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
    const Pedigree{T} = OrderedDict{Int64, T} where T <: AbstractIndividual

A particular case of an `OrderedDict` containing individuals accessed by ID.
"""
const Pedigree{T} = OrderedDict{Int64, T} where T <: AbstractIndividual

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
    pedigree = Pedigree{MutableIndividual}()
    for (rank, row) in enumerate(eachrow(dataframe))
        pedigree[row.ind] = MutableIndividual(
            row.ind,
            nothing,
            nothing,
            Int64[],
            row.sex,
            rank
        )
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
    if sort
        _order_pedigree!(pedigree)
    end
    _immutable_struct!(pedigree)
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
    pedigree = Pedigree{MutableIndividual}()
    open(filename) do file
        rank = 0
        while !eof(file)
            line = readline(file)
            rank += 1
            if rank == 1
                continue
            end
            (ind, father, mother, sex) = split(line)
            ind = parse(Int64, ind)
            father = parse(Int64, father)
            mother = parse(Int64, mother)
            sex = parse(Int64, sex)
            pedigree[ind] = MutableIndividual(
            ind,
            nothing,
            nothing,
            Int64[],
            sex,
            rank
        )
        end
    end
    open(filename) do file
        rank = 0
        while !eof(file)
            line = readline(file)
            rank += 1
            if rank == 1
                continue
            end
            (ind, father, mother, _) = split(line)
            ind = parse(Int64, ind)
            father = parse(Int64, father)
            if father != 0
                pedigree[ind].father = pedigree[father]
                push!(pedigree[father].children, pedigree[ind])
            end
            mother = parse(Int64, mother)
            if mother != 0
                pedigree[ind].mother = pedigree[mother]
                push!(pedigree[mother].children, pedigree[ind])
            end
        end
    end
    if sort
        _order_pedigree!(pedigree)
    end
    pedigree
    _immutable_struct!(pedigree)
end

"""
    _order_pedigree(pedigree::Pedigree)

Return a reordered pedigree where the individuals are in chronological order,
i.e. any individual's parents appear before them.
"""
function _order_pedigree!(pedigree::Pedigree)
    individuals = [individual for (_, individual) in pedigree]
    depths = [max_depth(individual) for individual in individuals]
    order = sortperm(depths)
    sorted_individuals = individuals[order]
    empty!(pedigree)
    for (rank, individual) in enumerate(sorted_individuals)
        individual.rank = rank
        pedigree[individual.ID] = individual
    end
    pedigree
end

"""
"""
function _immutable_struct!(pedigree::Pedigree)
    temporary_pedigree = copy(pedigree)
    pedigree = Pedigree{Individual}()
    for individual in collect(values(temporary_pedigree))
        pedigree[individual.ID] = Individual(
            individual.ID,
            isnothing(individual.father) ? nothing : pedigree[individual.father.ID],
            isnothing(individual.mother) ? nothing : pedigree[individual.mother.ID],
            Individual[],
            individual.sex,
            individual.rank
        )
        if !isnothing(individual.father)
            push!(pedigree[individual.father.ID].children, pedigree[individual.ID])
        end
        if !isnothing(individual.mother)
            push!(pedigree[individual.mother.ID].children, pedigree[individual.ID])
        end
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
function save_genealogy(pedigree::Pedigree, path::String; sorted::Bool = false)
    df = genout(pedigree, sorted = sorted)
    open(path, "w") do file
        firstline = "ind\tfather\tmother\tsex"
        println(file, firstline)
        for row in eachrow(df)
            line = join(row, "\t")
            println(file, line)
        end
    end
end