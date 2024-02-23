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
    end

The temporary unit structure of a [`Pedigree`](@ref).
"""
struct IntIndividual <: AbstractIndividual
    ID::Int64
    father::Int64
    mother::Int64
    sex::Int64
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
    pedigree = Pedigree{IntIndividual}()
    for row in eachrow(dataframe)
        pedigree[row.ind] = IntIndividual(
            row.ind,
            row.father,
            row.mother,
            row.sex,
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
                sex
            )
        end
    end
    if sort
        pedigree = _ordered_pedigree(pedigree)
    end
    _immutable_struct(pedigree)
end

"""
    _max_depth(pedigree::Pedigree{IntIndividual}, individual::Int64)

Return the maximum depth of an individual's pedigree.
"""
function _max_depth(pedigree::Pedigree{IntIndividual}, individual::IntIndividual)
    depth = 1
    father = individual.father
    mother = individual.mother
    if father != 0 && mother != 0
        depth += max(_max_depth(pedigree, pedigree[father]), _max_depth(pedigree, pedigree[mother]))
    elseif father != 0
        depth += _max_depth(pedigree, pedigree[father])
    elseif mother != 0
        depth += _max_depth(pedigree, pedigree[mother])
    end
    depth
end

"""
    _ordered_pedigree(pedigree::Pedigree)

Return a reordered pedigree where the individuals are in chronological order,
i.e. any individual's parents appear before them.
"""
function _ordered_pedigree(pedigree::Pedigree{IntIndividual})
    IDs = collect(keys(pedigree))
    depths = [_max_depth(pedigree, pedigree[ID]) for ID in IDs]
    order = sortperm(depths)
    sortedIDs = IDs[order]
    ordered_pedigree = Pedigree{IntIndividual}()
    for ID in sortedIDs
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
    for (rank, individual) in enumerate(collect(values(temporary_pedigree)))
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