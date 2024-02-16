"""
    phi(individualᵢ::Individual, individualⱼ::Individual, Ψ::Union{Nothing, Matrix{Float64}} = nothing)

Return the kinship coefficient between two individuals.

A matrix of the founders' kinships may optionally be provided.

Adapted from [Karigl, 1981](@ref), and [Kirkpatrick et al., 2019](@ref).

# Example

```julia
import GenLib as gen
genea140 = gen.genea140
ped = gen.genealogy(genea140)
pro1 = ped[10033]
pro2 = ped[113470]
gen.phi(pro1, pro2)
```
"""
function phi(individualᵢ::Individual, individualⱼ::Individual, Ψ::Union{Nothing, Matrix{Float64}} = nothing)
    if !isempty(individualᵢ.stats) && !isempty(individualⱼ.stats) # They are both founders
        return Ψ[individualᵢ.stats["index"], individualⱼ.stats["index"]]
    else # At least one of the individuals is not a founder
        value = 0.
        if individualᵢ.index > individualⱼ.index # From the genealogical order, i cannot be an ancestor of j
            # Φᵢⱼ = (Φₚⱼ + Φₘⱼ) / 2, if i is not an ancestor of j (Karigl, 1981)
            if !isnothing(individualᵢ.father)
                value += phi(individualᵢ.father, individualⱼ, Ψ) / 2
            end
            if !isnothing(individualᵢ.mother)
                value += phi(individualᵢ.mother, individualⱼ, Ψ) / 2
            end
        elseif individualⱼ.index > individualᵢ.index # Reverse the order since a > b
            # Φⱼᵢ = (Φₚⱼ + Φₘⱼ) / 2, if j is not an ancestor of i (Karigl, 1981)
            if !isnothing(individualⱼ.father)
                value += phi(individualⱼ.father, individualᵢ, Ψ) / 2
            end
            if !isnothing(individualⱼ.mother)
                value += phi(individualⱼ.mother, individualᵢ, Ψ) / 2
            end
        elseif individualᵢ.index == individualⱼ.index # Same individual
            # Φₐₐ = (1 + Φₚₘ) / 2 (Karigl, 1981)
            value += 1/2
            if !isnothing(individualᵢ.father) & !isnothing(individualᵢ.mother)
                value += phi(individualᵢ.father, individualᵢ.mother, Ψ) / 2
            end
        end
        return value
    end
end

"""
    _cut_vertex(individual::Individual, candidateID::Int64)

Return whether an individual can be used as a cut vertex.

A cut vertex is an individual that "when removed,
disrupt every path from any source [founder]
to any sink [proband]" ([Kirkpatrick et al., 2019](@ref)).
"""
function _cut_vertex(individual::Individual, candidateID::Int64, probandIDs::Vector{Int64})
    value = true
    for child in individual.children
        # Check if going down the pedigree
        # while avoiding the candidate ID
        # reaches a proband (sink) anyway
        if child.ID ∈ probandIDs
            return false
        elseif child.ID != candidateID
            value = value && _cut_vertex(child, candidateID, probandIDs)
        end
    end
    value
end

"""
    _cut_vertices(pedigree::Pedigree)

Return the IDs of the cut vertices as defined in [Kirkpatrick et al., 2019](@ref).

A cut vertex is an individual that "when removed,
disrupt every path from any source [founder]
to any sink [proband]" ([Kirkpatrick et al., 2019](@ref)).
"""
function _cut_vertices(pedigree::Pedigree)
    vertices = Int64[]
    probandIDs = pro(pedigree)
    probands = [pedigree[ID] for ID in probandIDs]
    stack = Stack{Individual}()
    push!(stack, probands...)
    while !isempty(stack)
        candidate = pop!(stack)
        # Check if avoiding paths from a "source" (ancestor)
        # down the pedigree through the candidate ID
        # never reaches a "sink" (proband) individual
        sourceIDs = ancestor(pedigree, candidate.ID)
        is_candidate = true
        for sourceID in sourceIDs
            if !_cut_vertex(pedigree[sourceID], candidate.ID, probandIDs)
                is_candidate = false
                break
            end
        end
        if is_candidate
            push!(vertices, candidate.ID)
        else
            if !isnothing(candidate.father)
                push!(stack, candidate.father)
            end
            if !isnothing(candidate.mother)
                push!(stack, candidate.mother)
            end
        end
    end
    sort(collect(Set(vertices)))
end

"""
    phi(pedigree::Pedigree, Ψ::Matrix{Float64})

Return a square matrix of pairwise kinship coefficients
between all probands given the founders' kinships.

An implementation of the recursive-cut algorithm presented in [Kirkpatrick et al., 2019](@ref).
"""
function phi(pedigree::Pedigree, Ψ::Matrix{Float64})
    probandIDs = filter(x -> isempty(pedigree[x].children), collect(keys(pedigree)))
    probands = [pedigree[ID] for ID in probandIDs]
    founderIDs = filter(x -> isnothing(pedigree[x].father) && isnothing(pedigree[x].mother), collect(keys(pedigree)))
    founders = [pedigree[ID] for ID in founderIDs]
    if probandIDs == founderIDs
        return Ψ
    end
    Φ = zeros(length(pedigree), length(pedigree))
    for (index₁, founder₁) in enumerate(founders)
        for (index₂, founder₂) in enumerate(founders)
            if index₁ ≤ index₂
                Φ[founder₁.index, founder₂.index] = Φ[founder₂.index, founder₁.index] = Ψ[index₁, index₂]
            end
        end
    end
    for individualᵢ in values(pedigree)
        i = individualᵢ.index
        for individualⱼ in values(pedigree)
            j = individualⱼ.index
            if Φ[i, j] > 0
                continue
            elseif i > j # i cannot be an ancestor of j
                father = individualᵢ.father
                mother = individualᵢ.mother
                if !isnothing(father)
                    Φ[i, j] += Φ[father.index, individualⱼ.index] / 2
                    Φ[j, i] += Φ[father.index, individualⱼ.index] / 2
                end
                if !isnothing(mother)
                    Φ[i, j] += Φ[mother.index, individualⱼ.index] / 2
                    Φ[j, i] += Φ[mother.index, individualⱼ.index] / 2
                end
            elseif j > i # j cannot be an ancestor of i
                father = individualⱼ.father
                mother = individualⱼ.mother
                if !isnothing(father)
                    Φ[i, j] += Φ[father.index, individualᵢ.index] / 2
                    Φ[j, i] += Φ[father.index, individualᵢ.index] / 2
                end
                if !isnothing(mother)
                    Φ[i, j] += Φ[mother.index, individualᵢ.index] / 2
                    Φ[j, i] += Φ[mother.index, individualᵢ.index] / 2
                end
            else # i == j, same individual
                Φ[i, i] += 0.5
                father = individualᵢ.father
                mother = individualᵢ.mother
                if !isnothing(father) && !isnothing(mother)
                    Φ[i, i] += Φ[mother.index, father.index] / 2
                end
            end
        end
    end
    indices = [proband.index for proband in probands]
    Φ[indices, indices]
end

"""
    phi(pedigree::Pedigree, probandIDs::Vector{Int64} = pro(pedigree); verbose::Bool = false)

Return the square matrix of the pairwise kinship coefficients of a set of probands.

If no probands are given, return the square matrix for all probands in the pedigree.

An implementation of the recursive-cut algorithm presented in [Kirkpatrick et al., 2019](@ref).

# Example

```julia
import GenLib as gen
geneaJi = gen.geneaJi
ped = gen.genealogy(geneaJi)
gen.phi(ped)
```
"""
function phi(pedigree::Pedigree, probandIDs::Vector{Int64} = pro(pedigree); verbose::Bool = false, estimate::Bool = false)
    founderIDs = founder(pedigree)
    Ψ = zeros(length(founderIDs), length(founderIDs))
    for f in eachindex(founderIDs)
        Ψ[f, f] = 0.5
    end
    isolated_pedigree = branching(pedigree, pro = probandIDs)
    upperIDs = _cut_vertices(isolated_pedigree)
    lowerIDs = probandIDs
    C = [upperIDs, lowerIDs]
    while upperIDs != lowerIDs
        isolated_pedigree = branching(pedigree, pro = upperIDs)
        lowerIDs = copy(upperIDs)
        upperIDs = _cut_vertices(isolated_pedigree)
        pushfirst!(C, upperIDs)
    end
    if verbose || estimate
        for i in 1:length(C)-1
            upperIDs = C[i]
            lowerIDs = C[i+1]
            Vᵢ = branching(pedigree, pro = lowerIDs, ancestors = upperIDs)
            println("Segment ", i, "/", length(C)-1, " contains ", length(Vᵢ),
                    " individuals and ", length(upperIDs), " founders.")
        end
        if estimate
            return
        end
    end
    for i in 1:length(C)-1
        upperIDs = C[i]
        lowerIDs = C[i+1]
        Vᵢ = branching(pedigree, pro = lowerIDs, ancestors = upperIDs)
        Ψ = phi(Vᵢ, Ψ)
        if verbose
            println("Kinships for segment ", i, "/", length(C)-1, " completed.")
        end
    end
    probandIDs = filter(x -> x ∈ probandIDs, collect(keys(pedigree)))
    order = sortperm(probandIDs)
    Ψ[order, order]
end