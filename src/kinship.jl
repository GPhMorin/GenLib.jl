"""
    struct PossibleFounder <: AbstractIndividual
        ID::Int64
        father::Union{Nothing, PossibleFounder}
        mother::Union{Nothing, PossibleFounder}
        children::Vector{PossibleFounder}
        sex::Int64
        rank::Int64
        index::Int64
    end

An individual with an index to access the founder's kinships in the Ψ matrix.
"""
struct PossibleFounder <: AbstractIndividual
    ID::Int64
    father::Union{Nothing, PossibleFounder}
    mother::Union{Nothing, PossibleFounder}
    children::Vector{PossibleFounder}
    sex::Int64
    rank::Int64
    index::Int64
end

"""
    phi(individualᵢ::T, individualⱼ::T) where T <: AbstractIndividual

Return the kinship coefficient between two individuals.

Adapted from [Karigl, 1981](@ref).

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
function phi(individualᵢ::T, individualⱼ::T) where T <: AbstractIndividual
    value = 0.
    if individualᵢ.rank > individualⱼ.rank # From the genealogical order, i cannot be an ancestor of j
        # Φᵢⱼ = (Φₚⱼ + Φₘⱼ) / 2, if i is not an ancestor of j (Karigl, 1981)
        if !isnothing(individualᵢ.father)
            value += phi(individualᵢ.father, individualⱼ) / 2
        end
        if !isnothing(individualᵢ.mother)
            value += phi(individualᵢ.mother, individualⱼ) / 2
        end
    elseif individualⱼ.rank > individualᵢ.rank # Reverse the order since a > b
        # Φⱼᵢ = (Φₚⱼ + Φₘⱼ) / 2, if j is not an ancestor of i (Karigl, 1981)
        if !isnothing(individualⱼ.father)
            value += phi(individualⱼ.father, individualᵢ) / 2
        end
        if !isnothing(individualⱼ.mother)
            value += phi(individualⱼ.mother, individualᵢ) / 2
        end
    elseif individualᵢ.rank == individualⱼ.rank # Same individual
        # Φₐₐ = (1 + Φₚₘ) / 2 (Karigl, 1981)
        value += 1/2
        if !isnothing(individualᵢ.father) & !isnothing(individualᵢ.mother)
            value += phi(individualᵢ.father, individualᵢ.mother) / 2
        end
    end
    return value
end

"""
    phi(individualᵢ::PossibleFounder, individualⱼ::PossibleFounder, Ψ::Matrix{Float64})

Return the kinship coefficient between two individuals given a matrix of the founders' kinships.

Adapted from [Karigl, 1981](@ref), and [Kirkpatrick et al., 2019](@ref).
"""
function phi(individualᵢ::PossibleFounder, individualⱼ::PossibleFounder, Ψ::Matrix{Float64})
    if individualᵢ.index != 0 && individualⱼ.index != 0 # They are both founders
        return Ψ[individualᵢ.index, individualⱼ.index]
    end
    value = 0.
    if individualᵢ.rank > individualⱼ.rank # From the genealogical order, i cannot be an ancestor of j
        # Φᵢⱼ = (Φₚⱼ + Φₘⱼ) / 2, if i is not an ancestor of j (Karigl, 1981)
        if !isnothing(individualᵢ.father)
            value += phi(individualᵢ.father, individualⱼ, Ψ) / 2
        end
        if !isnothing(individualᵢ.mother)
            value += phi(individualᵢ.mother, individualⱼ, Ψ) / 2
        end
    elseif individualⱼ.rank > individualᵢ.rank # Reverse the order since a > b
        # Φⱼᵢ = (Φₚⱼ + Φₘⱼ) / 2, if j is not an ancestor of i (Karigl, 1981)
        if !isnothing(individualⱼ.father)
            value += phi(individualⱼ.father, individualᵢ, Ψ) / 2
        end
        if !isnothing(individualⱼ.mother)
            value += phi(individualⱼ.mother, individualᵢ, Ψ) / 2
        end
    elseif individualᵢ.rank == individualⱼ.rank # Same individual
        # Φₐₐ = (1 + Φₚₘ) / 2 (Karigl, 1981)
        value += 1/2
        if !isnothing(individualᵢ.father) & !isnothing(individualᵢ.mother)
            value += phi(individualᵢ.father, individualᵢ.mother, Ψ) / 2
        end
    end
    return value
end

"""
    _cut_vertex(individual::T, candidateID::Int64) where T <: AbstractIndividual

Return whether an individual can be used as a cut vertex.

A cut vertex is an individual that "when removed,
disrupt every path from any source [founder]
to any sink [proband]" ([Kirkpatrick et al., 2019](@ref)).
"""
function _cut_vertex(individual::T, candidateID::Int64) where T <: AbstractIndividual
    value = true
    for child in individual.children
        # Check if going down the pedigree
        # while avoiding the candidate ID
        # reaches a proband (sink) anyway
        if isempty(child.children) # The child is a proband
            return false
        elseif child.ID != candidateID
            value = value && _cut_vertex(child, candidateID)
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
function _cut_vertices(pedigree::Pedigree{T}) where T <: AbstractIndividual
    vertices = Int64[]
    probandIDs = pro(pedigree)
    probands = [pedigree[ID] for ID in probandIDs]
    stack = T[]
    push!(stack, probands...)
    while !isempty(stack)
        candidate = pop!(stack)
        # Check if avoiding paths from a "source" (founder)
        # down the pedigree through the candidate ID
        # never reaches a "sink" (proband) individual
        sourceIDs = findFounders(pedigree, [candidate.ID])
        if all(_cut_vertex(pedigree[sourceID], candidate.ID) for sourceID in sourceIDs)
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
                coefficient = Ψ[index₁, index₂]
                Φ[founder₁.rank, founder₂.rank] = coefficient
                Φ[founder₂.rank, founder₁.rank] = coefficient
            end
        end
    end
    for individualᵢ in values(pedigree)
        i = individualᵢ.rank
        for individualⱼ in values(pedigree)
            j = individualⱼ.rank
            if Φ[i, j] > 0
                continue
            elseif i > j # i cannot be an ancestor of j
                father = individualᵢ.father
                mother = individualᵢ.mother
                if !isnothing(father)
                    coefficient = Φ[father.rank, individualⱼ.rank] / 2
                    Φ[i, j] += coefficient
                    Φ[j, i] += coefficient
                end
                if !isnothing(mother)
                    coefficient = Φ[mother.rank, individualⱼ.rank] / 2
                    Φ[i, j] += coefficient
                    Φ[j, i] += coefficient
                end
            elseif j > i # j cannot be an ancestor of i
                father = individualⱼ.father
                mother = individualⱼ.mother
                if !isnothing(father)
                    coefficient = Φ[father.rank, individualᵢ.rank] / 2
                    Φ[i, j] += coefficient
                    Φ[j, i] += coefficient
                end
                if !isnothing(mother)
                    coefficient = Φ[mother.rank, individualᵢ.rank] / 2
                    Φ[i, j] += coefficient
                    Φ[j, i] += coefficient
                end
            else # i == j, same individual
                father = individualᵢ.father
                mother = individualᵢ.mother
                if !isnothing(father) && !isnothing(mother)
                    Φ[i, i] += (1 + Φ[mother.rank, father.rank]) / 2
                else
                    Φ[i, i] += 0.5
                end
            end
        end
    end
    indices = [proband.rank for proband in probands]
    Φ[indices, indices]
end

"""
    phi(pedigree::Pedigree, probandIDs::Vector{Int64} = pro(pedigree); verbose::Bool = false)

Return a square matrix of pairwise kinship coefficients between probands.

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
function phi(pedigree::Pedigree, probandIDs::Vector{Int64} = pro(pedigree);
    verbose::Bool = false, estimate::Bool = false, MT::Bool = false)
    founderIDs = founder(pedigree)
    Ψ = zeros(length(founderIDs), length(founderIDs))
    for f in eachindex(founderIDs)
        Ψ[f, f] = 0.5
    end
    if verbose || estimate
        println("Isolating the pedigree...")
    end
    if verbose || estimate
        println("Cutting vertices...")
    end
    isolated_pedigree = branching(pedigree, pro = probandIDs)
    upperIDs = _cut_vertices(isolated_pedigree)
    lowerIDs = probandIDs
    C = [upperIDs, lowerIDs]
    segment = 1
    while upperIDs != lowerIDs
        segment += 1
        if verbose || estimate
            println("Vertices cut for segment ", segment, ".")
        end
        isolated_pedigree = branching(isolated_pedigree, pro = upperIDs)
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
        if MT
            index_pedigree = Pedigree{PossibleFounder}()
            index = 0
            for individual in collect(values(Vᵢ))
                father = individual.father
                mother = individual.mother
                if isnothing(father) && isnothing(mother)
                    index += 1
                    individual_index = copy(index)
                else
                    individual_index = 0
                end
                index_pedigree[individual.ID] = PossibleFounder(
                    individual.ID,
                    isnothing(father) ? nothing : index_pedigree[father.ID],
                    isnothing(mother) ? nothing : index_pedigree[mother.ID],
                    PossibleFounder[],
                    individual.sex,
                    individual.rank,
                    individual_index
                )
                if !isnothing(father)
                    push!(index_pedigree[father.ID].children, index_pedigree[individual.ID])
                end
                if !isnothing(mother)
                    push!(index_pedigree[mother.ID].children, index_pedigree[individual.ID])
                end
            end
            probands = filter(x -> isempty(x.children), collect(values(index_pedigree)))
            Φ = Matrix{Float64}(undef, length(probands), length(probands))
            Threads.@threads for i in eachindex(probands)
                Threads.@threads for j in eachindex(probands)
                    if i ≤ j
                        Φ[i, j] = Φ[j, i] = phi(probands[i], probands[j], Ψ)
                    end
                end
            end
            Ψ = copy(Φ)
        else
            Ψ = phi(Vᵢ, Ψ)
        end
        if verbose
            println("Kinships for segment ", i, "/", length(C)-1, " completed.")
        end
    end
    probandIDs = filter(x -> x ∈ probandIDs, collect(keys(pedigree)))
    order = sortperm(probandIDs)
    Ψ[order, order]
end

"""
    phi(pedigree::Pedigree, rowIDs::Vector{Int64}, columnIDs::Vector{Int64})

Return a rectangle matrix of kinship coefficients,
as defined by a list or row IDs and column IDs.
"""
function phi(pedigree::Pedigree, rowIDs::Vector{Int64}, columnIDs::Vector{Int64})
    Φ = Matrix{Float64}(undef, length(rowIDs), length(columnIDs))
    probandIDs = union(rowIDs, columnIDs)
    isolated_pedigree = branching(pedigree, pro = probandIDs)
    cut_vertices = _cut_vertices(isolated_pedigree)
    isolated_pedigree = branching(isolated_pedigree, ancestors = cut_vertices)
    Ψ = phi(pedigree, cut_vertices, MT = true)
    phi_pedigree = Pedigree{PossibleFounder}()
    founder_index = 1
    for (rank, individual) in enumerate(collect(values(isolated_pedigree)))
        father = individual.father
        mother = individual.mother
        if isnothing(father) && isnothing(mother)
            index = copy(founder_index)
            founder_index += 1
        else
            index = 0
        end
        phi_pedigree[individual.ID] = PossibleFounder(
            individual.ID,
            isnothing(father) ? nothing : phi_pedigree[father.ID],
            isnothing(mother) ? nothing : phi_pedigree[mother.ID],
            PossibleFounder[],
            individual.sex,
            rank,
            index
        )
        if !isnothing(father)
            push!(phi_pedigree[father.ID].children, phi_pedigree[individual.ID])
        end
        if !isnothing(mother)
            push!(phi_pedigree[mother.ID].children, phi_pedigree[individual.ID])
        end
    end
    Threads.@threads for i in eachindex(rowIDs)
        individualᵢ = phi_pedigree[rowIDs[i]]
        Threads.@threads for j in eachindex(columnIDs)
            individualⱼ = phi_pedigree[columnIDs[j]]
            Φ[i, j] = phi(individualᵢ, individualⱼ, Ψ)
        end
    end
    Φ
end