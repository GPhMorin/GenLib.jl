"""
    phi(individualᵢ::Individual, individualⱼ::Individual, Ψ::Union{Nothing, Matrix{Float64}} = nothing)

Return the kinship coefficient between two individuals.

A matrix of the founders' kinships may optionally be provided.
"""
function phi(individualᵢ::Individual, individualⱼ::Individual, Ψ::Union{Nothing, Matrix{Float64}} = nothing)
    # Adapted from GENLIB and Kirkpatrick et al., 2019.
    if (individualᵢ.sort != 0) && (individualⱼ.sort != 0) # They are both founders
        return Ψ[individualᵢ.sort, individualⱼ.sort]
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
    phi(pedigree::OrderedDict{Int64, Individual}, Ψ::Matrix{Float64})

Return a square matrix of the pairwise kinship coefficients between all probands
provided a matrix of the founders' kinships.
"""
function phi(pedigree::OrderedDict{Int64, Individual}, Ψ::Matrix{Float64})
    probandIDs = filter(x -> isempty(pedigree[x].children), collect(keys(pedigree)))
    probands = [pedigree[ID] for ID in probandIDs]
    founderIDs = filter(x -> isnothing(pedigree[x].father) && isnothing(pedigree[x].mother), collect(keys(pedigree)))
    founders = [pedigree[ID] for ID in founderIDs]
    if probandIDs == founderIDs
        return Ψ # The founders' kinships
    else
        Φ = Matrix{Float64}(undef, length(probandIDs), length(probandIDs))
        for f in eachindex(founders)
            founders[f].sort = f # To later access values in the Ψ matrix
        end
        Threads.@threads for i in eachindex(probandIDs)
            Threads.@threads for j in eachindex(probandIDs)
                if i == j # Calculate only once
                    individualᵢ = probands[i]
                    Φ[i, i] = phi(individualᵢ, individualᵢ, Ψ)
                elseif i < j # Calculate only once
                    individualᵢ = probands[i]
                    individualⱼ = probands[j]
                    Φ[i, j] = Φ[j, i] = phi(individualᵢ, individualⱼ, Ψ)
                end
            end
        end
        return Φ
    end
end

"""
    phi(pedigree::OrderedDict{Int64, Individual}, pro::Vector{Int64} = pro(pedigree); verbose::Bool = false)

Return a square matrix of the pairwise kinship coefficients
between the provided probands.

If no probands are provided, kinships for all of the pedigree's probands are computed.
"""
function phi(pedigree::OrderedDict{Int64, Individual}, probandIDs::Vector{Int64} = pro(pedigree); verbose::Bool = false)
    founderIDs = founder(pedigree)
    Ψ = zeros(length(founderIDs), length(founderIDs)) # Initialize the top founders' kinships
    for f in eachindex(founderIDs)
        Ψ[f, f] = 0.5
    end
    isolated_pedigree = branching(pedigree, pro = probandIDs)
    upperIDs = cut_vertices(isolated_pedigree)
    lowerIDs = probandIDs
    levels = [upperIDs, lowerIDs]
    while upperIDs != lowerIDs
        # Cut the pedigree into several sub-pedigrees
        isolated_pedigree = branching(isolated_pedigree, pro = upperIDs)
        lowerIDs = copy(upperIDs)
        upperIDs = cut_vertices(isolated_pedigree)
        pushfirst!(levels, upperIDs)
    end
    for i in 1:length(levels)-1
        # For each sub-pedigree, calculate the kinships using the previous founders' kinships
        upperIDs = levels[i]
        lowerIDs = levels[i+1]
        Vᵢ = branching(pedigree, pro = lowerIDs, ancestors = upperIDs)
        Ψ = phi(Vᵢ, Ψ)
        if verbose
            println("Kinships for segment ", i, "/", length(levels)-1, " completed.")
        end
    end
    # In GENLIB, the kinship matrix is in alphabetical ID order
    probandIDs = filter(x -> x ∈ probandIDs, collect(keys(pedigree)))
    order = sortperm(probandIDs)
    Ψ[order, order]
end

"""
    phi(pedigree::OrderedDict{Int64, Individual}, rowIDs::Vector{Int64}, columnIDs::Vector{Int64})

Return a rectangle matrix of kinship coefficients between row IDs and column IDs.
The kinship of someone with themself is replaced with their inbreeding.
"""
function phi(pedigree::OrderedDict{Int64, Individual}, rowIDs::Vector{Int64}, columnIDs::Vector{Int64})
    Φ = zeros(length(rowIDs), length(columnIDs)) # Initialize the kinship matrix
    Threads.@threads for i in eachindex(rowIDs)
        Threads.@threads for j in eachindex(columnIDs)
            IDᵢ = rowIDs[i]
            IDⱼ = columnIDs[j]
            individualᵢ = pedigree[IDᵢ]
            individualⱼ = pedigree[IDⱼ]
            Φ[i, j] = phi(individualᵢ, individualⱼ)
        end
    end
    Φ
end

"""
    cut_vertex(individual::Individual, candidateID::Int64)

Return whether an individual can be used as a cut vertex.

A cut vertex is an individual that "when removed,
disrupt every path from any source [founder]
to any sink [proband]" (Kirkpatrick et al., 2019).
"""
function cut_vertex(individual::Individual, candidateID::Int64)
    value = true
    for child in individual.children
        # Check if going down the pedigree
        # while avoiding the candidate ID
        # reaches a proband (sink) anyway
        if child.state == PROBAND
            return false
        elseif child.ID != candidateID
            value = value && cut_vertex(child, candidateID)
        end
    end
    value
end

"""
    cut_vertices(pedigree::OrderedDict{Int64, Individual})

Return the IDs of the cut vertices as defined in Kirkpatrick et al., 2019.

A cut vertex is an individual that "when removed,
disrupt every path from any source [founder]
to any sink [proband]" (Kirkpatrick et al., 2019).
"""
function cut_vertices(pedigree::OrderedDict{Int64, Individual})
    vertices = Int64[]
    probandIDs = pro(pedigree)
    probands = [pedigree[ID] for ID in probandIDs]
    for proband in probands # Mark the probands
        proband.state = PROBAND
    end
    founderIDs = founder(pedigree)
    founders = [pedigree[ID] for ID in founderIDs]
    for founder in founders # Mark the founders
        founder.state = FOUNDER
    end
    stack = probands
    while !isempty(stack)
        candidate = pop!(stack)
        # Check if avoiding paths from a "source" (ancestor)
        # down the pedigree through the candidate ID
        # never reaches a "sink" (proband) individual
        sourceIDs = ancestor(pedigree, candidate.ID)
        is_candidate = true
        for sourceID in sourceIDs
            if !cut_vertex(pedigree[sourceID], candidate.ID)
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
    for proband in probands # Unmark the probands
        proband.state = UNEXPLORED
    end
    for founder in founders # Unmark the founders
        founder.state = UNEXPLORED
    end
    sort(collect(Set(vertices)))
end

"""
    set_ancestors(pedigree::OrderedDict{Int64, Individual})

Return a lookup table of the individuals' ancestors.

As defined in Kirkpatrick et al., 2019.
"""
function set_ancestors(pedigree::OrderedDict{Int64, Individual})
    ancestors = Dict()
    for (ID, individual) in pedigree
        ancestors[ID] = [ID]
        mother = individual.mother
        father = individual.father
        if !isnothing(mother)
            ancestors[ID] = union(ancestors[mother.ID], ancestors[ID])
        end
        if !isnothing(father)
            ancestors[ID] = union(ancestors[father.ID], ancestors[ID])
        end
    end
    ancestors
end

"""
    ϕ(pedigree::OrderedDict{Int64, Individual}, Ψ::Matrix{Float64})

Return a square matrix of pairwise kinship coefficients
between all probands given the founders' kinships.

An implementation of the recursive-cut algorithm presented in Kirkpatrick et al., 2019.

The diagonal corresponds to inbreedings.
"""
function ϕ(pedigree::OrderedDict{Int64, Individual}, Ψ::Matrix{Float64})
    probandIDs = filter(x -> isempty(pedigree[x].children), collect(keys(pedigree)))
    probands = [pedigree[ID] for ID in probandIDs]
    founderIDs = filter(x -> isnothing(pedigree[x].father) && isnothing(pedigree[x].mother), collect(keys(pedigree)))
    founders = [pedigree[ID] for ID in founderIDs]
    if probandIDs == founderIDs
        return Ψ
    end
    n = length(pedigree)
    Φ = Matrix{Float64}(undef, n, n)
    for (f, founder) in enumerate(founders)
        Φ[founder.index, founder.index] = (1 + Ψ[f, f]) / 2
    end
    for (f, founder₁) in enumerate(founders), (g, founder₂) in enumerate(founders)
        if founder₁.ID < founder₂.ID
            Φ[founder₁.index, founder₂.index] = Φ[founder₂.index, founder₁.index] = Ψ[f, g]
        end
    end
    for (_, individualᵢ) in pedigree
        i = individualᵢ.index
        for (_, individualⱼ) in pedigree
            j = individualⱼ.index
            coefficient = 0.
            if i > j # i cannot be an ancestor of j
                father = individualᵢ.father
                if !isnothing(father)
                    coefficient += Φ[father.index, individualⱼ.index] / 2
                end
                mother = individualᵢ.mother
                if !isnothing(mother)
                    coefficient += Φ[mother.index, individualⱼ.index] / 2
                end
            elseif j > i # j cannot be an ancestor of i
                father = individualⱼ.father
                if !isnothing(father)
                    coefficient += Φ[father.index, individualᵢ.index] / 2
                end
                mother = individualⱼ.mother
                if !isnothing(mother)
                    coefficient += Φ[mother.index, individualᵢ.index] / 2
                end
            else # i.index == j.index, same individual
                coefficient += 0.5
                father = individualᵢ.father
                mother = individualᵢ.mother
                if !isnothing(father) && !isnothing(mother)
                    coefficient += Φ[mother.index, father.index] / 2
                end
            end
            Φ[i, j] = Φ[j, i] = coefficient
        end
    end
    for i in 1:n
        Φ[i, i] = (2 * Φ[i, i]) - 1
    end
    indices = [proband.index for proband in probands]
    Φ[indices, indices]
end

"""
    ϕ(pedigree::OrderedDict{Int64, Individual}, probandIDs::Vector{Int64} = pro(pedigree); verbose::Bool = false)

Return the square matrix of the pairwise kinship coefficients of a set of probands.

If no probands are given, return the square matrix for all probands in the pedigree.

An implementation of the recursive-cut algorithm presented in Kirkpatrick et al., 2019.

The diagonal corresponds to inbreedings.
"""
function ϕ(pedigree::OrderedDict{Int64, Individual}, probandIDs::Vector{Int64} = pro(pedigree); verbose::Bool = false)
    founderIDs = founder(pedigree)
    Ψ = zeros(length(founderIDs), length(founderIDs))
    isolated_pedigree = branching(pedigree, pro = probandIDs)
    upperIDs = cut_vertices(isolated_pedigree)
    lowerIDs = probandIDs
    C = [upperIDs, lowerIDs]
    while upperIDs != lowerIDs
        isolated_pedigree = branching(pedigree, pro = upperIDs)
        lowerIDs = copy(upperIDs)
        upperIDs = cut_vertices(isolated_pedigree)
        pushfirst!(C, upperIDs)
    end
    if verbose
        for i in 1:length(C)-1
            upperIDs = C[i]
            lowerIDs = C[i+1]
            Vᵢ = branching(pedigree, pro = lowerIDs, ancestors = upperIDs)
            println("Segment ", i, "/", length(C)-1, " contains ", length(Vᵢ),
                    " individuals and ", length(upperIDs), " founders.")
        end
    end
    for i in 1:length(C)-1
        upperIDs = C[i]
        lowerIDs = C[i+1]
        Vᵢ = branching(pedigree, pro = lowerIDs, ancestors = upperIDs)
        Ψ = ϕ(Vᵢ, Ψ)
        if verbose
            println("Kinships for segment ", i, "/", length(C)-1, " completed.")
        end
    end
    probandIDs = filter(x -> x ∈ probandIDs, collect(keys(pedigree)))
    order = sortperm(probandIDs)
    Ψ[order, order]
end