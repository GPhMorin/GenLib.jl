"""
    phi(individualᵢ::Individual, individualⱼ::Individual, Ψ::Union{Nothing, Matrix{Float64}} = nothing)

Return the kinship coefficient between two individuals.
The diagonal corresponds to inbreedings.

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
The diagonal corresponds to inbreedings.
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
The diagonal corresponds to inbreedings.

If no probands are provided, kinships for all of the pedigree's probands are computed.
"""
function phi(pedigree::OrderedDict{Int64, Individual}, probandIDs::Vector{Int64} = pro(pedigree); verbose::Bool = false)
    founderIDs = founder(pedigree)
    Ψ = zeros(length(founderIDs), length(founderIDs)) # Initialize the top founders' kinships
    for f in eachindex(founderIDs)
        Ψ[f, f] = 0.5
    end
    upperIDs = cut_vertices(pedigree)
    lowerIDs = probandIDs
    levels = [upperIDs, lowerIDs]
    while upperIDs != lowerIDs
        # Cut the pedigree into several sub-pedigrees
        isolated_pedigree = branching(pedigree, pro = upperIDs)
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
    for ID in probandIDs # Mark the probands
        pedigree[ID].state = PROBAND
    end
    founderIDs = founder(pedigree)
    candidateIDs = [ID for ID in collect(keys(pedigree))]
    for candidateID in candidateIDs
        # Check if avoiding paths from a "source" (founder)
        # down the pedigree through the candidate ID
        # never reaches a "sink" (proband) individual
        ancestorIDs = ancestor(pedigree, candidateID)
        sourceIDs = filter(x -> x ∈ founderIDs, ancestorIDs)
        candidate = true
        for sourceID in sourceIDs
            if !cut_vertex(pedigree[sourceID], candidateID)
                candidate = false
                break
            end
        end
        if candidate
            push!(vertices, candidateID)
        end
    end
    for ID in probandIDs # Unmark the probands
        pedigree[ID].state = UNEXPLORED
    end
    # Keep only the lowest candidates
    ancestorIDs = union([ancestor(pedigree, vertex) for vertex in vertices]...)
    vertices = filter(x -> (x ∈ vertices && x ∉ ancestorIDs) || (x ∈ founderIDs && x ∉ ancestorIDs), candidateIDs)
    vertices
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

An implementation of the recursive-cut algorithm presented in Kirkpatrick et al., 2019.

Return a square matrix of pairwise kinship coefficients
between all probands given the founders' kinships.
"""
function ϕ(pedigree::OrderedDict{Int64, Individual}, Ψ::Matrix{Float64})
    probandIDs = filter(x -> isempty(pedigree[x].children), collect(keys(pedigree)))
    probands = [pedigree[ID] for ID in probandIDs]
    founderIDs = filter(x -> isnothing(pedigree[x].father) && isnothing(pedigree[x].mother), collect(keys(pedigree)))
    founders = [pedigree[ID] for ID in founderIDs]
    if probandIDs == founderIDs
        return Matrix{Float64}(undef, length(founderIDs), length(founderIDs))
    end
    n = length(pedigree)
    Φ = zeros(n, n) * -1
    for (f, founder) in enumerate(founders)
        Φ[founder.index, founder.index] = (1 + Ψ[f, f]) / 2
    end
    for (f, founder₁) in enumerate(founders), (g, founder₂) in enumerate(founders)
        if founder₁.ID < founder₂.ID
            Φ[founder₁.index, founder₂.index] = Φ[founder₂.index, founder₁.index] = Ψ[f, g]
        end
    end
    ancestors = set_ancestors(pedigree)
    for (ID, individual) in pedigree
        i = individual.index
        for founder in founders
            f = founder.index
            if (ID != founder.ID) && (founder.ID ∉ ancestors[ID])
                Φ[i, f] = Φ[f, i] = 0
            end
        end
    end
    for (IDᵢ, individualᵢ) in pedigree, (IDⱼ, individualⱼ) in pedigree
        i = individualᵢ.index
        j = individualⱼ.index
        coefficient = 0.
        if i == j
            coefficient += 0.5
            father = individualᵢ.father
            mother = individualᵢ.mother
            if !isnothing(father)
                if !isnothing(mother)
                    m = mother.index
                    f = father.index
                    coefficient += Φ[m, f] / 2
                end
            end
        elseif IDᵢ ∉ ancestors[IDⱼ]
            father = individualᵢ.father
            if !isnothing(father)
                f = father.index
                coefficient += Φ[f, j] / 2
            end
            mother = individualᵢ.mother
            if !isnothing(mother)
                m = mother.index
                coefficient += Φ[m, j] / 2
            end
        end
        Φ[i, j] = Φ[j, i] = coefficient
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
"""
function ϕ(pedigree::OrderedDict{Int64, Individual}, probandIDs::Vector{Int64} = pro(pedigree); verbose::Bool = false)
    founderIDs = founder(pedigree)
    Ψ = zeros(length(founderIDs), length(founderIDs))
    upperIDs = cut_vertices(pedigree)
    lowerIDs = probandIDs
    C = [upperIDs, lowerIDs]
    while upperIDs != lowerIDs
        isolated_pedigree = branching(pedigree, pro = upperIDs)
        lowerIDs = copy(upperIDs)
        upperIDs = cut_vertices(isolated_pedigree)
        pushfirst!(C, upperIDs)
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