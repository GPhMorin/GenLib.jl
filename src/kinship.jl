"""
phi(individual₁::ReferenceIndividual, individual₂::ReferenceIndividual)

Computes the kinship coefficient between two individuals using references.
"""
function phi(individual₁::ReferenceIndividual, individual₂::ReferenceIndividual)
    # Ported from GENLIB's Kinship
    value = 0.
    if individual₂.index > individual₁.index
        if !isnothing(individual₂.father)
            value += phi(individual₂.father, individual₁) / 2
        end
        if !isnothing(individual₂.mother)
            value += phi(individual₂.mother, individual₁) / 2
        end
    elseif individual₁.index == individual₂.index
        value += 1/2
        if !isnothing(individual₁.father) & !isnothing(individual₁.mother)
            value += phi(individual₁.father, individual₁.mother) / 2
        end
    else
        if !isnothing(individual₁.father)
            value += phi(individual₁.father, individual₂) / 2
        end
        if !isnothing(individual₁.mother)
            value += phi(individual₁.mother, individual₂) / 2
        end
    end
    value
end

"""
phi(genealogy::Dict{Int64, Individual}, IDs::Vector{Int64} = pro(genealogy))

Takes a `genealogy` dictionary, computes the kinship coefficient
between all probands using a vector of `IDs` and returns a matrix.

For faster processing, use `ϕ` instead if a child cannot have a single unknown parent.
"""

function phi(genealogy::OrderedDict{Int64, Individual}, IDs::Vector{Int64} = pro(genealogy))
    matrix = zeros(length(IDs), length(IDs))
    reference = refer(genealogy)
    individuals = [reference[ID] for ID in IDs]
    Threads.@threads for j in eachindex(individuals)
        Threads.@threads for i in eachindex(individuals)
            individual₁ = individuals[i]
            individual₂ = individuals[j]
            if individual₂.ID > individual₁.ID
                matrix[i, j] = matrix[j, i] = phi(individual₁, individual₂)
            elseif individual₁.ID == individual₂.ID
                matrix[i, j] = phi(individual₁, individual₁)
            end
        end
    end
    matrix
end

"""
ϕ(individual₁::ReferenceIndividual, individual₂::ReferenceIndividual)

Computes the kinship coefficient between two individuals using references.
"""
function ϕ(individual₁::ReferenceIndividual, individual₂::ReferenceIndividual)
    # Ported from GENLIB's Kinship
    value = 0.
    if individual₂.index > individual₁.index
        if !isnothing(individual₂.father)
            value += ϕ(individual₂.father, individual₁) / 2
            value += ϕ(individual₂.mother, individual₁) / 2
        end
    elseif individual₁.index == individual₂.index
        value += 1/2
        if !isnothing(individual₁.father)
            value += ϕ(individual₁.father, individual₁.mother) / 2
        end
    else
        if !isnothing(individual₁.father)
            value += ϕ(individual₁.father, individual₂) / 2
            value += ϕ(individual₁.mother, individual₂) / 2
        end
    end
    value
end

"""
"""
function initialize_ancestry(genealogy::OrderedDict{Int64, Individual})
    A = Bool.(zeros(length(genealogy), length(genealogy)))
    for founderID in founder(genealogy)
        i = genealogy[founderID].index
        A[i, i] = true
    end
    for (ID₁, individual₁) in genealogy
        i = individual₁.index
        if !A[i, i]
            father = individual₁.father
            mother = individual₁.mother
            for (ID₂, individual₂) in genealogy
                v = individual₂.index
                if father != 0
                    f = genealogy[father].index
                    if A[f, v]
                        A[i, v] = true
                    end
                end
                if mother != 0
                    m = genealogy[mother].index
                    if A[m, v]
                        A[i, v] = true
                    end
                end
                if i == v
                    A[i, i] = true
                end
            end
        end
    end
    A
end

function cut_vertex(individual::ReferenceIndividual, candidateID::Int64)
    value = true
    for child in individual.children
        if child.state == PROBAND
            return false
        elseif child.ID != candidateID
            value = value && cut_vertex(child, candidateID)
        end
    end
    value
end

function cut_vertices(genealogy::OrderedDict{Int64, Individual})
    vertices = Int64[]
    reference = refer(genealogy)
    probandIDs = pro(genealogy)
    for ID in probandIDs
        reference[ID].state = PROBAND
    end
    founderIDs = founder(genealogy)
    unorderedIDs = [ID for ID in collect(keys(genealogy))]
    indices = [genealogy[ID].index for ID in unorderedIDs]
    order = sortperm(indices)
    candidateIDs = unorderedIDs[order]
    for candidateID in candidateIDs
        ancestorIDs = ancestor(genealogy, candidateID)
        sourceIDs = filter(x -> x ∈ founderIDs, ancestorIDs)
        if all(cut_vertex(reference[sourceID], candidateID) for sourceID in sourceIDs)
            push!(vertices, candidateID)
        end
    end
    ancestorIDs = union([ancestor(genealogy, vertex) for vertex in vertices]...)
    vertices = filter(x -> (x ∈ vertices && x ∉ ancestorIDs) || (x ∈ founderIDs && x ∉ ancestorIDs), candidateIDs)
    vertices
end

"""
"""

function ϕ(genealogy::OrderedDict{Int64, Individual}, Ψ::Matrix{Float64})
    probands = filter(x -> isempty(genealogy[x].children), collect(keys(genealogy)))
    founders = filter(x -> (genealogy[x].father == 0) && (genealogy[x].mother == 0), collect(keys(genealogy)))
    if probands == founders
        return zeros(length(founders), length(founders))
    end
    n = length(genealogy)
    Φ = zeros(n, n) * -1
    for (f, founder) in enumerate(founders)
        i = genealogy[founder].index
        Φ[i, i] = (1 + Ψ[f, f]) / 2
    end
    for (f, ID₁) in enumerate(founders), (g, ID₂) in enumerate(founders)
        if ID₁ != ID₂
            Φ[genealogy[ID₁].index, genealogy[ID₂].index] = Ψ[f, g]
        end
    end
    ancestors = set_ancestors(genealogy)
    for (i, (ID, individual)) in enumerate(genealogy)
        for founder in founders
            if (ID != founder) && (founder ∉ ancestors[ID])
                f = genealogy[founder].index
                Φ[i, f] = Φ[f, i] = 0
            end
        end
    end
    for (i, (ID₁, individual₁)) in enumerate(genealogy), (j, (ID₂, individual₂)) in enumerate(genealogy)
        coefficient = 0.
        if i > j # i cannot be an ancestor of j
            father = individual₁.father
            if (father != 0)
                p = genealogy[father].index
                coefficient += Φ[p, j] / 2
            end
            mother = individual₁.mother
            if mother != 0
                m = genealogy[mother].index
                coefficient += Φ[m, j] / 2
            end
        elseif i < j # i can be an ancestor of j
            if ID₁ ∉ ancestors[ID₂]
                father = individual₁.father
                if father != 0
                    p = genealogy[father].index
                    coefficient += Φ[p, j] / 2
                end
                mother = individual₁.mother
                if mother != 0
                    m = genealogy[mother].index
                    coefficient += Φ[m, j] / 2
                end
            end
        else # i == j
            coefficient += 0.5
            father = individual₁.father
            mother = individual₁.mother
            if father != 0
                if mother != 0
                    p = genealogy[father].index
                    m = genealogy[mother].index
                    coefficient += Φ[m, p] / 2
                end
            end
        end
        Φ[i, j] = Φ[j, i] = coefficient
    end
    for i in 1:n
        Φ[i, i] = (2 * Φ[i, i]) - 1
    end
    indices = [genealogy[ID].index for ID in probands]
    Φ[indices, indices]
end

function ϕ(genealogy::OrderedDict{Int64, Individual}; pro = pro(genealogy), verbose = false)
    founders = founder(genealogy)
    Ψ = zeros(length(founders), length(founders))
    upperIDs = cut_vertices(genealogy)
    lowerIDs = pro
    C = [upperIDs, lowerIDs]
    while upperIDs != lowerIDs
        isolated_genealogy = branching(genealogy; pro = upperIDs)
        lowerIDs = copy(upperIDs)
        upperIDs = cut_vertices(isolated_genealogy)
        pushfirst!(C, upperIDs)
    end
    for i in 1:length(C)-1
        upperIDs = C[i]
        lowerIDs = C[i+1]
        Vᵢ = branching(genealogy; pro = lowerIDs, ancestors = upperIDs)
        Ψ = ϕ(Vᵢ, Ψ)
        if verbose
            println("Kinships for generation ", i, "/", length(C)-1, " completed.")
        end
    end
    probands = filter(x -> x ∈ pro, collect(keys(genealogy)))
    order = sortperm(probands)
    Ψ[order, order]
end

function set_ancestors(genealogy::OrderedDict{Int64, Individual})
    ancestors = Dict()
    for (ID, individual) in genealogy
        ancestors[ID] = [ID]
        mother = individual.mother
        father = individual.father
        if (mother != 0)
            ancestors[ID] = union(ancestors[mother], ancestors[ID])
        end
        if (father != 0)
            ancestors[ID] = union(ancestors[father], ancestors[ID])
        end
    end
    ancestors
end