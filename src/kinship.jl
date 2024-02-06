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
function phi(genealogy::Dict{Int64, Individual}, IDs::Vector{Int64} = pro(genealogy))
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
function initialize_ancestry(genealogy)
    A = SparseMatrixCSC{Bool, Int64}(undef, length(genealogy), length(genealogy))
    unorderedIDs = [ID for ID in collect(keys(genealogy))]
    indices = [genealogy[ID].index for ID in unorderedIDs]
    order = sortperm(indices)
    IDs = unorderedIDs[order]
    for founderID in founder(genealogy)
        i = findfirst(founderID .== IDs)
        A[i, i] = true
    end
    for ID in IDs
        i = findfirst(ID .== IDs)
        if !A[i, i]
            f = findfirst(IDs .== genealogy[ID].father)
            m = findfirst(IDs .== genealogy[ID].mother)
            for possible_ancestor in IDs
                v = findfirst(possible_ancestor .== IDs)
                if !isnothing(f)
                    if A[f, v]
                        A[i, v] = true
                    end
                end
                if !isnothing(m)
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

"""
"""
function cut_vertices(genealogy, probandIDs, founderIDs)
    vertex_cuts = Int64[]
    for candidateID in collect(keys(genealogy))
        if (candidateID ∉ probandIDs) && (candidateID ∉ founderIDs)
            stack = ancestor(genealogy, candidateID)
            is_candidate = true
            while !isempty(stack)
                ID = pop!(stack)
                if ID ∈ probandIDs
                    is_candidate = false
                    break
                end
                if ID != candidateID
                    push!(stack, genealogy[ID].children...)
                end
            end
            if is_candidate
                push!(vertex_cuts, candidateID)
            end
        end
    end
    vertex_cuts = [vertex₁ for vertex₁ in vertex_cuts if !any(vertex₁ ∈ ancestor(genealogy, vertex₂) for vertex₂ in vertex_cuts)]
    vertex_cuts
end

"""
function cut_vertices(A)
    indices = Int64[]
    number = size(A, 1)
    for candidate in 1:number
        candidate_ancestors = Array(A[candidate, :])
        ancestors_descendants = [any(row) for row in eachrow(A[:, candidate_ancestors])]
        candidate_descendants = Array(A[:, candidate])
        leaks = ancestors_descendants .- (candidate_ancestors .|| candidate_descendants) .> 0
        if !any(leaks)
            push!(indices, candidate)
        end
    end
    indices
end
"""

function cut_vertex(individual, candidateID)
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

function cut_vertices(genealogy)
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
ϕ(genealogy, probandIDs = pro(genealogy), founderIDs = founder(genealogy))

Takes a `genealogy` dictionary, computes the kinship coefficient
between all probands using a vector of `probandIDs`
and a vector of `founderIDs` and returns a matrix.

An implementation of algorithm 1 found in Kirkpatrick et al., 2019.
"""
"""
function ϕ(genealogy, probandIDs = pro(genealogy))
    founderIDs = founder(genealogy)
    if probandIDs == founderIDs
        return zeros(length(founderIDs), length(founderIDs))
    end

    isolated_genealogy = branching(genealogy; probandIDs = probandIDs)
    vertex_cuts = cut_vertices(isolated_genealogy)
    Ψ = ϕ(genealogy, vertex_cuts)

    unorderedIDs = [ID for ID in collect(keys(isolated_genealogy))]
    indices = [genealogy[ID].index for ID in unorderedIDs]
    order = sortperm(indices)
    IDs = unorderedIDs[order]

    Φ = ones(length(IDs), length(IDs)) .* -1

    for founderID in vertex_cuts
        f₁ = findfirst(founderID .== IDs)
        f₂ = findfirst(founderID .== vertex_cuts)
        Φ[f₁, f₁] = (1 + Ψ[f₂, f₂]) / 2
    end
    for founderID₁ in vertex_cuts, founderID₂ in vertex_cuts
        f₁ = findfirst(founderID₁ .== IDs)
        g₁ = findfirst(founderID₂ .== IDs)
        if f₁ != g₁
            f₂ = findfirst(founderID₁ .== vertex_cuts)
            g₂ = findfirst(founderID₂ .== vertex_cuts)
            Φ[f₁, g₁] = Ψ[f₂, g₂]
        end
    end

    A = initialize_ancestry(isolated_genealogy)

    for founderID in founderIDs, ID in IDs
        f = findfirst(founderID .== IDs)
        j = findfirst(ID .== IDs)
        if !A[j, f]
            Φ[j, f] = 0
        end
    end
    for ID₁ in IDs
        for ID₂ in IDs
            p = findfirst(IDs .== genealogy[ID₁].father)
            m = findfirst(IDs .== genealogy[ID₁].mother)
            if !isnothing(p) && !isnothing(m)
                i = findfirst(IDs .== ID₁)
                j = findfirst(IDs .== ID₂)
                if i == j
                    Φ[i, i] = (1 + Φ[m, p]) / 2
                elseif !A[j, i]
                    Φ[i, j] = Φ[j, i] =  (Φ[m, j] + Φ[p, j]) / 2
                end
            elseif !isnothing(p)
                i = findfirst(IDs .== ID₁)
                j = findfirst(IDs .== ID₂)
                if i == j
                    Φ[i, i] = 0.5
                elseif !A[j, i]
                    Φ[i, j] = Φ[j, i] =  Φ[p, j] / 2
                end
            elseif !isnothing(m)
                i = findfirst(IDs .== ID₁)
                j = findfirst(IDs .== ID₂)
                if i == j
                    Φ[i, i] = 0.5
                elseif !A[j, i]
                    Φ[i, j] = Φ[j, i] =  Φ[m, j] / 2
                end
            end
        end
    end
    for i in eachindex(IDs)
        Φ[i, i] = (2 * Φ[i, i]) - 1
    end
    indices = [findfirst(ID .== IDs) for ID in IDs if ID ∈ probandIDs]
    Φ[indices, indices]
end
"""

"""
"""

function ϕ(genealogy, Ψ, probandIDs = pro(genealogy))
    founderIDs = founder(genealogy)
    if probandIDs == founderIDs
        return zeros(length(founderIDs), length(founderIDs))
    end

    unorderedIDs = [ID for ID in collect(keys(isolated_genealogy))]
    indices = [genealogy[ID].index for ID in unorderedIDs]
    order = sortperm(indices)
    IDs = unorderedIDs[order]

    Φ = ones(length(IDs), length(IDs)) .* -1

    for founderID in vertex_cuts
        f₁ = findfirst(founderID .== IDs)
        f₂ = findfirst(founderID .== vertex_cuts)
        Φ[f₁, f₁] = (1 + Ψ[f₂, f₂]) / 2
    end
    for founderID₁ in vertex_cuts, founderID₂ in vertex_cuts
        f₁ = findfirst(founderID₁ .== IDs)
        g₁ = findfirst(founderID₂ .== IDs)
        if f₁ != g₁
            f₂ = findfirst(founderID₁ .== vertex_cuts)
            g₂ = findfirst(founderID₂ .== vertex_cuts)
            Φ[f₁, g₁] = Ψ[f₂, g₂]
        end
    end

    A = initialize_ancestry(isolated_genealogy)

    for founderID in founderIDs, ID in IDs
        f = findfirst(founderID .== IDs)
        j = findfirst(ID .== IDs)
        if !A[j, f]
            Φ[j, f] = 0
        end
    end
    for ID₁ in IDs
        for ID₂ in IDs
            p = findfirst(IDs .== genealogy[ID₁].father)
            m = findfirst(IDs .== genealogy[ID₁].mother)
            if !isnothing(p) && !isnothing(m)
                i = findfirst(IDs .== ID₁)
                j = findfirst(IDs .== ID₂)
                if i == j
                    Φ[i, i] = (1 + Φ[m, p]) / 2
                elseif !A[j, i]
                    Φ[i, j] = Φ[j, i] =  (Φ[m, j] + Φ[p, j]) / 2
                end
            elseif !isnothing(p)
                i = findfirst(IDs .== ID₁)
                j = findfirst(IDs .== ID₂)
                if i == j
                    Φ[i, i] = 0.5
                elseif !A[j, i]
                    Φ[i, j] = Φ[j, i] =  Φ[p, j] / 2
                end
            elseif !isnothing(m)
                i = findfirst(IDs .== ID₁)
                j = findfirst(IDs .== ID₂)
                if i == j
                    Φ[i, i] = 0.5
                elseif !A[j, i]
                    Φ[i, j] = Φ[j, i] =  Φ[m, j] / 2
                end
            end
        end
    end
    for i in eachindex(IDs)
        Φ[i, i] = (2 * Φ[i, i]) - 1
    end
    indices = [findfirst(ID .== IDs) for ID in IDs if ID ∈ probandIDs]
    Φ[indices, indices]
end

function ϕ(genealogy)
    founderIDs = founder(genealogy)
    probandIDs = pro(genealogy)
    Ψ = zeros(length(founderIDs), length(founderIDs))
    ablated_genealogy = genealogy
    while !(probandIDs == founderIDs)
        ablated_genealogy = ablate(ablated_genealogy, founderIDs)
        founderIDs = founder(ablated_genealogy)
        println(founderIDs)
        Ψ = ϕ(genealogy, Ψ, founderIDs)
    end
    ϕ(genealogy, Ψ, probandIDs)
end