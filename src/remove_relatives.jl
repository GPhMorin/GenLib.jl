"""
remove_relativeIDs!(probandIDs::Vector{Int64}, genealogy::Dict{Int64, Individual})::Vector{Int64}

Takes a list of `probandIDs` and, according to a given `genealogy`,
removes IDs of individuals who are first cousins or closer in the genealogy.
"""
function remove_relatives!(probandIDs::Vector{Int64}, genealogy::Dict{Int64, Individual})::Vector{Int64}
    pointer = point(genealogy)
    candidateIDs = copy(probandIDs)
    empty!(probandIDs)
    for probandID in candidateIDs
        proband = pointer[probandID]
        relativeIDs = Vector{Int64}()
        add_relatives!(relativeIDs, proband, Int8(0))
        if isempty(relativeIDs ∩ probandIDs)
            push!(probandIDs, proband.ID)
        end
    end
    probandIDs
end

function add_relatives!(relativeIDs::Vector{Int64}, individual::PointerIndividual, depth::Int8)
    push!(relativeIDs, individual.ID)
    if depth < 4
        if !isnothing(individual.father)
            add_relatives!(relativeIDs, individual.father, depth += Int8(1))
        end
        if !isnothing(individual.mother)
            add_relatives!(relativeIDs, individual.mother, depth += Int8(1))
        end
        for child in individual.children
            add_relatives!(relativeIDs, child, depth += Int8(1))
        end
    end
end

function remove_relatives(genealogy::Dict{Int64, Individual}, probandIDs::Vector{Int64}, threshold::Float64 = 1/16)::Vector{Int64}
    indices = ones(length(probandIDs)) .> 0
    pointer = point(genealogy)
    Threads.@threads for j in eachindex(probandIDs)
        Threads.@threads for i in eachindex(probandIDs)
            if i < j
                if indices[j]
                    if indices[i]
                        proband₁ = pointer[probandIDs[i]]
                        proband₂ = pointer[probandIDs[j]]
                        if ϕ(proband₁, proband₂) > threshold
                            indices[j] = false
                        end
                    end
                end
            end
        end
    end
    probandIDs[indices]
end

function remove_relatives2(genealogy::Dict{Int64, Individual}, probandIDs::Vector{Int64}, threshold::Float64 = 1/16)::Vector{Int64}
    matrix = zeros(length(probandIDs), length(probandIDs))
    pointer = point(genealogy)
    Threads.@threads for i in eachindex(probandIDs)
        proband₁ = pointer[probandIDs[i]]
        Threads.@threads for j in eachindex(probandIDs)
            if i < j
                proband₂ = pointer[probandIDs[j]]
                matrix[i, j] = ϕ(proband₁, proband₂)
            end
        end
    end
    matrix = matrix .≥ threshold
    indices = [!any(row[:]) for row in eachrow(matrix)]
    probandIDs[indices]
end