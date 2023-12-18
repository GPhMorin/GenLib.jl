"""
phi(individual₁::PointerIndividual, individual₂::PointerIndividual)::Float64

Computes the kinship coefficient between two individuals using pointers.
"""
function phi(individual₁::PointerIndividual, individual₂::PointerIndividual)::Float64
    # Ported from GENLIB's Kinship
    value = 0.
    if individual₁ == individual₂
        value += 1/2
        if !isnothing(individual₁.father) & !isnothing(individual₁.mother)
            value += phi(individual₁.father, individual₁.mother) / 2
        end
        return value
    elseif individual₂.index > individual₁.index
        if !isnothing(individual₂.father)
            value += phi(individual₂.father, individual₁) / 2
        end
        if !isnothing(individual₂.mother)
            value += phi(individual₂.mother, individual₁) / 2
        end
        return value
    else
        if !isnothing(individual₁.father)
            value += phi(individual₁.father, individual₂) / 2
        end
        if !isnothing(individual₁.mother)
            value += phi(individual₁.mother, individual₂) / 2
        end
        return value
    end
    
end

"""
phi(genealogy::Dict{Int64, Individual}, IDs::Vector{Int64} = pro(genealogy))::Matrix{Float64}

Takes a `genealogy` dictionary, computes the kinship coefficient between all probands using a vector of `IDs` and returns a matrix.
"""
function phi(genealogy::Dict{Int64, Individual}, IDs::Vector{Int64} = pro(genealogy))::Matrix{Float64}
    matrix = zeros(length(IDs), length(IDs))
    pointer = point(genealogy)
    individuals = [pointer[ID] for ID in IDs]
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