"""
phi(individual₁::PointerIndividual, individual₂::PointerIndividual)::Float64

Computes the kinship coefficient between two individuals using pointers.
"""
function phi(individual₁::PointerIndividual, individual₂::PointerIndividual)::Float64
    # Ported from GENLIB's Kinship
    if individual₁ === individual₂
        if !isnothing(individual₁.father) & !isnothing(individual.mother)
            value = phi(individual₁.father, individual₁.mother)
        else
            value = 0.
        end
        return (1. + value) / 2.
    end
    if individual₂.index > individual₁.index
        individual₁, individual₂ = individual₂, individual₁
    end
    if isnothing(individual₁.father) & isnothing(individual.mother)
        return 0.
    end
    mother_value = 0.
    father_value = 0.
    if !isnothing(individual₁.father)
        father_value = phi(individual₁.father, individual₂)
    end
    if !isnothing(individual.mother)
        mother_value = phi(individual₁.mother, individual₂)
    end
    return (father_value + mother_value) / 2.
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