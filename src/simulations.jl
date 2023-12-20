function simulate_allele(genealogy::Dict{Int64, Individual}, ancestorID::Int64, probandIDs::Vector{Int64}, frequency::Float64, max_simulations::Int64)::Vector{Int64}
    number_of_simulations = 1
    current_frequency = 0
    affected = Vector{Int64}()
    pointer = point(genealogy)
    ancestor = pointer[ancestorID]
    ancestor.genotype = [:m, :M]
    while (current_frequency < 0.9 * frequency) | (current_frequency > 1.1 * frequency) & (number_of_simulations ≤ max_simulations)
        transmit!(ancestor)
        affected = [probandID for probandID in probandIDs if :m ∈ pointer[probandID].genotype]
        current_frequency = length(affected) / length(probandIDs)
        number_of_simulations += 1
    end
    affected
end

function transmit!(individual::PointerIndividual)::Nothing
    for child in individual.children
        transmitted_allele = rand(individual.genotype)
        if individual.sex == MALE
            child.genotype[1] = transmitted_allele
            transmit!(child)
        else
            child.genotype[2] = transmitted_allele
            transmit!(child)
        end
    end
end