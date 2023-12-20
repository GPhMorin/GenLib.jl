function simulate_allele(genealogy::Dict{Int64, Individual}, ancestorID::Int64, probandIDs::Vector{Int64}, frequency::Float64)::Vector{Int64}
    current_frequency = 0
    while 0.9 * frequency < current_frequency < 1.1 * frequency
        pointer = point(genealogy)
        ancestor = pointer[ancestorID]
        ancestor.genotype = [:m, :M]
        transmit!(ancestor)
        affected = [probandID for probandID in probandIDs if :m ∈ pointer[probandID].genotype]
        current_frequency = length(affected) / length(probandIDs)
    end
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