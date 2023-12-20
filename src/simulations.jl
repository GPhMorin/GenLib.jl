"""
simulate_allele(genealogy::Dict{Int64, Individual}, ancestorID::Int64, probandIDs::Vector{Int64}, frequency::Float64, patience::Int64)::Vector{Int64}

Takes a given `genealogy`, an `ancestorID`, a vector of `probandIDs`,
an allele `frequency` and a maximum number of simulations (`patience`)
and simulates an allele's dropping from that ancestor.
Returns the IDs of probands who carry the minor allele (`:m`).
If patience is exhausted, the affected probands from the last simulation are given.
"""
function simulate_allele(genealogy::Dict{Int64, Individual}, ancestorID::Int64, probandIDs::Vector{Int64}, frequency::Float64, patience::Int64)::Vector{Int64}
    number_of_simulations = 1
    current_frequency = 0
    carriers = Vector{Int64}()
    pointer = point(genealogy)
    ancestor = pointer[ancestorID]
    ancestor.genotype = [:M, :m] # :M is for the major allele, :m is for the minor allele
    while (current_frequency < 0.9 * frequency) | (current_frequency > 1.1 * frequency) & (number_of_simulations ≤ patience)
        transmit!(ancestor)
        carriers = [probandID for probandID in probandIDs if :m ∈ pointer[probandID].genotype]
        current_frequency = length(affected) / length(probandIDs)
        number_of_simulations += 1
    end
    carriers
end

"""
transmit!(individual::PointerIndividual)::Nothing

A recursive function that randomly transmits one of an
`individual`'s two alleles down do their children.
"""
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