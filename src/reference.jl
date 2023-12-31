@enum STATE begin
    PROBAND
    ANCESTOR
    EXPLORED
    EXPLOREDPROBAND
    UNEXPLORED
    NODE
    START
end

mutable struct ReferenceIndividual
    ID::Int64
    father::Union{Nothing, ReferenceIndividual}
    mother::Union{Nothing, ReferenceIndividual}
    index::Int64
    children::Vector{ReferenceIndividual}
    sex::SEX
    state::STATE
    probability::Float64
    allele::Int64
    sort::Int64
    haplotype₁::Int64
    haplotype₂::Int64
    ancestor::Bool
    descendant::Bool
    occurrence::Int64
end

"""
refer(genealogy::Dict{Int64, Individual})

Takes a `genealogy` dictionary and returns a dictionary of references to individuals.
In REPL: to avoid crash, end function call with `;`.
"""
function refer(genealogy::Dict{Int64, Individual})
    reference::Dict{Int64, ReferenceIndividual} = Dict()
    for (ID, individual) in genealogy
        reference[ID] = ReferenceIndividual(
            ID, # ID (ind in the ASC file)
            nothing, # father
            nothing, # mother
            individual.index, # index in the original ASC file
            [], # children
            individual.sex, # sex (1 for male, 2 for female)
            UNEXPLORED, # status
            0., # probability
            0, # allele
            0, # sort
            0, # haplotype₁
            0, # haplotype₂
            false, # whether the individual is an ancestor
            false, # whether the individual is a descendant
            0) # occurrence
    end
    for (ID, individual) in genealogy
        reference_individual = reference[ID]
        if individual.father != 0
            reference_individual.father = reference[individual.father]
            push!(reference_individual.father.children, reference_individual)
        end
        if individual.mother != 0
            reference_individual.mother = reference[individual.mother]
            push!(reference_individual.mother.children, reference_individual)
        end
    end
    reference
end