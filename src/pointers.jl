@enum STATE begin
    PROBAND
    ANCESTOR
    EXPLORED
    EXPLOREDPROBAND
    UNEXPLORED
    NODE
    START
end

mutable struct PointerIndividual
    ID::Int64
    father::Union{Nothing, PointerIndividual}
    mother::Union{Nothing, PointerIndividual}
    index::Int64
    children::Vector{PointerIndividual}
    sex::SEX
    state::STATE
    probability::Float64
    allele::Int64
    sort::Int64
    haplotype₁::Int64
    haplotype₂::Int64
    ancestor::Bool
    descendant::Bool
end

"""
point(genealogy::Dict{Int64, Individual})

Takes a `genealogy` dictionary and returns a dictionary of pointers to individuals.
In REPL: to avoid crash, end function call with `;`.
"""
function point(genealogy::Dict{Int64, Individual})
    pointer::Dict{Int64, PointerIndividual} = Dict()
    for (ID, individual) in genealogy
        pointer[ID] = PointerIndividual(ID, # ID (ind in the ASC file)
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
        false) # whether the individual is a descendant
    end
    for (ID, individual) in genealogy
        pointer_individual = pointer[ID]
        if individual.father != 0
            pointer_individual.father = pointer[individual.father]
            push!(pointer_individual.father.children, pointer_individual)
        end
        if individual.mother != 0
            pointer_individual.mother = pointer[individual.mother]
            push!(pointer_individual.mother.children, pointer_individual)
        end
    end
    pointer
end