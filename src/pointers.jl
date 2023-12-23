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
    allele::Int64
    sort::Int64
    haplotype₁::Int64
    haplotype₂::Int64
end

"""
point(genealogy::Dict{Int64, Individual})::Dict{Int64, PointerIndividual}

Takes a `genealogy` dictionary and returns a dictionary of pointers to individuals.
In REPL: to avoid crash, end function call with `;`.
"""
function point(genealogy::Dict{Int64, Individual})::Dict{Int64, PointerIndividual}
    pointer::Dict{Int64, PointerIndividual} = Dict()
    for (ID, individual) in genealogy
        pointer[ID] = PointerIndividual(ID, nothing, nothing, individual.index, [], individual.sex, UNEXPLORED, 0, 0, 0, 0)
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