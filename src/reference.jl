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
    sex::Int64
    state::STATE
    probability::Float64
    sort::Int64
    ancestor::Bool
    descendant::Bool
    occurrence::Int64
end

function Base.show(io::IO, individual::ReferenceIndividual)
    println(io, "ind: ", individual.ID)
    println(io, "father: ", !isnothing(individual.father) ? individual.father.ID : 0)
    println(io, "mother: ", !isnothing(individual.mother) ? individual.mother.ID : 0)
    print(io, "sex: ", individual.sex)
end

"""
    refer(genealogy::OrderedDict{Int64, Individual})

Return a dictionary of references to individuals.

In REPL: to avoid crash, end function call with `;`.
"""
function refer(genealogy::OrderedDict{Int64, Individual})
    reference = OrderedDict{Int64, ReferenceIndividual}()
    for (ID, individual) in genealogy
        reference[ID] = ReferenceIndividual(
            ID, # ID (ind in the ASC file)
            nothing, # father
            nothing, # mother
            individual.index, # index in the genealogy
            [], # children
            individual.sex, # sex (1 for male, 2 for female)
            UNEXPLORED, # status
            0., # probability
            0, # sort
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