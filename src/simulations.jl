mutable struct SimulationIndividual
    ID::Int64
    father::Union{Nothing, SimulationIndividual}
    mother::Union{Nothing, SimulationIndividual}
    children::Vector{SimulationIndividual}
    sex::SEX
    state::STATE
end

"""
simulate(genealogy::Dict{Int64, Individual})::Dict{Int64, SimulationIndividual}

Takes a `genealogy` dictionary and returns a dictionary of pointers to simulated individuals.
In REPL: to avoid crash, end function call with `;`.
"""
function simulate(genealogy::Dict{Int64, Individual})::Dict{Int64, SimulationIndividual}
    simulation::Dict{Int64, SimulationIndividual} = Dict()
    for (ID, individual) in genealogy
        simulation[ID] = SimulationIndividual(ID, nothing, nothing, [], individual.sex, UNEXPLORED)
    end
    for (ID, individual) in genealogy
        simulation_individual = simulation[ID]
        if individual.father != 0
            simulation_individual.father = simulation[individual.father]
            push!(simulation_individual.father.children, simulation_individual)
        end
        if individual.mother != 0
            simulation_individual.mother = simulation[individual.mother]
            push!(simulation_individual.mother.children, simulation_individual)
        end
    end
    simulation
end