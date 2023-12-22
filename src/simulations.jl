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

"""
"""
function simuHaplo(
    genealogy::Dict{Int64, Individual};
    probands::Vector{Int64} = pro(genealogy),
    ancestors::Vector{Int64} = founder(genealogy),
    iterations::Int64 = 1,
    model::Symbol = :poisson,
    parameters::Vector{Float64} = [1, 1],
    cM_length::Vector{Int64} = [100, 100],
    BP_length::Vector{Int64} = 100000000,
    physical_map_mother::Union{Nothing, DataFrame} = nothing,
    physical_map_father::Union{Nothing, DataFrame} = nothing,
    seed::Int64 = 42,
    all_nodes::Bool = false,
    output_directory::String = "."
    )::Nothing


end

function simuHaplo_convert(directory::String = ".")::Nothing
end

function simuHaplo_IBD_compare(
    probandID₁::Int64,
    probandID₂::Int64,
    BP_length::Int64,
    proband_haplotypes_path::String
    )::DataFrame
end

function simuHaplo_traceback(
    genealogy::Dict{Int64, Individual},
    probandID::Int64,
    ancestorID::Int64,
    all_nodes_path::String,
    proband_haplotypes_path::String
    )::DataFrame
end