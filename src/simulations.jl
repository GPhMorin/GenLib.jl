struct Haplotype
    haplotype::String
    position::Int64
    fixed::Bool
end

function simuHaplo(
    genealogy::Dict{Int64, Individual};
    probandIDs::Vector{Int64} = pro(genealogy),
    ancestorIDs::Vector{Int64} = founder(genealogy),
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

    # Ported from GENLIB's simuHaplo

    # Initialize all the nodes
    pointer = point(genealogy)
    haplotype = Dict{Int64, Haplotype}()

    # Label the nodes that are probands
    for ID in probandIDs
        pointer[ID].state = EXPLOREDPROBAND
    end

    # Label the starting points and ancestor haplotypes
    key = 1
    for ID in ancestorIDs
        ancestor = pointer[ID]
        ancestor.state = START
        ancestor.haplotype₁ = key
        haplotype₁ = Haplotype("$(ancestor.ID).1", -1, true)
        haplotype[key] = haplotype₁
        key += 1
        ancestor.haplotype₂ = key
        haplotype₂ = Haplotype("$(ancestor.ID).2", -1, true)
        haplotype[key] = haplotype₂
        key += 1
    end

    # Label the nodes' relevance
    for ID in ancestorIDs
        ancestor = pointer[ID]
        explore_tree(ancestor)
    end

    prepare_priority_sort(pointer)

    
end

function simuHaplo_convert(directory::String = ".")::Nothing
end

function simuHaplo_IBD_compare(
    probandID₁::Int64,
    probandID₂::Int64,
    BP_length::Int64,
    proband_haplotypes_path::String
    )::DataFrame

    # Ported from GENLIB's simuHaplo_IBD_compare
end

function simuHaplo_traceback(
    genealogy::Dict{Int64, Individual},
    probandID::Int64,
    ancestorID::Int64,
    all_nodes_path::String,
    proband_haplotypes_path::String
    )::DataFrame

    # Ported from GENLIB's simuHaplo_traceback
end