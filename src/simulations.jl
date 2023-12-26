struct Haplotype
    haplotype::String
    position::Int64
    fixed::Bool
end

function poisson_crossover(sex::SEX,
                           parameters::Vector{Float64},
                           cM_length::Vector{Float64},
                           recombinations::Int64)::Vector{Int64}
    crossovers = Int64[]

    uniform_distribution = Uniform(0, 1)
    poisson_distribution₁ = Poisson(parameters[1])
    poisson_distribution₂ = Poisson(parameters[2])

    if sex == MALE
        recombinations = rand(poisson_distribution₁)
        for _ in 1:recombinations
            position = rand(uniform_distribution)
            push!(crossovers, position)
        end
    end
    if sex == FEMALE
        recombinations = rand(poisson_distribution₂)
        for _ in 1:recombinations
            position = rand(uniform_distribution)
            push!(crossovers, position)
        end
    end
    sort!(crossovers)
end

function ztpoisson_crossover()

end

function gamma_crossover()

end

function simuHaplo(
    genealogy::Dict{Int64, Individual};
    probandIDs::Vector{Int64} = pro(genealogy),
    ancestorIDs::Vector{Int64} = founder(genealogy),
    iterations::Int64 = 1,
    model::Symbol = :poisson,
    parameters::Vector{Float64} = [1, 1],
    cM_length::Vector{Float64} = [100, 100],
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

    nodes = sort(collect(keys(pointer)))
    prepare_priority_sort(nodes)

    index = 1
    order = Dict{Int64, PointerIndividual}()
    jumps = Dict{Int64, Int64}()
    for ID in ancestorIDs
        ancestor = pointer(ID)
        start_priority_sort(ancestor, order, index, jumps)
    end

    if model == :poisson
        crossover = poisson_crossover
    elseif model == :ztpoisson
        crossover = zero_truncated_crossover
    elseif model == :gamma
        crossover = gamma_crossover
    end
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

function parse_output(filename::String, founder_haplotype::String)::Matrix{Int64}
    descent = DefaultDict{Int64, Vector{Tuple{Int64, Int64}}}([])
    file = open(filename)
    lines = readlines(file)
    close(file)
    proband_indices = parse(Int64, split(lines[1], ';')[2])
    Threads.@threads for proband_index in 1:proband_indices
        line = lines[proband_index+1]
        information, chromosome₁, chromosome₂ = filter(!isempty, split(line, ['{', '}']))
        proband = parse(Int64, split(information, ';')[2])
        chromosome₁ = split(chromosome₁, ';')
        chromosome₂ = split(chromosome₂, ';')
        current_index = 2
        while current_index < length(chromosome₁)
            current_start = parse(Int64, chromosome₁[current_index-1]) + 1
            current_haplotype = chromosome₁[current_index]
            current_end = parse(Int64, chromosome₁[current_index+1])
            if current_haplotype == founder_haplotype
                println("YES")
                push!(descent[proband], (current_start, current_end))
            end
        end
        current_index = 2
        while current_index < length(chromosome₂)
            current_start = parse(Int64, chromosome₂[current_index-1]) + 1
            current_haplotype = chromosome₂[current_index]
            current_end = parse(Int64, chromosome₂[current_index+1])
            if current_haplotype == founder_haplotype
                println("YA")
                push!(descent[proband], (current_start, current_end))
            end
            current_index += 2
        end
    end
    probands = sort(collect(keys(descent)))
    ibd = Matrix{Int64}(undef, length(probands), length(probands))
    for (j, proband₁) in enumerate(probands)
        for (i, proband₂) in enumerate(probands)
            ibd[i, j] = ibd[j, i] = length(findall(descent[proband₁] .& descent[proband₂]))
        end
    end
    ibd
end