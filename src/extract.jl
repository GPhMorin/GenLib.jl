"""
    mutable struct Candidate <: AbstractIndividual
        ID::Int
        father::Union{Nothing, Candidate}
        mother::Union{Nothing, Candidate}
        children::Vector{Candidate}
        sex::Int
        rank::Int
        is_ancestor::Bool
        is_descendant::Bool
    end

A mutable structure used internally to mark whether
an individual is an ancestor and/or a descendant.
"""
mutable struct Candidate <: AbstractIndividual
    ID::Int
    father::Union{Nothing, Candidate}
    mother::Union{Nothing, Candidate}
    children::Vector{Candidate}
    sex::Int
    rank::Int
    is_ancestor::Bool
    is_descendant::Bool
end

"""
    _mark_ancestors!(individual::Candidate)

Recursively mark the ancestors of an `individual`.
"""
function _mark_ancestors!(individual::Candidate)
    if !individual.is_ancestor
        individual.is_ancestor = true
        if !isnothing(individual.father)
            _mark_ancestors!(individual.father)
        end
        if !isnothing(individual.mother)
            _mark_ancestors!(individual.mother)
        end
    end
end

"""
    _mark_descendants!(individual::Candidate)

Recursively mark the descendants of an `individual`.
"""
function _mark_descendants!(individual::Candidate)
    if !individual.is_descendant
        individual.is_descendant = true
        for child ∈ individual.children
            _mark_descendants!(child)
        end
    end
end

"""
    branching(pedigree::Pedigree; pro::Union{Vector{Int}, Nothing} = nothing,
    ancestors::Union{Vector{Int}, Nothing} = nothing)

Return a pedigree that filters individuals who are in the paths
between select probands and ancestors.
"""
function branching(pedigree::Pedigree;
    pro::Union{Vector{Int}, Nothing} = nothing,
    ancestors::Union{Vector{Int}, Nothing} = nothing)
    isolated_pedigree = Pedigree{Individual}()
    marking_pedigree = Pedigree{Candidate}()
    for individual ∈ values(pedigree)
        father = isnothing(individual.father) ?
            nothing : marking_pedigree[individual.father.ID]
        mother = isnothing(individual.mother) ?
            nothing : marking_pedigree[individual.mother.ID]
        marking_pedigree[individual.ID] = Candidate(
            individual.ID,
            isnothing(father) ? nothing : marking_pedigree[father.ID],
            isnothing(mother) ? nothing : marking_pedigree[mother.ID],
            Candidate[],
            individual.sex,
            individual.rank,
            false,
            false)
        if !isnothing(father)
            push!(marking_pedigree[father.ID].children, marking_pedigree[individual.ID])
        end
        if !isnothing(mother)
            push!(marking_pedigree[mother.ID].children, marking_pedigree[individual.ID])
        end
    end
    if !isnothing(pro)
        Threads.@threads for i ∈ eachindex(pro)
            ID = pro[i]
            proband = marking_pedigree[ID]
            _mark_ancestors!(proband)
        end
    end
    if !isnothing(ancestors)
        Threads.@threads for j ∈ eachindex(ancestors)
            ID = ancestors[j]
            ancestor = marking_pedigree[ID]
            _mark_descendants!(ancestor)
        end
    end
    rank = 0
    if !isnothing(pro) && !isnothing(ancestors)
        for individual ∈ values(marking_pedigree)
            if individual.is_ancestor && individual.is_descendant
                rank += 1
                father = individual.father
                if !isnothing(father)
                    father = father.is_ancestor && father.is_descendant ? father : nothing
                end
                mother = individual.mother
                if !isnothing(mother)
                    mother = mother.is_ancestor && mother.is_descendant ? mother : nothing
                end
                isolated_pedigree[individual.ID] = Individual(individual.ID,
                isnothing(father) ? nothing : isolated_pedigree[father.ID],
                isnothing(mother) ? nothing : isolated_pedigree[mother.ID],
                Int[],
                individual.sex,
                rank)
                if !isnothing(father)
                    push!(isolated_pedigree[father.ID].children,
                        isolated_pedigree[individual.ID])
                end
                if !isnothing(mother)
                    push!(isolated_pedigree[mother.ID].children,
                        isolated_pedigree[individual.ID])
                end
            end
        end
    elseif !isnothing(pro)
        for individual ∈ values(marking_pedigree)
            if individual.is_ancestor
                rank += 1
                father = individual.father
                mother = individual.mother
                isolated_pedigree[individual.ID] = Individual(individual.ID,
                isnothing(father) ? nothing : isolated_pedigree[father.ID],
                isnothing(mother) ? nothing : isolated_pedigree[mother.ID],
                Int[],
                individual.sex,
                rank)
                if !isnothing(father)
                    push!(isolated_pedigree[father.ID].children,
                        isolated_pedigree[individual.ID])
                end
                if !isnothing(mother)
                    push!(isolated_pedigree[mother.ID].children,
                        isolated_pedigree[individual.ID])
                end
            end
        end
    elseif !isnothing(ancestors)
        for individual ∈ values(marking_pedigree)
            if individual.is_descendant
                rank += 1
                father = individual.father
                if !isnothing(father)
                    father = father.is_descendant ? father : nothing
                end
                mother = individual.mother
                if !isnothing(mother)
                    mother = mother.is_descendant ? mother : nothing
                end
                isolated_pedigree[individual.ID] = Individual(individual.ID,
                isnothing(father) ? nothing : isolated_pedigree[father.ID],
                isnothing(mother) ? nothing : isolated_pedigree[mother.ID],
                Int[],
                individual.sex,
                rank)
                if !isnothing(father)
                    push!(isolated_pedigree[father.ID].children,
                        isolated_pedigree[individual.ID])
                end
                if !isnothing(mother)
                    push!(isolated_pedigree[mother.ID].children,
                        isolated_pedigree[individual.ID])
                end
            end
        end
    end
    isolated_pedigree
end