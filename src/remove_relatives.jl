"""
remove_relativeIDs!(probandIDs::Vector{Int64}, genealogy::Dict{Int64, Individual})::Vector{Int64}

Takes a list of `probandIDs` and, according to a given `genealogy`,
removes IDs of individuals who are first cousins or closer in the genealogy.
"""
function remove_relativeIDs!(probandIDs::Vector{Int64}, genealogy::Dict{Int64, Individual})::Vector{Int64}
    pointer = point(genealogy)
    candidateIDs = copy(probandIDs)
    empty!(probandIDs)
    for probandID in candidateIDs
        proband = pointer[probandID]
        relativeIDs = Vector{Int64}()
        for child in proband.children
            push!(relativeIDs, child.ID)
        end
        if !isnothing(proband.father)
            push!(relativeIDs, proband.father.ID)
            for sibling in proband.father.children
                push!(relativeIDs, sibling.ID)
            end
            if !isnothing(proband.father.father)
                push!(relativeIDs, proband.father.father.ID)
                for avuncular in proband.father.father.children
                    push!(relativeIDs, avuncular.ID)
                    for cousin in avuncular.children
                        push!(relativeIDs, cousin.ID)
                    end
                end
            end
            if !isnothing(proband.father.mother)
                push!(relativeIDs, proband.father.mother.ID)
                for avuncular in proband.father.mother.children
                    push!(relativeIDs, avuncular.ID)
                    for cousin in avuncular.children
                        push!(relativeIDs, cousin.ID)
                    end
                end
            end
        end
        if !isnothing(proband.mother)
            push!(relativeIDs, proband.mother.ID)
            for sibling in proband.mother.children
                push!(relativeIDs, sibling.ID)
            end
            if !isnothing(proband.mother.father)
                push!(relativeIDs, proband.mother.father.ID)
                for avuncular in proband.mother.father.children
                    push!(relativeIDs, avuncular.ID)
                    for cousin in avuncular.children
                        push!(relativeIDs, cousin.ID)
                    end
                end
            end
            if !isnothing(proband.mother.mother)
                push!(relativeIDs, proband.mother.mother.ID)
                for avuncular in proband.mother.mother.children
                    push!(relativeIDs, avuncular.ID)
                    for cousin in avuncular.children
                        push!(relativeIDs, cousin.ID)
                    end
                end
            end
        end
        if isempty(relativeIDs ∩ probandIDs)
            push!(probandIDs, proband.ID)
        end
    end
    probandIDs
end