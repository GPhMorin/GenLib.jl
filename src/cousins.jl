"""
remove_relatives!(probandIDs::Vector{Int64}, genealogy::Dict{Int64, Individual})::Vector{Int64}

Takes a list of `probandIDs` and, according to a given `genealogy`,
removes IDs of individuals who are first cousins or closer in the genealogy.
"""
function remove_relatives!(probandIDs::Vector{Int64}, genealogy::Dict{Int64, Individual})::Vector{Int64}
    pointer = point(genealogy)
    preserved_probands = Vector{Int64}()
    for probandID in probandIDs
        proband = pointer[probandID]
        relatives = Vector{Int64}()
        for child in proband.children
            push!(relatives, child.ID)
        end
        if !isnothing(proband.father)
            push!(relatives, proband.father.ID)
            for sibling in proband.father.children
                push!(relatives, sibling.ID)
            end
            if !isnothing(proband.father.father)
                push!(relatives, proband.father.father.ID)
                for avuncular in proband.father.father.children
                    push!(relatives, avuncular.ID)
                    for cousin in avuncular.children
                        push!(relatives, cousin.ID)
                    end
                end
            end
            if !isnothing(proband.father.mother)
                push!(relatives, proband.father.mother.ID)
                for avuncular in proband.father.mother.children
                    push!(relatives, avuncular.ID)
                    for cousin in avuncular.children
                        push!(relatives, cousin.ID)
                    end
                end
            end
        end
        if !isnothing(proband.mother)
            push!(relatives, proband.mother.ID)
            for sibling in proband.mother.children
                push!(relatives, sibling.ID)
            end
            if !isnothing(proband.mother.father)
                push!(relatives, proband.mother.father.ID)
                for avuncular in proband.mother.father.children
                    push!(relatives, avuncular.ID)
                    for cousin in avuncular.children
                        push!(relatives, cousin.ID)
                    end
                end
            end
            if !isnothing(proband.mother.mother)
                push!(relatives, proband.mother.mother.ID)
                for avuncular in proband.mother.mother.children
                    push!(relatives, avuncular.ID)
                    for cousin in avuncular.children
                        push!(relatives, cousin.ID)
                    end
                end
            end
        end
        if !isnothing(relatives ∩ preserved_probands)
            push!(preserved_probands, proband.ID)
        end
    end
    probandIDs = preserved_probands
end

"""
remove_cousins!(probandIDs::Vector{Int64}, genealogy::Dict{Int64, Individual})::Vector{Int64}

Takes a list of `probandIDs` and, given a `genealogy`,
removes IDs of individuals with a coefficient of kinship ≥ 1/16
which corresponds to genealogical first cousins or more related.
"""
function remove_cousins!(probandIDs::Vector{Int64}, genealogy::Dict{Int64, Individual})::Vector{Int64}
    remove_relatives!(probandIDs, genealogy)
    kinship_matrix = phi(genealogy, probandIDs)
    preserved_probands = Vector{Int64}()
    for (index₁, pro₁) in enumerate(probandIDs)
        for pro₂ in preserved_probands
            index₂ = findfirst(preserved_probands .== pro₂)
            kinship = kinship_matrix[index₁, index₂]
            if kinship ≥ 1/16
                break
            end
        end
        push!(preserved_probands, pro₁)
    end
    probandIDs = preserved_probands
end