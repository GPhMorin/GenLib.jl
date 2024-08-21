"""
    genout(pedigree::Pedigree, sorted::Bool = false)

Return a pedigree as a `DataFrame`.

If `sorted` is `false` (the default), then the individuals will appear in the same order as
in the pedigree.

If `sorted` is `true`, then the individuals will appear in alphabetical ID order.

```julia
import GenLib as gen
geneaJi = gen.geneaJi
ped = gen.genealogy(geneaJi)
gen.genout(ped)
```
"""
function genout(pedigree::Pedigree; sorted::Bool = false)
    inds = Int64[]
    fathers = Int64[]
    mothers = Int64[]
    sexes = Int64[]
    if !sorted
        for (ID, individual) ∈ pedigree
            push!(inds, ID)
            push!(fathers, !isnothing(individual.father) ? individual.father.ID : 0)
            push!(mothers, !isnothing(individual.mother) ? individual.mother.ID : 0)
            push!(sexes, individual.sex)
        end
    else # if sorted
        inds = sort(collect(keys(pedigree)))
        for ID ∈ inds
            individual = pedigree[ID]
            push!(fathers, !isnothing(individual.father) ? individual.father.ID : 0)
            push!(mothers, !isnothing(individual.mother) ? individual.mother.ID : 0)
            push!(sexes, individual.sex)
        end
    end
    df = DataFrame([inds, fathers, mothers, sexes], ["ind", "father", "mother", "sex"])
    df
end