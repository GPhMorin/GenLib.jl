import GenLib as gen
using DataFrames
using Test

@testset "GenLib.jl" begin
    inds = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    fathers = [0, 0, 0, 1, 1, 0, 3, 3, 6, 6]
    mothers = [0, 0, 0, 2, 2, 0, 4, 4, 5, 5]
    sexes = [1, 2, 1, 2, 2, 1, 2, 1, 1, 2]
    df = DataFrame([inds, fathers, mothers, sexes], [:ind, :father, :mother, :sex])
    ped = gen.genealogy(df)
    @test ped[9].mother.sex == 2

    genea140 = gen.genea140
    ped = gen.genealogy(genea140)
    @test ped[33724].children[2].father.ID == 33724
    @test sum(gen.gc(ped)) == 140

    geneaJi = gen.geneaJi
    ped = gen.genealogy(geneaJi)
    @test gen.pro(ped) == [1, 2, 29]
    @test gen.founder(ped) == [17, 19, 20, 23, 25, 26]
    @test gen.findMRCA(ped, [1, 2, 29]) == ([14, 20], [4 4; 4 4; 3 3])
    @test gen.phi(ped[1], ped[2]) == 0.37109375
    @test gen.phi(ped) == [0.591796875 0.37109375 0.072265625;
                           0.37109375 0.591796875 0.072265625;
                           0.072265625 0.072265625 0.53515625]
end
