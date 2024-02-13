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
    @test repr(MIME("text/plain"), ped) == "A pedigree with:\n41523 individuals;\n68248 parent-child relations;\n20773 men;\n20750 women;\n140 subjects;\n18 generations."
    @test repr(MIME("text/plain"), ped[33724]) == "ind: 33724\nfather: 10086\nmother: 10087\nsex: 1"
    @test gen.children(ped, 33724) == [10033, 113470]
    @test ped[33724].children[2].father.ID == 33724
    @test sum(gen.gc(ped)) == 140
    @test gen.descendant(ped, 10086) == [10009, 10018, 10033, 33724, 105379, 113470,
                                         408065, 408069, 408075, 408375, 408937, 409808,
                                         623919, 712249, 712256, 860834, 860838, 868738,
                                         868740, 868743]

    geneaJi = gen.geneaJi
    ped = gen.genealogy(geneaJi)
    @test gen.pro(ped) == [1, 2, 29]
    @test gen.founder(ped) == [17, 19, 20, 23, 25, 26]
    @test gen.findMRCA(ped, [1, 2, 29]) == ([14, 20], [4 4; 4 4; 3 3])
    @test gen.phi(ped[1], ped[2]) == 0.37109375
    @test gen.phi(ped, verbose=true) == [0.591796875 0.37109375 0.072265625;
                           0.37109375 0.591796875 0.072265625;
                           0.072265625 0.072265625 0.53515625]
    @test gen.phi(ped, [1, 2], [1, 2, 29]) == [0.591796875 0.37109375 0.072265625; 
                                               0.37109375 0.591796875 0.072265625]
    founder1 = ped[17]
    founder1.sort = 1
    founder2 = ped[19]
    founder2.sort = 2
    Ψ = [0.5 0; 0 0.5]
    @test gen.phi(founder1, founder2, Ψ) == 0
    probands = gen.pro(ped)
    @test all(ped[proband].state == gen.UNEXPLORED for proband in probands)
    @test gen.findFounders(ped, [1, 2, 29]) == [17, 19, 20, 25, 26]
    @test gen.rec(ped) == [3, 3, 3, 1, 3, 3]
    @test gen.findDistance(ped, [1, 2], 25) == 12

    gen.save_genealogy(ped, "test.asc")
    ped = gen.genealogy("test.asc")
    @test length(ped) == 29
    rm("test.asc")
end
