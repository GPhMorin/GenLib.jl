import GenLib as gen
using DataFrames
using Test

@testset "Custom pedigree" begin
    inds = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    fathers = [0, 0, 0, 1, 1, 0, 3, 3, 6, 6]
    mothers = [0, 0, 0, 2, 2, 0, 4, 4, 5, 5]
    sexes = [1, 2, 1, 2, 2, 1, 2, 1, 1, 2]
    df = DataFrame([inds, fathers, mothers, sexes], [:ind, :father, :mother, :sex])
    ped = gen.genealogy(df)
    @test ped[9].mother.sex == 2
end

@testset "genea140" begin
    genea140 = gen.genea140
    ped = gen.genealogy(genea140)
    @test repr(MIME("text/plain"), ped) ==
        "A pedigree with:\n41523 individuals;\n68248 parent-child relations;\n20773 men;" *
        "\n20750 women;\n140 subjects;\n18 generations."
    @test gen.nomen(ped) == 20773
    @test gen.nowomen(ped) == 20750
    @test gen.noind(ped) == 41523
    @test repr(MIME("text/plain"), ped[33724]) == "ind: 33724\nfather: 10086\n" *
        "mother: 10087\nsex: 1"
    @test gen.children(ped, 33724) == [10033, 113470]
    @test ped[33724].children[2].father.ID == 33724
    @test sum(gen.gc(ped)) == 140
    @test gen.descendant(ped, 10086) == [10009, 10018, 10033, 33724, 105379, 113470,
        408065, 408069, 408075, 408375, 408937, 409808, 623919, 712249, 712256, 860834,
        860838, 868738, 868740, 868743]
    @test gen.depth(ped) == 18
    @test gen.sibship(ped, [11520]) == [15397, 39369, 49658]
    @test gen.sibship(ped, [11520], halfSibling=false) == []
end

@testset "geneaJi" begin
    geneaJi = gen.geneaJi
    ped = gen.genealogy(geneaJi)
    @test gen._lowest_founders(ped) == [17, 19, 20, 25, 26, 29]
    @test gen.pro(ped) == [1, 2, 29]
    @test gen.founder(ped) == [17, 19, 20, 23, 25, 26]
    genMatrix = gen.findMRCA(ped, [1, 2, 29])
    @test genMatrix.individuals == [1, 2, 29]
    @test genMatrix.ancestors == [14, 20]
    @test genMatrix.meioses == [4 4; 4 4; 3 3]
    @test gen.f(ped, [1]) == [0.18359375]
    @test gen.f(ped, [17]) == [0.]
    @test gen.phi(ped[1], ped[2]) == 0.37109375
    phi = gen.phi(ped, verbose = true)
    @test phi == [0.591796875 0.37109375 0.072265625; 0.37109375 0.591796875 0.072265625;
        0.072265625 0.072265625 0.53515625]
    @test gen.phiMean(phi) == 0.171875
    phi = gen.vector_phi(ped, verbose = true)
    @test gen.phiMean(phi) == 0.171875
    @test repr(MIME("text/plain"), phi) == "3×3 KinshipMatrix with 6 stored entries."
    @test phi[1, 2] == 0.37109375
    phi = gen.dict_phi(ped, verbose = true)
    @test gen.phiMean(phi) == 0.171875
    @test repr(MIME("text/plain"), phi) == "3×3 KinshipMatrix with 6 stored entries."
    @test sort(collect(keys(phi))) == gen.pro(ped)
    founder1 = ped[17]
    founder2 = ped[19]
    @test gen.phi(founder1, founder2) == 0
    probands = gen.pro(ped)
    @test gen.findFounders(ped, [1, 2, 29]) == [17, 19, 20, 25, 26]
    @test gen.rec(ped) == [3, 3, 3, 1, 3, 3]
    @test gen.findDistance(ped, [1, 2], 25) == 12
    @test gen.occ(ped) == [6 6 2; 8 8 2; 1 1 2; 0 0 1; 8 8 3; 8 8 3]
    @test gen.occ(ped, typeOcc="TOTAL") == [14; 18; 4; 1; 19; 19;;]
    @test gen._findMinDistanceMRCA(ped, [2, 29]) == 7
    @test gen.completeness(ped, type="IND")[8, 1] == 3.125
    @test gen.completeness(ped, genNo=[0, 4, 6]) == [100.; 62.5; 18.75;;]
    iso_ped = gen.branching(ped, pro=[1])
    @test gen.founder(iso_ped) == [17, 19, 20, 25, 26]
    iso_ped = gen.branching(ped, ancestors=[13])
    @test gen.pro(iso_ped) == [1, 2]
    iso_ped = gen.branching(ped, pro=[1], ancestors=[13])
    @test collect(keys(iso_ped)) == [13, 8, 4, 1]

    gen._save("test.asc", ped)
    ped = gen.genealogy("test.asc")
    @test length(ped) == 29
    rm("test.asc")

    gen._save("test.asc", ped, sorted = true)
    ped = gen.genealogy("test.asc")
    IDs = [ID for ID ∈ collect(keys(ped))]
    order = sort(IDs)
    @test IDs == IDs[order]
    rm("test.asc")
end

@testset "pop140" begin
    pop140 = gen.pop140
    pop = gen._pop(pop140)
    @test pop[217891] == "Saguenay"
end
