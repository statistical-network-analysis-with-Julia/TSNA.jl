using TSNA
using Test

@testset "TSNA.jl" begin
    @testset "Module loading" begin
        @test @isdefined(TSNA)
    end

    @testset "tPath construction" begin
        p = tPath([1, 2, 3], [0.5, 1.0], [(1, 2), (2, 3)])
        @test p isa tPath{Int, Float64}
        @test length(p) == 2

        # Mismatched lengths should error
        @test_throws ArgumentError tPath([1, 2], [0.5, 1.0], [(1, 2), (2, 3)])
        @test_throws ArgumentError tPath([1, 2, 3], [0.5], [(1, 2), (2, 3)])
    end

    @testset "Contact types" begin
        c = Contact(1, 2, 0.0, 1.0)
        @test c isa Contact{Int, Float64}
        @test c.source == 1
        @test c.target == 2

        cs = ContactSequence([c], 5)
        @test cs isa ContactSequence{Int, Float64}
        @test length(cs) == 1
        @test cs.n_vertices == 5
    end

    @testset "Temporal centrality API" begin
        @test isdefined(TSNA, :tDegree)
        @test isdefined(TSNA, :tBetweenness)
        @test isdefined(TSNA, :tCloseness)
        @test isdefined(TSNA, :tEigenvector)
        @test isdefined(TSNA, :tPagerank)
    end

    @testset "Temporal network measures API" begin
        @test isdefined(TSNA, :tDensity)
        @test isdefined(TSNA, :tReciprocity)
        @test isdefined(TSNA, :tTransitivity)
        @test isdefined(TSNA, :tEdgeDuration)
        @test isdefined(TSNA, :tVertexDuration)
        @test isdefined(TSNA, :tEdgePersistence)
        @test isdefined(TSNA, :tTurnover)
    end

    @testset "Temporal path API" begin
        @test isdefined(TSNA, :temporalDistance)
        @test isdefined(TSNA, :forwardReachableSet)
        @test isdefined(TSNA, :backwardReachableSet)
        @test isdefined(TSNA, :shortestTemporalPath)
    end

    @testset "Aggregation API" begin
        @test isdefined(TSNA, :tSnaStats)
        @test isdefined(TSNA, :tAggregate)
        @test isdefined(TSNA, :windowSnaStats)
        @test isdefined(TSNA, :as_contact_sequence)
    end

    @testset "Contact sequence utilities" begin
        @test isdefined(TSNA, :tieDecay)
    end
end
