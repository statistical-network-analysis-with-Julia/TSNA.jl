using TSNA
using NetworkDynamic
using Network
using Test

# Temporal chain: 1→2 active [0,10), 2→3 active [0,10) (equal onsets),
# 3→4 active [20,30), 4→5 active [5,8) (closes before walker can arrive)
function chain_fixture()
    dnet = DynamicNetwork(5; observation_start=0.0, observation_end=40.0)
    for v in 1:5
        activate!(dnet, 0.0, 40.0; vertex=v)
    end
    activate!(dnet, 0.0, 10.0; edge=(1, 2))
    activate!(dnet, 0.0, 10.0; edge=(2, 3))
    activate!(dnet, 20.0, 30.0; edge=(3, 4))
    activate!(dnet, 5.0, 8.0; edge=(4, 5))
    return dnet
end

@testset "TSNA.jl" begin
    @testset "Earliest arrival: interval semantics" begin
        dnet = chain_fixture()
        arrival, _ = earliestArrival(dnet, 1, 0.0)

        # Equal-onset chain traversed instantly at t=0 (the old
        # onset-ordered pass missed this)
        @test arrival[1] == 0.0
        @test arrival[2] == 0.0
        @test arrival[3] == 0.0
        # Must wait for the [20,30) spell (mid-window boarding at onset)
        @test arrival[4] == 20.0
        # 4→5 closed at 8 < 20: unreachable
        @test !haskey(arrival, 5)

        # Mid-spell boarding: starting at t=4, edge [0,10) is boarded at 4,
        # not rejected (the old code required arrival ≤ onset)
        arr4, _ = earliestArrival(dnet, 1, 4.0)
        @test arr4[2] == 4.0
        @test arr4[3] == 4.0

        # Spells active before start_time still count (old filter dropped
        # them entirely)
        arr9, _ = earliestArrival(dnet, 1, 9.9)
        @test arr9[2] == 9.9

        # After the spell ends: no path
        arr11, _ = earliestArrival(dnet, 1, 11.0)
        @test !haskey(arr11, 2)
    end

    @testset "temporalDistance and paths" begin
        dnet = chain_fixture()

        @test temporalDistance(dnet, 1, 3, 0.0) == 0.0
        @test temporalDistance(dnet, 1, 4, 0.0) == 20.0
        @test temporalDistance(dnet, 1, 5, 0.0) === nothing   # unreachable
        @test temporalDistance(dnet, 1, 4, 25.0) === nothing  # window too late

        p = temporalPath(dnet, 1, 4, 0.0)
        @test p isa tPath
        @test p.vertices == [1, 2, 3, 4]
        @test p.times == [0.0, 0.0, 20.0]
        @test issorted(p.times)
        @test path_duration(p) == 20.0
        @test temporalPath(dnet, 1, 5, 0.0) === nothing
        @test shortestTemporalPath === temporalPath

        # Single-vertex path has zero duration
        p1 = temporalPath(dnet, 1, 1, 0.0)
        @test length(p1) == 0
        @test path_duration(p1) == 0.0
    end

    @testset "Reachable sets and duality" begin
        dnet = chain_fixture()

        @test forwardReachableSet(dnet, 1, 0.0) == [1, 2, 3, 4]
        @test forwardReachableSet(dnet, 1, 15.0) == [1]
        @test forwardReachableSet(dnet, 3, 0.0) == [3, 4]

        # Backward reachability is the exact dual
        back = backwardReachableSet(dnet, 4, 40.0)
        @test back == [1, 2, 3, 4]
        for v in back
            @test 4 in forwardReachableSet(dnet, v, 0.0)
        end
        @test backwardReachableSet(dnet, 5, 40.0) == [4, 5]  # via [5,8) spell
    end

    @testset "Point spells as instantaneous contacts" begin
        dnet = DynamicNetwork(3; observation_start=0.0, observation_end=10.0)
        for v in 1:3
            activate!(dnet, 0.0, 10.0; vertex=v)
        end
        activate!(dnet, 5.0, 5.0; edge=(1, 2))  # point contact at t=5

        arr, _ = earliestArrival(dnet, 1, 0.0)
        @test arr[2] == 5.0
        arr6, _ = earliestArrival(dnet, 1, 6.0)
        @test !haskey(arr6, 2)
    end

    @testset "Point measures are identity-stable" begin
        # Vertex 1 inactive at t=7 — vectors must still be full length,
        # indexed by original IDs
        dnet = DynamicNetwork(4; observation_start=0.0, observation_end=10.0)
        activate!(dnet, 0.0, 5.0; vertex=1)
        for v in 2:4
            activate!(dnet, 0.0, 10.0; vertex=v)
        end
        activate!(dnet, 0.0, 10.0; edge=(2, 3))
        activate!(dnet, 0.0, 10.0; edge=(3, 4))
        activate!(dnet, 0.0, 5.0; edge=(1, 2))

        deg = tDegree(dnet, 7.0)
        @test length(deg) == 4
        @test deg[1] == 0.0            # inactive vertex scores 0
        @test deg[3] == 2.0            # in + out over 2→3, 3→4

        bc = tBetweenness(dnet, 7.0)
        @test length(bc) == 4
        @test bc[3] > 0

        @test length(tCloseness(dnet, 7.0)) == 4
        @test length(tEigenvector(dnet, 7.0)) == 4
        @test length(tPagerank(dnet, 7.0)) == 4

        # Density over the full vertex set
        @test tDensity(dnet, 7.0) ≈ 2 / 12
        @test tTransitivity(dnet, 7.0) >= 0.0
    end

    @testset "Reciprocity" begin
        dnet = DynamicNetwork(3; observation_start=0.0, observation_end=10.0)
        for v in 1:3
            activate!(dnet, 0.0, 10.0; vertex=v)
        end
        activate!(dnet, 0.0, 10.0; edge=(1, 2))
        activate!(dnet, 0.0, 10.0; edge=(2, 1))
        activate!(dnet, 0.0, 10.0; edge=(2, 3))
        @test tReciprocity(dnet, 5.0) ≈ 2 / 3
    end

    @testset "Durations" begin
        dnet = DynamicNetwork(3; observation_start=0.0, observation_end=20.0)
        activate!(dnet, 0.0, 5.0; edge=(1, 2))
        activate!(dnet, 10.0, 12.0; edge=(1, 2))  # second spell, same edge
        activate!(dnet, 0.0, 10.0; edge=(2, 3))

        # Per-spell (tsna::edgeDuration semantics): (5 + 2 + 10)/3
        @test tEdgeDuration(dnet) ≈ 17 / 3
        # Per-edge totals: (7 + 10)/2
        @test tEdgeDuration(dnet; mode=:total) ≈ 8.5
        @test tEdgeDuration(dnet; aggregate=:total) ≈ 17.0
        @test length(tEdgeDuration(dnet; aggregate=:all)) == 3

        activate!(dnet, 0.0, 20.0; vertex=1)
        @test tVertexDuration(dnet) ≈ 20.0
    end

    @testset "Formation/dissolution events and turnover" begin
        dnet = DynamicNetwork(4; observation_start=0.0, observation_end=30.0)
        activate!(dnet, 0.0, 8.0; edge=(1, 2))     # forms at 0, dissolves at 8
        activate!(dnet, 12.0, 18.0; edge=(1, 2))   # re-forms at 12, dissolves 18
        activate!(dnet, 5.0, 25.0; edge=(2, 3))    # forms at 5, dissolves 25
        add_spell!(dnet, Spell(20.0, 30.0; terminus_censored=true); edge=(3, 4))

        @test tEdgeFormation(dnet, 0.0, 10.0) == 2   # onsets at 0 and 5
        @test tEdgeFormation(dnet, 10.0, 30.0) == 2  # 12 and 20
        @test tEdgeDissolution(dnet, 0.0, 10.0) == 1   # terminus 8
        @test tEdgeDissolution(dnet, 10.0, 30.0) == 2  # 18 and 25
        # Right-censored spell terminus is not a dissolution event
        @test tEdgeDissolution(dnet, 0.0, 31.0) == 3

        # A formation+dissolution INSIDE one window is still counted
        # (point sampling would have missed it)
        @test tEdgeFormation(dnet, 10.0, 20.0) == 1
        @test tEdgeDissolution(dnet, 10.0, 20.0) == 1

        windows = tTurnover(dnet, 10.0)
        @test length(windows) == 3
        # Consistent shape on every element
        @test all(haskey(pairs(w), :n_formations) for w in windows)
        @test windows[1].n_formations == 2
        @test windows[1].formation_rate ≈ 0.2
    end

    @testset "tieDecay" begin
        dnet = DynamicNetwork(3; observation_start=0.0, observation_end=10.0)
        activate!(dnet, 0.0, 10.0; edge=(1, 2))   # active at the end
        activate!(dnet, 0.0, 5.0; edge=(2, 3))    # ended 5 before the end

        w = tieDecay(dnet; rate=0.1)
        @test w[(1, 2)] ≈ 1.0
        @test w[(2, 3)] ≈ exp(-0.5)

        wl = tieDecay(dnet; method=:linear, rate=0.1)
        @test wl[(2, 3)] ≈ 0.5
        @test_throws ArgumentError tieDecay(dnet; method=:bogus)
    end

    @testset "Contact sequences" begin
        dnet = chain_fixture()
        cs = as_contact_sequence(dnet)
        @test length(cs) == 4
        @test issorted([c.time for c in cs])
        @test cs.n_vertices == 5
    end

    @testset "tSnaStats and windowSnaStats" begin
        dnet = chain_fixture()
        stats = tSnaStats(dnet, [1.0, 25.0])
        @test length(stats) == 2
        @test stats[1].density > stats[2].density  # 2 edges vs 1
        @test all(isfinite, (stats[1].mean_degree, stats[1].transitivity))

        ws = windowSnaStats(dnet, 20.0)
        @test length(ws) == 2
        @test_throws ArgumentError tSnaStats(dnet, [1.0]; measures=[:bogus])
    end

    @testset "tAggregate" begin
        dnet = DynamicNetwork(3; observation_start=0.0, observation_end=10.0)
        activate!(dnet, 0.0, 10.0; edge=(1, 2))
        activate!(dnet, 2.0, 6.0; edge=(2, 3))

        u = tAggregate(dnet)
        @test ne(u) == 2

        inter = tAggregate(dnet; method=:intersection)
        @test has_edge(inter, 1, 2)
        @test !has_edge(inter, 2, 3)

        # This used to throw a MethodError (wrong argument order)
        w = tAggregate(dnet; method=:weighted)
        @test get_edge_attribute(w, :weight, 1, 2) ≈ 10.0
        @test get_edge_attribute(w, :weight, 2, 3) ≈ 4.0

        @test_throws ArgumentError tAggregate(dnet; method=:bogus)
    end
end
