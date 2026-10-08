using TSNA
using DynamicNetworks
using NetworkCore
using Random
using Test
using Dates
using Aqua
import SNA

# Reference implementation of earliest arrival (the pre-heap linear
# min-scan version with a per-call contact index), kept verbatim for
# exact-equivalence testing of the optimized search.
function _reference_earliest_arrival(dnet::DynamicNetwork{T, Time}, source::T,
                                     start_time;
                                     end_time=TSNA._path_end(dnet)) where {T, Time}
    start_time = convert(Time, start_time)
    end_time = convert(Time, end_time)
    out = TSNA._out_contacts(dnet)

    arrival = Dict{T, Time}(source => start_time)
    parent = Dict{T, Tuple{T, Time}}()
    settled = Set{T}()

    while true
        v = nothing
        best = nothing
        for (u, t) in arrival
            u in settled && continue
            if isnothing(best) || t < best
                v, best = u, t
            end
        end
        isnothing(v) && break
        push!(settled, v)
        t = arrival[v]

        for (w, spell) in get(out, v, Tuple{T, Spell{Time}}[])
            w in settled && continue
            depart = TSNA._board_time(spell, t, end_time)
            isnothing(depart) && continue
            if !haskey(arrival, w) || depart < arrival[w]
                arrival[w] = depart
                parent[w] = (v, depart)
            end
        end
    end

    return arrival, parent
end

# Random temporal network: n vertices, ~m spells with random windows,
# including occasional point spells and repeated edges
function random_temporal_network(rng, n, m; directed=true)
    dnet = DynamicNetwork(n; observation_start=0.0, observation_end=100.0,
                          directed=directed)
    for v in 1:n
        activate!(dnet, 0.0, 100.0; vertex=v)
    end
    for _ in 1:m
        i, j = rand(rng, 1:n), rand(rng, 1:n)
        i == j && continue
        onset = round(100 * rand(rng); digits=2)
        if rand(rng) < 0.15
            activate!(dnet, onset, onset; edge=(i, j))  # point contact
        else
            terminus = min(100.0, onset + round(30 * rand(rng); digits=2))
            activate!(dnet, onset, terminus; edge=(i, j))
        end
    end
    return dnet
end

# Temporal chain: 1→2 active [0,10), 2→3 active [0,10) (equal onsets),
# 3→4 active [20,30), 4→5 active [5,8) (closes before walker can arrive)
# `T` is the vertex id type; the path testsets run with Int and Int32 ids and
# call every function with literal (Int) vertex ids.
function chain_fixture(::Type{T}=Int) where T
    dnet = DynamicNetwork{T, Float64}(5; observation_start=0.0, observation_end=40.0)
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
    @testset "Missing dyads: all statistical and path entry points" begin
        dnet = chain_fixture()
        queries = [
            (d; kw...) -> t_degree(d, 5.0; kw...),
            (d; kw...) -> t_betweenness(d, 5.0; kw...),
            (d; kw...) -> t_closeness(d, 5.0; kw...),
            (d; kw...) -> t_eigenvector(d, 5.0; kw...),
            (d; kw...) -> t_pagerank(d, 5.0; kw...),
            (d; kw...) -> t_density(d, 5.0; kw...),
            (d; kw...) -> t_reciprocity(d, 5.0; kw...),
            (d; kw...) -> t_transitivity(d, 5.0; kw...),
            (d; kw...) -> t_edge_duration(d; kw...),
            (d; kw...) -> t_vertex_duration(d; kw...),
            (d; kw...) -> t_edge_formation(d, 0.0, 10.0; kw...),
            (d; kw...) -> t_edge_dissolution(d, 0.0, 10.0; kw...),
            (d; kw...) -> t_edge_persistence(d, 10.0; kw...),
            (d; kw...) -> t_turnover(d, 10.0; kw...),
            (d; kw...) -> tie_decay(d; kw...),
            (d; kw...) -> t_sna_stats(d, [0.0, 5.0]; kw...),
            (d; kw...) -> window_sna_stats(d, 10.0; kw...),
            (d; kw...) -> earliest_arrival(d, 1, 0.0; kw...),
            (d; kw...) -> earliest_arrival_all(d, 0.0; kw...),
            (d; kw...) -> temporal_distance(d, 1, 3, 0.0; kw...),
            (d; kw...) -> temporal_distance_matrix(d, 0.0; kw...),
            (d; kw...) -> reachability_matrix(d, 0.0; kw...),
            (d; kw...) -> forward_reachable_set(d, 1, 0.0; kw...),
            (d; kw...) -> backward_reachable_set(d, 3, 40.0; kw...),
            (d; kw...) -> length(as_contact_sequence(d; kw...)),
        ]
        for (i, j) in ((1, 2), (1, 5))  # present and absent face values
            for query in queries
                baseline = query(dnet)
                set_missing_dyad!(dnet.network, i, j)
                @test_throws ArgumentError query(dnet)
                @test isequal(query(dnet; missing=:face), baseline)
                @test_throws ArgumentError query(dnet; missing=:invent)
                clear_missing_dyads!(dnet.network)
            end
        end
        set_missing_dyad!(dnet.network, 1, 2)
        @test_throws ArgumentError temporal_path(dnet, 1, 3, 0.0)
        @test temporal_path(dnet, 1, 3, 0.0; missing=:face).vertices == [1, 2, 3]
        @test_throws ArgumentError earliest_arrival!(TemporalPathWorkspace{Int,Float64}(), dnet, 1, 0.0)
        @test n_missing_dyads(t_aggregate(dnet)) == 1
        @test NetworkCore.missing_policies(t_density) == (:error, :face)
    end

    @testset "Active vertex denominators and explicit legacy policies" begin
        dnet = DynamicNetwork(5; observation_start=0.0, observation_end=10.0)
        for v in 1:4
            activate!(dnet, 0.0, 10.0; vertex=v)
        end
        deactivate!(dnet, -Inf, Inf; vertex=5)   # never present (no spells = active)
        activate!(dnet, 0.0, 10.0; edge=(1, 2))
        activate!(dnet, 0.0, 10.0; edge=(2, 1))
        activate!(dnet, 0.0, 10.0; edge=(2, 3))
        activate!(dnet, 0.0, 10.0; edge=(3, 4))
        @test t_density(dnet, 1.0) == 4/12
        @test t_density(dnet, 1.0; active_only=false) == 4/20
        @test t_reciprocity(dnet, 1.0) == 4/6 # one mutual and three null dyads
        @test t_reciprocity(dnet, 1.0; method=:edgewise) == 2/4
        @test length(t_degree(dnet, 1.0)) == 5
        @test isnan(t_degree(dnet, 1.0)[5])
        @test t_degree(dnet, 1.0; active_only=false)[5] == 0.0
        @test t_sna_stats(dnet, [1.0])[1].density == t_density(dnet, 1.0)
        # Even a masked dyad involving an inactive vertex must be acknowledged.
        set_missing_dyad!(dnet.network, 1, 5)
        @test_throws ArgumentError t_density(dnet, 1.0)
        @test t_density(dnet, 1.0; missing=:face) == 4/12
    end

    @testset "DateTime durations and window validation" begin
        t0 = DateTime(2026, 1, 1)
        dnet = DynamicNetwork{Int,DateTime}(3; observation_start=t0,
                                          observation_end=t0 + Hour(3))
        for v in 1:3
            activate!(dnet, t0, t0 + Hour(3); vertex=v)
        end
        activate!(dnet, t0, t0 + Hour(1); edge=(1, 2))
        activate!(dnet, t0 + Minute(30), t0 + Hour(2); edge=(2, 3))
        @test t_edge_duration(dnet; aggregate=:total) == 9000
        @test t_vertex_duration(dnet; aggregate=:total) == 32400
        D = temporal_distance_matrix(dnet, t0)
        @test D[1, 3] == Minute(30)
        @test D[3, 1] === nothing
        @test D[1, 1] == Millisecond(0)
        @test t_turnover(dnet, Hour(1))[1].formation_rate == 2/3600
        @test get_edge_attribute(t_aggregate(dnet; method=:weighted), :weight, 1, 2) == 3600
        for query in (t_edge_persistence, t_turnover, window_sna_stats)
            @test_throws ArgumentError query(dnet, Hour(0))
            @test_throws ArgumentError query(dnet, Hour(-1))
        end
        @test_throws ArgumentError earliest_arrival(dnet, 0, t0)
        @test_throws ArgumentError earliest_arrival(dnet, 1, t0; target=4)
        @test_throws ArgumentError earliest_arrival(dnet, 1, t0; end_time=t0 - Hour(1))
    end

    @testset "Earliest arrival: interval semantics ($T ids)" for T in (Int, Int32)
        dnet = chain_fixture(T)
        arrival, _ = earliest_arrival(dnet, 1, 0.0)
        @test arrival isa Dict{T, Float64}

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
        arr4, _ = earliest_arrival(dnet, 1, 4.0)
        @test arr4[2] == 4.0
        @test arr4[3] == 4.0

        # Spells active before start_time still count (old filter dropped
        # them entirely)
        arr9, _ = earliest_arrival(dnet, 1, 9.9)
        @test arr9[2] == 9.9

        # After the spell ends: no path
        arr11, _ = earliest_arrival(dnet, 1, 11.0)
        @test !haskey(arr11, 2)
    end

    @testset "Heap search: exact equivalence with reference" begin
        rng = Xoshiro(42)
        for trial in 1:20
            directed = isodd(trial)
            n = rand(rng, 3:12)
            dnet = random_temporal_network(rng, n, rand(rng, 5:40);
                                           directed=directed)
            source = rand(rng, 1:n)
            start_time = round(100 * rand(rng) * rand(rng); digits=2)

            arr_new, par_new = earliest_arrival(dnet, source, start_time)
            arr_ref, _ = _reference_earliest_arrival(dnet, source, start_time)

            # Arrival times must match the reference exactly
            @test arr_new == arr_ref

            # Parent maps may break ties differently, but every parent
            # chain must reconstruct a valid time-respecting path with the
            # reported arrival times
            for w in keys(arr_new)
                w == source && continue
                v = w
                while v != source
                    u, t = par_new[v]
                    @test t == arr_new[v]        # boarding = arrival label
                    @test arr_new[u] <= t        # time-respecting order
                    v = u
                end
            end

            # Windowed query equivalence
            end_time = start_time + 30.0
            arr_w, _ = earliest_arrival(dnet, source, start_time;
                                        end_time=end_time)
            arr_wr, _ = _reference_earliest_arrival(dnet, source, start_time;
                                                    end_time=end_time)
            @test arr_w == arr_wr

            # Early-exit variants agree with the full search
            for target in 1:n
                d = temporal_distance(dnet, source, target, start_time)
                if haskey(arr_ref, target)
                    @test d == arr_ref[target] - start_time
                else
                    @test d === nothing
                end
            end
        end
    end

    @testset "Contact index memoization" begin
        dnet = chain_fixture()
        idx1 = TSNA._contact_index(dnet)
        # Unchanged network: same index object is reused
        @test TSNA._contact_index(dnet) === idx1

        # Mutation invalidates the cache and results stay correct
        arr_before, _ = earliest_arrival(dnet, 1, 0.0)
        @test !haskey(arr_before, 5)
        activate!(dnet, 25.0, 35.0; edge=(4, 5))
        @test TSNA._contact_index(dnet) !== idx1
        arr_after, _ = earliest_arrival(dnet, 1, 0.0)
        @test arr_after[5] == 25.0

        # deactivate! invalidates too
        deactivate!(dnet, 0.0, 40.0; edge=(4, 5))
        arr_gone, _ = earliest_arrival(dnet, 1, 0.0)
        @test !haskey(arr_gone, 5)
    end

    @testset "camelCase names of earlier versions are gone" begin
        # They were advertised as "R tsna-style", but most are not tsna
        # functions and the five that share a tsna name (tPath, tDegree,
        # tEdgeFormation, tEdgeDissolution, tSnaStats) mean something else
        # there. None was ever released; the documentation's rename table maps
        # each to its snake_case name.
        for nm in (:tPath, :temporalDistance, :earliestArrival, :earliestArrivalAll,
                   :forwardReachableSet, :backwardReachableSet, :temporalPath,
                   :shortest_temporal_path, :shortestTemporalPath, :tDegree,
                   :tBetweenness, :tCloseness, :tEigenvector, :tPagerank, :tDensity,
                   :tReciprocity, :tTransitivity, :tEdgeDuration, :tVertexDuration,
                   :tEdgeFormation, :tEdgeDissolution, :tEdgePersistence, :tTurnover,
                   :tieDecay, :tSnaStats, :windowSnaStats, :tAggregate, :tEdgeFormationAt)
            @test !isdefined(TSNA, nm)
        end
        @test isempty([nm for nm in names(TSNA) if Base.isdeprecated(TSNA, nm)])
    end

    @testset "temporal_distance and paths ($T ids)" for T in (Int, Int32)
        dnet = chain_fixture(T)

        @test temporal_distance(dnet, 1, 3, 0.0) == 0.0
        @test temporal_distance(dnet, 1, 4, 0.0) == 20.0
        @test temporal_distance(dnet, 1, 5, 0.0) === nothing   # unreachable
        @test temporal_distance(dnet, 1, 4, 25.0) === nothing  # window too late

        p = temporal_path(dnet, 1, 4, 0.0)
        @test p isa TemporalPath{T, Float64}
        @test p.vertices == [1, 2, 3, 4]
        @test p.times == [0.0, 0.0, 20.0]
        @test issorted(p.times)
        @test path_duration(p) == 20.0
        @test temporal_path(dnet, 1, 5, 0.0) === nothing

        # Single-vertex path has zero duration
        p1 = temporal_path(dnet, 1, 1, 0.0)
        @test length(p1) == 0
        @test path_duration(p1) == 0.0

        # Any Integer vertex id is accepted (a DynamicNetwork{Int32} used to
        # need Int32(1)); out-of-range ids are refused by name, never by an
        # InexactError from the conversion.
        @test temporal_distance(dnet, Int8(1), UInt(4), 0.0) == 20.0
        @test temporal_path(dnet, big(1), 4, 0.0).vertices == [1, 2, 3, 4]
        for bad in (0, 6, typemax(Int64))
            @test_throws ArgumentError temporal_distance(dnet, bad, 1, 0.0)
            @test_throws ArgumentError temporal_distance(dnet, 1, bad, 0.0)
            @test_throws ArgumentError temporal_path(dnet, bad, 1, 0.0)
            @test_throws ArgumentError earliest_arrival(dnet, bad, 0.0)
            @test_throws ArgumentError earliest_arrival(dnet, 1, 0.0; target=bad)
            @test_throws ArgumentError forward_reachable_set(dnet, bad, 0.0)
            @test_throws ArgumentError backward_reachable_set(dnet, bad, 40.0)
        end
        ws = TemporalPathWorkspace{T, Float64}()
        @test earliest_arrival!(ws, dnet, 1, 0.0; target=4)[1][4] == 20.0
        @test earliest_arrival_all(dnet, 0.0; sources=[1, 3])[3][4] == 20.0
    end

    @testset "Reachable sets and duality ($T ids)" for T in (Int, Int32)
        dnet = chain_fixture(T)

        @test forward_reachable_set(dnet, 1, 0.0) == [1, 2, 3, 4]
        @test forward_reachable_set(dnet, 1, 15.0) == [1]
        @test forward_reachable_set(dnet, 3, 0.0) == [3, 4]

        # Backward reachability is the exact dual
        back = backward_reachable_set(dnet, 4, 40.0)
        @test back == [1, 2, 3, 4]
        for v in back
            @test 4 in forward_reachable_set(dnet, v, 0.0)
        end
        @test backward_reachable_set(dnet, 5, 40.0) == [4, 5]  # via [5,8) spell
    end

    @testset "Point spells as instantaneous contacts ($T ids)" for T in (Int, Int32)
        dnet = DynamicNetwork{T, Float64}(3; observation_start=0.0, observation_end=10.0)
        for v in 1:3
            activate!(dnet, 0.0, 10.0; vertex=v)
        end
        activate!(dnet, 5.0, 5.0; edge=(1, 2))  # point contact at t=5

        arr, _ = earliest_arrival(dnet, 1, 0.0)
        @test arr[2] == 5.0
        arr6, _ = earliest_arrival(dnet, 1, 6.0)
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

        deg = t_degree(dnet, 7.0)
        @test length(deg) == 4
        @test isnan(deg[1])            # absent actor: NaN, not a score
        @test t_degree(dnet, 7.0; active_only=false)[1] == 0.0
        @test deg[3] == 2.0            # in + out over 2→3, 3→4

        bc = t_betweenness(dnet, 7.0)
        @test length(bc) == 4
        @test bc[3] > 0

        @test length(t_closeness(dnet, 7.0)) == 4
        @test all(isnan(x[1]) for x in (t_closeness(dnet, 7.0), t_eigenvector(dnet, 7.0),
                                         t_pagerank(dnet, 7.0), bc))
        @test length(t_eigenvector(dnet, 7.0)) == 4
        @test length(t_pagerank(dnet, 7.0)) == 4

        # Default denominator uses active vertices; the legacy universe is explicit.
        @test t_density(dnet, 7.0) ≈ 2 / 6
        @test t_density(dnet, 7.0; active_only=false) ≈ 2 / 12
        @test t_transitivity(dnet, 7.0) >= 0.0
    end

    @testset "Reciprocity" begin
        dnet = DynamicNetwork(3; observation_start=0.0, observation_end=10.0)
        for v in 1:3
            activate!(dnet, 0.0, 10.0; vertex=v)
        end
        activate!(dnet, 0.0, 10.0; edge=(1, 2))
        activate!(dnet, 0.0, 10.0; edge=(2, 1))
        activate!(dnet, 0.0, 10.0; edge=(2, 3))
        @test t_reciprocity(dnet, 5.0) ≈ 2 / 3
    end

    @testset "Durations" begin
        dnet = DynamicNetwork(3; observation_start=0.0, observation_end=20.0)
        activate!(dnet, 0.0, 5.0; edge=(1, 2))
        activate!(dnet, 10.0, 12.0; edge=(1, 2))  # second spell, same edge
        activate!(dnet, 0.0, 10.0; edge=(2, 3))

        # Default = tsna::edgeDuration(nd): per-edge totals, edges(net) order
        @test t_edge_duration(dnet) == [7.0, 10.0]
        @test t_edge_duration(dnet; aggregate=:mean) ≈ 8.5
        # Per spell (tsna subject = "spells"): (5 + 2 + 10)/3
        @test t_edge_duration(dnet; mode=:spell, aggregate=:mean) ≈ 17 / 3
        @test t_edge_duration(dnet; aggregate=:total) ≈ 17.0
        @test length(t_edge_duration(dnet; mode=:spell)) == 3
        @test_throws ArgumentError t_edge_duration(dnet; mode=:bogus)
        @test_throws ArgumentError t_edge_duration(dnet; aggregate=:bogus)

        activate!(dnet, 0.0, 10.0; vertex=1)
        # vertices 2 and 3 have no spells: active over the whole window
        @test t_vertex_duration(dnet) == [10.0, 20.0, 20.0]
    end

    @testset "Formation/dissolution events and turnover" begin
        dnet = DynamicNetwork(4; observation_start=0.0, observation_end=30.0)
        activate!(dnet, 0.0, 8.0; edge=(1, 2))     # forms at 0, dissolves at 8
        activate!(dnet, 12.0, 18.0; edge=(1, 2))   # re-forms at 12, dissolves 18
        activate!(dnet, 5.0, 25.0; edge=(2, 3))    # forms at 5, dissolves 25
        add_spell!(dnet, Spell(20.0, 30.0; terminus_censored=true); edge=(3, 4))

        @test t_edge_formation(dnet, 0.0, 10.0) == 2   # onsets at 0 and 5
        @test t_edge_formation(dnet, 10.0, 30.0) == 2  # 12 and 20
        @test t_edge_dissolution(dnet, 0.0, 10.0) == 1   # terminus 8
        @test t_edge_dissolution(dnet, 10.0, 30.0) == 2  # 18 and 25
        # Right-censored spell terminus is not a dissolution event
        @test t_edge_dissolution(dnet, 0.0, 31.0) == 3

        # A formation+dissolution INSIDE one window is still counted
        # (point sampling would have missed it)
        @test t_edge_formation(dnet, 10.0, 20.0) == 1
        @test t_edge_dissolution(dnet, 10.0, 20.0) == 1

        windows = t_turnover(dnet, 10.0)
        @test length(windows) == 3
        # Consistent shape on every element
        @test all(haskey(pairs(w), :n_formations) for w in windows)
        @test windows[1].n_formations == 2
        @test windows[1].formation_rate ≈ 0.2
    end

    @testset "t_edge_persistence" begin
        dnet = DynamicNetwork(4; observation_start=0.0, observation_end=40.0)
        activate!(dnet, 0.0, 40.0; edge=(1, 2))   # active at every window start
        activate!(dnet, 0.0, 15.0; edge=(2, 3))   # active at 0 and 10 only
        activate!(dnet, 25.0, 40.0; edge=(3, 4))  # appears at 30 (never in prev)

        # Window starts 0,10,20,30. Pairs: (0→10): 2/2 persist,
        # (10→20): 1/2, (20→30): 1/1 → pooled 4/5
        @test t_edge_persistence(dnet, 10.0) ≈ 4 / 5

        # One big window pair: edges at 0 are {12,23}; at 20 only 12 → 1/2
        @test t_edge_persistence(dnet, 20.0) ≈ 1 / 2

        # Fewer than two windows: undefined
        @test isnan(t_edge_persistence(dnet, 40.0))
        @test isnan(t_edge_persistence(dnet, 100.0))

        # No edges to track: undefined
        empty_net = DynamicNetwork(3; observation_start=0.0, observation_end=40.0)
        @test isnan(t_edge_persistence(empty_net, 10.0))

        # Perfectly stable network
        stable = DynamicNetwork(3; observation_start=0.0, observation_end=40.0)
        activate!(stable, 0.0, 40.0; edge=(1, 2))
        activate!(stable, 0.0, 40.0; edge=(2, 3))
        @test t_edge_persistence(stable, 10.0) == 1.0
    end

    @testset "tie_decay" begin
        dnet = DynamicNetwork(3; observation_start=0.0, observation_end=10.0)
        activate!(dnet, 0.0, 10.0; edge=(1, 2))   # active at the end
        activate!(dnet, 0.0, 5.0; edge=(2, 3))    # ended 5 before the end

        w = tie_decay(dnet; rate=0.1)
        @test w[(1, 2)] ≈ 1.0
        @test w[(2, 3)] ≈ exp(-0.5)

        wl = tie_decay(dnet; method=:linear, rate=0.1)
        @test wl[(2, 3)] ≈ 0.5
        @test_throws ArgumentError tie_decay(dnet; method=:bogus)
    end

    @testset "Contact sequences" begin
        dnet = chain_fixture()
        cs = as_contact_sequence(dnet)
        @test length(cs) == 4
        @test issorted([c.time for c in cs])
        @test cs.n_vertices == 5
    end

    @testset "Contact sequences from unbounded spells" begin
        # A no-record edge cut by deactivate! keeps the axis extremes as its
        # outer bounds. Subtracting them gave a contact at -Inf lasting Inf on
        # a float axis, a negative duration on an Int axis and a meaningless
        # one on a DateTime axis.
        t0 = DateTime(2020)
        axes = ((Float64, 3.0, 4.0, Inf),
                (Int, 3, 4, typemax(Int)),
                (Int32, Int32(3), Int32(4), typemax(Int32)),
                (DateTime, t0 + Day(3), t0 + Day(4), Millisecond(typemax(Int64))),
                (Date, Date(2020, 1, 4), Date(2020, 1, 5), Day(typemax(Int64))))
        for (Time, a, b, open_dur) in axes, T in (Int, Int32)
            d = DynamicNetwork{T, Time}(3)
            add_edge!(d.network, 1, 2)
            deactivate!(d, a, b; edge=(1, 2))      # leaves [min, a) and [b, max)
            cs, rep = as_contact_sequence(d; report=true)
            @test [(c.source, c.target, c.time, c.duration) for c in cs] == [(1, 2, b, open_dur)]
            @test :unbounded_onsets in dropped_fields(rep)
            @test :unbounded_termini in dropped_fields(rep)
            # A finite spell keeps its measured duration, and nothing is reported
            f = DynamicNetwork{T, Time}(3)
            activate!(f, a, b; edge=(2, 3))
            cs, rep = as_contact_sequence(f; report=true)
            @test only(cs).duration == spell_duration(Spell(a, b))
            @test !(:unbounded_onsets in dropped_fields(rep))
            @test !(:unbounded_termini in dropped_fields(rep))
        end
        # A point spell at the end of the axis is an instant, not an open contact
        p = DynamicNetwork(2)
        activate!(p, Inf, Inf; edge=(1, 2))
        cs, rep = as_contact_sequence(p; report=true)
        @test only(cs).duration == 0.0 && !(:unbounded_termini in dropped_fields(rep))
    end

    @testset "A window reaching the axis extremes is unbounded" begin
        # On integer and calendar axes typemax stands for Inf; a window
        # reaching it was tiled as if finite (here into four huge windows;
        # with a unit window size, practically forever).
        t0 = DateTime(2020)
        for (TA, lo, hi, w) in ((Int, 0, typemax(Int), typemax(Int) ÷ 4),
                                (Int32, Int32(0), typemax(Int32), typemax(Int32) ÷ Int32(4)),
                                (DateTime, t0, typemax(DateTime), (typemax(DateTime) - t0) ÷ 4),
                                (Float64, 0.0, Inf, 1.0))
            q = DynamicNetwork{Int, TA}(3; observation_start=lo, observation_end=hi)
            activate!(q, lo, hi; edge=(1, 2))
            err = try t_turnover(q, w); nothing catch e; e end
            @test err isa ArgumentError && occursin("finite observation period", err.msg)
            @test_throws ArgumentError window_sna_stats(q, w)
            @test_throws ArgumentError t_edge_persistence(q, w)
        end
    end

    @testset "t_sna_stats and window_sna_stats" begin
        dnet = chain_fixture()
        stats = t_sna_stats(dnet, [1.0, 25.0])
        @test length(stats) == 2
        @test stats[1].density > stats[2].density  # 2 edges vs 1
        @test all(isfinite, (stats[1].mean_degree, stats[1].transitivity))

        ws = window_sna_stats(dnet, 20.0)
        @test length(ws) == 2
        @test_throws ArgumentError t_sna_stats(dnet, [1.0]; measures=[:bogus])
    end

    @testset "t_aggregate" begin
        dnet = DynamicNetwork(3; observation_start=0.0, observation_end=10.0)
        activate!(dnet, 0.0, 10.0; edge=(1, 2))
        activate!(dnet, 2.0, 6.0; edge=(2, 3))

        u = t_aggregate(dnet)
        @test ne(u) == 2

        inter = t_aggregate(dnet; method=:intersection)
        @test has_edge(inter, 1, 2)
        @test !has_edge(inter, 2, 3)

        # This used to throw a MethodError (wrong argument order)
        w = t_aggregate(dnet; method=:weighted)
        @test get_edge_attribute(w, :weight, 1, 2) ≈ 10.0
        @test get_edge_attribute(w, :weight, 2, 3) ≈ 4.0

        @test_throws ArgumentError t_aggregate(dnet; method=:bogus)
    end

    # =========================================================================
    # Conversion invariants (see docs/src/guide/conversion_invariants.md)
    #
    # A ContactSequence has no slot for an unobserved dyad, so a masked
    # DynamicNetwork is REJECTED rather than flattened into contacts that read
    # as observed. t_aggregate inherits network_collapse's invariants, so the
    # mask survives it.
    # =========================================================================
    @testset "Conversion invariants: as_contact_sequence" begin
        for directed in (true, false)
            dnet = DynamicNetwork(4; observation_start=0.0, observation_end=10.0,
                                  directed=directed)
            activate_vertices!(dnet, [1, 2, 3, 4], 0.0, 10.0)
            activate!(dnet, 0.0, 4.0; edge=(1, 2))
            activate!(dnet, 3.0, 6.0; edge=(1, 2))     # overlapping: one activity
            activate!(dnet, 7.0, 7.0; edge=(2, 3))     # point spell
            add_spell!(dnet, Spell(8.0, 9.0); edge=(3, 4), merge=false)
            add_spell!(dnet, Spell(9.0, 9.5); edge=(3, 4), merge=false)   # adjacent

            cs, rep = as_contact_sequence(dnet; report=true)

            # Preserved: one contact per contiguous activity (spells merged
            # per edge, as R's activate.edges stores them), the vertex count,
            # directedness, onsets and durations.
            @test length(cs) == 3
            @test cs.n_vertices == 4
            @test cs.directed == directed
            @test [c.time for c in cs] == [0.0, 7.0, 8.0]
            @test [c.duration for c in cs] == [6.0, 0.0, 1.5]   # point spell = 0
            @test !(:default_active_edges in dropped_fields(rep))

            # A base edge with no spell record has no finite onset: reported
            add_edge!(dnet.network, 1, 4)
            cs2, rep2 = as_contact_sequence(dnet; report=true)
            @test length(cs2) == 3
            @test :default_active_edges in dropped_fields(rep2)

            # Dropped by nature, and named.
            @test !is_lossless(rep)
            @test :spell_censoring in dropped_fields(rep)
            @test :vertex_spells in dropped_fields(rep)
            @test :observation_period in dropped_fields(rep)
            @test !(:missing_dyads in dropped_fields(rep))   # no mask here
        end
    end

    @testset "Conversion invariants: as_contact_sequence rejects a masked network" begin
        dnet = DynamicNetwork(4; observation_start=0.0, observation_end=10.0)
        activate_vertices!(dnet, [1, 2, 3, 4], 0.0, 10.0)
        activate!(dnet, 0.0, 5.0; edge=(1, 2))
        activate!(dnet, 0.0, 5.0; edge=(2, 3))
        set_missing_dyad!(dnet.network, 2, 3)     # masked, PRESENT face value
        set_missing_dyad!(dnet.network, 3, 4)     # masked, ABSENT face value

        # Default policy: refuse. A contact cannot say "unobserved".
        @test_throws ArgumentError as_contact_sequence(dnet)
        @test_throws ArgumentError as_contact_sequence(dnet; missing=:error)
        @test_throws ArgumentError as_contact_sequence(dnet; missing=:bogus)

        # Explicit opt-in, and the report says what that cost.
        cs, rep = as_contact_sequence(dnet; missing=:face, report=true)
        @test length(cs) == 2
        @test :missing_dyads in dropped_fields(rep)

        # Clearing the mask makes the network observed, and it converts.
        clear_missing_dyads!(dnet.network)
        @test length(as_contact_sequence(dnet)) == 2
    end

    @testset "Conversion invariants: t_aggregate carries the mask" begin
        for directed in (true, false)
            net = network(5; directed=directed, loops=true)
            add_edge!(net, 1, 2)
            add_edge!(net, 2, 3)
            add_edge!(net, 3, 3)
            set_vertex_attribute!(net, :grp, 1, "a")
            set_network_attribute!(net, :title, "demo")
            set_missing_dyad!(net, 2, 3)          # PRESENT face value
            set_missing_dyad!(net, 4, 5)          # ABSENT face value
            dnet = as_dynamic_network(net; onset=0.0, terminus=10.0)

            for method in (:union, :intersection, :weighted)
                agg, rep = t_aggregate(dnet; method=method, report=true)
                @test is_directed(agg) == directed
                @test agg.loops
                @test has_edge(agg, 3, 3)                    # self-loop survives
                @test get_vertex_attribute(agg, :grp, 1) == "a"
                @test get_network_attribute(agg, :title) == "demo"
                # THE regression: the mask must not aggregate away to zero.
                @test n_missing_dyads(agg) == 2
                @test is_missing_dyad(agg, 2, 3)
                @test is_missing_dyad(agg, 4, 5)
                @test has_edge(agg, 2, 3) && !has_edge(agg, 4, 5)
                @test :spells in dropped_fields(rep)
            end

            @test get_edge_attribute(t_aggregate(dnet; method=:weighted),
                                     :weight, 1, 2) == 10.0
        end
    end

    @testset "Batch temporal paths (TSNA.jl#1)" begin
        # An all-source analysis runs one search per vertex. Each search used to
        # allocate its own arrival/parent/settled/heap; a workspace reuses them.
        # The results must be IDENTICAL — this is a pure allocation change.

        rng = Random.Xoshiro(11)
        n = 25
        dn = DynamicNetwork(n; observation_start=0.0, observation_end=100.0,
                            directed=true)
        for v in 1:n
            activate!(dn, 0.0, 100.0; vertex=v)
        end
        for _ in 1:120
            i, j = rand(rng, 1:n), rand(rng, 1:n)
            i == j && continue
            onset = round(80 * rand(rng); digits=2)
            activate!(dn, onset, min(100.0, onset + round(20 * rand(rng); digits=2));
                      edge=(i, j))
        end

        @testset "batch == per-source loop, exactly" begin
            batch = earliest_arrival_all(dn, 0.0)
            for v in 1:n
                single, _ = earliest_arrival(dn, v, 0.0)
                @test batch[v] == single          # same keys AND same times
            end
            @test length(batch) == n
        end

        @testset "the workspace is genuinely reused" begin
            ws = TemporalPathWorkspace{Int, Float64}()
            a1, _ = earliest_arrival!(ws, dn, 1, 0.0)
            snapshot = copy(a1)
            # The returned dict ALIASES the workspace: the next search overwrites
            # it. That is the documented contract, and the reason the batch entry
            # point copies before moving on.
            a2, _ = earliest_arrival!(ws, dn, 2, 0.0)
            @test a2 === a1                        # same container, reused
            ref2, _ = earliest_arrival(dn, 2, 0.0) # fresh search agrees
            @test a2 == ref2
            # ...and the reuse did not corrupt the second search with the first's
            # state (the bug a workspace invites)
            @test snapshot == earliest_arrival(dn, 1, 0.0)[1]
        end

        @testset "an early target break leaves a clean workspace" begin
            # The heap is NOT drained when a search stops early at `target`, so
            # the reset has to clear it or the next search inherits stale labels.
            ws = TemporalPathWorkspace{Int, Float64}()
            earliest_arrival!(ws, dn, 1, 0.0; target=3)     # may break early
            a, _ = earliest_arrival!(ws, dn, 5, 0.0)        # full search after it
            @test a == earliest_arrival(dn, 5, 0.0)[1]
        end

        @testset "distance and reachability matrices" begin
            D = temporal_distance_matrix(dn, 0.0)
            R = reachability_matrix(dn, 0.0)
            @test size(D) == (n, n) && size(R) == (n, n)
            for v in 1:n
                @test D[v, v] == 0.0
                @test R[v, v]
            end
            # Agree with the single-source functions they batch
            for i in 1:n, j in 1:n
                td = temporal_distance(dn, i, j, 0.0)
                @test D[i, j] == td
                @test R[i, j] == !isnothing(td) || (i == j)
            end
        end

        @testset "fewer allocations than the naive loop" begin
            # The point of the exercise. Warm up first: a cold call measures
            # compilation, not the algorithm.
            earliest_arrival_all(dn, 0.0)
            [earliest_arrival(dn, v, 0.0) for v in 1:n]

            batched = @allocated earliest_arrival_all(dn, 0.0)
            naive = @allocated [earliest_arrival(dn, v, 0.0) for v in 1:n]
            @test batched < naive
        end
    end
end

@testset "Empty temporal batches still validate time windows" begin
    empty = DynamicNetwork(0; observation_start=0.0, observation_end=10.0)
    @test_throws ArgumentError earliest_arrival_all(empty, 3.0; end_time=2.0)
    @test_throws ArgumentError temporal_distance_matrix(empty, 3.0; end_time=2.0)
    network = DynamicNetwork(2; observation_start=0.0, observation_end=10.0)
    @test_throws ArgumentError earliest_arrival_all(network, 3.0; sources=Int[], end_time=2.0)
    @test isempty(earliest_arrival_all(network, 2.0; sources=Int[], end_time=3.0))
end

@testset "R semantics" begin
    @testset "Per-vertex measures use the active network, NaN for absent actors" begin
        # 1 - 2 - 3 tied, actor 4 joins at t = 6.
        d = DynamicNetwork(4; directed=false, observation_start=0.0, observation_end=10.0)
        activate!(d, 0.0, 10.0; edge=(1, 2))
        activate!(d, 0.0, 10.0; edge=(2, 3))
        activate!(d, 6.0, 10.0; vertex=4)
        cl = t_closeness(d, 2.0)
        @test cl[1:3] ≈ [2/3, 1.0, 2/3]                      # was [0, 0, 0, 0]
        @test isnan(cl[4])
        @test t_betweenness(d, 2.0; normalized=true)[2] ≈ 1.0   # was 1/3
        @test t_degree(d, 2.0; normalized=true)[2] ≈ 1.0        # was 2/3
        @test isnan(t_pagerank(d, 2.0)[4]) && sum(t_pagerank(d, 2.0)[1:3]) ≈ 1
        # Once actor 4 is present it is scored (an isolate here)
        @test t_closeness(d, 7.0)[4] == 0.0
        # The legacy whole-universe computation is explicit
        @test all(iszero, t_closeness(d, 2.0; active_only=false))
        # Exactly the static measure on the active extract, scattered back
        snap = network_extract(d, 2.0)
        pid = [get_vertex_attribute(snap, :vertex_pid, v) for v in 1:nv(snap)]
        @test t_betweenness(d, 2.0)[pid] == SNA.betweenness(snap)
        # Nobody present: all NaN, no error
        empty = DynamicNetwork(3; observation_start=0.0, observation_end=10.0)
        for v in 1:3
            deactivate!(empty, -Inf, Inf; vertex=v)
        end
        @test all(isnan, t_degree(empty, 1.0))
    end

    @testset "Spells are merged per edge before lifetimes are counted" begin
        # [0,5) + [5,10) and [0,6) + [4,10): two ties, each active [0,10).
        for merge in (true, false)
            d = DynamicNetwork(3; observation_start=0.0, observation_end=20.0)
            add_spell!(d, Spell(0.0, 5.0); edge=(1, 2), merge=merge)
            add_spell!(d, Spell(5.0, 10.0); edge=(1, 2), merge=merge)
            add_spell!(d, Spell(0.0, 6.0); edge=(2, 3), merge=merge)
            add_spell!(d, Spell(4.0, 10.0); edge=(2, 3), merge=merge)
            @test t_edge_formation(d, 0.0, 20.0) == 2           # was 4
            @test t_edge_dissolution(d, 0.0, 20.0) == 2         # was 4
            @test t_edge_duration(d) == [10.0, 10.0]            # was [10, 12]
            @test t_edge_duration(d; mode=:spell) == [10.0, 10.0]
            @test sum(w.n_formations for w in t_turnover(d, 5.0)) == 2
            @test length(as_contact_sequence(d)) == 2
            @test get_edge_attribute(t_aggregate(d; method=:weighted), :weight, 2, 3) == 10.0
        end
        # Panel data, one activate! per wave: a persisting tie forms once
        d = DynamicNetwork(2; observation_start=0.0, observation_end=4.0)
        for t in 0:3
            activate!(d, t, t + 1; edge=(1, 2))
        end
        @test t_edge_formation(d, 0.0, 4.0) == 1
        @test t_edge_dissolution(d, 0.0, 5.0) == 1
    end

    @testset "Durations are clipped to the observation window" begin
        # R tsna::edgeDuration gives 5 8 4.
        d = DynamicNetwork(3; observation_start=0.0, observation_end=10.0)
        activate!(d, 0.0, 5.0; edge=(1, 2))
        activate!(d, 2.0, Inf; edge=(1, 3))
        activate!(d, 1.0, 3.0; edge=(2, 3))
        activate!(d, 6.0, 8.0; edge=(2, 3))
        @test t_edge_duration(d) == [5.0, 8.0, 4.0]          # was mean = Inf
        @test all(isfinite, t_edge_duration(d; mode=:spell))
        @test t_edge_duration(d; aggregate=:median) == 5.0
        # Spells before the window are truncated at its start
        activate!(d, -5.0, 1.0; edge=(3, 1))
        @test t_edge_duration(d)[end] == 1.0
        # No window given: nothing is clipped, as in tsna (an open spell lasts
        # Inf; a window derived from the spell bounds used to give [4, 2])
        d2 = DynamicNetwork(3)
        activate!(d2, 2.0, 6.0; edge=(1, 2))
        activate!(d2, 4.0, Inf; edge=(2, 3))
        @test t_edge_duration(d2) == [4.0, Inf]
    end

    @testset "Censored bounds are not formation/dissolution events" begin
        d = DynamicNetwork(3; observation_start=0.0, observation_end=10.0)
        add_spell!(d, Spell(0.0, 4.0; onset_censored=true); edge=(1, 2))
        add_spell!(d, Spell(2.0, 10.0; terminus_censored=true); edge=(2, 3))
        @test t_edge_formation(d, 0.0, 10.0) == 1             # was 2
        @test t_edge_dissolution(d, 0.0, 11.0) == 1
        @test t_edge_formation(d, 0.0, 10.0; include_censored=true) == 2
        @test t_edge_dissolution(d, 0.0, 11.0; include_censored=true) == 2
        w = t_turnover(d, 5.0)
        @test w[1].n_formations == 1 && w[1].formation_rate ≈ 0.2   # was 0.4
        # A spell starting before the window is left-censored by position
        activate!(d, -3.0, 6.0; edge=(3, 1))
        @test t_edge_formation(d, -5.0, 10.0) == 1
        @test t_edge_dissolution(d, 0.0, 10.0) == 2           # 4 and 6
        # Edges with no spell record are active throughout: censored both ways
        add_edge!(d.network, 1, 3)
        @test t_edge_formation(d, -5.0, 10.0) == 1
    end

    @testset "Edges with no spell record are active (R active.default)" begin
        d = DynamicNetwork(3; observation_start=0.0, observation_end=10.0)
        add_edge!(d.network, 1, 2)                 # no spells: always active
        activate!(d, 5.0, 6.0; edge=(2, 3))
        @test forward_reachable_set(d, 1, 0.0) == [1, 2, 3]   # via the default edge
        @test temporal_distance(d, 1, 3, 0.0) == 5.0
        @test t_edge_duration(d) == [10.0, 1.0]
        @test tie_decay(d)[(1, 2)] == 1.0
        @test t_degree(d, 1.0) == [1.0, 1.0, 0.0]
        # Adding a base edge directly invalidates the memoized contact index
        add_edge!(d.network, 3, 1)
        @test forward_reachable_set(d, 3, 0.0) == [1, 2, 3]
    end

    @testset "Keywords forwarded to SNA.jl's sna-named functions" begin
        # Every keyword TSNA documents must be one the SNA function accepts
        # (t_transitivity used to forward `type` to SNA.gtrans, which has none).
        d = DynamicNetwork(5; directed=false, observation_start=0.0, observation_end=10.0)
        activate_edges!(d, [(1, 2), (2, 3), (1, 3), (3, 4)], 0.0, 10.0)
        activate!(d, 6.0, 10.0; vertex=5)
        snap = network_extract(d, 5.0)
        @test t_transitivity(d, 5.0) == SNA.gtrans(snap) ≈ 0.6
        @test t_transitivity(d, 5.0; measure=:strong) == SNA.gtrans(snap; measure=:strong)
        @test t_transitivity(d, 5.0; measure=:weakcensus) == SNA.gtrans(snap; measure=:weakcensus)
        @test t_transitivity(d, 5.0; type=:global) == t_transitivity(d, 5.0)
        @test t_transitivity(d, 5.0; type=:average) ≈
              SNA.transitivity(snap; type=:average) ≈ (1 + 1 + 1/3) / 3
        @test t_transitivity(d, 5.0; type=:average, cmode=:weak) ≈
              SNA.transitivity(snap; type=:average, cmode=:weak)
        @test_throws ArgumentError t_transitivity(d, 5.0; type=:local)
        @test_throws ArgumentError t_transitivity(d, 5.0; measure=:bogus)
        @test t_transitivity(d, 5.0; active_only=false) isa Float64
        # per-vertex measures: each documented pass-through keyword
        @test t_degree(d, 5.0; rescale=true)[1:4] ≈ SNA.degreecent(snap; rescale=true)
        @test t_degree(d, 5.0; normalized=true, diag=false)[3] ≈ 1.0
        @test t_degree(d, 5.0; ignore_eval=true, attr=:weight)[3] == 3.0
        @test_throws ArgumentError t_degree(d, 5.0; mode=:bogus)
        @test t_betweenness(d, 5.0; cmode=:undirected, rescale=true)[1:4] ≈
              SNA.betweenness(snap; cmode=:undirected, rescale=true)
        @test t_closeness(d, 5.0; cmode=:suminvundir)[1:4] ≈ SNA.closeness(snap; cmode=:suminvundir)
        @test t_closeness(d, 5.0; rescale=true)[1:4] ≈ SNA.closeness(snap; rescale=true)
        @test t_eigenvector(d, 5.0; rescale=true)[1:4] ≈ SNA.evcent(snap; rescale=true)
        @test t_density(d, 5.0; diag=false) == SNA.gden(snap)
        for m in (:dyadic, :dyadic_nonnull, :edgewise)
            @test isequal(t_reciprocity(d, 5.0; method=m), SNA.grecip(snap; measure=m))
        end
        @test t_sna_stats(d, [5.0]; reciprocity_method=:edgewise)[1].transitivity ≈ 0.6
        # Graphs-style or retired SNA keywords are a MethodError, not ignored
        @test_throws MethodError t_closeness(d, 5.0; normalize=true)
        @test_throws MethodError t_transitivity(d, 5.0; kind=:local)
    end

    @testset "tied_duration and t_edge_density (tsna)" begin
        d = DynamicNetwork(3; observation_start=0.0, observation_end=10.0)
        activate!(d, 0.0, 4.0; edge=(1, 2))
        activate!(d, 2.0, Inf; edge=(1, 3))
        @test tied_duration(d) == [12.0, 0.0, 0.0]
        @test tied_duration(d; neighborhood=:in) == [0.0, 4.0, 8.0]
        @test tied_duration(d; neighborhood=:combined, mode=:counts) == [2.0, 1.0, 1.0]
        @test t_edge_density(d) ≈ 12 / 20
        @test t_edge_density(d; agg_unit=:dyad) ≈ 12 / 60
        @test t_edge_density(d; mode=:event) ≈ 2 / 20
        @test t_edge_density(DynamicNetwork(3)) == 0.0
        @test_throws ArgumentError t_edge_density(d; mode=:event, agg_unit=:dyad)
        @test_throws ArgumentError tied_duration(d; neighborhood=:bogus)
        set_missing_dyad!(d.network, 1, 2)
        @test_throws ArgumentError tied_duration(d)
        @test_throws ArgumentError t_edge_density(d)
        @test tied_duration(d; missing=:face) == [12.0, 0.0, 0.0]
    end

    @testset "No observation window: tsna's rules, no placeholder range" begin
        # Networks d, e and p are the golden fixture's fixed no-window cases
        # (test/fixtures/r/tsna_semantics.R, nowin_41..43), so their event,
        # duration, tiedDuration, tEdgeDensity and tPath values are R tsna
        # 0.3.6 output. The remaining values (tie_decay, t_aggregate,
        # window_sna_stats, t_edge_persistence, reachability, the refusals and
        # the integer-axis network q) are derived by hand from tsna's rules.
        # A tie that forms at the last change time and stays open: the
        # derived window [2, 6) used to drop its formation and duration.
        d = DynamicNetwork(3)
        activate!(d, 2.0, 6.0; edge=(1, 2))
        activate!(d, 6.0, Inf; edge=(2, 3))
        @test get_observation_period(d) === nothing
        @test [t_edge_formation(d, t, t + 1) for t in 0.0:8.0] == [0, 0, 1, 0, 0, 0, 1, 0, 0]
        @test [t_edge_dissolution(d, t, t + 1) for t in 0.0:8.0] == [0, 0, 0, 0, 0, 0, 1, 0, 0]
        @test [r.n_formations for r in t_turnover(d, 1.0)] == [1, 0, 0, 0, 1]   # tEdgeFormation(d)
        @test [r.window_start for r in t_turnover(d, 1.0)] == 2.0:6.0
        @test t_edge_duration(d) == [4.0, Inf]
        @test t_edge_duration(d; mode=:spell) == [4.0, Inf]
        @test t_vertex_duration(d) == [Inf, Inf, Inf]
        @test tied_duration(d) == [4.0, 0.0, 0.0]          # tsna clips to [2, 6]
        @test t_edge_density(d) ≈ 0.5
        @test t_edge_density(d; agg_unit=:dyad) ≈ 4 / 24
        # Paths run without an end, as tsna::tPath's default end = Inf
        @test forward_reachable_set(d, 1, 0.0) == [1, 2, 3]   # was [1, 2]
        @test earliest_arrival(d, 1, 0.0)[1] == Dict(1 => 0.0, 2 => 2.0, 3 => 6.0)
        @test temporal_distance(d, 1, 3, 0.0) == 6.0
        @test reachability_matrix(d, 0.0)[1, 3]
        @test backward_reachable_set(d, 3, 10.0) == [1, 2, 3]
        @test temporal_path(d, 1, 3, 0.0).times == [2.0, 6.0]
        @test tie_decay(d)[(2, 3)] == 1.0                   # at the last change time
        @test ne(t_aggregate(d)) == 2                        # the whole axis
        @test get_edge_attribute(t_aggregate(d; method=:weighted), :weight, 2, 3) == Inf
        @test_throws ArgumentError t_aggregate(d; onset=0.0)
        @test [r.n_edges for r in window_sna_stats(d, 2.0; measures=[:n_edges])] == [1.0, 1.0, 1.0]
        @test t_edge_persistence(d, 4.0) == 0.0               # starts 2 and 6 (closed range)
        # A left-censored tie that dissolves at the first change time: the
        # dissolution used to be lost, and the duration was 2 instead of Inf.
        e = DynamicNetwork(2)
        activate!(e, -Inf, 5.0; edge=(1, 2))
        activate!(e, 8.0, 10.0; edge=(1, 2))
        @test [t_edge_dissolution(e, t, t + 1) for t in 0.0:10.0] ==
              [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]
        @test [r.n_dissolutions for r in t_turnover(e, 1.0)] == [1, 0, 0, 0, 0, 1]
        @test [r.n_formations for r in t_turnover(e, 1.0; include_censored=true)] ==
              [0, 0, 0, 1, 0, 0]                             # -Inf is no time step
        @test t_edge_duration(e) == [Inf]
        @test t_edge_duration(e; mode=:spell) == [Inf, 2.0]
        @test tied_duration(e) == [2.0, 0.0]
        @test t_edge_density(e) ≈ 0.4
        # No finite spell bound at all: tsna falls back to the range (0, 1);
        # TSNA refuses wherever that range would enter the answer.
        p = DynamicNetwork(3)
        add_edge!(p.network, 1, 2)                           # no record: always active
        activate!(p, -Inf, Inf; edge=(2, 3))
        @test t_edge_duration(p) == [Inf, Inf]
        @test t_vertex_duration(p) == [Inf, Inf, Inf]
        for f in (tied_duration, t_edge_density, tie_decay)
            err = try f(p); nothing catch x; x end
            @test err isa ArgumentError && occursin("set_observation_period!", err.msg)
        end
        for f in (t_turnover, window_sna_stats, t_edge_persistence)
            @test_throws ArgumentError f(p, 1.0)
        end
        @test earliest_arrival(p, 2, 3.0)[1] == Dict(2 => 3.0, 3 => 3.0)   # threw before
        @test temporal_distance(p, 1, 3, 3.0) == 0.0
        @test tie_decay(p; at=4.0)[(1, 2)] == 1.0
        @test t_edge_formation(p, -Inf, Inf; include_censored=true) == 2     # at -Inf
        @test t_edge_formation(p, 0.0, 10.0; include_censored=true) == 0
        # The same rules on an integer time axis, where the axis extremes
        # stand for -Inf/Inf: unbounded durations are Inf, never an overflow.
        q = DynamicNetwork{Int, Int}(3)
        add_edge!(q.network, 1, 2)
        activate!(q, 2, 5; edge=(2, 3))
        deactivate!(q, 3, 4; edge=(1, 2))                    # cuts (typemin, typemax)
        @test t_edge_duration(q) == [Inf, 3.0]
        @test t_edge_duration(q; mode=:spell) == [Inf, Inf, 3.0]
        @test t_edge_density(q) ≈ (1 + 1 + 3) / (2 * 3)     # clipped to [2, 5]
        @test [r.n_formations for r in t_turnover(q, 1)] == [1, 0, 1, 0]
        @test earliest_arrival(q, 1, 0)[1] == Dict(1 => 0, 2 => 0, 3 => 2)
        # With a window nothing changes: clip to it
        set_observation_period!(d, 0.0, 6.0)
        @test t_edge_duration(d) == [4.0]
        @test forward_reachable_set(d, 1, 0.0) == [1, 2]
    end

    @testset "t_edge_density counts the base network's dyads" begin
        # tsna's network.dyadcount: n1*n2 on a two-mode network, n^2 with
        # loops (R tsna 0.3.6 gives 0.25 and 0.2222; n(n-1) gave 0.15 and 0.333).
        d = as_dynamic_network(Network(5; directed=false, bipartite=2); onset=0.0, terminus=10.0)
        activate!(d, 0.0, 10.0; edge=(1, 3)); activate!(d, 0.0, 5.0; edge=(2, 4))
        @test t_edge_density(d; agg_unit=:dyad) ≈ 0.25
        @test TSNA._dyad_count(d.network) == 6
        l = as_dynamic_network(Network(3; directed=true, loops=true); onset=0.0, terminus=10.0)
        activate!(l, 0.0, 10.0; edge=(1, 1)); activate!(l, 0.0, 10.0; edge=(1, 2))
        @test t_edge_density(l; agg_unit=:dyad) ≈ 2 / 9
        @test TSNA._dyad_count(network(4; directed=false, loops=true)) == 10
        @test TSNA._dyad_count(network(4; directed=true, bipartite=1)) == 6
        @test TSNA._dyad_count(network(4)) == 12
        @test TSNA._dyad_count(network(4; directed=false)) == 6
        # A zero-length range has no density (NaN, as tsna's 0/0); an
        # infinite window has no length to divide by
        z = DynamicNetwork(2; observation_start=1.0, observation_end=1.0)
        activate!(z, 1.0, 1.0; edge=(1, 2))
        @test isnan(t_edge_density(z))
        z2 = DynamicNetwork(2)                         # one change time: [5, 5]
        activate!(z2, 0.0, 5.0; edge=(1, 2)); deactivate!(z2, 0.0, 5.0; edge=(1, 2))
        activate!(z2, -Inf, 3.0; edge=(2, 1)); deactivate!(z2, -Inf, 3.0; edge=(2, 1))
        activate!(z2, 5.0, Inf; edge=(1, 2))
        @test get_change_times(z2) == [5.0] && isnan(t_edge_density(z2))
        set_observation_period!(z2, -Inf, Inf)
        @test_throws ArgumentError t_edge_density(z2)
    end

    @testset "Int32 vertex ids give the Int answers" begin
        # SNA's geodesic and max-flow kernels used to take n::Int, so
        # t_closeness threw a MethodError on DynamicNetwork{Int32, ...}.
        build(T) = begin
            d = DynamicNetwork{T, Float64}(5; directed=false, observation_start=0.0,
                                           observation_end=10.0)
            for (i, j, a, b) in ((1, 2, 0.0, 6.0), (2, 3, 0.0, 10.0), (3, 4, 2.0, 9.0), (4, 5, 5.0, 10.0))
                activate!(d, a, b; edge=(T(i), T(j)))
            end
            activate!(d, 4.0, 10.0; vertex=T(5))
            d
        end
        d64, d32 = build(Int), build(Int32)
        for t in (1.0, 3.0, 7.0), f in (t_closeness, t_betweenness, t_degree, t_eigenvector)
            @test isequal(f(d32, t), f(d64, t))
        end
        @test t_closeness(d32, 7.0)[1:5] ≈ t_closeness(d64, 7.0)[1:5]
        @test temporal_distance_matrix(d32, 0.0) == temporal_distance_matrix(d64, 0.0)
        @test t_edge_duration(d32) == t_edge_duration(d64)
        @test length(ContactSequence(Contact{Int32, Float64}[], Int32(3))) == 0
    end

    @testset "Concretely typed result vectors" begin
        d = DynamicNetwork(3; observation_start=0.0, observation_end=10.0)
        activate!(d, 0.0, 5.0; edge=(1, 2))
        @test isconcretetype(eltype(t_sna_stats(d, [1.0, 2.0])))
        @test isconcretetype(eltype(window_sna_stats(d, 5.0)))
        @test isconcretetype(eltype(t_turnover(d, 5.0)))
        @test_throws ArgumentError t_sna_stats(d, [1.0]; measures=[:bogus])
    end

    @testset "Golden fixture: durations, events, paths and snapshot measures vs R tsna/sna" begin
        fx = NetworkCore.load_golden(joinpath(@__DIR__, "fixtures", "tsna_semantics.toml"))
        V = fx.values
        tol = fx.tolerance["centrality"]
        num(x) = parse(Float64, x)                        # "Inf"/"-Inf" included
        ints(x) = parse.(Int, split(x))
        # R's labelled values, put in the order TSNA documents for its plain
        # vectors (edges by (i, j) -- (min, max) on undirected networks --
        # spells by edge then onset, vertices by id). The expected vectors are
        # built from the fixture alone, not from TSNA's internals.
        by_edge(xs) = [num(split(x)[3]) for x in sort(xs; by=x -> ints(join(split(x)[1:2], " ")))]
        by_spell(xs) = [num(split(x)[4]) for x in
                        sort(xs; by=x -> (ints(join(split(x)[1:2], " ")), num(split(x)[3])))]
        by_vertex(xs) = Dict(parse(Int, split(x)[1]) => num(split(x)[2]) for x in xs)
        function build(c, n; window::Bool)
            directed, bip, loops = c["directed"], c["bipartite"], c["loops"]
            d = window ? DynamicNetwork{Int, Float64, directed}(n; observation_start=0.0,
                                                                observation_end=12.0) :
                         DynamicNetwork{Int, Float64, directed}(n)
            d.network = Network{Int}(; n=n, directed=directed, loops=loops,
                                     bipartite=bip > 0 ? bip : nothing)
            for e in c["edges"]
                @test add_edge!(d.network, ints(e)...)
            end
            for op in c["ops"]
                k, a, b, on, te = split(op)
                a, b = parse(Int, a), parse(Int, b)
                k == "av" ? activate!(d, num(on), num(te); vertex=a) :
                            activate!(d, num(on), num(te); edge=(a, b))
            end
            return d
        end
        function check_tied(d, c)
            for (k, nb) in enumerate(c["tied_neighborhoods"])
                @test tied_duration(d; neighborhood=Symbol(nb)) == num.(split(c["tied_duration"][k], ","))
                @test tied_duration(d; mode=:counts, neighborhood=Symbol(nb)) ==
                      num.(split(c["tied_counts"][k], ","))
            end
        end

        # Groups 1 and 2: observed over [0, 12] (one-mode; two-mode and looped)
        n_checked = 0
        for g in 1:V["n_cases"]
            c = V["case_$g"]
            d = build(c, 6; window=true)
            # tsna::edgeDuration(nd), subject "edges" and "spells"
            @test t_edge_duration(d) == by_edge(c["edge_duration"])
            @test t_edge_duration(d; mode=:spell) == by_spell(c["spell_duration"])
            # TSNA clips vertex spells to the window (computed in R from
            # networkDynamic's stored spells); tsna::vertexDuration replaces
            # only infinite bounds, so the two agree on the vertices whose
            # finite bounds lie in the window.
            clipped = by_vertex(c["vertex_duration_clipped"])
            jl_vd = Dict(zip(sort(collect(keys(clipped))), t_vertex_duration(d)))
            @test jl_vd == clipped
            r_vd = by_vertex(c["vertex_duration"])
            for x in c["vertex_inside"]
                v, inside = parse(Int, split(x)[1]), split(x)[2] == "true"
                inside && @test jl_vd[v] == r_vd[v]
            end
            # tsna::tEdgeFormation / tEdgeDissolution(nd, start = 0, end = 12)
            @test [t_edge_formation(d, t, t + 1) for t in 0:12] == c["formation"]
            @test [t_edge_dissolution(d, t, t + 1) for t in 0:12] == c["dissolution"]
            @test [t_edge_formation(d, t, t + 1; include_censored=true) for t in 0:12] ==
                  c["formation_censored"]
            @test [t_edge_dissolution(d, t, t + 1; include_censored=true) for t in 0:12] ==
                  c["dissolution_censored"]
            # sna on network.extract(nd, at = t), scattered back by vertex id
            for snapline in c["snapshots"]
                f = split(snapline, "|")
                t = num(f[1])
                ids = isempty(f[2]) ? Int[] : parse.(Int, split(f[2], ","))
                absent = setdiff(1:6, ids)
                for (col, fn) in ((3, t_degree), (4, t_betweenness), (5, t_closeness))
                    jl = fn(d, t)
                    @test all(isnan, jl[absent])
                    isempty(ids) && continue
                    @test isapprox(jl[ids], num.(split(f[col], ",")); atol=tol, nans=true)
                end
                isempty(ids) || @test isapprox(t_density(d, t), num(f[6]); atol=tol, nans=true)
                n_checked += 1
            end
            # tsna::tiedDuration and tsna::tEdgeDensity (network.dyadcount dyads)
            check_tied(d, c)
            @test TSNA._dyad_count(d.network) == c["dyad_count"]
            @test isapprox([t_edge_density(d), t_edge_density(d; agg_unit=:dyad),
                            t_edge_density(d; mode=:event)], c["edge_density"]; atol=tol)
        end
        @test n_checked == 4 * V["n_onemode"]
        @test count(g -> V["case_$g"]["bipartite"] > 0, 1:V["n_cases"]) >= 8
        @test count(g -> V["case_$g"]["loops"], 1:V["n_cases"]) >= 8

        # Group 3: no observation window. tsna reads lifetimes and events
        # unclipped (Inf for open spells), clips tiedDuration/tEdgeDensity to
        # the range of the change times, and runs default event series over
        # the closed range of the change times; with no change time at all it
        # uses the placeholder (0, 1), which TSNA refuses.
        n_paths = 0
        for g in 1:V["n_nowin"]
            c = V["nowin_$g"]
            n = c["n"]
            d = build(c, n; window=false)
            @test get_observation_period(d) === nothing
            @test get_change_times(d) == Float64.(c["change_times"])
            @test isequal(t_edge_duration(d), by_edge(c["edge_duration"]))
            @test isequal(t_edge_duration(d; mode=:spell), by_spell(c["spell_duration"]))
            r_vd = by_vertex(c["vertex_duration"])
            @test Dict(zip(sort(collect(keys(r_vd))), t_vertex_duration(d))) == r_vd
            if isempty(c["change_times"])
                @test_throws ArgumentError tied_duration(d)
                @test_throws ArgumentError t_edge_density(d)
                @test_throws ArgumentError t_turnover(d, 1.0)
                @test_throws ArgumentError window_sna_stats(d, 1.0)
            else
                lo, hi = extrema(c["change_times"])
                rows = t_turnover(d, 1.0)
                @test first(rows).window_start == c["series_start"] == lo
                @test last(rows).window_start == hi            # the closed range
                @test [r.n_formations for r in rows] == c["formation"]
                @test [r.n_dissolutions for r in rows] == c["dissolution"]
                rows_c = t_turnover(d, 1.0; include_censored=true)
                @test [r.n_formations for r in rows_c] == c["formation_censored"]
                @test [r.n_dissolutions for r in rows_c] == c["dissolution_censored"]
                @test [t_edge_formation(d, t, t + 1) for t in lo:hi] == c["formation"]
                @test [t_edge_dissolution(d, t, t + 1; include_censored=true) for t in lo:hi] ==
                      c["dissolution_censored"]
                check_tied(d, c)
                @test isapprox([t_edge_density(d), t_edge_density(d; agg_unit=:dyad)],
                               c["edge_density"][1:2]; atol=tol)
                # tsna 0.3.6's event density divides by (n_edges · end − start)
                lo == 0 && @test isapprox(t_edge_density(d; mode=:event), c["edge_density"][3]; atol=tol)
            end
            # tsna::tPath(direction = "fwd") from the default start, end = Inf
            for p in c["paths"]
                v, tdist = parse(Int, split(p, "|")[1]), num.(split(split(p, "|")[2], ","))
                arr, _ = earliest_arrival(d, v, c["path_start"])
                @test [haskey(arr, u) ? arr[u] - c["path_start"] : Inf for u in 1:n] == tdist
                n_paths += 1
            end
        end
        @test n_paths >= 200

        # Group 4: samplk-style panels, networkDynamic(network.list = ...)
        for g in 1:V["n_panels"]
            c = V["panel_$g"]
            n, directed = c["n"], c["directed"]
            nets = map(1:3) do k
                w = network(n; directed=directed)
                for e in c["wave_$k"]
                    add_edge!(w, ints(e)...)
                end
                w
            end
            d = DynamicNetwork(nets)
            @test t_edge_duration(d) == by_edge(c["edge_duration"])
            @test [t_edge_formation(d, t, t + 1) for t in 0:3] == c["formation"]
            @test [t_edge_dissolution(d, t, t + 1) for t in 0:3] == c["dissolution"]
            check_tied(d, c)
            @test isapprox([t_edge_density(d), t_edge_density(d; agg_unit=:dyad)],
                           c["edge_density"]; atol=tol)
            @test isapprox([t_density(d, t) for t in 0:2], c["gden"][1:3]; atol=tol)
            @test isapprox([t_reciprocity(d, t) for t in 0:2], c["grecip"][1:3]; atol=tol)
            @test isnan(c["gden"][4])                         # nobody is active at t = 3
            arr, _ = earliest_arrival(d, 1, 0.0)
            @test [haskey(arr, u) ? arr[u] : Inf for u in 1:n] == c["tpath_1"]
        end
    end
end

@testset "Every exported docstring carries a runnable example" begin
    # Mirrors ERGM.jl's testset: every TSNA-owned docstring of an exported,
    # non-deprecated binding contains a ```julia block, and every block runs in
    # a fresh module (an example that needs another package says so itself).
    meta = Base.Docs.meta(TSNA)
    undocumented = String[]; missing_example = String[]
    blocks = Tuple{String,String}[]
    for nm in names(TSNA)
        (nm === :TSNA || Base.isdeprecated(TSNA, nm)) && continue
        b = Base.Docs.Binding(TSNA, nm)
        if !haskey(meta, b)
            push!(undocumented, string(nm))
            continue
        end
        has_example = false
        for (_, ds) in meta[b].docs
            txt = ds.text isa AbstractString ? ds.text : join(string.(ds.text), "\n")
            for m in eachmatch(r"```julia\n(.*?)```"s, txt)
                has_example = true
                push!(blocks, (string(nm), String(m.captures[1])))
            end
        end
        has_example || push!(missing_example, string(nm))
    end
    @test isempty(undocumented)
    @test isempty(missing_example)
    @test length(blocks) >= 30
    for (nm, code) in blocks
        m = Module(Symbol("DocExample_", nm))
        ok = try
            Core.eval(m, :(using TSNA))
            Core.eval(m, Meta.parseall(code; filename="docstring:$nm"))
            true
        catch err
            println(stderr, "docstring example of $nm failed: ", sprint(showerror, err))
            false
        end
        @test ok
    end
end

@testset "Aqua" begin
    # Deprecated bindings are exported for one release; Aqua's undefined-export
    # check is satisfied by them.
    Aqua.test_all(TSNA)
    @test isempty(Test.detect_ambiguities(TSNA))
end
