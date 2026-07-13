using TSNA
using NetworkDynamic
using Networks
using Random
using Test

# Reference implementation of earliest arrival (the pre-heap linear
# min-scan version with a per-call contact index), kept verbatim for
# exact-equivalence testing of the optimized search.
function _reference_earliest_arrival(dnet::DynamicNetwork{T, Time}, source::T,
                                     start_time;
                                     end_time=dnet.observation_period[2]) where {T, Time}
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

    @testset "snake_case / camelCase aliases" begin
        @test tPath === TemporalPath
        @test earliestArrival === earliest_arrival
        @test temporalDistance === temporal_distance
        @test forwardReachableSet === forward_reachable_set
        @test backwardReachableSet === backward_reachable_set
        @test temporalPath === temporal_path
        @test shortestTemporalPath === temporal_path
        @test shortest_temporal_path === temporal_path
        @test tDegree === t_degree
        @test tBetweenness === t_betweenness
        @test tCloseness === t_closeness
        @test tEigenvector === t_eigenvector
        @test tPagerank === t_pagerank
        @test tDensity === t_density
        @test tReciprocity === t_reciprocity
        @test tTransitivity === t_transitivity
        @test tEdgeDuration === t_edge_duration
        @test tVertexDuration === t_vertex_duration
        @test tEdgeFormation === t_edge_formation
        @test tEdgeDissolution === t_edge_dissolution
        @test tEdgePersistence === t_edge_persistence
        @test tTurnover === t_turnover
        @test tieDecay === tie_decay
        @test tSnaStats === t_sna_stats
        @test windowSnaStats === window_sna_stats
        @test tAggregate === t_aggregate

        # Both spellings give identical results
        dnet = chain_fixture()
        @test t_degree(dnet, 5.0) == tDegree(dnet, 5.0)
        @test t_density(dnet, 5.0) == tDensity(dnet, 5.0)
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
            activate!(dnet, 3.0, 6.0; edge=(1, 2))     # overlapping spell
            activate!(dnet, 7.0, 7.0; edge=(2, 3))     # point spell

            cs, rep = as_contact_sequence(dnet; report=true)

            # Preserved: one contact per spell (overlaps stay separate), the
            # vertex count, directedness, onsets and durations.
            @test length(cs) == 3
            @test cs.n_vertices == 4
            @test cs.directed == directed
            @test [c.time for c in cs] == [0.0, 3.0, 7.0]
            @test [c.duration for c in cs] == [4.0, 3.0, 0.0]   # point spell = 0

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
