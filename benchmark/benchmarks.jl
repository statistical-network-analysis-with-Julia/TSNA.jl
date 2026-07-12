#!/usr/bin/env julia
# benchmark/benchmarks.jl — BenchmarkTools suite for TSNA.jl's hot loops.
#
# Covers the optimized `earliest_arrival` temporal search (heap-based
# frontier over the memoized per-network contact index, replacing the
# linear min-scan with a per-call index). Benchmarked from a fixed sample
# of sources on random temporal networks at n = 500 and n = 2000 with the
# spell count scaled proportionally.
#
# Defines the standard `SUITE::BenchmarkGroup`. Run standalone with
#     julia --project=benchmark benchmark/benchmarks.jl
# which tunes + runs the suite and prints one tab-separated `BENCHJL` line
# per benchmark (consumed by the site repo's tools/run_benchmarks.jl).

using BenchmarkTools
using NetworkDynamic
using Random
using TSNA

# ---------------------------------------------------------------------------
# Fixtures: random temporal networks, ~10 spells per vertex over [0, 100]
# ---------------------------------------------------------------------------

const SIZES = (500, 2000)
const SPELLS_PER_VERTEX = 10
const N_SOURCES = 10

function random_temporal_network(rng::AbstractRNG, n::Int, m::Int)
    dnet = DynamicNetwork(n; observation_start=0.0, observation_end=100.0,
                          directed=true)
    for v in 1:n
        activate!(dnet, 0.0, 100.0; vertex=v)
    end
    for _ in 1:m
        i, j = rand(rng, 1:n), rand(rng, 1:n)
        i == j && continue
        onset = round(100 * rand(rng); digits=2)
        if rand(rng) < 0.15
            activate!(dnet, onset, onset; edge=(i, j))      # point contact
        else
            terminus = min(100.0, onset + round(30 * rand(rng); digits=2))
            activate!(dnet, onset, terminus; edge=(i, j))
        end
    end
    return dnet
end

const DNETS = Dict(n => random_temporal_network(Random.Xoshiro(n), n,
                                                SPELLS_PER_VERTEX * n)
                   for n in SIZES)
const SOURCES = Dict(n => rand(Random.Xoshiro(n + 1), 1:n, N_SOURCES)
                     for n in SIZES)

"Earliest-arrival searches from a fixed sample of sources (index warm)."
function arrival_sweep(dnet, sources)
    total = 0
    for s in sources
        arrival, _ = earliest_arrival(dnet, s, 0.0)
        total += length(arrival)
    end
    return total
end

# Warm the memoized contact indexes so the benchmarks measure the search,
# not the one-off index build.
for n in SIZES
    arrival_sweep(DNETS[n], SOURCES[n])
end

# ---------------------------------------------------------------------------
# Suite
# ---------------------------------------------------------------------------

const SUITE = BenchmarkGroup()

let g = addgroup!(SUITE, "earliest_arrival")
    for n in SIZES
        g["sweep_n$(n)"] =
            @benchmarkable arrival_sweep($(DNETS[n]), $(SOURCES[n]))
    end
end

# ---------------------------------------------------------------------------
# Standalone entry point
# ---------------------------------------------------------------------------

function print_benchjl(results::BenchmarkGroup)
    for (path, trial) in BenchmarkTools.leaves(results)
        est = median(trial)
        println("BENCHJL\t", join(path, "/"), "\t",
                BenchmarkTools.time(est), "\t",
                BenchmarkTools.allocs(est), "\t",
                BenchmarkTools.memory(est))
    end
end

function main()
    tune!(SUITE)
    results = run(SUITE; verbose=false, seconds=1)
    print_benchjl(results)
    return results
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
