# TSNA.jl


[![Network Analysis](https://img.shields.io/badge/Network-Analysis-orange.svg)](https://github.com/statistical-network-analysis-with-Julia/TSNA.jl)
[![Build Status](https://github.com/statistical-network-analysis-with-Julia/TSNA.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/statistical-network-analysis-with-Julia/TSNA.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://statistical-network-analysis-with-Julia.github.io/TSNA.jl/stable/)
[![Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://statistical-network-analysis-with-Julia.github.io/TSNA.jl/dev/)
[![Julia](https://img.shields.io/badge/Julia-1.9+-purple.svg)](https://julialang.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

<p align="center">
  <img src="docs/src/assets/logo.svg" alt="TSNA.jl icon" width="160">
</p>

Temporal Social Network Analysis for Julia.

## Overview

TSNA.jl provides descriptive analysis tools for dynamic networks, including temporal centrality measures, temporal path analysis, reachability analysis, and duration metrics.

This package is a Julia port of the R `tsna` package from the StatNet collection.

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/statistical-network-analysis-with-Julia/TSNA.jl")
```

## Features

- **Temporal centrality**: Time-varying degree, betweenness, closeness
- **Temporal paths**: Time-respecting paths and reachability
- **Duration metrics**: Edge/vertex activity duration, persistence, turnover
- **Aggregation**: Compute statistics over time windows

## Quick Start

```julia
using NetworkDynamic
using TSNA

# Create dynamic network
dnet = DynamicNetwork(10; observation_start=0.0, observation_end=100.0)
# ... add activity spells ...

# Temporal centrality at time 50
deg = tDegree(dnet, 50.0)
bet = tBetweenness(dnet, 50.0)

# Temporal path finding
dist = temporalDistance(dnet, source, target, start_time)
path = shortestTemporalPath(dnet, source, target, start_time)

# Reachability
reachable = forwardReachableSet(dnet, source, start_time)
```

## Temporal Centrality

```julia
# At a specific time point
tDegree(dnet, at; mode=:total)      # :in, :out, or :total
tBetweenness(dnet, at)
tCloseness(dnet, at)
tEigenvector(dnet, at)
tPagerank(dnet, at; damping=0.85)
```

## Temporal Network Measures

```julia
# At a time point
tDensity(dnet, at)
tReciprocity(dnet, at)
tTransitivity(dnet, at)

# Over time series
times = 0.0:10.0:100.0
stats = tSnaStats(dnet, times; stats=[:density, :reciprocity])
```

## Temporal Paths

A temporal path must have non-decreasing times - you can only traverse edges when they're active.

```julia
# Temporal path structure
struct tPath{T, Time}
    vertices::Vector{T}     # Sequence of vertices
    times::Vector{Time}     # Time at each transition
    edges::Vector{Tuple}    # Edges traversed
end

# Find earliest arrival time
arrival = temporalDistance(dnet, source, target, start_time)

# Find actual path
path = shortestTemporalPath(dnet, source, target, start_time)
```

## Reachability

```julia
# Who can source reach starting at start_time?
forward = forwardReachableSet(dnet, source, start_time)

# Who can reach target by end_time?
backward = backwardReachableSet(dnet, target, end_time)
```

## Duration Metrics

```julia
# Edge duration statistics
mean_dur = tEdgeDuration(dnet; aggregate=:mean)
all_durs = tEdgeDuration(dnet; aggregate=:all)  # Per-edge Dict

# Vertex duration
tVertexDuration(dnet; aggregate=:mean)

# Edge persistence across time windows
persistence = tEdgePersistence(dnet, window_size)

# Turnover rates
turnover = tTurnover(dnet, window_size)
# Returns: (formation_rate, dissolution_rate, n_formations, n_dissolutions)

# Tie decay rate
decay_rate = tieDecay(dnet; method=:exponential)
```

## Aggregation

```julia
# Statistics at multiple time points
times = collect(0.0:5.0:100.0)
stats = tSnaStats(dnet, times; stats=[:density, :n_edges, :mean_degree])

# Statistics in sliding windows
stats = windowSnaStats(dnet, window_size;
                       stats=[:density],
                       step=window_size/2)

# Aggregate to static network
static = tAggregate(dnet; method=:union)       # Ever active
static = tAggregate(dnet; method=:intersection) # Always active
static = tAggregate(dnet; method=:weighted)     # Weight = total time
```

## Contact Sequences

```julia
# Convert to contact sequence
cs = as_contact_sequence(dnet)

# Contact structure
struct Contact{T, Time}
    source::T
    target::T
    time::Time
    duration::Time
end
```

## Example: Information Spread

```julia
# Who could receive information from source within time window?
start_time = 0.0
reachable = forwardReachableSet(dnet, source, start_time)
println("$(length(reachable)) nodes reachable from source")

# Track how reachability grows over time
for t in 10.0:10.0:100.0
    r = forwardReachableSet(dnet, source, start_time)
    r_by_t = filter(v -> temporalDistance(dnet, source, v, start_time) <= t, r)
    println("t=$t: $(length(r_by_t)) nodes reachable")
end
```

## Documentation

For more detailed documentation, see:

- [Stable Documentation](https://statistical-network-analysis-with-Julia.github.io/TSNA.jl/stable/)
- [Development Documentation](https://statistical-network-analysis-with-Julia.github.io/TSNA.jl/dev/)

## References

1. Bender-deMoll, S., Morris, M. (2019). tsna: Tools for Temporal Social Network Analysis. R package, StatNet collection.

2. Moody, J. (2002). The importance of relationship timing for diffusion. *Social Forces*, 81(1), 25-56.

3. Holme, P., Saramaki, J. (2012). Temporal networks. *Physics Reports*, 519(3), 97-125.

## License

MIT License - see [LICENSE](LICENSE) for details.
