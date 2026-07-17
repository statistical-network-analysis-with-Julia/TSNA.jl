# TSNA.jl


[![Network Analysis](https://img.shields.io/badge/Network-Analysis-orange.svg)](https://github.com/statistical-network-analysis-with-Julia/TSNA.jl)
[![Build Status](https://github.com/statistical-network-analysis-with-Julia/TSNA.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/statistical-network-analysis-with-Julia/TSNA.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://statistical-network-analysis-with-Julia.github.io/TSNA.jl/stable/)
[![Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://statistical-network-analysis-with-Julia.github.io/TSNA.jl/dev/)
[![Julia](https://img.shields.io/badge/Julia-1.12+-purple.svg)](https://julialang.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

<p align="center">
  <img src="docs/src/assets/logo.svg" alt="TSNA.jl icon" width="160">
</p>

Temporal Social Network Analysis for Julia.

## Overview

TSNA.jl provides descriptive analysis tools for dynamic networks, including temporal centrality measures, temporal path analysis, reachability analysis, and duration metrics.

This package is a Julia port of the R `tsna` package from the StatNet collection.

## Installation

Requires Julia 1.12+. TSNA.jl depends on the unregistered
[Networks.jl](https://github.com/statistical-network-analysis-with-Julia/Networks.jl), [NetworkDynamic.jl](https://github.com/statistical-network-analysis-with-Julia/NetworkDynamic.jl), and [SNA.jl](https://github.com/statistical-network-analysis-with-Julia/SNA.jl) packages, which must be added first (in this order):

```julia
using Pkg
Pkg.add(url="https://github.com/statistical-network-analysis-with-Julia/Networks.jl")
Pkg.add(url="https://github.com/statistical-network-analysis-with-Julia/NetworkDynamic.jl")
Pkg.add(url="https://github.com/statistical-network-analysis-with-Julia/SNA.jl")
Pkg.add(url="https://github.com/statistical-network-analysis-with-Julia/TSNA.jl")
```

For development, you can instead clone all ecosystem repositories side by
side (the monorepo layout) and start Julia with the root workspace project
(`julia --project=.` in the clone root): the `[sources]` path dependencies
then wire the packages together with no ordered installs needed.

## Features

- **Temporal centrality**: Time-varying degree, betweenness, closeness
- **Temporal paths**: Time-respecting paths and reachability
- **Duration metrics**: Edge/vertex activity duration, formation/dissolution counts, persistence, turnover
- **Aggregation**: Compute statistics over time windows

## Naming

All public functions have **snake_case primary names** (the Julia
convention) with the R tsna-style **camelCase aliases** kept for
migration: `earliest_arrival`/`earliestArrival`, `t_degree`/`tDegree`,
`t_sna_stats`/`tSnaStats`, ... Both are exported and interchangeable.

## Quick Start

```julia
using NetworkDynamic
using TSNA

# Create dynamic network
dnet = DynamicNetwork(10; observation_start=0.0, observation_end=100.0)
for i in 1:10
    activate!(dnet, 0.0, 100.0; vertex=i)
end
activate!(dnet, 0.0, 60.0; edge=(1, 2))
activate!(dnet, 10.0, 80.0; edge=(2, 5))

# Temporal centrality at time 50
deg = t_degree(dnet, 50.0)
bet = t_betweenness(dnet, 50.0)

# Temporal path finding (from vertex 1 to vertex 5, starting at t=0)
dist = temporal_distance(dnet, 1, 5, 0.0)
path = shortest_temporal_path(dnet, 1, 5, 0.0)

# Reachability
reachable = forward_reachable_set(dnet, 1, 0.0)
```

## Temporal Centrality

```julia
at = 50.0   # evaluation time used below

# At a specific time point
t_degree(dnet, at; mode=:total)      # :in, :out, or :total
t_betweenness(dnet, at)
t_closeness(dnet, at)
t_eigenvector(dnet, at)
t_pagerank(dnet, at; damping=0.85)
```

## Temporal Network Measures

```julia
# At a time point
t_density(dnet, at)
t_reciprocity(dnet, at)
t_transitivity(dnet, at)

# Over time series
times = 0.0:10.0:100.0
stats = t_sna_stats(dnet, times; measures=[:density, :reciprocity])
```

## Temporal Paths

A temporal path must have non-decreasing times - you can only traverse edges when they're active.

<!-- skip-check -->
```julia
# Temporal path structure
struct tPath{T, Time}
    vertices::Vector{T}     # Sequence of vertices
    times::Vector{Time}     # Time at each transition
    edges::Vector{Tuple}    # Edges traversed
end

# Find earliest arrival time
arrival = temporal_distance(dnet, source, target, start_time)

# Find actual path
path = shortest_temporal_path(dnet, source, target, start_time)
```

## Reachability

<!-- skip-check -->
```julia
# Who can source reach starting at start_time?
forward = forward_reachable_set(dnet, source, start_time)

# Who can reach target by end_time?
backward = backward_reachable_set(dnet, target, end_time)
```

## Duration Metrics

```julia
# Edge duration statistics
mean_dur = t_edge_duration(dnet; aggregate=:mean)
all_durs = t_edge_duration(dnet; aggregate=:all)  # Raw per-spell vector

# Vertex duration
t_vertex_duration(dnet; aggregate=:mean)

# Formation / dissolution events in [onset, terminus), counted as spell
# onset/terminus events (like tsna::tEdgeFormationAt)
n_form = t_edge_formation(dnet, 0.0, 50.0)
n_diss = t_edge_dissolution(dnet, 0.0, 50.0)

# Proportion of edges surviving across consecutive 10-unit windows
# (temporal-correlation measure of Nicosia et al.)
persistence = t_edge_persistence(dnet, 10.0)

# Turnover rates in 10-unit windows
turnover = t_turnover(dnet, 10.0)
# Returns: (formation_rate, dissolution_rate, n_formations, n_dissolutions)

# Tie decay rate
decay_rate = tie_decay(dnet; method=:exponential)
```

## Aggregation

```julia
# Statistics at multiple time points
times = collect(0.0:5.0:100.0)
stats = t_sna_stats(dnet, times; measures=[:density, :n_edges, :mean_degree])

# Snapshot statistics at the start of each window of width 20
stats = window_sna_stats(dnet, 20.0; measures=[:density])

# Aggregate to static network
static = t_aggregate(dnet; method=:union)        # Ever active
static = t_aggregate(dnet; method=:intersection)  # Always active
static = t_aggregate(dnet; method=:weighted)      # Weight = total time
```

## Contact Sequences

<!-- skip-check -->
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

<!-- skip-check -->
```julia
# Who could receive information from source within time window?
start_time = 0.0
reachable = forward_reachable_set(dnet, source, start_time)
println("$(length(reachable)) nodes reachable from source")

# Track how reachability grows over time
for t in 10.0:10.0:100.0
    r = forward_reachable_set(dnet, source, start_time)
    r_by_t = filter(v -> temporal_distance(dnet, source, v, start_time) <= t, r)
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

## Citation

If you use TSNA.jl in your work, please cite it using the entry in
[`CITATION.bib`](CITATION.bib):

```biblatex
@misc{SNWJTSNAJL,
  author = {{Statistical Network Analysis with Julia}},
  title = {TSNA.jl: Temporal Social Network Analysis for Julia},
  year = {2026},
  url = {https://github.com/statistical-network-analysis-with-Julia/TSNA.jl},
  note = {Homepage: https://statistical-network-analysis-with-Julia.github.io/TSNA.jl; GitHub: https://github.com/statistical-network-analysis-with-Julia}
}
```

## License

MIT License - see [LICENSE](LICENSE) for details.
