# TSNA.jl

*Temporal Social Network Analysis for Julia*

A Julia package for descriptive analysis of dynamic networks, including temporal centrality measures, time-respecting path analysis, reachability, and duration metrics.

## Overview

Social network analysis traditionally operates on static network snapshots. TSNA.jl extends these methods to dynamic networks where vertices and edges appear and disappear over time. It provides tools for computing temporal versions of classical SNA measures, finding time-respecting paths, analyzing reachability, and computing duration and turnover metrics.

TSNA.jl is a port of the R [tsna](https://cran.r-project.org/package=tsna) package from the [StatNet](https://statnet.org/) collection.

### What is Temporal SNA?

Temporal SNA applies descriptive network analysis to time-varying networks. Unlike static SNA, temporal analysis respects the ordering of time -- an actor can only influence another through a sequence of contacts that occur in non-decreasing time order.

```text
Static view:    A -- B -- C     (A can reach C through B)
Temporal view:  A --(t=5)--> B --(t=3)--> C    (A cannot reach C: edge B->C occurs before A->B)
```

### Key Concepts

| Concept | Description |
|---------|-------------|
| **Temporal Centrality** | Centrality measures computed on network snapshots at specific times |
| **Temporal Path** | A path where edge traversals occur in non-decreasing time order |
| **Temporal Distance** | Earliest arrival time from source to target |
| **Forward Reachability** | Set of vertices reachable from a source after a start time |
| **Backward Reachability** | Set of vertices that can reach a target before an end time |
| **Edge Duration** | Total time an edge is active |
| **Turnover** | Rate of edge formation and dissolution |

### Applications

Temporal SNA is used in:

- **Epidemiology**: Modeling disease spread through time-varying contact networks
- **Information diffusion**: How quickly can information reach all actors?
- **Organizational dynamics**: Measuring evolving influence and centrality
- **Communication networks**: Analyzing temporal patterns in email or messaging data
- **Transportation**: Time-constrained routing through dynamic networks
- **Animal behavior**: Studying temporal interaction patterns

## Features

- **Temporal centrality**: Time-varying degree, betweenness, closeness, eigenvector, and PageRank
- **Temporal paths**: Find shortest time-respecting paths and earliest arrival times
- **Reachability analysis**: Forward and backward reachable sets from any vertex
- **Duration metrics**: Edge and vertex activity duration, persistence, and turnover rates
- **Contact sequences**: Convert dynamic networks to temporal contact sequences
- **Aggregation**: Compute statistics at multiple time points or in sliding windows
- **Network aggregation**: Collapse dynamic networks using union, intersection, or weighted methods

### Naming

All public functions have snake_case primary names (the Julia convention)
with the R tsna-style camelCase aliases kept for migration:
[`earliest_arrival`](@ref)/`earliestArrival`, [`t_degree`](@ref)/`tDegree`,
[`t_sna_stats`](@ref)/`tSnaStats`, ... Both are exported and
interchangeable. The earliest-arrival search is a heap-based Dijkstra
label-setting over a memoized per-network contact index, so repeated path
and reachability queries on the same network are fast.

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/statistical-network-analysis-with-Julia/Network.jl")
Pkg.add(url="https://github.com/statistical-network-analysis-with-Julia/NetworkDynamic.jl")
Pkg.add(url="https://github.com/statistical-network-analysis-with-Julia/SNA.jl")
Pkg.add(url="https://github.com/statistical-network-analysis-with-Julia/TSNA.jl")
```

Or for development:

```julia
using Pkg
Pkg.develop(path="/path/to/TSNA.jl")
```

## Quick Start

```julia
using NetworkDynamic
using TSNA

# Create a dynamic network
dnet = DynamicNetwork(10; observation_start=0.0, observation_end=100.0)
activate_vertices!(dnet, collect(1:10), 0.0, 100.0)
activate!(dnet, 0.0, 50.0; edge=(1, 2))
activate!(dnet, 20.0, 80.0; edge=(2, 3))
activate!(dnet, 40.0, 100.0; edge=(3, 4))

# Temporal centrality at time 50
deg = t_degree(dnet, 50.0)
bet = t_betweenness(dnet, 50.0)

# Temporal path finding
dist = temporal_distance(dnet, 1, 4, 0.0)
path = shortest_temporal_path(dnet, 1, 4, 0.0)

# Duration metrics
mean_dur = t_edge_duration(dnet; aggregate=:mean)
turnover = t_turnover(dnet, 20.0)
```

## Choosing Analyses

| Question | Function |
|----------|----------|
| Who is most central at time t? | [`t_degree`](@ref), [`t_betweenness`](@ref), [`t_closeness`](@ref) |
| How fast can information spread? | [`temporal_distance`](@ref), [`forward_reachable_set`](@ref) |
| Who can reach whom? | [`forward_reachable_set`](@ref), [`backward_reachable_set`](@ref) |
| How stable are ties? | [`t_edge_duration`](@ref), [`t_edge_persistence`](@ref) |
| How much turnover is there? | [`t_turnover`](@ref), [`tie_decay`](@ref) |
| How does the network evolve? | [`t_sna_stats`](@ref), [`window_sna_stats`](@ref) |

## Documentation

```@contents
Pages = [
    "getting_started.md",
    "guide/centrality.md",
    "guide/paths.md",
    "guide/metrics.md",
    "api/centrality.md",
    "api/paths.md",
    "api/metrics.md",
]
Depth = 2
```

## Theoretical Background

### Time-Respecting Paths

In a temporal network, a valid path from $s$ to $r$ starting at time $t_0$ must traverse edges in non-decreasing time order:

$$s = v_0 \xrightarrow{t_1} v_1 \xrightarrow{t_2} v_2 \xrightarrow{t_3} \ldots \xrightarrow{t_k} v_k = r$$

Where $t_0 \leq t_1 \leq t_2 \leq \ldots \leq t_k$ and each edge $(v_{i-1}, v_i)$ is active at time $t_i$.

This constraint means that temporal reachability is **not symmetric** and **not transitive** in general.

### Temporal Distance

The temporal distance from $s$ to $r$ starting at $t_0$ is the earliest time at which $r$ can be reached:

$$d_T(s, r, t_0) = \min\{t_k : \exists \text{ time-respecting path } s \to r \text{ starting at } t_0 \text{ arriving at } t_k\}$$

## References

1. Bender-deMoll, S., Morris, M. (2012). `tsna`: Tools for Temporal Social Network Analysis. R package.

2. Holme, P. (2015). Modern temporal network theory: a colloquium. *European Physical Journal B*, 88(9), 1-30.

3. Holme, P., Saramaki, J. (2012). Temporal networks. *Physics Reports*, 519(3), 97-125.

4. Nicosia, V., Tang, J., Mascolo, C., Musolesi, M., Russo, G., Latora, V. (2013). Graph metrics for temporal networks. In *Temporal Networks* (pp. 15-40). Springer.
