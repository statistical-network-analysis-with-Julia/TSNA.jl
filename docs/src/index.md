# TSNA.jl

Describe a dynamic network while respecting the order and duration of contacts. TSNA.jl measures snapshot structure, time-respecting reachability, earliest arrival, and the duration and turnover of observed ties.

| First analysis | Learn the model or data | Reference and detail |
|:--|:--|:--|
| [Analyse a dynamic network](getting_started.md) | [Understand temporal paths](guide/paths.md) | [Choose temporal measures](guide/metrics.md) |

!!! note "Supported scope"

    Temporal paths describe possible reachability under the implemented contact and waiting rules; they do not predict transmission or causal influence. Analyses reject missing dyads by default (`missing=:error`); `missing=:face` explicitly uses stored values. Density uses active actors and reciprocity is the proportion of symmetric dyads (mutual or null) by default; legacy conventions require explicit keywords.

## Installation

```@raw html
<p>Use Julia <strong>1.12 or newer</strong> and the <a href="/getting-started/">shared workspace installation guide</a>. These development packages are not yet registered; the guide prepares the required sibling checkouts and a Julia environment for the examples.</p>
```

## Quick Start

A contact at time 1 can be followed by a later contact at time 5, provided the intermediate actor can wait:

```julia
using NetworkDynamic, TSNA

dnet = DynamicNetwork(4; observation_start=0.0, observation_end=10.0)
activate_vertices!(dnet, collect(1:4), 0.0, 10.0)
activate!(dnet, 1.0, 4.0; edge=(1, 2))
activate!(dnet, 5.0, 8.0; edge=(2, 3))
elapsed = temporal_distance(dnet, 1, 3, 0.0)
println((elapsed_to_actor3=elapsed, density_at_time2=t_density(dnet, 2.0)))
@assert elapsed == 5.0
```

The temporal distance is elapsed time from the requested start, not a hop count. Actor 2 remains active while waiting. Aggregating both contacts into a static graph would discard the timing condition that makes this route possible.

## Choosing Analyses

| Question | Function |
|----------|----------|
| Who is most central at time t? | [`t_degree`](@ref), [`t_betweenness`](@ref), [`t_closeness`](@ref) |
| What is the earliest reachable time? | [`temporal_distance`](@ref), [`forward_reachable_set`](@ref) |
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

[`temporal_distance`](@ref) reports this as the *elapsed* time
$d_T(s, r, t_0) - t_0$, and returns `nothing` when no such path exists.

## Module Reference

```@docs
TSNA
```

## References

1. Bender-deMoll, S., Morris, M. (2012). `tsna`: Tools for Temporal Social Network Analysis. R package.

2. Holme, P. (2015). Modern temporal network theory: a colloquium. *European Physical Journal B*, 88(9), 1-30.

3. Holme, P., Saramaki, J. (2012). Temporal networks. *Physics Reports*, 519(3), 97-125.

4. Nicosia, V., Tang, J., Mascolo, C., Musolesi, M., Russo, G., Latora, V. (2013). Graph metrics for temporal networks. In *Temporal Networks* (pp. 15-40). Springer.

## Citation

If you use TSNA.jl in your work, please cite it using the entry in
[`CITATION.bib`](https://github.com/statistical-network-analysis-with-Julia/TSNA.jl/blob/main/CITATION.bib):

```biblatex
@misc{SNWJTSNAJL,
  author = {{Statistical Network Analysis with Julia}},
  title = {TSNA.jl: Temporal Social Network Analysis for Julia},
  year = {2026},
  url = {https://github.com/statistical-network-analysis-with-Julia/TSNA.jl},
  note = {Homepage: https://statistical-network-analysis-with-Julia.github.io/TSNA.jl; GitHub: https://github.com/statistical-network-analysis-with-Julia}
}
```
