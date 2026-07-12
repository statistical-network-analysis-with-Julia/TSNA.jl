# Temporal Centrality

This guide covers the temporal centrality measures available in TSNA.jl. These measures compute standard centrality metrics on network snapshots extracted at specific time points.

## Overview

Temporal centrality measures answer the question: "Who is most central **at a given time**?" Unlike static centrality, temporal centrality can change as the network evolves -- a vertex that is highly central at one time may be peripheral at another.

All temporal centrality functions in TSNA.jl follow the same pattern:

<!-- skip-check -->
```julia
centrality_values = t_centrality(dnet, at_time; options...)
```

They work by:

1. Extracting a static network snapshot at the specified time using
   `network_extract` with `retain_all_vertices=true`
2. Computing the standard centrality measure on the snapshot using SNA.jl
3. Returning a vector of centrality values of length `nv(dnet)`, indexed
   by the network's own vertex IDs (inactive vertices score 0)

## Degree Centrality

### t_degree

Degree centrality counts the number of edges incident to each vertex:

```julia
using NetworkDynamic
using TSNA

# Create dynamic network
dnet = DynamicNetwork(10; observation_start=0.0, observation_end=100.0)
activate_vertices!(dnet, collect(1:10), 0.0, 100.0)
activate!(dnet, 0.0, 50.0; edge=(1, 2))
activate!(dnet, 0.0, 50.0; edge=(1, 3))
activate!(dnet, 0.0, 100.0; edge=(2, 3))
activate!(dnet, 50.0, 100.0; edge=(4, 5))

# Total degree at t=25
deg = t_degree(dnet, 25.0)
println("Degree at t=25: ", deg)
```

### Directed Degree Modes

For directed networks, specify the degree mode:

```julia
# In-degree (number of incoming edges)
in_deg = t_degree(dnet, 25.0; mode=:in)

# Out-degree (number of outgoing edges)
out_deg = t_degree(dnet, 25.0; mode=:out)

# Total degree (in + out, default)
total_deg = t_degree(dnet, 25.0; mode=:total)
```

| Mode | Counts | Interpretation |
|------|--------|----------------|
| `:in` | Incoming edges | Popularity, prestige |
| `:out` | Outgoing edges | Activity, influence |
| `:total` | Both directions | Overall connectivity |

### Degree Over Time

Track how a vertex's degree changes:

```julia
times = collect(0.0:10.0:100.0)
v = 1  # Track vertex 1

println("Vertex $v degree over time:")
for t in times
    deg = t_degree(dnet, t)
    if v <= length(deg)
        println("  t=$t: degree=$(deg[v])")
    end
end
```

## Betweenness Centrality

### t_betweenness

Betweenness centrality measures how often a vertex lies on shortest paths between other vertices:

```julia
bet = t_betweenness(dnet, 50.0)
println("Betweenness at t=50: ", round.(bet, digits=3))

# Normalized by (n-1)(n-2)
bet_norm = t_betweenness(dnet, 50.0; normalized=true)
```

| Parameter | Description | Default |
|-----------|-------------|---------|
| `normalized` | Normalize by $(n-1)(n-2)$ | `false` (raw scores, as in SNA.jl) |

**Interpretation**: A vertex with high betweenness acts as a bridge or broker between different parts of the network. Temporal betweenness reveals when these brokerage positions exist.

### Temporal Brokerage Analysis

```julia
using Network   # for nv

# Track betweenness over time to find when vertices serve as brokers
times = collect(0.0:5.0:100.0)
n = nv(dnet)

println("Peak brokerage times:")
for v in 1:n
    max_bet = 0.0
    max_t = 0.0
    for t in times
        bet = t_betweenness(dnet, t)
        if v <= length(bet) && bet[v] > max_bet
            max_bet = bet[v]
            max_t = t
        end
    end
    if max_bet > 0
        println("  Vertex $v: peak betweenness=$(round(max_bet, digits=3)) at t=$max_t")
    end
end
```

## Closeness Centrality

### t_closeness

Closeness centrality measures how close a vertex is to all other vertices (inverse of average shortest path length):

```julia
clo = t_closeness(dnet, 50.0)
println("Closeness at t=50: ", round.(clo, digits=3))
```

**Interpretation**: High closeness indicates a vertex can quickly reach (or be reached by) all others. In temporal networks, closeness can change dramatically as edges activate and deactivate.

### Comparing Centrality Measures

```julia
t = 50.0
deg = t_degree(dnet, t)
bet = t_betweenness(dnet, t)
clo = t_closeness(dnet, t)

println("Vertex\t_degree\tBetween\t_closeness")
for v in 1:length(deg)
    println("$v\t$(deg[v])\t$(round(bet[v], digits=3))\t$(round(clo[v], digits=3))")
end
```

## Eigenvector Centrality

### t_eigenvector

Eigenvector centrality measures a vertex's influence based on the influence of its neighbors:

```julia
eig = t_eigenvector(dnet, 50.0)
println("Eigenvector centrality at t=50: ", round.(eig, digits=3))
```

**Interpretation**: A vertex has high eigenvector centrality when it is connected to other well-connected vertices. This captures "influence" in a recursive sense.

**Note**: Eigenvector centrality may not converge for disconnected networks. If the snapshot at time $t$ is disconnected, results should be interpreted with caution.

## PageRank

### t_pagerank

PageRank is a variant of eigenvector centrality with a damping factor, originally designed for ranking web pages:

```julia
pr = t_pagerank(dnet, 50.0)
println("PageRank at t=50: ", round.(pr, digits=3))

# With custom damping factor
pr_low = t_pagerank(dnet, 50.0; damping=0.5)
pr_high = t_pagerank(dnet, 50.0; damping=0.95)
```

| Parameter | Description | Default |
|-----------|-------------|---------|
| `damping` | Probability of following an edge (vs. random jump) | `0.85` |

**Interpretation**: PageRank handles disconnected networks better than eigenvector centrality due to the random jump component. Lower damping means more random jumps, leading to more uniform values.

## Network-Level Temporal Measures

### t_density

Network density at a specific time:

```julia
d = t_density(dnet, 50.0)
println("Density at t=50: $(round(d, digits=3))")
```

Density equals the number of active edges divided by the maximum possible edges.

### t_reciprocity

Proportion of edges that are reciprocated (directed networks):

```julia
r = t_reciprocity(dnet, 50.0)
println("Reciprocity at t=50: $(round(r, digits=3))")
```

For undirected networks, reciprocity is always 1.0.

### t_transitivity

Global clustering coefficient (transitivity):

```julia
tr = t_transitivity(dnet, 50.0)
println("Transitivity at t=50: $(round(tr, digits=3))")
```

Transitivity measures the proportion of connected triples that form triangles.

## Computing Statistics Over Time

### Point-in-Time Statistics

Use `t_sna_stats` to compute multiple statistics at a series of time points:

```julia
times = collect(0.0:10.0:100.0)
stats = t_sna_stats(dnet, times;
    measures=[:density, :reciprocity, :transitivity, :n_edges, :mean_degree]
)

# Access results (one NamedTuple row per time point)
densities = [row.density for row in stats]
reciprocities = [row.reciprocity for row in stats]
edge_counts = [row.n_edges for row in stats]

# Print table
println("Time\t_density\tRecip\tEdges\tMean Deg")
for row in stats
    println("$(row.time)\t$(round(row.density, digits=3))\t",
            "$(round(row.reciprocity, digits=3))\t",
            "$(Int(row.n_edges))\t",
            "$(round(row.mean_degree, digits=2))")
end
```

### Available Statistics

| Symbol | Description |
|--------|-------------|
| `:density` | Network density |
| `:reciprocity` | Proportion of reciprocated edges |
| `:transitivity` | Global clustering coefficient |
| `:n_edges` | Number of active edges |
| `:mean_degree` | Mean degree (`2 * edges / vertices` undirected, `edges / vertices` directed) |

### Window-Based Statistics

Use `window_sna_stats` to sample statistics at the start of consecutive
windows spanning the observation period (a convenience wrapper around
`t_sna_stats` — no need to build the time grid yourself):

```julia
# Statistics sampled every 10 time units
stats = window_sna_stats(dnet, 10.0;
    measures=[:density, :n_edges]
)

println("Window densities: ", round.([row.density for row in stats], digits=3))
println("Window edge counts: ", [Int(row.n_edges) for row in stats])
```

Each window of length `window_size` contributes one snapshot taken at the
window's start time.

## Practical Considerations

### Empty Snapshots

At some time points, no vertices or edges may be active:

```julia
deg = t_degree(dnet, 200.0)  # After observation period
# Returns a length-nv(dnet) vector of zeros
```

Check that the snapshot has active edges before interpreting centrality
scores.

### Vertex Indexing

Snapshots retain **every** vertex (`retain_all_vertices=true`), so the
returned vectors always have length `nv(dnet)` and are indexed by the
dynamic network's own vertex IDs — even when some vertices are inactive
at the query time (they simply score 0):

```julia
# If only vertices 3, 5, 7 are active at t=50
deg = t_degree(dnet, 50.0)
# length(deg) == nv(dnet); deg[3], deg[5], deg[7] carry the activity,
# all other entries are 0
```

### Disconnected Components

Some centrality measures (closeness, eigenvector) may behave unexpectedly on disconnected networks. Consider:

- Using PageRank instead of eigenvector centrality for disconnected networks
- Restricting analysis to the largest connected component
- Checking the number of components before computing centrality

## Best Practices

1. **Check snapshot size**: Verify the extracted network has enough vertices for meaningful centrality
2. **Use appropriate measures**: Degree for local connectivity, betweenness for brokerage, PageRank for global influence
3. **Track over time**: Compute centrality at multiple time points to reveal dynamics
4. **Compare measures**: Different centrality measures highlight different aspects of temporal importance
5. **Consider direction**: Use `:in` and `:out` modes for directed networks to distinguish popularity from activity
6. **Grid vs. custom times**: Use `t_sna_stats` with your own time points, `window_sna_stats` for a regular grid spanning the observation period
