# Getting Started

This tutorial walks through common use cases for TSNA.jl, from computing temporal centrality to analyzing reachability and network dynamics.

## Installation

Install TSNA.jl from GitHub:

```julia
using Pkg
Pkg.add(url="https://github.com/Statistical-network-analysis-with-Julia/TSNA.jl")
```

TSNA.jl depends on NetworkDynamic.jl, Network.jl, and SNA.jl, which will be installed automatically.

## Basic Workflow

The typical TSNA.jl workflow consists of four steps:

1. **Create or load a dynamic network** - Using NetworkDynamic.jl
2. **Compute temporal measures** - Centrality, density, and other metrics at specific times
3. **Analyze temporal paths** - Find reachable sets and shortest paths
4. **Compute duration metrics** - Edge persistence, turnover, and decay rates

## Step 1: Create a Dynamic Network

TSNA.jl operates on `DynamicNetwork` objects from NetworkDynamic.jl:

```julia
using Network   # nv, ne on snapshots
using NetworkDynamic
using TSNA

# Create a dynamic network with 10 vertices
dnet = DynamicNetwork(10;
    observation_start=0.0,
    observation_end=100.0
)

# Activate all vertices for the entire period
activate_vertices!(dnet, collect(1:10), 0.0, 100.0)

# Add edges with different activity periods
activate!(dnet, 0.0, 40.0; edge=(1, 2))    # Active early
activate!(dnet, 10.0, 60.0; edge=(2, 3))   # Overlaps with (1,2)
activate!(dnet, 30.0, 80.0; edge=(3, 4))   # Overlaps with (2,3)
activate!(dnet, 50.0, 100.0; edge=(4, 5))  # Active late
activate!(dnet, 20.0, 70.0; edge=(1, 5))   # Shortcut
activate!(dnet, 0.0, 100.0; edge=(5, 6))   # Always active
activate!(dnet, 5.0, 95.0; edge=(6, 7))    # Nearly always active
activate!(dnet, 60.0, 90.0; edge=(7, 8))   # Active mid-late
activate!(dnet, 0.0, 50.0; edge=(8, 9))    # Active early
activate!(dnet, 25.0, 75.0; edge=(9, 10))  # Active middle
```

## Step 2: Compute Temporal Measures

### Centrality at a Point in Time

Compute standard centrality measures on the network snapshot at a specific time:

```julia
# Degree centrality at t=30
deg = tDegree(dnet, 30.0)
println("Degree at t=30: ", deg)

# Directed degree
in_deg = tDegree(dnet, 30.0; mode=:in)
out_deg = tDegree(dnet, 30.0; mode=:out)

# Betweenness centrality at t=30
bet = tBetweenness(dnet, 30.0)
println("Betweenness at t=30: ", round.(bet, digits=3))

# Closeness centrality at t=30
clo = tCloseness(dnet, 30.0)
println("Closeness at t=30: ", round.(clo, digits=3))
```

### Network-Level Measures

```julia
# Density at different times
for t in 0.0:20.0:100.0
    d = tDensity(dnet, t)
    println("t=$t: density = $(round(d, digits=3))")
end

# Reciprocity (directed networks)
r = tReciprocity(dnet, 50.0)
println("Reciprocity at t=50: $(round(r, digits=3))")

# Transitivity (clustering coefficient)
tr = tTransitivity(dnet, 50.0)
println("Transitivity at t=50: $(round(tr, digits=3))")
```

### Time Series of Statistics

```julia
# Compute density and reciprocity at regular intervals
times = collect(0.0:10.0:100.0)
stats = tSnaStats(dnet, times; measures=[:density, :reciprocity, :n_edges])

println("Time\tDensity\tReciprocity\tEdges")
for row in stats   # one NamedTuple per time point
    d = round(row.density, digits=3)
    r = round(row.reciprocity, digits=3)
    e = Int(row.n_edges)
    println("$(row.time)\t$d\t$r\t\t$e")
end
```

## Step 3: Analyze Temporal Paths

### Temporal Distance

Find the earliest arrival time from one vertex to another:

```julia
# How quickly can information travel from vertex 1 to vertex 5?
dist = temporalDistance(dnet, 1, 5, 0.0)
println("Earliest arrival at v5 from v1 starting at t=0: $dist")

# Starting at a later time
dist_late = temporalDistance(dnet, 1, 5, 50.0)
println("Earliest arrival at v5 from v1 starting at t=50: $dist_late")
```

### Shortest Temporal Path

Find the actual path, not just the arrival time:

```julia
path = shortestTemporalPath(dnet, 1, 5, 0.0)

if !isnothing(path)
    println("Path found: ", path)
    println("Path length: ", length(path), " edges")
    println("Vertices: ", path.vertices)
    println("Times: ", path.times)
else
    println("No temporal path exists")
end
```

### Forward Reachability

Find all vertices reachable from a source:

```julia
# Who can vertex 1 reach starting at t=0?
reachable = forwardReachableSet(dnet, 1, 0.0)
println("Vertices reachable from v1: ", reachable)
println("Reachability: $(length(reachable))/$(nv(dnet))")
```

### Backward Reachability

Find all vertices that can reach a target:

```julia
# Who can reach vertex 10 by t=100?
sources = backwardReachableSet(dnet, 10, 100.0)
println("Vertices that can reach v10: ", sources)
```

## Step 4: Compute Duration Metrics

### Edge Duration

```julia
# Mean edge duration
mean_dur = tEdgeDuration(dnet; aggregate=:mean)
println("Mean edge duration: $(round(mean_dur, digits=1))")

# Median edge duration
med_dur = tEdgeDuration(dnet; aggregate=:median)
println("Median edge duration: $(round(med_dur, digits=1))")

# Total across all edges
total_dur = tEdgeDuration(dnet; aggregate=:total)
println("Total edge-time: $(round(total_dur, digits=1))")

# Per-edge durations
all_durs = tEdgeDuration(dnet; aggregate=:all)
for (edge, dur) in all_durs
    println("Edge $(edge[1])->$(edge[2]): duration = $dur")
end
```

### Vertex Duration

```julia
# Mean vertex activity duration
v_dur = tVertexDuration(dnet; aggregate=:mean)
println("Mean vertex duration: $(round(v_dur, digits=1))")
```

### Edge Persistence

Measure how many edges persist across time windows:

```julia
# What fraction of edges persist across 20-unit windows?
persistence = tEdgePersistence(dnet, 20.0)
println("Edge persistence (window=20): $(round(persistence, digits=3))")

# Try different window sizes
for w in [10.0, 20.0, 30.0, 50.0]
    p = tEdgePersistence(dnet, w)
    println("Window $w: persistence = $(round(p, digits=3))")
end
```

### Turnover

Compute edge formation and dissolution rates:

```julia
turnover = tTurnover(dnet, 20.0)
println("Formation rate: $(round(turnover.formation_rate, digits=4))")
println("Dissolution rate: $(round(turnover.dissolution_rate, digits=4))")
println("New edges: $(turnover.n_formations)")
println("Lost edges: $(turnover.n_dissolutions)")
```

### Tie Decay

Estimate the rate at which ties decay:

```julia
# Exponential decay rate
decay = tieDecay(dnet; method=:exponential)
println("Decay rate: $(round(decay, digits=4))")

# Halflife method
halflife_rate = tieDecay(dnet; method=:halflife)
println("Halflife decay rate: $(round(halflife_rate, digits=4))")
```

## Complete Example

```julia
using NetworkDynamic
using TSNA

# Build a dynamic network representing office communication
n = 8
dnet = DynamicNetwork(n; observation_start=0.0, observation_end=50.0)
activate_vertices!(dnet, collect(1:n), 0.0, 50.0)

# Morning communications (0-20)
activate!(dnet, 0.0, 15.0; edge=(1, 2))
activate!(dnet, 5.0, 20.0; edge=(2, 3))
activate!(dnet, 3.0, 18.0; edge=(1, 4))
activate!(dnet, 8.0, 20.0; edge=(3, 4))

# Afternoon communications (15-40)
activate!(dnet, 15.0, 35.0; edge=(4, 5))
activate!(dnet, 20.0, 40.0; edge=(5, 6))
activate!(dnet, 25.0, 40.0; edge=(6, 7))
activate!(dnet, 18.0, 38.0; edge=(5, 8))

# Late communications (30-50)
activate!(dnet, 30.0, 50.0; edge=(7, 8))
activate!(dnet, 35.0, 50.0; edge=(8, 1))  # Feedback loop
activate!(dnet, 32.0, 48.0; edge=(2, 6))  # Cross-team link

# === Analysis ===

# 1. Track centrality over time
println("=== Degree Centrality Over Time ===")
for t in [5.0, 15.0, 25.0, 35.0, 45.0]
    deg = tDegree(dnet, t)
    top_v = argmax(deg)
    println("t=$t: most central vertex = $top_v (degree=$(deg[top_v]))")
end

# 2. Reachability analysis
println("\n=== Reachability from Vertex 1 ===")
for t_start in [0.0, 10.0, 20.0, 30.0]
    reach = forwardReachableSet(dnet, 1, t_start)
    println("Starting at t=$t_start: $(length(reach)) vertices reachable")
end

# 3. Temporal paths
println("\n=== Shortest Paths from v1 ===")
for target in [4, 6, 8]
    path = shortestTemporalPath(dnet, 1, target, 0.0)
    if !isnothing(path)
        println("v1 -> v$target: $(length(path)) edges, arrival at t=$(path.times[end])")
    else
        println("v1 -> v$target: no path")
    end
end

# 4. Duration metrics
println("\n=== Duration Metrics ===")
println("Mean edge duration: $(round(tEdgeDuration(dnet; aggregate=:mean), digits=1))")
println("Edge persistence (window=10): $(round(tEdgePersistence(dnet, 10.0), digits=3))")

turnover = tTurnover(dnet, 10.0)
println("Formation rate: $(round(turnover.formation_rate, digits=4))")
println("Dissolution rate: $(round(turnover.dissolution_rate, digits=4))")

# 5. Network evolution
println("\n=== Network Evolution ===")
times = collect(0.0:5.0:50.0)
stats = tSnaStats(dnet, times; measures=[:density, :n_edges, :mean_degree])
for row in stats
    println("t=$(row.time): edges=$(Int(row.n_edges)), ",
            "density=$(round(row.density, digits=3)), ",
            "mean_deg=$(round(row.mean_degree, digits=2))")
end
```

## Working with Contact Sequences

Convert a dynamic network to a sequence of contacts:

```julia
cs = as_contact_sequence(dnet)
println("Number of contacts: ", length(cs))

for contact in cs
    println("$(contact.source) -> $(contact.target) at t=$(contact.time), ",
            "duration=$(contact.duration)")
end
```

## Aggregating to Static Networks

Collapse a dynamic network to static using different methods:

```julia
# Union: include any edge that was ever active
static_union = tAggregate(dnet; method=:union)
println("Union: $(ne(static_union)) edges")

# Intersection: include only edges active throughout
static_inter = tAggregate(dnet; method=:intersection)
println("Intersection: $(ne(static_inter)) edges")

# Weighted: weight by total activity time
static_weighted = tAggregate(dnet; method=:weighted)
println("Weighted: $(ne(static_weighted)) edges")
```

## Best Practices

1. **Define observation period**: Always set `observation_start` and `observation_end` on the dynamic network
2. **Activate vertices**: Ensure all vertices are active before querying centrality at specific times
3. **Check for empty snapshots**: Some time points may have no active vertices or edges
4. **Use appropriate window sizes**: For persistence and turnover, window size should match your analysis timescale
5. **Compare with static analysis**: Temporal measures should reveal dynamics invisible to static SNA
6. **Consider edge direction**: Temporal paths in directed networks can only follow edge direction

## Next Steps

- Learn about [Temporal Centrality](guide/centrality.md) measures in detail
- Understand [Temporal Paths](guide/paths.md) and reachability analysis
- Explore [Duration Metrics](guide/metrics.md) for measuring network dynamics
