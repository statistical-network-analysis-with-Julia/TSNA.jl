# Temporal Paths

This guide covers temporal path analysis in TSNA.jl, including time-respecting paths, temporal distance, and reachability analysis.

## What are Temporal Paths?

In a static network, a path from $A$ to $C$ through $B$ is valid whenever edges $A \to B$ and $B \to C$ exist. In a temporal network, a path must also respect time -- you can only traverse an edge while it is active, and edges must be traversed in non-decreasing time order.

```text
Static:     A → B → C    (always valid if both edges exist)

Temporal:   A --(t=5)--> B --(t=10)--> C    (valid: 5 <= 10)
            A --(t=5)--> B --(t=3)--> C     (INVALID: 3 < 5, going back in time)
```

This constraint fundamentally changes reachability. In temporal networks:

- Reachability depends on the **start time**
- Reachability is **not symmetric**: A can reach B does not imply B can reach A
- Reachability is **not transitive**: A reaches B and B reaches C does not guarantee A reaches C

## The tPath Type

Temporal paths are represented by the `tPath{T, Time}` type:

```julia
struct tPath{T, Time}
    vertices::Vector{T}       # Sequence of vertices visited
    times::Vector{Time}       # Time at each edge traversal
    edges::Vector{Tuple{T,T}} # Edges traversed
end
```

### Properties

- `length(path)` returns the number of edges (hops)
- `path.vertices[1]` is the source
- `path.vertices[end]` is the target
- `path.times[i]` is the time of the $i$-th edge traversal
- `path.edges[i]` is the $i$-th edge

### Example

```julia
# A path: 1 --(t=5)--> 2 --(t=10)--> 3 --(t=15)--> 4
path = tPath(
    [1, 2, 3, 4],                    # vertices
    [5.0, 10.0, 15.0],               # times
    [(1, 2), (2, 3), (3, 4)]         # edges
)

println(path)
# tPath: 1 --(5.0)--> 2 --(10.0)--> 3 --(15.0)--> 4

println("Length: ", length(path.edges))  # 3 edges
println("Source: ", path.vertices[1])  # 1
println("Target: ", path.vertices[end])  # 4
println("Departure: ", path.times[1])  # 5.0
println("Arrival: ", path.times[end])  # 15.0
```

## Temporal Distance

### Finding Earliest Arrival

The `temporal_distance` function finds the earliest time at which a target vertex can be reached from a source, starting at a given time:

```julia
using NetworkDynamic
using TSNA

# Create network with sequential edges
dnet = DynamicNetwork(5; observation_start=0.0, observation_end=100.0)
activate_vertices!(dnet, collect(1:5), 0.0, 100.0)

activate!(dnet, 0.0, 20.0; edge=(1, 2))   # Active early
activate!(dnet, 10.0, 40.0; edge=(2, 3))  # Overlaps
activate!(dnet, 30.0, 60.0; edge=(3, 4))  # Later
activate!(dnet, 50.0, 80.0; edge=(4, 5))  # Even later

# Earliest arrival from 1 to 5, starting at t=0
dist = temporal_distance(dnet, 1, 5, 0.0)
println("Earliest arrival at v5: $dist")
# Path: 1→2 at t=0, 2→3 at t=10, 3→4 at t=30, 4→5 at t=50
# Arrival: t=50

# Starting later
dist2 = temporal_distance(dnet, 1, 5, 25.0)
println("Starting at t=25: $dist2")
# Edge (1,2) is still active at t=25, but we need subsequent edges
```

### When No Path Exists

If no temporal path exists, `temporal_distance` returns `nothing`:

```julia
# No path exists if we start too late
dist = temporal_distance(dnet, 1, 5, 90.0)
if dist === nothing
    println("No temporal path exists starting at t=90")
end
```

### Self-Distance

The distance from a vertex to itself is the start time:

```julia
dist = temporal_distance(dnet, 1, 1, 5.0)
println(dist)  # 5.0
```

## Shortest Temporal Path

### Finding the Path

The `shortest_temporal_path` function returns the actual path, not just the arrival time:

```julia
path = shortest_temporal_path(dnet, 1, 5, 0.0)

if !isnothing(path)
    println("Path found!")
    println("Vertices: ", path.vertices)
    println("Times: ", path.times)
    println("Edges: ", path.edges)
    println("Total hops: ", length(path))
else
    println("No path exists")
end
```

### Path to Self

```julia
path = shortest_temporal_path(dnet, 1, 1, 0.0)
println(path)
# tPath: 1  (trivial path, no edges)
println(length(path))  # 0
```

### Algorithm

The shortest temporal path algorithm works as follows:

1. Initialize earliest arrival times: $\text{arrival}[s] = t_{\text{start}}$, all others $= \infty$
2. Sort all edge spells by onset time (ascending)
3. For each edge activation $(i, j)$ at time $t$:
   - If vertex $i$ has been reached by time $t$ (i.e., $\text{arrival}[i] \leq t$):
     - If $t < \text{arrival}[j]$: update $\text{arrival}[j] = t$ and record predecessor
4. Reconstruct path from predecessor chain

This is a temporal adaptation of Dijkstra's algorithm, taking advantage of the fact that edge activations can be processed in time order.

## Reachability Analysis

### Forward Reachability

The forward reachable set contains all vertices that can be reached from a source starting at a given time:

```julia
using Network   # for nv

# Who can vertex 1 reach starting at t=0?
reachable = forward_reachable_set(dnet, 1, 0.0)
println("Forward reachable from v1 at t=0: ", reachable)
println("Fraction reachable: $(length(reachable))/$(nv(dnet))")
```

### Backward Reachability

The backward reachable set contains all vertices that can reach a target by a given time:

```julia
# Who can reach vertex 5 by t=100?
sources = backward_reachable_set(dnet, 5, 100.0)
println("Backward reachable to v5 by t=100: ", sources)
```

### Reachability Properties

```julia
# Forward reachability always includes the source
reach = forward_reachable_set(dnet, 1, 0.0)
@assert 1 in reach

# Backward reachability always includes the target
sources = backward_reachable_set(dnet, 5, 100.0)
@assert 5 in sources
```

### Reachability Over Time

Track how reachability changes with the start time:

```julia
source = 1
println("Forward reachability from v$source:")
for t_start in 0.0:10.0:90.0
    reach = forward_reachable_set(dnet, source, t_start)
    println("  Start at t=$t_start: $(length(reach)) vertices reachable")
end
```

As the start time increases, fewer future edge activations are available, so reachability typically decreases.

## Practical Examples

### Information Spread Analysis

Analyze how quickly information can spread through a network:

```julia
dnet = DynamicNetwork(20; observation_start=0.0, observation_end=100.0)
activate_vertices!(dnet, collect(1:20), 0.0, 100.0)

# Add edges (communication events)
# ... (add many edges with various activity spells)

# How many vertices can be reached from vertex 1 at different times?
source = 1
println("Information spread from vertex $source:")
for t in 0.0:10.0:100.0
    reach = forward_reachable_set(dnet, source, t)
    println("  Starting at t=$t: $(length(reach)) vertices reachable")
end

# What fraction of all pairs are temporally connected?
n = nv(dnet)
connected_pairs = 0
total_pairs = n * (n - 1)

for s in 1:n
    reach = forward_reachable_set(dnet, s, 0.0)
    connected_pairs += length(reach) - 1  # Exclude self
end

println("Temporal connectivity: $(connected_pairs)/$(total_pairs) = ",
        "$(round(connected_pairs/total_pairs, digits=3))")
```

### Earliest Delivery Analysis

Find the earliest time each vertex can receive information from a source:

```julia
source = 1
start_time = 0.0

println("Earliest arrival from vertex $source:")
for target in 1:nv(dnet)
    target == source && continue
    dist = temporal_distance(dnet, source, target, start_time)
    if dist === nothing
        println("  v$target: unreachable")
    else
        println("  v$target: arrives at t=$dist")
    end
end
```

### Path Analysis

Compare paths from different starting times:

```julia
source = 1
target = 10

println("Paths from v$source to v$target:")
for t_start in 0.0:20.0:80.0
    path = shortest_temporal_path(dnet, source, target, t_start)
    if !isnothing(path)
        println("  Start t=$t_start: $(length(path)) hops, ",
                "arrival at t=$(path.times[end])")
    else
        println("  Start t=$t_start: no path")
    end
end
```

### Bidirectional Reachability

Find vertices that can both reach and be reached by a focal vertex:

```julia
focal = 1

# Forward: who can focal reach?
forward = forward_reachable_set(dnet, focal, 0.0)

# Backward: who can reach focal?
backward = backward_reachable_set(dnet, focal, 100.0)

# Bidirectional: both forward and backward reachable
bidirectional = intersect(forward, backward)

println("Forward reachable: $(length(forward))")
println("Backward reachable: $(length(backward))")
println("Bidirectional: $(length(bidirectional))")
```

## Computational Complexity

| Function | Complexity | Notes |
|----------|-----------|-------|
| `temporal_distance` | $O(E \log E)$ | $E$ = total edge spells |
| `shortest_temporal_path` | $O(E \log E + P)$ | $P$ = path reconstruction |
| `forward_reachable_set` | $O(E \log E)$ | Same as temporal_distance |
| `backward_reachable_set` | $O(E \log E)$ | Reverse time processing |

All algorithms sort edge spells by onset time, then process them in order. The dominant cost is the sort, making them efficient for networks with many vertices but moderate numbers of edge spells.

## Best Practices

1. **Choose appropriate start times**: Reachability depends strongly on the start time
2. **Consider direction**: In directed networks, forward and backward reachability differ
3. **Check for Inf**: Always check if `temporal_distance` returns `Inf` before using the result
4. **Use forward reachability for spread**: Forward reachable sets model information or disease spread
5. **Use backward reachability for influence**: Backward reachable sets identify potential sources of influence
6. **Track over time**: Compute reachability at multiple start times to understand dynamics
