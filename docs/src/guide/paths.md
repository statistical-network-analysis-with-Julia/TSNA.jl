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

## The TemporalPath Type

Temporal paths are represented by the `TemporalPath{T, Time}` type (also
exported under the R-style alias `tPath`):

```julia
struct TemporalPath{T, Time}
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
- `path_duration(path)` is the elapsed time between the first and last
  traversal (a proper time difference, also for `Date`/`DateTime` axes)

### Example

```julia
# A path: 1 --(t=5)--> 2 --(t=10)--> 3 --(t=15)--> 4
path = TemporalPath(
    [1, 2, 3, 4],                    # vertices
    [5.0, 10.0, 15.0],               # times
    [(1, 2), (2, 3), (3, 4)]         # edges
)

println(path)
# TemporalPath: 1 --(5.0)--> 2 --(10.0)--> 3 --(15.0)--> 4

println("Length: ", length(path))  # 3 edges
println("Source: ", path.vertices[1])  # 1
println("Target: ", path.vertices[end])  # 4
println("Departure: ", path.times[1])  # 5.0
println("Arrival: ", path.times[end])  # 15.0
println("Duration: ", path_duration(path))  # 10.0
```

## Temporal Distance

### Finding the Fastest Route

The `temporal_distance` function returns the **elapsed time** of the
earliest time-respecting path from a source to a target — the earliest
possible arrival time minus the start time:

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

# Elapsed time from 1 to 5, starting at t=0
dist = temporal_distance(dnet, 1, 5, 0.0)
println("Fastest route to v5 takes: $dist")
# Path: 1→2 at t=0, 2→3 at t=10, 3→4 at t=30, 4→5 at t=50
# Arrival at t=50, so elapsed time = 50.0 - 0.0 = 50.0

# Starting later can be faster: edges can be boarded mid-spell,
# so a late start skips the waiting
dist2 = temporal_distance(dnet, 1, 5, 15.0)
println("Starting at t=15: $dist2")
# Path: 1→2 at t=15, 2→3 at t=15, 3→4 at t=30, 4→5 at t=50
# Arrival at t=50, so elapsed time = 50.0 - 15.0 = 35.0
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

The distance from a vertex to itself is zero:

```julia
dist = temporal_distance(dnet, 1, 1, 5.0)
println(dist)  # 0.0
```

## Earliest-Arrival Temporal Path

### Finding the Path

The `temporal_path` function returns the actual earliest-arrival path,
not just the elapsed time (`shortest_temporal_path` is an alias kept for
API compatibility — note it is the *earliest-arrival* path, not the
fewest-hops path):

```julia
path = temporal_path(dnet, 1, 5, 0.0)

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
path = temporal_path(dnet, 1, 1, 0.0)
println(path)
# TemporalPath: 1  (trivial path, no edges)
println(length(path))  # 0
```

### Algorithm

All path and reachability functions delegate to `earliest_arrival`, a
heap-based Dijkstra label-setting search over a memoized per-vertex
contact index:

1. Initialize earliest arrival times: $\text{arrival}[s] = t_{\text{start}}$, all others unknown
2. Push $(t_{\text{start}}, s)$ onto a binary min-heap keyed by arrival time
3. Pop the vertex $v$ with the smallest arrival time $t$ (its label is now final);
   for each outgoing contact $(v, w)$ with spell $[\text{onset}, \text{terminus})$:
   - The spell can be boarded at $\max(t, \text{onset})$ provided $t < \text{terminus}$
     (interval semantics: boarding mid-spell is allowed, and spells that
     began before the start time but are still active count)
   - If that boarding time improves $\text{arrival}[w]$: update it, record
     the predecessor, and push $w$ onto the heap
4. Reconstruct the path from the predecessor chain

Because boarding times never precede the label being settled, labels are
final once popped — this is Dijkstra's algorithm with temporal edge
availability. The per-vertex contact index is memoized per network (and
invalidated on mutation), so repeated queries on the same network skip
rebuilding it. With a target given, the search stops as soon as the
target is settled.

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
    path = temporal_path(dnet, source, target, t_start)
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
| `earliest_arrival` | $O(E \log V)$ | $E$ = total edge spells, $V$ = vertices; stops early with a `target` |
| `temporal_distance` | $O(E \log V)$ | One targeted `earliest_arrival` search |
| `temporal_path` | $O(E \log V + P)$ | $P$ = path reconstruction |
| `forward_reachable_set` | $O(E \log V)$ | One full `earliest_arrival` search |
| `backward_reachable_set` | $O(V \cdot E \log V)$ | One targeted forward search per vertex |

All functions delegate to the heap-based `earliest_arrival` search over a
memoized per-vertex contact index. Building the index costs $O(E)$ but is
cached per network (invalidated on mutation), so repeated queries — such
as the per-vertex searches inside `backward_reachable_set` — do not
rebuild it.

## Best Practices

1. **Choose appropriate start times**: Reachability depends strongly on the start time
2. **Consider direction**: In directed networks, forward and backward reachability differ
3. **Check for `nothing`**: Always check if `temporal_distance` (or `temporal_path`) returned `nothing` before using the result
4. **Use forward reachability for spread**: Forward reachable sets model information or disease spread
5. **Use backward reachability for influence**: Backward reachable sets identify potential sources of influence
6. **Track over time**: Compute reachability at multiple start times to understand dynamics
