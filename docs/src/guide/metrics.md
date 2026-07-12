# Duration Metrics

This guide covers the duration, persistence, turnover, and aggregation metrics available in TSNA.jl for characterizing network dynamics.

## Overview

While centrality and path analysis describe the network at specific times, duration metrics characterize **how the network changes over time**. They answer questions like:

- How long do ties last?
- How stable is the network structure?
- At what rate are new ties forming and old ties dissolving?
- Is the network becoming denser or sparser over time?

## Edge Duration

### tEdgeDuration

Compute the total duration of edge activity:

```julia
using NetworkDynamic
using TSNA

dnet = DynamicNetwork(5; observation_start=0.0, observation_end=100.0)
activate_vertices!(dnet, collect(1:5), 0.0, 100.0)

activate!(dnet, 0.0, 50.0; edge=(1, 2))    # Duration: 50
activate!(dnet, 10.0, 30.0; edge=(2, 3))   # Duration: 20
activate!(dnet, 40.0, 80.0; edge=(3, 4))   # Duration: 40
activate!(dnet, 0.0, 100.0; edge=(4, 5))   # Duration: 100
activate!(dnet, 20.0, 40.0; edge=(1, 3))   # Duration: 20

# Mean edge duration
mean_dur = tEdgeDuration(dnet; aggregate=:mean)
println("Mean edge duration: $mean_dur")  # (50+20+40+100+20)/5 = 46.0

# Median edge duration
med_dur = tEdgeDuration(dnet; aggregate=:median)
println("Median edge duration: $med_dur")  # 40.0

# Total edge-time
total_dur = tEdgeDuration(dnet; aggregate=:total)
println("Total edge-time: $total_dur")  # 230.0

# Per-edge durations
all_durs = tEdgeDuration(dnet; aggregate=:all)
for (edge, dur) in all_durs
    println("  Edge $(edge[1])->$(edge[2]): $dur")
end
```

### Aggregation Options

| Option | Description |
|--------|-------------|
| `:mean` | Mean duration across all edges |
| `:median` | Median duration |
| `:total` | Sum of all durations |
| `:all` | Dictionary mapping each edge to its duration |

### Multiple Spells

When an edge has multiple activity spells, `tEdgeDuration` sums them:

```julia
# Edge (1,2) active in two periods
activate!(dnet, 0.0, 20.0; edge=(1, 2))
activate!(dnet, 40.0, 60.0; edge=(1, 2))

# Duration = 20 + 20 = 40
dur = tEdgeDuration(dnet; aggregate=:all)
println("Edge (1,2) total duration: ", dur[(1, 2)])  # 40.0
```

## Vertex Duration

### tVertexDuration

Compute the total duration of vertex activity:

```julia
# Mean vertex activity duration
v_dur = tVertexDuration(dnet; aggregate=:mean)
println("Mean vertex duration: $v_dur")

# All vertex durations
all_v_dur = tVertexDuration(dnet; aggregate=:all)
for (i, dur) in enumerate(all_v_dur)
    println("  Vertex $i: $dur")
end
```

The same aggregation options (`:mean`, `:median`, `:total`, `:all`) apply.

## Edge Persistence

### tEdgePersistence

Measures the proportion of edges that persist (survive) across time windows:

```julia
# Edge persistence across 20-unit windows
persistence = tEdgePersistence(dnet, 20.0)
println("Persistence (window=20): $(round(persistence, digits=3))")
```

**How it works:**

1. Divide the observation period into windows of the specified size
2. For each consecutive pair of windows $(w_i, w_{i+1})$:
   - Count edges active at the start of $w_i$
   - Count how many are also active at the start of $w_{i+1}$
3. Return the overall proportion: persisted / total

### Interpreting Persistence

| Persistence | Interpretation |
|-------------|----------------|
| ~1.0 | Very stable network (almost no edge turnover) |
| ~0.5 | Moderate turnover (half of edges change per window) |
| ~0.0 | High turnover (almost complete edge replacement) |

### Window Size Sensitivity

```julia
# Persistence at different window sizes
for w in [5.0, 10.0, 20.0, 30.0, 50.0]
    p = tEdgePersistence(dnet, w)
    println("Window=$w: persistence=$(round(p, digits=3))")
end
```

Larger windows allow more time for changes, so persistence generally decreases with window size.

## Turnover

### tTurnover

Compute edge formation and dissolution rates:

```julia
turnover = tTurnover(dnet, 20.0)

println("Formation rate: $(round(turnover.formation_rate, digits=4))")
println("Dissolution rate: $(round(turnover.dissolution_rate, digits=4))")
println("Number of formations: $(turnover.n_formations)")
println("Number of dissolutions: $(turnover.n_dissolutions)")
```

### Returned Fields

| Field | Description |
|-------|-------------|
| `formation_rate` | Formations / at-risk non-edges |
| `dissolution_rate` | Dissolutions / at-risk edges |
| `n_formations` | Total number of new edges across all windows |
| `n_dissolutions` | Total number of dissolved edges across all windows |

### Understanding Rates

The **formation rate** is the probability that a non-edge becomes an edge in the next window:

$$\text{formation rate} = \frac{\text{new edges at } t+1}{\text{non-edges at } t}$$

The **dissolution rate** is the probability that an edge dissolves in the next window:

$$\text{dissolution rate} = \frac{\text{dissolved edges at } t+1}{\text{edges at } t}$$

### Turnover Analysis

```julia
# Compare turnover at different timescales
for w in [10.0, 20.0, 30.0]
    t = tTurnover(dnet, w)
    println("Window $w:")
    println("  Formation: $(round(t.formation_rate, digits=4)) ($(t.n_formations) new)")
    println("  Dissolution: $(round(t.dissolution_rate, digits=4)) ($(t.n_dissolutions) lost)")
    println("  Net change: $(t.n_formations - t.n_dissolutions)")
end
```

## Tie Decay

### tieDecay

Estimate the rate at which ties decay (end) from observed edge spell durations:

```julia
# Exponential decay rate (MLE for exponential distribution)
decay = tieDecay(dnet; method=:exponential)
println("Decay rate: $(round(decay, digits=4))")
println("Expected duration: $(round(1.0/decay, digits=1))")

# Halflife method
halflife_rate = tieDecay(dnet; method=:halflife)
println("Halflife decay rate: $(round(halflife_rate, digits=4))")
println("Halflife: $(round(log(2)/halflife_rate, digits=1))")
```

### Methods

| Method | Formula | Interpretation |
|--------|---------|----------------|
| `:exponential` | $\hat{\lambda} = 1/\bar{d}$ | Rate parameter of exponential distribution |
| `:halflife` | $\hat{\lambda}_h = \ln(2)/\bar{d}$ | Rate such that $P(\text{survive } > \text{halflife}) = 0.5$ |

Where $\bar{d}$ is the mean edge duration.

## Contact Sequences

### as_contact_sequence

Convert a dynamic network to a flat sequence of contacts:

```julia
cs = as_contact_sequence(dnet)
println("Number of contacts: $(length(cs))")

for contact in cs
    println("  $(contact.source) -> $(contact.target): ",
            "start=$(contact.time), duration=$(contact.duration)")
end
```

### The Contact Type

```julia
struct Contact{T, Time}
    source::T       # Source vertex
    target::T       # Target vertex
    time::Time      # Start time
    duration::Time  # Duration of contact
end
```

### The ContactSequence Type

```julia
struct ContactSequence{T, Time}
    contacts::Vector{Contact{T, Time}}  # Sorted by time
    n_vertices::Int                     # Number of vertices
    directed::Bool                      # Whether directed
end
```

Contact sequences are sorted by time and support iteration:

```julia
cs = as_contact_sequence(dnet)

# Iterate over contacts
for c in cs
    println("$(c.source) -> $(c.target) at t=$(c.time)")
end
```

## Network Aggregation

### tAggregate

Collapse a dynamic network to a static network:

```julia
# Union: include any edge that was ever active
static_union = tAggregate(dnet; method=:union)
println("Union: $(ne(static_union)) edges")

# Intersection: include edges active throughout the entire observation period
static_inter = tAggregate(dnet; method=:intersection)
println("Intersection: $(ne(static_inter)) edges")

# Weighted: weight by total activation time
static_weighted = tAggregate(dnet; method=:weighted)
println("Weighted: $(ne(static_weighted)) edges")
```

### Aggregation Methods

| Method | Include Edge If | Weight |
|--------|----------------|--------|
| `:union` | Edge was ever active | 1.0 |
| `:intersection` | Edge active throughout entire observation period | 1.0 |
| `:weighted` | Edge was ever active | Total activation time |

### Weighted Aggregation

The weighted method creates a static network where edge weights represent total activity time:

```julia
static = tAggregate(dnet; method=:weighted)

# Access edge weights
for e in edges(static)
    w = get_edge_attribute(static, src(e), dst(e), :weight)
    println("Edge $(src(e))->$(dst(e)): weight = $w")
end
```

## Time Series Analysis

### Regular Interval Statistics

Use `tSnaStats` to compute SNA statistics at regular time points:

```julia
times = collect(0.0:5.0:100.0)
stats = tSnaStats(dnet, times;
    measures=[:density, :reciprocity, :n_edges, :mean_degree]
)

# Plot-ready data (one NamedTuple row per time point)
for row in stats
    println("t=$(row.time): density=$(round(row.density, digits=3)), ",
            "edges=$(Int(row.n_edges))")
end
```

### Sliding Window Statistics

Use `windowSnaStats` to compute statistics in sliding windows:

```julia
stats = windowSnaStats(dnet, 20.0;
    measures=[:density, :n_edges]
)

println("Window densities: ", round.([row.density for row in stats], digits=3))
```

### Comparing Point-in-Time vs. Window

```julia
times = collect(0.0:10.0:100.0)

# Point-in-time: what is the density at each instant?
point_stats = tSnaStats(dnet, times; measures=[:density])

# Window: what is the density aggregated over each 10-unit window?
window_stats = windowSnaStats(dnet, 10.0; measures=[:density])

println("Point densities: ", round.([row.density for row in point_stats], digits=3))
println("Window densities: ", round.([row.density for row in window_stats], digits=3))
```

Point-in-time measures capture the instantaneous state. Window measures capture aggregate activity over a period and are typically higher (more edges are active at some point during the window than at any single instant).

## Complete Dynamics Example

```julia
using NetworkDynamic
using TSNA

# Create a network with known dynamics
n = 15
dnet = DynamicNetwork(n; observation_start=0.0, observation_end=100.0)
activate_vertices!(dnet, collect(1:n), 0.0, 100.0)

# Phase 1 (0-30): Initial connections form
for i in 1:5
    activate!(dnet, 0.0, 40.0; edge=(i, i+1))
end

# Phase 2 (20-60): Network densifies
for i in 1:10
    activate!(dnet, 20.0, 60.0; edge=(i, mod1(i+2, n)))
end

# Phase 3 (50-100): Some ties persist, others dissolve
for i in 1:8
    activate!(dnet, 50.0, 100.0; edge=(i, mod1(i+3, n)))
end

# === Full Dynamics Analysis ===

println("=== Edge Duration ===")
dur = tEdgeDuration(dnet; aggregate=:mean)
println("Mean edge duration: $(round(dur, digits=1))")

println("\n=== Persistence ===")
for w in [10.0, 20.0, 30.0]
    p = tEdgePersistence(dnet, w)
    println("Window $w: $(round(p, digits=3))")
end

println("\n=== Turnover ===")
t = tTurnover(dnet, 20.0)
println("Formation rate: $(round(t.formation_rate, digits=4))")
println("Dissolution rate: $(round(t.dissolution_rate, digits=4))")

println("\n=== Tie Decay ===")
decay = tieDecay(dnet; method=:exponential)
println("Decay rate: $(round(decay, digits=4))")
println("Expected tie lifetime: $(round(1/decay, digits=1))")

println("\n=== Network Evolution ===")
times = collect(0.0:10.0:100.0)
stats = tSnaStats(dnet, times; measures=[:density, :n_edges])
for row in stats
    println("t=$(row.time): $(Int(row.n_edges)) edges, ",
            "density=$(round(row.density, digits=3))")
end
```

## Best Practices

1. **Choose appropriate window sizes**: Window sizes should match the timescale of interest in your research question
2. **Report multiple metrics**: Use persistence, turnover, and duration together for a complete picture
3. **Consider censoring**: Edge spells at the boundaries of the observation period may be censored (truncated)
4. **Use aggregation wisely**: Union is most inclusive, intersection most conservative, weighted preserves duration information
5. **Check for empty windows**: Some time windows may have no edges, producing division-by-zero or empty results
6. **Compare with expectations**: Null models or benchmark values help interpret whether observed turnover is high or low
