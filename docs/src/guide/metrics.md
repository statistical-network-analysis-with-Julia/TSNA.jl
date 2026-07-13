# Duration Metrics

This guide covers the duration, persistence, turnover, and aggregation metrics available in TSNA.jl for characterizing network dynamics.

## Overview

While centrality and path analysis describe the network at specific times, duration metrics characterize **how the network changes over time**. They answer questions like:

- How long do ties last?
- How stable is the network structure?
- At what rate are new ties forming and old ties dissolving?
- Is the network becoming denser or sparser over time?

## Edge Duration

### t_edge_duration

Summarize edge activity durations:

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
mean_dur = t_edge_duration(dnet; aggregate=:mean)
println("Mean edge duration: $mean_dur")  # (50+20+40+100+20)/5 = 46.0

# Median edge duration
med_dur = t_edge_duration(dnet; aggregate=:median)
println("Median edge duration: $med_dur")  # 40.0

# Total edge-time
total_dur = t_edge_duration(dnet; aggregate=:total)
println("Total edge-time: $total_dur")  # 230.0

# Raw per-spell durations (a Vector{Float64}; use mode=:total to sum
# spells per edge first)
all_durs = t_edge_duration(dnet; aggregate=:all)
println("  Spell durations: ", all_durs)
```

### Aggregation Options

| Option | Description |
|--------|-------------|
| `:mean` | Mean duration across all edges |
| `:median` | Median duration |
| `:total` | Sum of all durations |
| `:all` | Raw vector of durations (one per spell, or per edge with `mode=:total`) |

### Multiple Spells

By default (`mode=:spell`, matching `tsna::edgeDuration`) each activity
spell contributes its own duration; pass `mode=:total` to sum the spells
of each edge first:

```julia
# Edge (1,2) active in two periods
activate!(dnet, 0.0, 20.0; edge=(1, 2))
activate!(dnet, 40.0, 60.0; edge=(1, 2))

# mode=:spell (default) keeps the two spells separate;
# mode=:total sums them per edge: 20 + 20 = 40
durs = t_edge_duration(dnet; mode=:total, aggregate=:all)
println("Per-edge total durations: ", durs)
```

## Vertex Duration

### t_vertex_duration

Summarize vertex activity durations:

```julia
# Mean vertex activity duration
v_dur = t_vertex_duration(dnet; aggregate=:mean)
println("Mean vertex duration: $v_dur")

# Raw durations (per spell; use mode=:total for per-vertex totals)
all_v_dur = t_vertex_duration(dnet; mode=:total, aggregate=:all)
println("Per-vertex total durations: ", all_v_dur)
```

The same aggregation options (`:mean`, `:median`, `:total`, `:all`) apply.

## Edge Persistence

### t_edge_persistence

Measures the proportion of edges that persist (survive) across time windows:

```julia
# Edge persistence across 20-unit windows
persistence = t_edge_persistence(dnet, 20.0)
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
| `NaN` | Fewer than two windows, or no active edges to track |

### Window Size Sensitivity

```julia
# Persistence at different window sizes
for w in [5.0, 10.0, 20.0, 30.0, 50.0]
    p = t_edge_persistence(dnet, w)
    println("Window=$w: persistence=$(round(p, digits=3))")
end
```

Larger windows allow more time for changes, so persistence generally decreases with window size.

## Formation and Dissolution Events

### t_edge_formation and t_edge_dissolution

Count edge-spell events in an arbitrary window `[onset, terminus)`:

```julia
# How many spells started in [0, 50)?
nf = t_edge_formation(dnet, 0.0, 50.0)
println("Formations in [0, 50): $nf")

# How many spells ended in [0, 50)?
nd = t_edge_dissolution(dnet, 0.0, 50.0)
println("Dissolutions in [0, 50): $nd")
```

`t_edge_formation` counts spell **onsets** and `t_edge_dissolution`
counts spell **termini** falling inside the window. Right-censored spells
(those still active when observation ended) are excluded from the
dissolution count — a censored terminus is an artifact of the observation
window, not a real dissolution event.

## Turnover

### t_turnover

Compute edge formation and dissolution rates:

```julia
# One NamedTuple per window of length 20
for w in t_turnover(dnet, 20.0)
    println("[$(w.window_start), $(w.window_end)): ",
            "formation rate = $(round(w.formation_rate, digits=4)), ",
            "dissolution rate = $(round(w.dissolution_rate, digits=4)), ",
            "+$(w.n_formations)/-$(w.n_dissolutions) edges")
end
```

### Returned Fields (per window)

| Field | Description |
|-------|-------------|
| `window_start`, `window_end` | The window boundaries |
| `n_formations` | Spell onset events in the window |
| `n_dissolutions` | Spell terminus events in the window |
| `formation_rate` | Formations per unit time |
| `dissolution_rate` | Dissolutions per unit time |

### Understanding Rates

Formations and dissolutions are counted as **spell events** (onsets and
termini falling inside the window, right-censored termini excluded), and
each rate is the event count divided by the window length:

$$\text{formation rate} = \frac{\text{spell onsets in window}}{\text{window length}}
\qquad
\text{dissolution rate} = \frac{\text{spell termini in window}}{\text{window length}}$$

### Turnover Analysis

```julia
# Compare turnover at different timescales (totals over all windows)
for w in [10.0, 20.0, 30.0]
    windows = t_turnover(dnet, w)
    total_form = sum(x.n_formations for x in windows)
    total_diss = sum(x.n_dissolutions for x in windows)
    println("Window $w: +$total_form/-$total_diss, ",
            "net change = $(total_form - total_diss)")
end
```

## Tie Decay

### tie_decay

Per-edge tie weights decayed by the time since each edge was last active
(1.0 for currently active ties):

```julia
# Exponential decay: exp(-rate·Δ), Δ = time since last activity
weights = tie_decay(dnet; method=:exponential, rate=0.1)
println("Edge decay weights: ", weights)

# Linear decay: max(0, 1 - rate·Δ)
weights_lin = tie_decay(dnet; method=:linear, rate=0.05)

# Evaluate at a specific time instead of the observation end
weights_25 = tie_decay(dnet; at=25.0)
```

### Methods

| Method | Formula | Interpretation |
|--------|---------|----------------|
| `:exponential` | $e^{-\text{rate}\,\Delta}$ | Smooth exponential forgetting |
| `:linear` | $\max(0,\ 1 - \text{rate}\,\Delta)$ | Weight hits zero after $1/\text{rate}$ time units |

Where $\Delta$ is the time from the end of the edge's most recent spell
to `at` (0 for currently active ties).

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
    source::T   # Source vertex
    target::T   # Target vertex
    time::Time  # Start time
    duration    # Duration of contact (terminus - onset)
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

### t_aggregate

Collapse a dynamic network to a static network:

```julia
using Networks   # for ne on the aggregated static networks

# Union: include any edge that was ever active
static_union = t_aggregate(dnet; method=:union)
println("Union: $(ne(static_union)) edges")

# Intersection: include edges active throughout the entire observation period
static_inter = t_aggregate(dnet; method=:intersection)
println("Intersection: $(ne(static_inter)) edges")

# Weighted: weight by total activation time
static_weighted = t_aggregate(dnet; method=:weighted)
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
using Graphs: src, dst

static = t_aggregate(dnet; method=:weighted)

# Access edge weights
for e in edges(static)
    w = get_edge_attribute(static, :weight, src(e), dst(e))
    println("Edge $(src(e))->$(dst(e)): weight = $w")
end
```

## Time Series Analysis

### Regular Interval Statistics

Use `t_sna_stats` to compute SNA statistics at regular time points:

```julia
times = collect(0.0:5.0:100.0)
stats = t_sna_stats(dnet, times;
    measures=[:density, :reciprocity, :n_edges, :mean_degree]
)

# Plot-ready data (one NamedTuple row per time point)
for row in stats
    println("t=$(row.time): density=$(round(row.density, digits=3)), ",
            "edges=$(Int(row.n_edges))")
end
```

### Window-Based Statistics

Use `window_sna_stats` to sample statistics at the start of consecutive
windows spanning the observation period:

```julia
stats = window_sna_stats(dnet, 20.0;
    measures=[:density, :n_edges]
)

println("Window densities: ", round.([row.density for row in stats], digits=3))
```

### Comparing Custom Times vs. Window Grid

```julia
times = collect(0.0:10.0:100.0)

# Custom time points (includes the endpoint t=100)
point_stats = t_sna_stats(dnet, times; measures=[:density])

# Regular grid: one snapshot at the start of each 10-unit window
window_stats = window_sna_stats(dnet, 10.0; measures=[:density])

println("Point densities: ", round.([row.density for row in point_stats], digits=3))
println("Window densities: ", round.([row.density for row in window_stats], digits=3))
```

Both compute the same point-in-time snapshot statistics;
`window_sna_stats` just builds the time grid for you (one snapshot at the
start of each window, so the observation end itself is not sampled).

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
dur = t_edge_duration(dnet; aggregate=:mean)
println("Mean edge duration: $(round(dur, digits=1))")

println("\n=== Persistence ===")
for w in [10.0, 20.0, 30.0]
    p = t_edge_persistence(dnet, w)
    println("Window $w: $(round(p, digits=3))")
end

println("\n=== Turnover ===")
for w in t_turnover(dnet, 20.0)
    println("[$(w.window_start), $(w.window_end)): ",
            "+$(w.n_formations)/-$(w.n_dissolutions)")
end

println("\n=== Tie Decay ===")
weights = tie_decay(dnet; method=:exponential)
println("Edge decay weights: ", weights)

println("\n=== Network Evolution ===")
times = collect(0.0:10.0:100.0)
stats = t_sna_stats(dnet, times; measures=[:density, :n_edges])
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
