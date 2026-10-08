# Duration Metrics

This guide covers the duration, persistence, turnover, and aggregation metrics available in TSNA.jl for characterizing network dynamics.

All statistical queries require observed dyads unless `missing=:face` is
explicitly supplied. Density and `t_sna_stats` use active vertices by default
(`active_only=false` keeps the full universe). Reciprocity uses the dyadic
definition, which includes null dyads; pass `method=:edgewise` for the
fraction of edges reciprocated. Temporal conversions retain masks where the
target can represent them.

Date and DateTime durations are measured in seconds; rates are per second.
Numeric time axes retain their native units. Windows must advance time and
require a finite observation period.

## Lifetimes, the observation window and censoring

All duration and event metrics follow R tsna's conventions (through
networkDynamic's `as.data.frame`):

- **Spells are merged per element first.** Adjacent or overlapping spells of
  one tie are one lifetime: a tie recorded wave by wave (`activate!(dnet, t,
  t + 1; edge=...)` for each wave) forms once and dissolves once.
  (`activate!` already merges as R's `activate.edges` does; spells stored with
  `add_spell!(...; merge=false)` are merged on the fly.)
- **Elements with no spells are active throughout** (R's `active.default =
  TRUE`), so a base-network edge without spells counts as one lifetime spanning
  the window, censored at both ends.
- **Spells are clipped to the observation window** (`get_observation_period`).
  An open-ended spell contributes the time up to the window end; a spell
  starting before the window is truncated at its start.
- **Without a window, tsna's rules apply.** Lifetimes and events read the
  spells unclipped, so an open-ended spell lasts `Inf` and only an infinite
  bound is censored; `tied_duration` and `t_edge_density` clip to the range of
  the change times; `t_turnover`, `window_sna_stats` and `t_edge_persistence`
  start windows at every step of the *closed* range `[first, last]` of the
  change times. A network with no finite spell bound has no such range, and
  those functions raise an `ArgumentError` asking for
  `set_observation_period!` (tsna would use the placeholder `(0, 1)`).
- **A bound cut by the window, or stored as censored, is censored.** A
  censored onset is not a formation and a censored terminus is not a
  dissolution, unless `include_censored=true` (tsna's `include.censored`).

## Overview

While centrality and path analysis describe the network at specific times, duration metrics characterize **how the network changes over time**. They answer questions like:

- How long do ties last?
- How stable is the network structure?
- At what rate are new ties forming and old ties dissolving?
- Is the network becoming denser or sparser over time?

## Edge Duration

### t_edge_duration

With the defaults (`mode=:total`, `aggregate=:all`), `t_edge_duration(dnet)`
is R's `tsna::edgeDuration(nd)`: one entry per edge (in the order of
`edges(dnet.network)`), the total time it is active inside the window. Edges
with no activity in the window are left out.

```julia
using DynamicNetworks
using TSNA

dnet = DynamicNetwork(5; observation_start=0.0, observation_end=100.0)
activate_vertices!(dnet, collect(1:5), 0.0, 100.0)

activate!(dnet, 0.0, 50.0; edge=(1, 2))    # Duration: 50
activate!(dnet, 10.0, 30.0; edge=(2, 3))   # Duration: 20
activate!(dnet, 40.0, 80.0; edge=(3, 4))   # Duration: 40
activate!(dnet, 0.0, 100.0; edge=(4, 5))   # Duration: 100
activate!(dnet, 20.0, 40.0; edge=(1, 3))   # Duration: 20

# Per-edge durations: (1,2), (1,3), (2,3), (3,4), (4,5)
per_edge = t_edge_duration(dnet)
println("Per-edge durations: ", per_edge)  # [50.0, 20.0, 20.0, 40.0, 100.0]

# Mean edge duration
mean_dur = t_edge_duration(dnet; aggregate=:mean)
println("Mean edge duration: $mean_dur")  # (50+20+40+100+20)/5 = 46.0

# Median edge duration
med_dur = t_edge_duration(dnet; aggregate=:median)
println("Median edge duration: $med_dur")  # 40.0

# Total edge-time
total_dur = t_edge_duration(dnet; aggregate=:total)
println("Total edge-time: $total_dur")  # 230.0

# One entry per merged spell instead (tsna subject = "spells")
spell_durs = t_edge_duration(dnet; mode=:spell)
println("  Spell durations: ", spell_durs)
```

### Aggregation Options

| Option | Description |
|--------|-------------|
| `:mean` | Mean duration across all edges |
| `:median` | Median duration |
| `:total` | Sum of all durations |
| `:all` (default) | The vector of durations (one per edge, or per spell with `mode=:spell`) |

### Multiple Spells

By default (`mode=:total`) the spells of each edge are summed; `mode=:spell`
keeps one entry per spell. Adjacent or overlapping spells are one spell:

```julia
d2 = DynamicNetwork(3; observation_start=0.0, observation_end=100.0)
activate!(d2, 0.0, 20.0; edge=(1, 2))
activate!(d2, 40.0, 60.0; edge=(1, 2))     # a second, separate spell
activate!(d2, 0.0, 30.0; edge=(2, 3))
activate!(d2, 30.0, 50.0; edge=(2, 3))     # adjacent: merged to [0, 50)

t_edge_duration(d2)                    # [40.0, 50.0]  per edge
t_edge_duration(d2; mode=:spell)       # [20.0, 20.0, 50.0]  per spell
```

### Censoring

A spell still active at the window end (or that began before the window
start) is right- (left-) censored: its clipped duration is shorter than the tie's
true lifetime. The plain mean, median or total of `t_edge_duration` therefore
underestimates the typical tie lifetime when many ties are censored.
TSNA.jl offers no censoring-aware estimator (such as Kaplan–Meier); count the
censored spells, or restrict to spells that both start and end inside the
window, before reading a mean as a lifetime.

```julia
d3 = DynamicNetwork(2; observation_start=0.0, observation_end=10.0)
activate!(d3, 2.0, Inf; edge=(1, 2))   # still active when observation ends
t_edge_duration(d3)                    # [8.0]: clipped at 10, never Inf
```

## Vertex Duration

### t_vertex_duration

Summarize vertex activity durations:

Vertex durations follow the same conventions: one entry per vertex active in
the window (in ID order), its spells merged and clipped to the window.
R's `tsna::vertexDuration` differs here: it replaces only *infinite* bounds by
the window and neither drops nor truncates finite spells outside it, so in R a
vertex can be active longer than the window lasts. TSNA.jl clips, as R does for
edges.

```julia
# Mean vertex activity duration
v_dur = t_vertex_duration(dnet; aggregate=:mean)
println("Mean vertex duration: $v_dur")

# Per-vertex totals (the default)
all_v_dur = t_vertex_duration(dnet)
println("Per-vertex total durations: ", all_v_dur)
```

The same `mode` and `aggregate` options apply.

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

`t_edge_formation` counts the **onsets** and `t_edge_dissolution` the
**termini** of merged, clipped spells falling inside `[onset, terminus)`. A
left-censored onset (a tie already present when observation starts) is not a
formation, and a right-censored terminus (a tie still present when observation
ends) is not a dissolution: both are artefacts of the observation window. Pass
`include_censored=true` to count them anyway, at their clipped position.

```julia
d4 = DynamicNetwork(3; observation_start=0.0, observation_end=10.0)
add_spell!(d4, Spell(0.0, 4.0; onset_censored=true); edge=(1, 2))
activate!(d4, 2.0, Inf; edge=(2, 3))
t_edge_formation(d4, 0.0, 10.0)                          # 1 (at t = 2)
t_edge_formation(d4, 0.0, 10.0; include_censored=true)   # 2
t_edge_dissolution(d4, 0.0, 11.0)                        # 1 (at t = 4)
```

R's `tsna::tEdgeFormation(nd, start, end, time.interval = 1)` returns a series:
the number of onsets equal to each `t` in `seq(start, end, 1)`. On
integer-valued spells that is `[t_edge_formation(dnet, t, t + 1) for t in
start:end]`, and likewise for dissolutions.

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

Formations and dissolutions are counted as **lifetime events** (onsets and
termini of merged spells falling inside the window, censored bounds excluded
unless `include_censored=true`), and each rate is the event count divided by
the window length:

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
to `at` (0 for currently active ties, including base-network edges with no
spells, which are active throughout).

## Contact Sequences

### as_contact_sequence

Convert a dynamic network to a flat sequence of contacts, one per contiguous
activity of an edge (its spells merged, as R's `activate.edges` stores them).
A base-network edge with no spells has no finite onset and is left out, and
so is a spell that is open to the left (onset `-Inf`, or the axis minimum on
integer and calendar axes); a contact open to the right lasts the axis's
unbounded duration (`Inf`, or `typemax`). The contacts are not clipped to an
observation window. With `report=true` the `ConversionReport` names each
(`:default_active_edges`, `:unbounded_onsets`, `:unbounded_termini`):

```julia
cs = as_contact_sequence(dnet)
println("Number of contacts: $(length(cs))")

for contact in cs
    println("  $(contact.source) -> $(contact.target): ",
            "start=$(contact.time), duration=$(contact.duration)")
end
```

### The Contact Type

```text
struct Contact{T, Time}
    source::T   # Source vertex
    target::T   # Target vertex
    time::Time  # Start time
    duration    # Duration of contact (terminus - onset)
end
```

### The ContactSequence Type

```text
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
using NetworkCore   # for ne on the aggregated static networks

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
| `:union` | Edge and both endpoints active at some point in the window | 1.0 |
| `:intersection` | Edge and both endpoints active throughout the window | 1.0 |
| `:weighted` | As `:union` | Total activation time in the window (spells merged and clipped) |

`t_aggregate` is `DynamicNetworks.network_collapse` (R's `network.collapse`)
over the window, except that every vertex is kept so IDs stay stable.

### Weighted Aggregation

The weighted method creates a static network where edge weights represent total activity time:

```julia
using NetworkCore   # edges, src, dst, get_edge_attribute

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
using DynamicNetworks
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
3. **Consider censoring**: Spells cut by the observation window are censored; their clipped durations understate tie lifetimes, and their cut bounds are not formation or dissolution events
4. **Use aggregation wisely**: Union is most inclusive, intersection most conservative, weighted preserves duration information
5. **Check for empty windows**: Some time windows may have no edges, producing division-by-zero or empty results
6. **Compare with expectations**: Null models or benchmark values help interpret whether observed turnover is high or low
