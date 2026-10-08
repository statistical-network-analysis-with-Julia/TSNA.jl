# Coming from R tsna

TSNA.jl ports the descriptive core of R's `tsna` 0.3.6, built on
DynamicNetworks.jl (R `networkDynamic`) and SNA.jl (R `sna`). The functions use
snake_case names. Earlier development versions also had camelCase names,
described as "R tsna-style"; most were not tsna functions, and the five that
share a tsna name (`tPath`, `tDegree`, `tEdgeFormation`, `tEdgeDissolution`,
`tSnaStats`) take other arguments and return something else in tsna. They
were never released and are gone; [Renamed functions](@ref) maps each to its
snake_case name.

## Conventions shared with R

These conventions are those of networkDynamic and tsna:

- Spells are half-open `[onset, terminus)`, and a point spell `[t, t)` is an
  instantaneous contact.
- An element with no spells is active throughout (`active.default = TRUE`).
- Spells of one element are merged on activation, so adjacent or overlapping
  spells form one lifetime.
- Durations and event counts are taken over the observation window
  (`net.obs.period`; `get_observation_period`). Spells are truncated to the
  window and flagged censored where they are cut, as in networkDynamic's
  `as.data.frame`.
- **Without a window, each function follows tsna's own rule** (tsna's
  functions do not agree on one):

  | What | tsna, no `net.obs.period` | TSNA.jl, no window |
  |:--|:--|:--|
  | `edgeDuration`, `vertexDuration` | spells read over `(-Inf, Inf)`: an open spell lasts `Inf` | the same: [`t_edge_duration`](@ref), [`t_vertex_duration`](@ref) return `Inf` |
  | `tEdgeFormation`/`tEdgeDissolution` events | read over `(-Inf, Inf)`: only an infinite bound is censored | the same: [`t_edge_formation`](@ref), [`t_edge_dissolution`](@ref) |
  | default series of `tEdgeFormation`, `tEdgeDissolution`, `tSnaStats` | `seq(start, end)` over the **closed** range of `get.change.times` | [`t_turnover`](@ref), [`window_sna_stats`](@ref), [`t_edge_persistence`](@ref) start windows at every step of `[first, last]`, the last change time included |
  | `tiedDuration`, `tEdgeDensity` | clipped to `get_bounds`: the range of `get.change.times` | the same: [`tied_duration`](@ref), [`t_edge_density`](@ref) |
  | `tPath` | `end = Inf` | `end_time` defaults to `Inf` (with a window: the window end) |
  | any of the above, no finite change time | the placeholder range `(0, 1)` | `ArgumentError` asking for `set_observation_period!` |

One thing differs:

- **Snapshot measures are scattered back by vertex ID.** Per-vertex snapshot
  measures return a vector indexed by vertex ID, with `NaN` for actors absent at
  that time.

## Function table

Every exported `tsna` function is listed. "≈" marks a counterpart whose
signature or result differs as described.

| R tsna | TSNA.jl | Differences |
|:--|:--|:--|
| `tPath(nd, v, direction = "fwd", type = "earliest.arrive", start, end)` | [`earliest_arrival`](@ref), [`temporal_path`](@ref), [`temporal_distance`](@ref), [`forward_reachable_set`](@ref) | ≈ R returns a `tPath` object (`tdist`, `previous`, `gsteps`) for every vertex. Julia returns the arrival and parent dictionaries; `temporal_path` reconstructs one path as a [`TemporalPath`](@ref). Arrival times agree exactly. `start_time` is required (tsna defaults to the first change time). `end_time` defaults to `Inf`, as tsna's `end`, except that a set observation window ends the search at its end (tsna ignores the window). `graph.step.time` (a traversal cost) is not implemented: traversal is instantaneous, as with tsna's default `0`. |
| `tPath(..., direction = "bkwd")` | [`backward_reachable_set`](@ref) | ≈ Returns the vertex set only. tsna's backward search does not chain same-instant point spells, although its forward search does. `backward_reachable_set` is the exact dual of the forward search, so 5 of 298 fuzzed backward queries differ from tsna. |
| `tPath(..., type = "latest.depart")` | — | Not implemented. |
| `forward.reachable(nd, v, start, end, per.step.depth)` | [`forward_reachable_set`](@ref) | ≈ `per.step.depth` is not implemented. |
| `tReach(nd, direction, sample, start, end)` | `sum(reachability_matrix(dnet, start; end_time); dims = 2)` | ≈ [`reachability_matrix`](@ref) gives every forward reachable set at once. Seed sampling (`sample`), the backward direction and `graph.step.time` are not implemented. |
| `is.tPath`, `as.network.tPath`, `plot.tPath`, `plotPaths` | — | Not implemented. Plotting is out of scope (see NDTV.jl). |
| `tDegree(nd, start, end, time.interval, cmode)` | [`t_degree`](@ref) at each time | ≈ R returns a time × vertex series over `seq(start, end, time.interval)`. `t_degree(dnet, t; mode)` gives one time point: `mode = :total/:in/:out` is `cmode = "freeman"/"indegree"/"outdegree"`. |
| `tSnaStats(nd, snafun, start, end, time.interval, aggregate.dur, rule)` | [`t_degree`](@ref), [`t_betweenness`](@ref), [`t_closeness`](@ref), [`t_eigenvector`](@ref), [`t_pagerank`](@ref), [`t_density`](@ref), [`t_reciprocity`](@ref), [`t_transitivity`](@ref); [`t_sna_stats`](@ref), [`window_sna_stats`](@ref) | ≈ R applies any of 23 sna functions to `network.collapse(nd, at = t)` (the active vertices) at each grid time. TSNA.jl has one function per measure, also on the active network. When the vertex set changes, tsna `rbind`s per-vertex vectors of different lengths and so misaligns their columns; TSNA.jl scatters each value back to its vertex ID, with `NaN` for absent actors. `t_sna_stats` returns network-level measures (density, reciprocity, transitivity, mean degree, edge count) at the given times. The static measures are SNA.jl's ports under sna's names (`SNA.degreecent` for `degree`, `SNA.betweenness`, `SNA.closeness`, `SNA.evcent`, `SNA.gden`, `SNA.grecip`, `SNA.gtrans`); `t_pagerank` uses Graphs.jl's `pagerank`, which is not an sna measure. `aggregate.dur`/`rule` (interval snapshots) and arbitrary `snafun` are not implemented. |
| `tErgmStats(nd, formula, start, end, ...)` | — | Not implemented. Use `ERGM.summary_stats` on `network_extract` snapshots. |
| `edgeDuration(nd, mode = "duration", subject = "edges", active.default)` | [`t_edge_duration`](@ref) | Same result with the defaults: one total per edge, spells merged and truncated to the window. `subject = "spells"` is `mode = :spell`. ≈ `subject = "dyads"` and `mode = "counts"` are not implemented. The order follows `edges(dnet.network)`, not R's edge ids. `aggregate = :mean/:median/:total` summarises the vector. |
| `vertexDuration(nd, mode, subject, v, active.default)` | [`t_vertex_duration`](@ref) | ≈ R replaces only infinite bounds by the window. It neither drops nor truncates finite spells outside the window, so an R vertex can be "active" longer than the window lasts. TSNA.jl truncates vertex spells as `edgeDuration` truncates edge spells. The two agree on vertices whose spells lie inside the window. |
| `tiedDuration(nd, mode, active.default, neighborhood)` | [`tied_duration`](@ref)`(dnet; mode, neighborhood, active_default)` | Same values (golden fixture). `mode` is `:duration`/`:counts`. |
| `tEdgeFormation(nd, start, end, time.interval, result.type, include.censored)` | [`t_edge_formation`](@ref), [`t_turnover`](@ref) | ≈ R returns a series: the number of uncensored onsets equal to each `t` in `seq(start, end, time.interval)`. `t_edge_formation(dnet, a, b; include_censored)` counts onsets in `[a, b)`, so on integer-valued spells R's series is `[t_edge_formation(dnet, t, t + 1) for t in start:end]`. On a network without a window, R's default series (the closed range of the change times) is the `n_formations` column of `t_turnover(dnet, 1)`. `result.type = "fraction"` is not implemented. |
| `tEdgeDissolution(...)` | [`t_edge_dissolution`](@ref) | As for `tEdgeFormation`, with termini. |
| `tEdgeDensity(nd, mode, agg.unit, active.default)` | [`t_edge_density`](@ref)`(dnet; mode, agg_unit, active_default)` | Same values (golden fixture), including `agg.unit = "dyad"` on two-mode and looped networks (`network.dyadcount`). For `mode = "event"` tsna 0.3.6 divides by `n_edges · end − start`; TSNA.jl uses the intended `n_edges · (end − start)`, so the two differ when the range does not start at 0. |
| `pShiftCount(nd, start, end, output)` | — | Not implemented. The relational-event participation shifts are in REM.jl/Revel.jl. |
| `timeProjectedNetwork(nd, start, end, ...)` | — | Not implemented. |

These TSNA.jl functions have no tsna counterpart:

- [`t_edge_persistence`](@ref): the pooled proportion of edges surviving between window starts.
- [`t_turnover`](@ref): formation and dissolution counts and rates per window.
- [`tie_decay`](@ref): decayed tie weights.
- [`t_aggregate`](@ref): networkDynamic's `network.collapse` with an aggregation rule.
- [`as_contact_sequence`](@ref).
- [`temporal_distance_matrix`](@ref) and [`earliest_arrival_all`](@ref).

## Validation against R

The golden fixture `test/fixtures/tsna_semantics.toml` is generated by
`test/fixtures/r/tsna_semantics.R` (R 4.6.1, tsna 0.3.6, sna 2.8,
networkDynamic 0.12.0). Its networks have adjacent, overlapping, point,
open-ended and pre-window spells, and vertices and edges that are never
activated. The test suite compares TSNA.jl's plain vectors with R's labelled
values, ordered by the documented rule (edges by `(i, j)`, spells by edge and
onset, vertices by id), so the comparison does not go through TSNA.jl's own
internals. It has four groups:

- **60 one-mode networks observed over `[0, 12]`.** Exactly: `edgeDuration`
  (subject edges and spells); `tEdgeFormation`/`tEdgeDissolution` with and
  without `include.censored`; `tiedDuration` and `tEdgeDensity`; the clipped
  vertex durations TSNA.jl documents (computed in R from networkDynamic's stored
  spells), which equal `vertexDuration` on vertices whose spells lie in the
  window. To 1e-10: sna's `degree`, `betweenness`, `closeness` and `gden` on
  `network.extract(nd, at = t)` at four instants, matched by vertex id, with
  `NaN` for absent actors.
- **20 two-mode (undirected and directed) and looped networks** over
  `[0, 12]`: the same durations, events, `tiedDuration`, and `tEdgeDensity`
  with `network.dyadcount`'s dyads.
- **44 networks without an observation window** (one-mode, two-mode, looped,
  and four fixed edge cases: a tie that forms at the last change time, a
  dissolution at the first change time, no finite spell bound at all, and an
  undirected open tie): `edgeDuration` with its `Inf` values,
  `vertexDuration`, the default `tEdgeFormation`/`tEdgeDissolution` series,
  `tiedDuration`, `tEdgeDensity` and `tPath(direction = "fwd")` from every
  vertex (251 paths, exact). Where tsna uses the placeholder `(0, 1)`, TSNA.jl
  must raise an `ArgumentError`. tsna's `tPath` warns that it is unreliable
  when spell times are negative, and indeed its answers change when every time
  is shifted by a constant; the fixture therefore runs it on a copy shifted to
  non-negative times (elapsed times do not depend on the shift), while
  TSNA.jl runs on the original times.
- **10 samplk-style panels** (three waves on 18 vertices) built with
  `networkDynamic(network.list = ...)` and `DynamicNetwork(networks)`:
  `edgeDuration`, the event series, `tiedDuration`, `tEdgeDensity`,
  `tSnaStats(gden, grecip)` and `tPath` from vertex 1.

tsna's backward search (`tPath(direction = "bkwd")`) does not chain
same-instant point spells although its forward search does, so a randomised
comparison finds a few backward queries (5 of 298) where tsna stops short of
`backward_reachable_set`, which is the exact dual of the forward search.

To regenerate the fixture:

```text
Rscript test/fixtures/r/tsna_semantics.R > test/fixtures/tsna_semantics.toml
```

## Renamed functions

Earlier development versions exported camelCase names next to the snake_case
ones. None was released, and none is defined now. The tsna column gives the R
function of the same name, where there is one; it is a different function.

| Earlier name | TSNA.jl name | tsna function of that name |
|:--|:--|:--|
| `tPath` | [`TemporalPath`](@ref) (a type) | `tPath`: the path search, here [`earliest_arrival`](@ref) and friends |
| `temporalDistance` | [`temporal_distance`](@ref) | — |
| `earliestArrival` | [`earliest_arrival`](@ref) | — |
| `earliestArrivalAll` | [`earliest_arrival_all`](@ref) | — |
| `forwardReachableSet` | [`forward_reachable_set`](@ref) | — (tsna has `forward.reachable`) |
| `backwardReachableSet` | [`backward_reachable_set`](@ref) | — |
| `temporalPath`, `shortestTemporalPath`, `shortest_temporal_path` | [`temporal_path`](@ref) (the earliest-arrival path, not a shortest one) | — |
| `tDegree` | [`t_degree`](@ref) | `tDegree`: a time × vertex series |
| `tBetweenness`, `tCloseness`, `tEigenvector`, `tPagerank` | [`t_betweenness`](@ref), [`t_closeness`](@ref), [`t_eigenvector`](@ref), [`t_pagerank`](@ref) | — |
| `tDensity`, `tReciprocity`, `tTransitivity` | [`t_density`](@ref), [`t_reciprocity`](@ref), [`t_transitivity`](@ref) | — |
| `tEdgeDuration`, `tVertexDuration` | [`t_edge_duration`](@ref), [`t_vertex_duration`](@ref) | — (tsna has `edgeDuration`, `vertexDuration`) |
| `tEdgeFormation`, `tEdgeDissolution` | [`t_edge_formation`](@ref), [`t_edge_dissolution`](@ref) | `tEdgeFormation`, `tEdgeDissolution`: event series |
| `tEdgePersistence`, `tTurnover`, `tieDecay` | [`t_edge_persistence`](@ref), [`t_turnover`](@ref), [`tie_decay`](@ref) | — |
| `tSnaStats`, `windowSnaStats` | [`t_sna_stats`](@ref), [`window_sna_stats`](@ref) | `tSnaStats`: one sna function over a time grid |
| `tAggregate` | [`t_aggregate`](@ref) | — |
