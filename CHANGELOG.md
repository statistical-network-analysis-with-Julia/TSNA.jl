# Changelog

All notable changes to TSNA.jl are documented in this file. The format is
based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the
package adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.2.0] - Unreleased

Release driven by the 2026-07 expert-panel review: snake_case primary names
(with the R camelCase spellings kept as aliases), a heap-based
earliest-arrival core over a memoized contact index, and R-`tsna`-aligned
return types for the temporal metrics (breaking where the old returns were
ad hoc).

### Breaking

Although every camelCase name survives as an alias (see Changed), several
functions changed *return types or semantics*, under both spellings:

- **`temporal_distance`/`temporalDistance`** returns the elapsed time
  `arrival − start` and **`nothing` when unreachable** (was a raw
  arrival-style `Time` with `Inf`/`typemax` sentinels). *Migration:* test
  `isnothing(d)` instead of `d == Inf`.
- **`forward_reachable_set`/`backward_reachable_set`** return a sorted
  `Vector` (was a `Set`); `backward_reachable_set` is now
  `(dnet, target, end_time; start_time=)`. *Migration:* treat the result as
  a vector.
- **`t_edge_duration`/`t_vertex_duration`** return a scalar aggregate
  (`aggregate=:mean|:median|:total`, or `:all` for the raw vector; new
  `mode=:spell|:total`) and `NaN` on empty input (was per-element
  collections and `0.0`). *Migration:* pass `aggregate=:all` for raw
  durations; expect `NaN` not `0.0`.
- **`t_turnover`** returns a `Vector` of per-window NamedTuples with counts
  and rates (was a single NamedTuple); `t_edge_persistence` returns `NaN`
  (was `1.0`) with fewer than two windows or no edges. *Migration:* iterate
  the vector.
- **`tie_decay`/`tieDecay` redefined:** returns a `Dict` of decayed tie
  weights per edge (`method=:exponential|:linear`, `rate`, `at`) instead of
  a scalar decay-rate estimate. *Migration:* this is a different quantity —
  review every caller.
- **`t_sna_stats`/`window_sna_stats`** return a `Vector{NamedTuple}` (one
  row per window) with a `measures=` keyword (was a `Dict{Symbol,Vector}`
  with `stats=`). *Migration:* rename the keyword and consume rows.
- **`shortest_temporal_path` is an alias of `temporal_path`**
  (earliest-arrival path, not fewest hops). *Migration:* if you needed true
  fewest-hops paths, compute them separately.
- **Point-in-time centralities (`t_degree`, `t_betweenness`, ...)** delegate
  to SNA.jl on a `retain_all_vertices=true` snapshot: results are always
  length `nv(dnet)` indexed by stable vertex IDs (inactive vertices score
  0), and `t_betweenness` defaults to `normalized=false` (SNA.jl 0.2
  convention). *Migration:* expect full-length vectors; pass
  `normalized=true` for the old scaling.
- **Minimum Julia raised to 1.12**; package UUID regenerated. *Migration:*
  upgrade Julia and re-resolve environments pinning the old UUID.

### Added

- `earliest_arrival(dnet, source, start_time; end_time, target)` — the new
  heap-based (Dijkstra-style label-setting) earliest-arrival core returning
  `(arrival, parent)` dictionaries; all path/reachability functions
  delegate to it.
- New temporal metrics `t_edge_formation` / `t_edge_dissolution`
  (`tEdgeFormation`/`tEdgeDissolution`) counting onset/terminus events per
  window (dissolution excludes right-censored spells).
- `TemporalPath` result type with `path_duration` returning a proper time
  difference (works for `Date`/`DateTime` axes).
- BenchmarkTools suite (`benchmark/`) exercising `earliest_arrival`.

### Changed

- Public API renamed to snake_case primaries with every R-style camelCase
  name kept as an exported `const` alias — both spellings work:
  `earliest_arrival`↤`earliestArrival`, `temporal_distance`,
  `forward_reachable_set`, `backward_reachable_set`, `temporal_path`,
  `shortest_temporal_path`, `t_degree`/`t_betweenness`/`t_closeness`/
  `t_eigenvector`/`t_pagerank`, `t_density`/`t_reciprocity`/
  `t_transitivity`, `t_edge_duration`/`t_vertex_duration`/
  `t_edge_persistence`/`t_turnover`, `tie_decay`, `t_sna_stats`/
  `window_sna_stats`/`t_aggregate`; the path type is `TemporalPath` with
  `tPath` as alias.

### Performance

- Memoized per-network contact index (invalidated via NetworkDynamic's
  `mutation_count`, lock-guarded), so repeated queries stop rebuilding
  per-vertex contact lists — `backward_reachable_set` runs one search per
  vertex and benefits most.
- `earliest_arrival` uses a binary min-heap label-setting search instead of
  a linear scan over sorted edge events.
- `t_turnover` sorts event times once and answers each window with two
  binary searches.

## [0.1.0] - 2026-02-09

Initial release: temporal paths, reachability, temporal centrality, and
windowed metrics over NetworkDynamic.jl networks (camelCase API).
