# Changelog

All notable changes to TSNA.jl are documented in this file. The format is
based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the
package adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.2.0] - Unreleased

First public release: temporal social network analysis on DynamicNetworks.jl
networks, following R `tsna`'s conventions where tsna has the function.

**Dependencies renamed:** the foundation package is now `NetworkCore` (developed as `Networks`) and the dynamic-network package is now `DynamicNetworks` (developed as `NetworkDynamic`); update `using` lines accordingly. Types and functions keep their names.

### Breaking

- **Per-vertex snapshot measures use the active network.**
  `t_degree`, `t_betweenness`, `t_closeness`, `t_eigenvector` and
  `t_pagerank` are computed on the actors active at the instant, as tsna's
  `tSnaStats` does. They return a length-`nv(dnet)` vector by vertex ID, with
  `NaN` for absent actors. `active_only=false` computes on the whole vertex
  universe instead.
- **`t_transitivity` uses the active actors**, like `t_density`,
  `t_reciprocity` and `t_sna_stats` (`active_only=false` to opt out). Its
  keywords are explicit: `measure` (`SNA.gtrans`'s `:weak`, `:strong`, …) and
  `type=:global`/`:average` (`cmode` for the average local coefficient);
  `type=:local` is refused.
- The snapshot measures call SNA.jl's sna-named functions (`degreecent`,
  `betweenness`, `closeness`, `evcent`, `gden`, `grecip`, `gtrans`) and pass
  their keywords (`cmode`, `rescale`, …); `t_pagerank` is Graphs.jl's PageRank.
- **`t_edge_duration` defaults to `tsna::edgeDuration`**: `mode=:total`,
  `aggregate=:all`, so it returns one total active time per edge. Use
  `mode=:spell` for one entry per spell and `aggregate=:mean/:median/:total`
  for a summary.
- **One lifetime per contiguous activity.** Each element's spells are merged
  before durations, formations, dissolutions, turnover, weighted aggregation
  and contact sequences are computed. A tie recorded wave by wave forms once.
- **Durations are clipped to the observation window.** An open-ended spell
  contributes the time up to the window end instead of `Inf`.
  `t_vertex_duration` clips too, unlike `tsna::vertexDuration`, which truncates
  only infinite bounds.
- **Without an observation window, each statistic follows tsna's rule for
  that case**, instead of clipping everything to a window derived from the
  first and last change times (which dropped a tie forming at the last change
  time, a dissolution at the first, and path reach beyond the last, and fell
  back to `(0, 1)` when there was no finite spell bound):
  - `t_edge_duration`, `t_vertex_duration`, `t_edge_formation` and
    `t_edge_dissolution` read the spells unclipped, as tsna's
    `edgeDuration`/`vertexDuration`/`tEdgeFormation` do: an open spell lasts
    `Inf`, and only an infinite bound is censored;
  - `tied_duration` and `t_edge_density` clip to the range of the change
    times, as tsna's `tiedDuration`/`tEdgeDensity`;
  - `t_turnover`, `window_sna_stats` and `t_edge_persistence` start windows at
    every step of the closed range `[first, last]` of the change times (tsna's
    `seq(start, end)`), so events at the last change time count;
    `t_turnover(dnet, 1)` is tsna's default `tEdgeFormation`/`tEdgeDissolution`
    series on integer spells;
  - path functions default to `end_time = Inf` (tsna's `tPath`) and
    `backward_reachable_set` to `start_time = -Inf`;
  - `tie_decay` defaults to `at` = the last change time, and `t_aggregate`
    collapses the whole time axis;
  - where tsna would use the placeholder `(0, 1)` (a network with no finite
    spell bound), these raise an `ArgumentError` asking for
    `set_observation_period!`.
  With a window, nothing changes. A golden-fixture group of 44 networks without
  a window pins all of this against tsna 0.3.6.
- **`t_edge_density(agg_unit=:dyad)` counts the base network's dyads** as R's
  `network.dyadcount` does: `n₁·n₂` on a two-mode network (`2n₁n₂` directed)
  and `n` more with loops, instead of `n(n−1)` (or `/2`). A zero-length range
  gives `NaN` (as tsna) and an infinite window an `ArgumentError`.
- **Censored bounds are not events.** A left-censored onset is not a formation,
  and a right-censored terminus is not a dissolution, unless
  `include_censored=true` (tsna's `include.censored = FALSE`). A bound is
  censored when it is flagged so or lies outside the window.
- **Elements with no spells are active** (R's `active.default`, through
  DynamicNetworks.jl). They are traversable in paths, count in durations and
  `tie_decay`, and are present in snapshots.
- **`as_contact_sequence`** emits one contact per merged spell. It drops edges
  with no spell record and names them as `:default_active_edges` in the report.
- **The camelCase names of earlier development versions are removed**
  (27 names, among them `tDegree`, `earliestArrival`, `tSnaStats` and
  `shortest_temporal_path`), without a deprecation period, since none was
  released. They were advertised as "R tsna-style", but most are not tsna
  functions, and `tPath`, `tDegree`, `tEdgeFormation`, `tEdgeDissolution` and
  `tSnaStats` mean something else in tsna; `shortest_temporal_path` returned
  the earliest-arrival path, not a shortest one. The documentation's
  *Coming from R tsna* page has a rename table.
- Statistical and path queries reject unobserved dyads by default
  (`missing=:face` opts in to the recorded values).

### Added

- `tied_duration` (`tsna::tiedDuration`) and `t_edge_density`
  (`tsna::tEdgeDensity`), checked against R by the golden fixture.
- A PrecompileTools workload: the time to the first result drops from
  about 10 s to 1.2 s.
- `include_censored` on `t_edge_formation`, `t_edge_dissolution` and
  `t_turnover`; `active_only` on the per-vertex measures; `active_default` on
  `t_edge_duration` and `t_vertex_duration`.
- Earliest-arrival search: a heap-based label-setting search over a memoized
  contact index.
- Batch entry points `earliest_arrival_all`, `temporal_distance_matrix` and
  `reachability_matrix`, with a reusable `TemporalPathWorkspace` and
  `earliest_arrival!`.
- `report=true` on `t_aggregate` and `as_contact_sequence`, returning a
  `NetworkCore.ConversionReport` of what the conversion dropped.
- Validation against R:
  - a golden fixture (`test/fixtures/tsna_semantics.toml`, from
    `test/fixtures/r/tsna_semantics.R`) pinning `edgeDuration`,
    `vertexDuration`, `tEdgeFormation`/`tEdgeDissolution`, `tiedDuration`,
    `tEdgeDensity`, `tPath` and sna snapshot measures on 60 one-mode networks
    with a window, 20 two-mode and looped ones, 44 networks without a window
    and 10 samplk-style panels, compared without going through TSNA.jl's
    internals;
  - an "R tsna" concordance page in the documentation.
- Every export has a runnable docstring example, executed by the test suite.
- Aqua checks in the test suite.

### Fixed

- `t_closeness` was zero for every actor whenever any actor was absent.
- Normalised `t_betweenness` and `t_degree` were deflated by absent actors.
- Adjacent or overlapping spells of one tie were counted as separate
  lifetimes: four formations instead of two, and a total duration of 12 for a
  tie active 10.
- Open-ended spells made duration means `Inf`.
- Left-censored onsets were counted as formations.
- `t_sna_stats`, `window_sna_stats` and `t_turnover` now return concretely
  typed vectors.
- `earliest_arrival`, `earliest_arrival!`, `temporal_distance`,
  `forward_reachable_set`, `backward_reachable_set` and `temporal_path` take
  vertex ids of any `Integer` type. On a `DynamicNetwork{Int32}` they threw a
  `MethodError` unless the id was written `Int32(1)`. An id outside the
  network raises an `ArgumentError` naming it.
- `as_contact_sequence` no longer subtracts unbounded spell bounds. A spell
  with no finite onset gave a contact at `-Inf` on a float axis, a negative
  duration on an integer axis and a meaningless one on a `DateTime` axis; it
  is now dropped and reported (`:unbounded_onsets`). A contact with no finite
  end lasts the axis's unbounded duration (`spell_duration`: `Inf`, or
  `typemax`), and the report counts it (`:unbounded_termini`).
- `t_turnover` and `window_sna_stats` refuse an observation window that
  reaches the extremes of an integer or calendar axis (DynamicNetworks'
  unbounded bounds), as they refuse `±Inf`, instead of tiling it step by step.

### Known limitations

These items are listed under "Not implemented" in the README:

- Not implemented from tsna:
  - `tReach` (seed sampling and the backward direction);
  - `tErgmStats`, `pShiftCount`, `timeProjectedNetwork`;
  - the `tPath` object helpers and plots (`is.tPath`, `as.network.tPath`,
    `plot.tPath`, `plotPaths`).
- Paths:
  - latest-departure paths (`type = "latest.depart"`), traversal costs
    (`graph.step.time`) and `forward.reachable`'s `per.step.depth` are not
    supported;
  - vertex activity does not restrict paths (as in tsna).
- Snapshot statistics:
  - there is no `tSnaStats`-style application of an arbitrary sna function over
    a time grid;
  - interval snapshots (`aggregate.dur`, `rule`) are not supported.
- Durations and events:
  - durations are observed, not censoring-corrected; there is no Kaplan–Meier
    or other censoring-aware estimator, so means are biased downward when many
    spells are censored;
  - `edgeDuration`'s `subject = "dyads"` and `mode = "counts"`, and
    `tEdgeFormation`'s `result.type = "fraction"`, are not supported.
