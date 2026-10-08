# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

TSNA.jl is a Julia port of the R `tsna` package (from the StatNet collection) that provides temporal social network analysis tools for dynamic networks, including temporal centrality measures, time-respecting path analysis, reachability, and duration/turnover metrics.

## Development Commands

- **Run tests:** `julia --project -e 'using Pkg; Pkg.test()'` (includes the golden-fixture, docstring-example and Aqua testsets)
- **Regenerate the golden fixture:** `Rscript test/fixtures/r/tsna_semantics.R > test/fixtures/tsna_semantics.toml` (needs R with tsna, sna, networkDynamic)
- **Load package in REPL:** `julia --project -e 'using TSNA'`
- **Build docs:** `julia --project=docs docs/make.jl`
- **Install local deps:** `julia --project -e 'using Pkg; Pkg.instantiate()'` (requires sibling directories for Network, DynamicNetworks, SNA)

## Architecture

The package is a single-file module at `src/TSNA.jl` with no internal submodules. It is organized into these sections:

- **Temporal Path Types** -- `TemporalPath` and `Contact`/`ContactSequence` structs
- **Temporal Measures at a Point** -- `t_degree`, `t_density`, `t_reciprocity`, `t_transitivity`, `t_betweenness`, `t_closeness`, `t_eigenvector`, `t_pagerank` (all take a `DynamicNetwork` and a time point)
- **Temporal Path Finding** -- `earliest_arrival` (heap-based Dijkstra label-setting over a memoized per-network contact index, with INTERVAL semantics: an edge spell [onset, terminus) is boardable at any instant in it, mid-spell boarding allowed, spells active before start_time count, point spells [t,t) usable exactly at t), `temporal_distance` (returns `nothing` when unreachable), `forward_reachable_set`, `backward_reachable_set` (exact dual, computed via forward searches), `temporal_path` (earliest-arrival path; NOT a shortest path, which is why the never-released name `shortest_temporal_path` was removed)

- **Vertex arguments** of the path functions (`earliest_arrival[!]`, `temporal_distance`, `forward_/backward_reachable_set`, `temporal_path`, `target=`) are `::Integer`, range-checked and converted by `_vertex_id(dnet, v, role)` — never type them `::T`, which made literal ids fail on `DynamicNetwork{Int32}`. The "… ($T ids)" path testsets run with `Int` and `Int32` ids.
- **Contact sequences** (`as_contact_sequence`): stored activity over the whole axis, never clipped to a window. A merged spell whose onset is the unbounded bound (`unbounded_spell(Time).onset`) has no contact instant: dropped and reported (`:unbounded_onsets`); an open terminus gives `spell_duration`'s unbounded duration and is counted (`:unbounded_termini`). Never subtract spell bounds directly.
- **Batch / all-source temporal paths** (TSNA.jl#1) -- a single-source search allocates four containers (arrival, parent, settled, heap); an *all-source* analysis runs one search per vertex and used to pay that `nv(dnet)` times over, even though the searches are independent and none of the scratch outlives its own search. `TemporalPathWorkspace{T,Time}` holds the containers and `earliest_arrival!(ws, ...)` reuses them. Batch entry points manage the workspace for you: `earliest_arrival_all` (all sources), `temporal_distance_matrix` (all-pairs elapsed times, `nothing` where unreachable), `reachability_matrix`. `backward_reachable_set` — which is inherently one search per vertex — now routes through the workspace too. Measured: **12.7× fewer allocations and ~20× faster at n=50** (2.3× / 1.3× at n=150) versus the per-source loop.

  Two contracts worth knowing: `earliest_arrival!` returns dictionaries that **alias the workspace** and are overwritten by the next search (`earliest_arrival_all` copies them out for you), and `_reset!` must **drain the heap**, because a search that stops early at `target` leaves labels in it that would otherwise leak into the next source. Both are pinned in the "Batch temporal paths" testset, along with the sharpest check available: the batch result must equal the per-source loop exactly, since this is a pure allocation change.
- **Spell tables (R semantics)** -- `_edge_spell_table(dnet; active_default)` lists every BASE edge (in `edges(dnet.network)` order) with its spells MERGED (`DynamicNetworks.merge_spell_vector`); an edge with no spell record gets the unbounded spell under R's `active.default = TRUE` (via `get_edge_activity`). `_vertex_spell_table` likewise. `_clip(spells, lo, hi)` is networkDynamic's `as.data.frame`: keep spells overlapping the window `[lo, hi)`, truncate them to it, and flag a bound censored when stored censored or cut by the window (incl. ±Inf). `_clip` also flags an infinite (axis-extreme) bound censored, as R does (`onset == -Inf`). Every duration/event metric goes through these, so adjacent/overlapping spells are ONE lifetime. Durations go through `_span(s)` (Inf for an unbounded bound, via `DynamicNetworks.unbounded_spell`; never subtract the Int/DateTime extremes).
- **Window semantics (tsna's per-function rules).** `get_observation_period(dnet)` is the window or `nothing`. The helpers next to `_clip` encode what each statistic reads, so do not add a new "default window": `_lifetime_window` (window, else `_axis` = the unbounded spell) for `t_edge_duration`/`t_vertex_duration`/`_edge_events` (tsna's `as.data.frame` with start=-Inf, end=Inf: open spells give `Inf`); `_bounds(dnet, context)` (window, else `extrema(get_change_times)`; `ArgumentError` naming `set_observation_period!` when there is no change time — tsna's `get_bounds` would use `(0, 1)`) for `tied_duration`, `t_edge_density`, `tie_decay`'s default `at`; `_series_range`/`_window_starts` (window half-open, else the CLOSED change-time range, tsna's `seq(start, end)`) for `t_turnover`, `window_sna_stats`, `t_edge_persistence`; `_path_end`/`_path_start` (window, else ±Inf, tsna's `tPath` end = Inf) for every path entry point. `t_aggregate` without a window and without bounds collapses the whole axis.
- **Duration and Turnover** -- `t_edge_duration` (defaults `mode=:total, aggregate=:all` = `tsna::edgeDuration(nd)`: per-edge clipped totals, edges with nothing in the window omitted; `mode=:spell` = tsna `subject="spells"`), `t_vertex_duration` (same, per vertex; CLIPS finite spells, unlike `tsna::vertexDuration`, which only replaces infinite bounds — the fixture test compares only in-window vertices), `t_edge_formation`/`t_edge_dissolution` (onset/terminus EVENT counts of merged clipped spells in a window; censored bounds excluded unless `include_censored=true`, as tsna's `include.censored=FALSE`; tsna's per-time series = `[t_edge_formation(d, t, t+1) for t in start:end]` on integer spells), `t_edge_persistence` (pooled proportion of edges surviving across consecutive windows), `t_turnover` (per-window event counts+rates, concretely typed NamedTuple rows), `tie_decay`
- **Aggregation and Time Series** -- `t_sna_stats`, `window_sna_stats`, `t_aggregate`

### Conversion invariants

TSNA's two conversions honour the **ecosystem conversion contract** (NetworkCore.jl `src/conversion.jl`); the per-path table for the whole ecosystem is `NetworkCore.jl/docs/src/guide/conversion_invariants.md`.

- **`as_contact_sequence` rejects a masked network.** A `Contact` is a tie that happened at a time — there is no contact meaning "we do not know whether this pair ever met" — so a `DynamicNetwork` whose base network carries a missing-dyad mask raises (`missing=:error`, the ecosystem default) instead of being flattened into contacts that read as observed. `missing=:face` is the auditable opt-in. Preserved: vertex count, directedness, and the onset and duration of every MERGED spell (one contact per contiguous activity, as R's `activate.edges` stores it; a point spell `[t,t)` becomes a zero-duration contact). Dropped by nature and named in the report: censoring flags, vertex spells (actor presence), attributes, the observation window, and base edges with no spell record (`:default_active_edges`: active throughout, no finite onset).
- **`t_aggregate` is `network_collapse` plus an aggregation rule** and inherits its invariants: vertex IDs are stable, so directedness, `loops`, two-mode metadata, static attributes and the **missing-dyad mask** all survive all three methods (`:union`, `:intersection`, `:weighted`).
- Both take `report=true`, returning `(result, ::NetworkCore.ConversionReport)` naming each dropped field. Pinned by the three "Conversion invariants: ..." testsets in `test/runtests.jl`.

Per-vertex measures (`t_degree`, `t_betweenness`, `t_closeness`,
`t_eigenvector`, `t_pagerank`) go through `_per_vertex`: compute on the ACTIVE
extract (`network_extract(dnet, at)`, renumbered with `:vertex_pid`) and
scatter back into a length-`nv(dnet)` vector with `NaN` for absent actors
(tsna's tSnaStats uses `network.collapse(at = t)`; absent actors as isolates
zeroed every sna closeness and deflated normalisations). Scalar
density/reciprocity/transitivity and `t_sna_stats` also use active vertices by
default; `active_only=false` opts into the whole universe everywhere. Static
measures are called through one-line `_sna_*` indirections onto SNA.jl's
R-sna names (`SNA.degreecent` = sna::degree, `SNA.betweenness`,
`SNA.closeness`, `SNA.evcent`, `SNA.gden`, `SNA.grecip`, `SNA.gtrans`;
PageRank is `Graphs.pagerank`, not an sna measure), so an SNA.jl rename is a
one-line change. `t_transitivity` has explicit keywords (`measure` for
`gtrans`; `type=:average` + `cmode` for `SNA.transitivity`) — never forward
`kwargs...` to a function whose keyword set differs; the "Keywords forwarded
to SNA.jl's sna-named functions" testset exercises every documented keyword.
`tied_duration`/`t_edge_density` port tsna's `tiedDuration`/`tEdgeDensity`
(fixture-pinned; the event density uses the intended denominator, not tsna
0.3.6's `n_edges · end − start`). Never call the old `SNA.degree_centrality`-style names.
Paths: `_out_contacts` uses `_edge_spell_table`, so base edges with no spells
are traversable at any time; vertex activity does not restrict paths (as
tsna). The memoized contact index's version is `(mutation_count, ne, nv)` of
the network, because a base edge added directly to `dnet.network` changes the
index without a spell mutation. Delegates are qualified
`SNA.` calls and receive the caller's missing-data policy. All statistics and
paths guard the original network before dropping vertices or using cached
spells; `missing=:face` is explicit. `t_aggregate` preserves the mask instead.
`DynamicNetworks.spell_active_at` and `elapsed_seconds` are imported public
helpers; Date/DateTime durations use seconds, numeric axes their native unit.

## Key Dependencies

- **PrecompileTools.jl** -- the `@compile_workload` at the end of the module (README path, both directedness flavours).

- **DynamicNetworks.jl** (local sibling) -- provides `DynamicNetwork`, `network_extract`, `network_collapse`, `active_edges`, edge/vertex spell storage
- **NetworkCore.jl** (local sibling) -- base network type
- **SNA.jl** (local sibling) -- static SNA measures under R sna's names (`degreecent`, `betweenness`, `closeness`, `evcent`, `gden`, `grecip`, `gtrans`; call sites qualify `SNA.` because Graphs exports colliding names)
- **Graphs.jl** -- graph algorithms (centrality, clustering)
- **DataStructures.jl** -- `BinaryMinHeap` for the earliest-arrival search
- **Statistics** -- statistical summaries

Local dependencies NetworkCore, DynamicNetworks and SNA use `[sources]` paths to
sibling directories (`../NetworkCore.jl`, `../DynamicNetworks.jl`, `../SNA.jl`).

## Conventions

- Public functions have snake_case names only. The camelCase names of earlier development versions (and `shortest_temporal_path`) were removed, never released; a testset pins that none is defined and that nothing is deprecated. Do not add aliases. They were never tsna functions, or (`tPath`, `tDegree`, `tEdgeFormation`, `tEdgeDissolution`, `tSnaStats`) meant something else in tsna. The tsna mapping and the rename table live in `docs/src/guide/r_concordance.md`.
- Every exported, non-deprecated name has a docstring with a runnable ```` ```julia ```` example (using only `TSNA`, `DynamicNetworks`, `NetworkCore`); the "Every exported docstring carries a runnable example" testset executes them. `Aqua.test_all(TSNA)` and `detect_ambiguities` run in the suite.
- Golden fixture `test/fixtures/tsna_semantics.toml` (R tsna 0.3.6 / sna 2.8 / networkDynamic 0.12): `case_1..60` one-mode networks over `[0, 12]` (edgeDuration edges/spells, vertexDuration on in-window vertices plus R-computed clipped vertex durations, tEdgeFormation/tEdgeDissolution with/without include.censored, tiedDuration, tEdgeDensity — all exact — and sna degree/betweenness/closeness/gden on `network.extract(at = t)` at 1e-10); `case_61..80` two-mode (directed/undirected) and looped networks (adds `dyad_count`); `nowin_1..44` with NO window (Inf durations, default event series = `t_turnover(d, 1)`, tiedDuration/tEdgeDensity on the change-time range, `tPath` from every vertex, refusals where tsna uses `(0, 1)`); `panel_1..10` samplk-style panels through `DynamicNetwork(networks)`. The test orders R's labelled values by the documented rule (edges by `(i,j)`, spells by edge then onset, vertices by id) and never calls `_edge_spell_table`/`_clip`. Adding a group: put its count in the `[values]` header (TOML tables close over later keys) and write special numbers as `inf`/`nan` (`tnum`).
- `t_edge_density(agg_unit=:dyad)` divides by `_dyad_count(net)`, R's `network.dyadcount` (two-mode `n1*n2`, doubled if directed; `+n` with loops).
- All temporal measure functions are parameterized on `{T, Time}` matching `DynamicNetwork{T, Time}`.
- Aggregate functions accept an `aggregate` keyword (`:mean`, `:median`, `:total`, `:all`).
- The module uses `where {T, Time}` type parameters consistently on all public functions.
- Docstrings follow Julia convention with a signature line, description, and `# Arguments` / `# Methods` sections.
- All public symbols are explicitly exported at the top of the module.
