# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

TSNA.jl is a Julia port of the R `tsna` package (from the StatNet collection) that provides temporal social network analysis tools for dynamic networks, including temporal centrality measures, time-respecting path analysis, reachability, and duration/turnover metrics.

## Development Commands

- **Run tests:** `julia --project -e 'using Pkg; Pkg.test()'`
- **Load package in REPL:** `julia --project -e 'using TSNA'`
- **Build docs:** `julia --project=docs docs/make.jl`
- **Install local deps:** `julia --project -e 'using Pkg; Pkg.instantiate()'` (requires sibling directories for Network, NetworkDynamic, SNA)

## Architecture

The package is a single-file module at `src/TSNA.jl` with no internal submodules. It is organized into these sections:

- **Temporal Path Types** -- `TemporalPath` (R-style alias `tPath`) and `Contact`/`ContactSequence` structs
- **Temporal Measures at a Point** -- `t_degree`, `t_density`, `t_reciprocity`, `t_transitivity`, `t_betweenness`, `t_closeness`, `t_eigenvector`, `t_pagerank` (all take a `DynamicNetwork` and a time point)
- **Temporal Path Finding** -- `earliest_arrival` (heap-based Dijkstra label-setting over a memoized per-network contact index, with INTERVAL semantics: an edge spell [onset, terminus) is boardable at any instant in it, mid-spell boarding allowed, spells active before start_time count, point spells [t,t) usable exactly at t), `temporal_distance` (returns `nothing` when unreachable), `forward_reachable_set`, `backward_reachable_set` (exact dual, computed via forward searches), `temporal_path` (earliest-arrival path; `shortest_temporal_path` is an alias — NOT fewest hops)
- **Duration and Turnover** -- `t_edge_duration`/`t_vertex_duration` (per-spell default matching tsna::edgeDuration, `mode=:total` for per-element sums), `t_edge_formation`/`t_edge_dissolution` (spell onset/terminus EVENT counts in a window; right-censored termini excluded), `t_edge_persistence` (pooled proportion of edges surviving across consecutive windows), `t_turnover` (per-window event counts+rates, consistent NamedTuple shape), `tie_decay`
- **Aggregation and Time Series** -- `t_sna_stats`, `window_sna_stats`, `t_aggregate`

All point-measure functions extract snapshots with `network_extract(dnet, at; retain_all_vertices=true)` — vectors are ALWAYS length nv(dnet), indexed by the original vertex IDs (inactive vertices score 0) — and delegate to SNA.jl (qualified `SNA.` calls; Graphs also exports several of these names). Path/reachability functions operate directly on `dnet.edge_spells`.

## Key Dependencies

- **NetworkDynamic.jl** (local sibling) -- provides `DynamicNetwork`, `network_extract`, `network_collapse`, `active_edges`, edge/vertex spell storage
- **Network.jl** (local sibling) -- base network type
- **SNA.jl** (local sibling) -- static SNA measures (centralities/transitivity delegate to it; call sites must qualify `SNA.` because Graphs exports colliding names)
- **Graphs.jl** -- graph algorithms (centrality, clustering)
- **DataStructures.jl** -- `BinaryMinHeap` for the earliest-arrival search
- **Statistics** -- statistical summaries

Local dependencies (Network, NetworkDynamic, SNA) are referenced via `[sources]` paths in Project.toml pointing to sibling directories (`../Network`, `../NetworkDynamic`, `../SNA`).

## Conventions

- Public functions have snake_case primary names (`t_degree`, `t_edge_duration`, `earliest_arrival`, ...) with R tsna-style camelCase aliases (`tDegree`, `tEdgeDuration`, `earliestArrival`, ...); both are exported and interchangeable. Alias constants carry a one-line "R-style alias" docstring.
- All temporal measure functions are parameterized on `{T, Time}` matching `DynamicNetwork{T, Time}`.
- Aggregate functions accept an `aggregate` keyword (`:mean`, `:median`, `:total`, `:all`).
- The module uses `where {T, Time}` type parameters consistently on all public functions.
- Docstrings follow Julia convention with a signature line, description, and `# Arguments` / `# Methods` sections.
- All public symbols are explicitly exported at the top of the module.
