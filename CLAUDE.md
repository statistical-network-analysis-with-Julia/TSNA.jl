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

- **Temporal Path Types** -- `tPath` and `Contact`/`ContactSequence` structs
- **Temporal Measures at a Point** -- `tDegree`, `tDensity`, `tReciprocity`, `tTransitivity`, `tBetweenness`, `tCloseness`, `tEigenvector`, `tPagerank` (all take a `DynamicNetwork` and a time point)
- **Temporal Path Finding** -- `earliestArrival` (Dijkstra label-setting with INTERVAL semantics: an edge spell [onset, terminus) is boardable at any instant in it, mid-spell boarding allowed, spells active before start_time count, point spells [t,t) usable exactly at t), `temporalDistance` (returns `nothing` when unreachable), `forwardReachableSet`, `backwardReachableSet` (exact dual, computed via forward searches), `temporalPath` (earliest-arrival path; `shortestTemporalPath` is an alias — NOT fewest hops)
- **Duration and Turnover** -- `tEdgeDuration`/`tVertexDuration` (per-spell default matching tsna::edgeDuration, `mode=:total` for per-element sums), `tEdgeFormation`/`tEdgeDissolution` (spell onset/terminus EVENT counts in a window; right-censored termini excluded), `tTurnover` (per-window event counts+rates, consistent NamedTuple shape), `tieDecay`
- **Aggregation and Time Series** -- `tSnaStats`, `windowSnaStats`, `tAggregate`

All point-measure functions extract snapshots with `network_extract(dnet, at; retain_all_vertices=true)` — vectors are ALWAYS length nv(dnet), indexed by the original vertex IDs (inactive vertices score 0) — and delegate to SNA.jl (qualified `SNA.` calls; Graphs also exports several of these names). Path/reachability functions operate directly on `dnet.edge_spells`.

## Key Dependencies

- **NetworkDynamic.jl** (local sibling) -- provides `DynamicNetwork`, `network_extract`, `network_collapse`, `active_edges`, edge/vertex spell storage
- **Network.jl** (local sibling) -- base network type
- **SNA.jl** (local sibling) -- static SNA measures (centralities/transitivity delegate to it; call sites must qualify `SNA.` because Graphs exports colliding names)
- **Graphs.jl** -- graph algorithms (centrality, clustering)
- **StatsBase.jl** / **Statistics** -- statistical summaries

Local dependencies (Network, NetworkDynamic, SNA) are referenced via `[sources]` paths in Project.toml pointing to sibling directories (`../Network`, `../NetworkDynamic`, `../SNA`).

## Conventions

- Function names use camelCase prefixed with `t` for temporal measures (e.g., `tDegree`, `tEdgeDuration`), matching the R tsna naming convention.
- All temporal measure functions are parameterized on `{T, Time}` matching `DynamicNetwork{T, Time}`.
- Aggregate functions accept an `aggregate` keyword (`:mean`, `:median`, `:total`, `:all`).
- The module uses `where {T, Time}` type parameters consistently on all public functions.
- Docstrings follow Julia convention with a signature line, description, and `# Arguments` / `# Methods` sections.
- All public symbols are explicitly exported at the top of the module.
