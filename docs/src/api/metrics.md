# Metrics API Reference

This page documents the duration, persistence, and aggregation functions in TSNA.jl.

## Duration Metrics

```@docs
t_edge_duration
t_vertex_duration
```

## Formation and Dissolution Events

```@docs
t_edge_formation
t_edge_dissolution
```

## Persistence and Turnover

```@docs
t_edge_persistence
t_turnover
tie_decay
```

## Aggregation

```@docs
t_aggregate
```

## R-Style Aliases

Every function above is also exported under its R `tsna`-style camelCase
name; both spellings are interchangeable.

```@docs
tEdgeDuration
tVertexDuration
tEdgeFormation
tEdgeDissolution
tEdgePersistence
tTurnover
tieDecay
tAggregate
```
