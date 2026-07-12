# Centrality API Reference

This page documents the temporal centrality functions in TSNA.jl.

## Vertex Centrality

Centrality measures computed on network snapshots at specific time points.

```@docs
t_degree
t_betweenness
t_closeness
t_eigenvector
t_pagerank
```

## Network-Level Measures

Global network statistics computed at specific time points.

```@docs
t_density
t_reciprocity
t_transitivity
```

## Time Series

Functions for computing statistics at multiple time points or in sliding windows.

```@docs
t_sna_stats
window_sna_stats
```
