# Paths API Reference

This page documents the temporal path and reachability functions in TSNA.jl.

## Path Types

```@docs
TemporalPath
path_duration
Contact
ContactSequence
```

## Earliest Arrival and Path Finding

```@docs
earliest_arrival
temporal_distance
temporal_path
shortest_temporal_path
```

## Reachability

```@docs
forward_reachable_set
backward_reachable_set
```

## Batch and In-Place Computation

All-pairs routines and the reusable workspace behind them. `earliest_arrival!`
writes into a caller-supplied [`TemporalPathWorkspace`](@ref) so that a sweep
over many sources allocates once rather than once per source.

```@docs
TemporalPathWorkspace
earliest_arrival!
earliest_arrival_all
temporal_distance_matrix
reachability_matrix
```

## Contact Sequence Conversion

```@docs
as_contact_sequence
```

## R-Style Aliases

Every function above is also exported under its R `tsna`-style camelCase
name; both spellings are interchangeable.

```@docs
tPath
earliestArrival
earliestArrivalAll
temporalDistance
temporalPath
shortestTemporalPath
forwardReachableSet
backwardReachableSet
```
