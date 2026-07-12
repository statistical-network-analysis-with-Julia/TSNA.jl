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
temporalDistance
temporalPath
shortestTemporalPath
forwardReachableSet
backwardReachableSet
```
