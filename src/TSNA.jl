"""
    TSNA.jl - Temporal Social Network Analysis

Tools for analyzing dynamic networks: time-respecting paths and
reachability with full interval-spell semantics, temporal centrality
measures, duration/turnover metrics, and aggregation.

Time-respecting paths use the interval model: an edge with spell
`[onset, terminus)` can be traversed at any instant `t` with
`onset ≤ t < terminus`, so a walker may board mid-spell and waiting at a
vertex is free. Point spells `[t, t)` are instantaneous contacts usable
exactly at `t`.

Activity follows R networkDynamic (via DynamicNetworks.jl): an element with no
spell record is active (`active.default = TRUE`) and spells are merged per
element. With an observation window (`set_observation_period!`), durations
and event counts are taken over it, with spells clipped to it. Without one,
each statistic follows tsna's own rule for a network with no
`net.obs.period`: lifetimes and events over the whole time axis (open spells
last `Inf`), `tiedDuration`/`tEdgeDensity` over the range of the change
times, event series over the closed range of the change times, and paths
without an end. Where tsna would fall back to the placeholder range (0, 1),
TSNA raises an `ArgumentError` asking for a window instead.

Port of the R tsna package from the StatNet collection, with snake_case
names; the documentation's R concordance maps each tsna function.

Statistical and path queries reject masked (unobserved) dyads by default.
Pass `missing=:face` explicitly to use the recorded network and spells.
Conversions that preserve the missing mask, such as `t_aggregate`, retain it.
"""
module TSNA

using DataStructures: BinaryMinHeap
using Dates
using Graphs
using NetworkCore
using DynamicNetworks
using DynamicNetworks: spell_active_at, elapsed_seconds, merge_spell_vector
using SNA
using Statistics
using PrecompileTools: @setup_workload, @compile_workload

# Temporal paths and reachability
export TemporalPath, path_duration
export temporal_distance, earliest_arrival

# Batch / all-source temporal paths: one reused workspace instead of one set of
# scratch containers per source (TSNA.jl#1)
export TemporalPathWorkspace, earliest_arrival!
export earliest_arrival_all
export temporal_distance_matrix, reachability_matrix
export forward_reachable_set, backward_reachable_set
export temporal_path

# Point-in-time measures
export t_degree, t_betweenness, t_closeness, t_eigenvector, t_pagerank
export t_density, t_reciprocity, t_transitivity

# Duration and turnover
export t_edge_duration, t_vertex_duration
export t_edge_formation, t_edge_dissolution
export t_edge_persistence, t_turnover, tie_decay
export tied_duration, t_edge_density

# Contact sequences
export Contact, ContactSequence, as_contact_sequence

# Aggregation and time series
export t_sna_stats, window_sna_stats, t_aggregate

# =============================================================================
# Temporal Path Types
# =============================================================================

"""
    TemporalPath{T, Time}

A time-respecting path through a dynamic network: `times[k]` is the
instant edge `edges[k]` is traversed, and times are non-decreasing.
Returned by [`temporal_path`](@ref).

# Example
```julia
using TSNA
p = TemporalPath([1, 2, 3], [0.0, 4.0], [(1, 2), (2, 3)])
length(p), path_duration(p)        # (2, 4.0)
```
"""
struct TemporalPath{T, Time}
    vertices::Vector{T}
    times::Vector{Time}
    edges::Vector{Tuple{T, T}}

    function TemporalPath{T, Time}(vertices::Vector{T}, times::Vector{Time},
                                   edges::Vector{Tuple{T, T}}) where {T, Time}
        length(times) == length(edges) ||
            throw(ArgumentError("times and edges must have same length"))
        length(vertices) == length(edges) + 1 ||
            throw(ArgumentError("vertices must have length edges + 1"))
        new{T, Time}(vertices, times, edges)
    end
end

TemporalPath(vertices::Vector{T}, times::Vector{Time},
             edges::Vector{Tuple{T,T}}) where {T, Time} =
    TemporalPath{T, Time}(vertices, times, edges)

Base.length(p::TemporalPath) = length(p.edges)

function Base.show(io::IO, p::TemporalPath)
    print(io, "TemporalPath: ")
    for (i, v) in enumerate(p.vertices)
        print(io, v)
        if i <= length(p.times)
            print(io, " --($(p.times[i]))--> ")
        end
    end
end

"""
    path_duration(p::TemporalPath) -> Time difference

Elapsed time between the first and last traversal of the path (zero for a
path with no edges).

# Example
```julia
using TSNA
path_duration(TemporalPath([1, 2, 3], [1.0, 6.5], [(1, 2), (2, 3)]))   # 5.5
```
"""
function path_duration(p::TemporalPath{T, Time}) where {T, Time}
    isempty(p.times) && return _zero_duration(Time)
    return p.times[end] - p.times[1]
end

_zero_duration(::Type{Time}) where Time<:Number = zero(Time)
_zero_duration(::Type{DateTime}) = Millisecond(0)
_zero_duration(::Type{Date}) = Day(0)

# =============================================================================
# Spell tables: merged per element, R active.default, clipped to the window
# =============================================================================

_edge_label(dnet, i, j) = is_directed(dnet) ? (i, j) : (min(i, j), max(i, j))

# Merged spells of every base edge, in `edges(dnet.network)` order. An edge
# with no spell record is active on the whole axis under `active_default`
# (R's active.default = TRUE). Spells are merged per element (R's activate.*
# merges on insertion), so adjacent or overlapping spells are one lifetime.
function _edge_spell_table(dnet::DynamicNetwork{T, Time};
                           active_default::Bool=true) where {T, Time}
    table = Tuple{Tuple{T, T}, Vector{Spell{Time}}}[]
    for e in edges(dnet.network)
        i, j = T(src(e)), T(dst(e))
        sp = get_edge_activity(dnet, i, j; active_default)
        push!(table, (_edge_label(dnet, i, j), merge_spell_vector(sp)))
    end
    return table
end

function _vertex_spell_table(dnet::DynamicNetwork{T, Time};
                             active_default::Bool=true) where {T, Time}
    return [(T(v), merge_spell_vector(get_vertex_activity(dnet, T(v); active_default)))
            for v in 1:nv(dnet)]
end

# R's as.data.frame.networkDynamic: keep the spells that overlap the window
# [lo, hi) (networkDynamic's spells.overlap), truncate them to it, and flag a
# bound censored when it was stored censored, lies outside the window, or is
# infinite (R: `onset < start | onset == -Inf`). With a finite window the
# truncated durations are finite; over the whole axis an open spell keeps its
# infinite bound.
function _clip(spells::Vector{Spell{Time}}, lo::Time, hi::Time) where Time
    win = Spell(lo, hi)
    u = DynamicNetworks.unbounded_spell(Time)
    out = Spell{Time}[]
    for s in spells
        (spell_overlap(s, win) || s == win) || continue
        push!(out, Spell(max(s.onset, lo), min(s.terminus, hi);
                         onset_censored=s.onset_censored || s.onset < lo ||
                                        s.onset == u.onset,
                         terminus_censored=s.terminus_censored || s.terminus > hi ||
                                           s.terminus == u.terminus))
    end
    return out
end

# --- What a statistic reads when no observation window is set -------------
#
# With a window, every statistic clips spells to it. Without one, tsna's rule
# differs by function, and TSNA follows it:
#
# - lifetimes (edgeDuration, vertexDuration) and events (tEdgeFormation,
#   tEdgeDissolution) read the spells unclipped, i.e. over (-Inf, Inf): R's
#   as.data.frame.networkDynamic with start = -Inf, end = Inf. An open spell
#   lasts Inf; only infinite bounds are censored.
# - tiedDuration and tEdgeDensity clip to tsna's get_bounds: the range
#   [first, last] of the change times.
# - series (tEdgeFormation's default seq(start, end)) cover that closed range.
# - paths run to end = Inf (tPath's default end).
#
# tsna's get_bounds falls back to (0, 1) when there is no change time at all;
# TSNA refuses instead (`_bounds`).

# The whole time axis (-Inf, Inf), or the axis extremes on integer and
# calendar axes (DynamicNetworks' convention for unbounded bounds).
function _axis(::DynamicNetwork{T, Time}) where {T, Time}
    u = DynamicNetworks.unbounded_spell(Time)
    return (u.onset, u.terminus)
end

# Lifetimes and events: the window, else the whole axis.
function _lifetime_window(dnet::DynamicNetwork)
    w = get_observation_period(dnet)
    return isnothing(w) ? _axis(dnet) : w
end

# tsna's get_bounds: the window, else the range of the change times. Never
# the (0, 1) placeholder.
function _bounds(dnet::DynamicNetwork, context::AbstractString)
    w = get_observation_period(dnet)
    isnothing(w) || return w
    times = get_change_times(dnet)
    isempty(times) && throw(ArgumentError(
        "$context needs a finite time range, but this network has no observation " *
        "window and no finite spell bound to take one from; give it a window with " *
        "set_observation_period!(dnet, start, stop) (tsna would silently use the " *
        "placeholder range (0, 1) here)"))
    return (first(times), last(times))
end

# Grid of a series: the window, half-open [start, stop); without one, the
# closed range [first, last] of the change times (tsna's seq(start, end)).
function _series_range(dnet::DynamicNetwork, context::AbstractString)
    w = get_observation_period(dnet)
    isnothing(w) || return (w[1], w[2], false)
    lo, hi = _bounds(dnet, context)
    return (lo, hi, true)
end

# Path horizons: the window, else unbounded (tsna's tPath default end = Inf).
_path_end(dnet::DynamicNetwork) =
    (w = get_observation_period(dnet); isnothing(w) ? _axis(dnet)[2] : w[2])
_path_start(dnet::DynamicNetwork) =
    (w = get_observation_period(dnet); isnothing(w) ? _axis(dnet)[1] : w[1])

# Duration of a (clipped) spell in Float64 time units: Inf when a bound is
# unbounded (never a subtraction of the integer or calendar extremes).
function _span(s::Spell{Time}) where Time
    s.onset == s.terminus && return 0.0
    u = DynamicNetworks.unbounded_spell(Time)
    (s.onset == u.onset || s.terminus == u.terminus) && return Inf
    return _dur(s.terminus - s.onset)
end

# =============================================================================
# Contact Sequence
# =============================================================================

"""
    Contact{T, Time}

A single contact (one contiguous edge activity) in a temporal network:
`source`, `target`, onset `time` and `duration`.

# Example
```julia
using TSNA
c = Contact{Int, Float64}(1, 2, 3.0, 2.5)
c.time, c.duration                 # (3.0, 2.5)
```
"""
struct Contact{T, Time}
    source::T
    target::T
    time::Time
    duration
end

"""
    ContactSequence{T, Time}

A sequence of contacts ordered by onset time, with the vertex count and
directedness of the network they came from. Built by
[`as_contact_sequence`](@ref); iterable.

# Example
```julia
using TSNA
cs = ContactSequence([Contact{Int, Float64}(2, 3, 5.0, 1.0),
                      Contact{Int, Float64}(1, 2, 1.0, 2.0)], 3)
[c.time for c in cs]               # [1.0, 5.0]
```
"""
struct ContactSequence{T, Time}
    contacts::Vector{Contact{T, Time}}
    n_vertices::Int
    directed::Bool

    function ContactSequence(contacts::Vector{Contact{T, Time}}, n::Integer;
                             directed::Bool=true) where {T, Time}
        sorted = sort(contacts, by=c -> c.time)
        new{T, Time}(sorted, n, directed)
    end
end

Base.length(cs::ContactSequence) = length(cs.contacts)
Base.iterate(cs::ContactSequence, state=1) =
    state > length(cs) ? nothing : (cs.contacts[state], state + 1)

"""
    as_contact_sequence(dnet::DynamicNetwork; missing=:error, report=false) -> ContactSequence

Convert a dynamic network's edge activity to a contact sequence: one
`Contact` per contiguous activity of an edge (its spells merged, so adjacent
or overlapping spells are one contact, as R's `activate.edges` stores them),
carrying its onset and duration (DynamicNetworks' `spell_duration`). A point spell
`[t,t)` becomes a contact of zero duration.

The contacts are the stored activity over the whole time axis; they are not
clipped to an observation window. A spell with no finite start (onset `-Inf`,
or the axis minimum on integer and calendar axes;
see `DynamicNetworks.unbounded_spell`) has no instant at which to place a
contact, so it is dropped and reported, like an edge with no spell record. A
spell with a finite start and no finite end becomes a contact whose duration
is the axis's unbounded duration (`Inf` on a floating-point axis, `typemax` on
an integer axis, `Millisecond(typemax(Int64))` on a `DateTime` axis), as
`spell_duration` gives it; the report counts these contacts too.

# Conversion invariants

Preserved: the vertex count, directedness, and the onset and duration of
every (merged) edge spell, so the activity is reconstructable.

A `ContactSequence` has no slot for the rest of the dynamic network, so the
conversion is lossy by nature: **spell censoring flags**, **vertex activity
spells** (actor presence/composition), static and time-varying attributes, and
the observation window are dropped. Base-network edges with **no spell
record** (active throughout, by default) and spells with **no finite onset**
are dropped too, and contacts with **no finite end** are counted. Pass
`report=true` for `(cs, ::NetworkCore.ConversionReport)` naming each.

A `Contact` cannot record that a dyad is *unobserved*, so a network with a
missing-dyad mask is **rejected** by default (`missing=:error`) rather than
being flattened into contacts that read as observed. Pass `missing=:face` to
convert the recorded face values anyway (the mask is then dropped, and the
report says so).

# Example
```julia
using TSNA, DynamicNetworks
dnet = DynamicNetwork(3; observation_start=0.0, observation_end=10.0)
activate!(dnet, 0.0, 4.0; edge=(1, 2))
activate!(dnet, 7.0, 7.0; edge=(2, 3))     # instantaneous contact
cs = as_contact_sequence(dnet)
[(c.source, c.target, c.time, c.duration) for c in cs]
# [(1, 2, 0.0, 4.0), (2, 3, 7.0, 0.0)]
```
"""
function as_contact_sequence(dnet::DynamicNetwork{T, Time};
                             missing::Symbol=:error,
                             report::Bool=false) where {T, Time}
    require_observed(dnet.network, missing; context="as_contact_sequence")

    contacts = Contact{T, Time}[]
    n_default = 0
    n_open_onset = 0      # spells with no finite start: no contact instant
    n_open_end = 0        # contacts with no finite end: unbounded duration
    u = DynamicNetworks.unbounded_spell(Time)
    for e in edges(dnet.network)
        i, j = T(src(e)), T(dst(e))
        key = _edge_label(dnet, i, j)
        rec = get(dnet.edge_spells, key, nothing)
        if isnothing(rec)
            n_default += 1
            continue
        end
        for spell in merge_spell_vector(rec)
            # Never subtract the axis extremes (that overflowed on integer and
            # calendar axes, and placed contacts at -Inf on float axes).
            if spell.onset == u.onset
                n_open_onset += 1
                continue
            end
            spell.terminus == u.terminus && spell.onset != spell.terminus &&
                (n_open_end += 1)
            push!(contacts, Contact{T, Time}(key[1], key[2], spell.onset,
                                             spell_duration(spell)))
        end
    end

    cs = ContactSequence(contacts, Int(nv(dnet)); directed=is_directed(dnet))

    rep = ConversionReport(:DynamicNetwork, :ContactSequence)
    record_drop!(rep, :spell_censoring,
                 "a Contact records onset and duration only; onset/terminus " *
                 "censoring flags have no slot")
    record_drop!(rep, :vertex_spells,
                 "vertex activity (actor presence/composition) is not " *
                 "representable in a contact sequence")
    record_drop!(rep, :attributes,
                 "static and time-varying vertex/edge/network attributes are " *
                 "not carried")
    window = get_observation_period(dnet)
    isnothing(window) || record_drop!(rep, :observation_period,
                                      "the observation window $window has no " *
                                      "contact-sequence counterpart")
    n_default > 0 && record_drop!(rep, :default_active_edges,
                 "$n_default base-network edge(s) have no spell record (active " *
                 "throughout by default) and no finite onset, so no contact")
    n_open_onset > 0 && record_drop!(rep, :unbounded_onsets,
                 "$n_open_onset spell(s) have no finite onset (active since the " *
                 "start of the time axis), so no contact instant; they are dropped")
    n_open_end > 0 && record_drop!(rep, :unbounded_termini,
                 "$n_open_end contact(s) have no finite end; their duration is " *
                 "the unbounded duration $(repr(spell_duration(u))), not a measured one")
    n_mask = n_missing_dyads(dnet.network)
    n_mask > 0 && record_drop!(rep, :missing_dyads,
                               "$n_mask masked dyad(s) converted at face value " *
                               "under missing=:face; the contacts do not record " *
                               "that they are unobserved")

    return report ? (cs, rep) : cs
end

# =============================================================================
# Temporal Path Finding (interval semantics)
# =============================================================================

# Per-vertex outgoing contacts: v -> [(neighbor, spell), ...]. For
# undirected networks each spell is listed from both endpoints. A base edge
# with no spell record is active throughout (R's active.default), so it
# contributes the unbounded spell. Vertex activity does not restrict paths
# (as in tsna::tPath).
function _out_contacts(dnet::DynamicNetwork{T, Time}) where {T, Time}
    out = Dict{T, Vector{Tuple{T, Spell{Time}}}}()
    directed = is_directed(dnet)
    for ((i, j), spells) in _edge_spell_table(dnet)
        for spell in spells
            push!(get!(out, i, Tuple{T, Spell{Time}}[]), (j, spell))
            if !directed
                push!(get!(out, j, Tuple{T, Spell{Time}}[]), (i, spell))
            end
        end
    end
    return out
end

# Memoized contact index: rebuilding the per-vertex contact lists on every
# path query is O(total spells), which dominated repeated-query workloads
# (e.g. backward_reachable_set runs one search per vertex). The cache is
# keyed weakly by network identity and invalidated via the network's
# `mutation_count`, which DynamicNetworks bumps on every spell mutation, and
# by the base network's edge and vertex counts (edges with no spell record
# are active, so adding one directly to `dnet.network` changes the index).
const _CONTACT_INDEX_LOCK = ReentrantLock()
const _CONTACT_INDEX_CACHE = WeakKeyDict{Any, Tuple{Tuple{Int,Int,Int}, Any}}()

_index_version(dnet) = (dnet.mutation_count, Int(ne(dnet.network)), Int(nv(dnet.network)))

function _contact_index(dnet::DynamicNetwork{T, Time}) where {T, Time}
    lock(_CONTACT_INDEX_LOCK) do
        version = _index_version(dnet)
        entry = get(_CONTACT_INDEX_CACHE, dnet, nothing)
        if !isnothing(entry) && entry[1] == version
            return entry[2]::Dict{T, Vector{Tuple{T, Spell{Time}}}}
        end
        index = _out_contacts(dnet)
        _CONTACT_INDEX_CACHE[dnet] = (version, index)
        return index
    end
end

# A vertex argument of any Integer type, checked against the network and
# converted to its id type (so literal ids work on a DynamicNetwork{Int32}).
function _vertex_id(dnet::DynamicNetwork{T}, v::Integer, role::AbstractString) where T
    1 <= v <= nv(dnet) || throw(ArgumentError(
        "$role must be a vertex of the network (1:$(nv(dnet))); got $v"))
    return T(v)
end

# Can an edge with `spell` be boarded by a walker present from time `t`,
# before `end_time`? Returns the boarding instant or nothing.
function _board_time(spell, t, end_time)
    if spell.onset == spell.terminus
        # Point contact: usable exactly at its instant
        (t <= spell.onset && spell.onset < end_time) && return spell.onset
        return nothing
    end
    t >= spell.terminus && return nothing
    depart = max(t, spell.onset)
    depart < end_time || return nothing
    return depart
end

"""
    earliest_arrival(dnet, source, start_time; end_time=window end or Inf,
                     target=nothing, missing=:error) -> (arrival::Dict, parent::Dict)

Earliest-arrival times from `source` to every vertex, starting at
`start_time`, under interval semantics: an edge spell `[onset, terminus)`
is traversable at any instant in it (boarding mid-spell is allowed;
spells that began before `start_time` but are still active count). Edges
with no spell record are traversable at any time (R's `active.default`);
vertex activity does not restrict paths, as in `tsna::tPath`.
Implemented as a heap-based Dijkstra label-setting search over a memoized
per-network contact index, so chains of simultaneous (equal-onset) spells
are handled correctly and repeated queries do not rebuild the index.
`end_time` (exclusive) defaults to the end of the observation window
(`get_observation_period`) and, when no window is set, to `Inf` (no
horizon), as `tsna::tPath`'s `end` does. (tsna's `tPath` ignores a window;
here a window bounds the observation, so paths stop at its end.)

With `target` set, the search stops as soon as that vertex is settled
(its arrival time is already final); the returned dictionaries then only
cover the explored part of the network.

Returns the arrival-time dictionary (vertices absent = unreachable) and
the parent map `(vertex => (predecessor, boarding time))` for path
reconstruction. Agrees with `tsna::tPath(direction = "fwd")` on randomised
networks with and without a window (the test suite's golden fixture).

# Example
```julia
using TSNA, DynamicNetworks
dnet = DynamicNetwork(3; observation_start=0.0, observation_end=10.0)
activate!(dnet, 1.0, 3.0; edge=(1, 2))
activate!(dnet, 5.0, 6.0; edge=(2, 3))
arrival, parent = earliest_arrival(dnet, 1, 0.0)
arrival[3]                          # 5.0
```
"""
function earliest_arrival(dnet::DynamicNetwork{T, Time}, source::Integer, start_time;
                          end_time=_path_end(dnet),
                          target::Union{Nothing, Integer}=nothing,
                          missing::Symbol=:error) where {T, Time}
    ws = TemporalPathWorkspace{T, Time}()
    return earliest_arrival!(ws, dnet, source, start_time;
                             end_time=end_time, target=target, missing=missing)
end

"""
    TemporalPathWorkspace{T, Time}()

Reusable scratch space for the earliest-arrival search: the arrival and parent
maps, the settled set, and the heap.

A single-source search allocates all four containers. An **all-source** analysis
(temporal closeness, betweenness, `backward_reachable_set`, a reachability
matrix) runs one search per vertex, so it pays that allocation `nv(dnet)` times
over — even though the searches are independent and none of the scratch outlives
its own search. A workspace is filled, read, and emptied once per source instead.

Pass one to [`earliest_arrival!`](@ref), or use the batch entry points
([`earliest_arrival_all`](@ref), [`temporal_distance_matrix`](@ref),
[`reachability_matrix`](@ref)), which manage it for you.

The memoized contact index is shared across searches regardless (it is cached on
the network); the workspace is about the *per-search* containers.

# Example
```julia
using TSNA, DynamicNetworks
dnet = DynamicNetwork(3; observation_start=0.0, observation_end=10.0)
activate!(dnet, 1.0, 3.0; edge=(1, 2))
ws = TemporalPathWorkspace{Int, Float64}()
arrival, _ = earliest_arrival!(ws, dnet, 1, 0.0)
arrival[2]                          # 1.0
```
"""
struct TemporalPathWorkspace{T, Time}
    arrival::Dict{T, Time}
    parent::Dict{T, Tuple{T, Time}}
    settled::Set{T}
    heap::BinaryMinHeap{Tuple{Time, T}}
end

TemporalPathWorkspace{T, Time}() where {T, Time} =
    TemporalPathWorkspace{T, Time}(Dict{T, Time}(), Dict{T, Tuple{T, Time}}(),
                                   Set{T}(), BinaryMinHeap{Tuple{Time, T}}())

function _reset!(ws::TemporalPathWorkspace)
    empty!(ws.arrival)
    empty!(ws.parent)
    empty!(ws.settled)
    # BinaryMinHeap has no `empty!`; drain it. It is empty already whenever the
    # previous search ran to exhaustion, so this only costs anything after an
    # early `target` break.
    while !isempty(ws.heap)
        pop!(ws.heap)
    end
    return ws
end

"""
    earliest_arrival!(ws::TemporalPathWorkspace, dnet, source, start_time;
                      end_time=window end or Inf, target=nothing, missing=:error)
        -> (arrival::Dict, parent::Dict)

In-place [`earliest_arrival`](@ref): runs the same search, reusing `ws`'s
containers instead of allocating fresh ones.

**The returned dictionaries alias `ws`** and are overwritten by the next search
on the same workspace. Copy what you need to keep, or use
[`earliest_arrival_all`](@ref), which does that for you.

# Example
```julia
using TSNA, DynamicNetworks
dnet = DynamicNetwork(3; observation_start=0.0, observation_end=10.0)
activate!(dnet, 2.0, 4.0; edge=(1, 2))
ws = TemporalPathWorkspace{Int, Float64}()
first_run = copy(earliest_arrival!(ws, dnet, 1, 0.0)[1])   # copy: ws is reused
first_run[2]                        # 2.0
```
"""
function earliest_arrival!(ws::TemporalPathWorkspace{T, Time},
                           dnet::DynamicNetwork{T, Time}, source::Integer, start_time;
                           end_time=_path_end(dnet),
                           target::Union{Nothing, Integer}=nothing,
                           missing::Symbol=:error) where {T, Time}
    require_observed(dnet.network, missing; context="earliest_arrival!")
    source = _vertex_id(dnet, source, "source")
    target = isnothing(target) ? nothing : _vertex_id(dnet, target, "target")
    start_time = convert(Time, start_time)
    end_time = convert(Time, end_time)
    start_time <= end_time || throw(ArgumentError("start_time must be <= end_time"))
    out = _contact_index(dnet)
    no_contacts = Tuple{T, Spell{Time}}[]

    _reset!(ws)
    arrival, parent, settled, heap = ws.arrival, ws.parent, ws.settled, ws.heap
    arrival[source] = start_time

    # Label-setting search (Dijkstra with a binary min-heap and lazy
    # deletion; boarding times never precede the label being settled, so
    # labels are final once popped)
    push!(heap, (start_time, source))

    while !isempty(heap)
        t, v = pop!(heap)
        v in settled && continue
        push!(settled, v)
        v === target && break

        for (w, spell) in get(out, v, no_contacts)
            w in settled && continue
            depart = _board_time(spell, t, end_time)
            isnothing(depart) && continue
            if !haskey(arrival, w) || depart < arrival[w]
                arrival[w] = depart
                parent[w] = (v, depart)
                push!(heap, (depart, w))
            end
        end
    end

    return arrival, parent
end

"""
    earliest_arrival_all(dnet, start_time; sources=all vertices,
                         end_time=window end or Inf, missing=:error) -> Dict{T, Dict{T, Time}}

Earliest-arrival times **from every source in one pass**, reusing a single
[`TemporalPathWorkspace`](@ref) across the searches.

The searches are independent, so the result is identical to calling
[`earliest_arrival`](@ref) per source — but the scratch containers are allocated
once rather than once per source, and the memoized contact index is built at most
once. This is the entry point for any all-source analysis (temporal closeness,
betweenness, reachability); running the single-source function in a loop is the
pattern this exists to replace.

Each source's arrival map is copied out of the workspace before the next search
overwrites it, so the returned dictionaries are independent and safe to keep.

# Example
```julia
using TSNA, DynamicNetworks
dnet = DynamicNetwork(3; observation_start=0.0, observation_end=10.0)
activate!(dnet, 1.0, 3.0; edge=(1, 2))
activate!(dnet, 5.0, 6.0; edge=(2, 3))
all_arr = earliest_arrival_all(dnet, 0.0)
all_arr[1][3], haskey(all_arr[3], 1)   # (5.0, false)
```
"""
function earliest_arrival_all(dnet::DynamicNetwork{T, Time}, start_time;
                              sources=T.(1:nv(dnet)),
                              end_time=_path_end(dnet),
                              missing::Symbol=:error) where {T, Time}
    require_observed(dnet.network, missing; context="earliest_arrival_all")
    start_time = convert(Time, start_time)
    end_time = convert(Time, end_time)
    start_time <= end_time || throw(ArgumentError("start_time must be <= end_time"))
    ws = TemporalPathWorkspace{T, Time}()
    result = Dict{T, Dict{T, Time}}()
    for s in sources
        arrival, _ = earliest_arrival!(ws, dnet, s, start_time; end_time=end_time, missing=missing)
        result[T(s)] = copy(arrival)      # detach from the workspace (s is checked)
    end
    return result
end

"""
    temporal_distance_matrix(dnet, start_time; end_time=window end or Inf, missing=:error)
        -> Matrix{Union{Duration, Nothing}}

All-pairs temporal distances (elapsed time of the earliest time-respecting path,
`arrival − start_time`), in one batched pass. `nothing` where no path exists
before `end_time`; zero on the diagonal.

Computed with a single reused workspace via [`earliest_arrival_all`](@ref)
rather than `nv(dnet)²` independent [`temporal_distance`](@ref) calls.

# Example
```julia
using TSNA, DynamicNetworks
dnet = DynamicNetwork(3; observation_start=0.0, observation_end=10.0)
activate!(dnet, 1.0, 3.0; edge=(1, 2))
activate!(dnet, 5.0, 6.0; edge=(2, 3))
D = temporal_distance_matrix(dnet, 0.0)
D[1, 3], D[3, 1]                    # (5.0, nothing)
```
"""
function temporal_distance_matrix(dnet::DynamicNetwork{T, Time}, start_time;
                                  end_time=_path_end(dnet),
                                  missing::Symbol=:error) where {T, Time}
    n = nv(dnet)
    t0 = convert(Time, start_time)
    all_arr = earliest_arrival_all(dnet, t0; end_time=end_time, missing=missing)
    Duration = typeof(t0 - t0)
    D = Matrix{Union{Duration, Nothing}}(nothing, n, n)
    for i in 1:n
        arr = all_arr[T(i)]
        for j in 1:n
            haskey(arr, T(j)) && (D[i, j] = arr[T(j)] - t0)
        end
    end
    return D
end

"""
    reachability_matrix(dnet, start_time; end_time=window end or Inf, missing=:error) -> BitMatrix

`R[i, j]` is `true` when a time-respecting path runs from `i` to `j` within the
window (the diagonal is `true`). One batched pass; see
[`earliest_arrival_all`](@ref).

# Example
```julia
using TSNA, DynamicNetworks
dnet = DynamicNetwork(3; observation_start=0.0, observation_end=10.0)
activate!(dnet, 5.0, 6.0; edge=(1, 2))
activate!(dnet, 1.0, 3.0; edge=(2, 3))     # closes before 1 reaches 2
R = reachability_matrix(dnet, 0.0)
R[1, 2], R[1, 3]                    # (true, false)
```
"""
function reachability_matrix(dnet::DynamicNetwork{T, Time}, start_time;
                             end_time=_path_end(dnet),
                             missing::Symbol=:error) where {T, Time}
    n = nv(dnet)
    all_arr = earliest_arrival_all(dnet, start_time; end_time=end_time, missing=missing)
    R = falses(n, n)
    for i in 1:n
        arr = all_arr[T(i)]
        for j in 1:n
            R[i, j] = haskey(arr, T(j))
        end
    end
    return R
end

"""
    temporal_distance(dnet, source, target, start_time; end_time=window end or Inf, missing=:error)
        -> elapsed time or nothing

Elapsed time of the earliest time-respecting path from `source` to
`target` departing at `start_time` (arrival − start). Returns `nothing`
when no path exists before `end_time`.

# Example
```julia
using TSNA, DynamicNetworks
dnet = DynamicNetwork(3; observation_start=0.0, observation_end=10.0)
activate!(dnet, 1.0, 3.0; edge=(1, 2))
activate!(dnet, 5.0, 6.0; edge=(2, 3))
temporal_distance(dnet, 1, 3, 0.0)   # 5.0
temporal_distance(dnet, 3, 1, 0.0)   # nothing
```
"""
function temporal_distance(dnet::DynamicNetwork{T, Time}, source::Integer, target::Integer,
                           start_time; end_time=_path_end(dnet),
                           missing::Symbol=:error) where {T, Time}
    target = _vertex_id(dnet, target, "target")
    arrival, _ = earliest_arrival(dnet, source, start_time;
                                  end_time=end_time, target=target, missing=missing)
    haskey(arrival, target) || return nothing
    return arrival[target] - convert(Time, start_time)
end

"""
    forward_reachable_set(dnet, source, start_time; end_time=window end or Inf, missing=:error) -> Vector

Vertices reachable from `source` by a time-respecting path departing at
or after `start_time` and arriving before `end_time` (`source` included).
The sorted vertex set of `tsna::tPath(nd, v, direction = "fwd")`.

# Example
```julia
using TSNA, DynamicNetworks
dnet = DynamicNetwork(4; observation_start=0.0, observation_end=10.0)
activate!(dnet, 1.0, 3.0; edge=(1, 2))
activate!(dnet, 5.0, 6.0; edge=(2, 3))
forward_reachable_set(dnet, 1, 0.0)  # [1, 2, 3]
```
"""
function forward_reachable_set(dnet::DynamicNetwork{T, Time}, source::Integer, start_time;
                               end_time=_path_end(dnet),
                               missing::Symbol=:error) where {T, Time}
    arrival, _ = earliest_arrival(dnet, source, start_time; end_time=end_time, missing=missing)
    return sort(collect(keys(arrival)))
end

"""
    backward_reachable_set(dnet, target, end_time; start_time=window start or -Inf, missing=:error) -> Vector

Vertices from which `target` can be reached by a time-respecting path in
`[start_time, end_time)` — the exact dual of
[`forward_reachable_set`](@ref) (computed by forward searches that stop
as soon as `target` is settled, so both use identical traversal
semantics). `tsna::tPath(direction = "bkwd")` does not chain same-instant
point spells although its forward search does; this function stays the
exact dual of the forward search.

# Example
```julia
using TSNA, DynamicNetworks
dnet = DynamicNetwork(4; observation_start=0.0, observation_end=10.0)
activate!(dnet, 1.0, 3.0; edge=(1, 2))
activate!(dnet, 5.0, 6.0; edge=(2, 3))
backward_reachable_set(dnet, 3, 10.0)   # [1, 2, 3]
```
"""
function backward_reachable_set(dnet::DynamicNetwork{T, Time}, target::Integer, end_time;
                                start_time=_path_start(dnet),
                                missing::Symbol=:error) where {T, Time}
    require_observed(dnet.network, missing; context="backward_reachable_set")
    target = _vertex_id(dnet, target, "target")
    start_time <= end_time || throw(ArgumentError("start_time must be <= end_time"))
    # One search per vertex — the all-source pattern. Reuse a single workspace
    # rather than allocating the arrival/parent/settled/heap containers nv(dnet)
    # times over. (`target` still short-circuits each search, so this keeps the
    # early exit; only the allocation is shared.)
    ws = TemporalPathWorkspace{T, Time}()
    reachable = T[]
    for v in 1:nv(dnet)
        vT = T(v)
        if vT == target
            push!(reachable, vT)
            continue
        end
        arrival, _ = earliest_arrival!(ws, dnet, vT, start_time;
                                       end_time=end_time, target=target, missing=missing)
        haskey(arrival, target) && push!(reachable, vT)
    end
    return reachable
end

"""
    temporal_path(dnet, source, target, start_time; end_time=window end or Inf, missing=:error)
        -> Union{TemporalPath, Nothing}

The **earliest-arrival** (foremost) time-respecting path from `source` to
`target` departing at `start_time` — not necessarily the path with the fewest
hops or the shortest elapsed time. Returns `nothing` when no path exists.

# Example
```julia
using TSNA, DynamicNetworks
dnet = DynamicNetwork(3; observation_start=0.0, observation_end=10.0)
activate!(dnet, 1.0, 3.0; edge=(1, 2))
activate!(dnet, 5.0, 6.0; edge=(2, 3))
p = temporal_path(dnet, 1, 3, 0.0)
p.vertices, p.times                 # ([1, 2, 3], [1.0, 5.0])
```
"""
function temporal_path(dnet::DynamicNetwork{T, Time}, source::Integer, target::Integer,
                       start_time; end_time=_path_end(dnet),
                       missing::Symbol=:error) where {T, Time}
    source = _vertex_id(dnet, source, "source")
    target = _vertex_id(dnet, target, "target")
    arrival, parent = earliest_arrival(dnet, source, start_time;
                                       end_time=end_time, target=target, missing=missing)
    haskey(arrival, target) || return nothing

    verts = T[target]
    times = Time[]
    path_edges = Tuple{T, T}[]
    v = target
    while v != source
        u, t = parent[v]
        pushfirst!(verts, u)
        pushfirst!(times, t)
        pushfirst!(path_edges, (u, v))
        v = u
    end

    return TemporalPath(verts, times, path_edges)
end

# =============================================================================
# Temporal Measures at a Point
# =============================================================================
#
# Per-vertex measures are computed on the network of the vertices ACTIVE at
# the instant (as tsna's tSnaStats does, via network.collapse(at = t)) and
# scattered back into a vector indexed by the dynamic network's vertex IDs,
# with NaN for absent actors. Inactive actors must not enter the computation:
# as isolates they zero every sna closeness score and inflate the n of
# normalised betweenness, degree and PageRank. `active_only=false` computes
# on the whole vertex universe (inactive actors as isolates) instead.

function _snapshot(dnet, at; missing::Symbol=:error, active_only::Bool=true)
    # Check before dropping inactive endpoints: their missing ties must not
    # disappear before the caller has chosen a policy.
    require_observed(dnet.network, missing; context="temporal snapshot statistic")
    return network_extract(dnet, at; retain_all_vertices=!active_only)
end

function _per_vertex(f, dnet::DynamicNetwork, at; missing::Symbol, active_only::Bool)
    snap = _snapshot(dnet, at; missing, active_only)
    active_only || return Vector{Float64}(f(snap))
    out = fill(NaN, Int(nv(dnet)))
    nv(snap) == 0 && return out
    vals = f(snap)
    for v in 1:nv(snap)
        out[Int(get_vertex_attribute(snap, :vertex_pid, v))] = vals[v]
    end
    return out
end

# The static measures are SNA.jl's R-sna-named functions (qualified calls).
# One indirection per measure, so a rename in SNA.jl is one line here.
# PageRank is not an sna measure; it is Graphs.jl's `pagerank` on the snapshot.
const _DEGREE_CMODE = Dict(:total => :freeman, :in => :indegree, :out => :outdegree)
function _sna_degree(net; mode::Symbol=:total, kwargs...)
    haskey(_DEGREE_CMODE, mode) ||
        throw(ArgumentError("mode must be :total, :in or :out (got :$mode)"))
    return SNA.degreecent(net; cmode=_DEGREE_CMODE[mode], kwargs...)
end
_sna_betweenness(net; kwargs...) = SNA.betweenness(net; kwargs...)
_sna_closeness(net; kwargs...) = SNA.closeness(net; kwargs...)
_sna_eigenvector(net; kwargs...) = SNA.evcent(net; kwargs...)
_sna_pagerank(net; α::Float64=0.85, missing::Symbol=:error) =
    nv(net) == 0 ? Float64[] : Graphs.pagerank(net, α)
_sna_density(net; kwargs...) = SNA.gden(net; kwargs...)
_sna_reciprocity(net; method::Symbol=:dyadic, kwargs...) = SNA.grecip(net; measure=method, kwargs...)
_sna_transitivity(net; measure::Symbol=:weak, missing::Symbol=:error) =
    SNA.gtrans(net; measure, missing)
_sna_average_clustering(net; cmode::Symbol=:total, missing::Symbol=:error) =
    SNA.transitivity(net; type=:average, cmode, missing)

"""
    t_degree(dnet::DynamicNetwork, at; mode=:total, active_only=true, missing=:error, kwargs...)
        -> Vector{Float64}

Degree of each actor at time `at` (`SNA.degreecent`, R's `sna::degree`:
`mode=:total` is sna's `cmode = "freeman"`, `:in`/`:out` its
`"indegree"`/`"outdegree"`), computed on the network of the vertices
active at `at` (as tsna's `tSnaStats(nd, "degree")` does) and indexed by the
dynamic network's vertex IDs: an actor absent at `at` gets `NaN`, and
`normalized=true` divides by the number of *active* actors. With
`active_only=false` the whole vertex universe is used (absent actors as
isolates scoring 0). Further keywords go to `SNA.degreecent`: `rescale`,
`normalized`, `diag`, `ignore_eval`, `attr`.

# Example
```julia
using TSNA, DynamicNetworks
dnet = DynamicNetwork(4; observation_start=0.0, observation_end=10.0)
activate!(dnet, 0.0, 10.0; edge=(1, 2))
activate!(dnet, 0.0, 10.0; edge=(2, 3))
activate!(dnet, 6.0, 10.0; vertex=4)       # actor 4 joins at t = 6
t_degree(dnet, 2.0)                 # [1.0, 2.0, 1.0, NaN]
```
"""
function t_degree(dnet::DynamicNetwork{T, Time}, at; mode::Symbol=:total,
                  active_only::Bool=true, missing::Symbol=:error,
                  kwargs...) where {T, Time}
    return _per_vertex(net -> _sna_degree(net; mode, missing, kwargs...),
                       dnet, at; missing, active_only)
end

"""
    t_betweenness(dnet::DynamicNetwork, at; normalized=false, active_only=true,
                  missing=:error, kwargs...) -> Vector{Float64}

Betweenness of each actor at time `at` (`SNA.betweenness`, raw scores by
default as in sna), computed on the active network and scattered back with `NaN`
for absent actors (see [`t_degree`](@ref)). Further keywords go to
`SNA.betweenness`: `cmode` (`:directed`/`:undirected`) and `rescale`.

# Example
```julia
using TSNA, DynamicNetworks
dnet = DynamicNetwork(4; directed=false, observation_start=0.0, observation_end=10.0)
activate!(dnet, 0.0, 10.0; edge=(1, 2))
activate!(dnet, 0.0, 10.0; edge=(2, 3))
activate!(dnet, 6.0, 10.0; vertex=4)
t_betweenness(dnet, 2.0)            # [0.0, 1.0, 0.0, NaN]
```
"""
t_betweenness(dnet::DynamicNetwork, at; normalized::Bool=false, active_only::Bool=true,
              missing::Symbol=:error, kwargs...) =
    _per_vertex(net -> _sna_betweenness(net; normalized, missing, kwargs...),
                dnet, at; missing, active_only)

"""
    t_closeness(dnet::DynamicNetwork, at; active_only=true, missing=:error, kwargs...)
        -> Vector{Float64}

Closeness of each actor at time `at` (`SNA.closeness`, R's `sna::closeness`:
`(n−1)/Σd`, and 0 for an actor that cannot reach every other active actor), computed on the
active network and scattered back with `NaN` for absent actors. Absent actors
used to enter as isolates and zero every score. Further keywords go to
`SNA.closeness`: `cmode` (`:directed`, `:undirected`, `:suminvdir`,
`:suminvundir`) and `rescale`.

# Example
```julia
using TSNA, DynamicNetworks
dnet = DynamicNetwork(4; directed=false, observation_start=0.0, observation_end=10.0)
activate!(dnet, 0.0, 10.0; edge=(1, 2))
activate!(dnet, 0.0, 10.0; edge=(2, 3))
activate!(dnet, 6.0, 10.0; vertex=4)
round.(t_closeness(dnet, 2.0); digits=3)   # [0.667, 1.0, 0.667, NaN]
```
"""
t_closeness(dnet::DynamicNetwork, at; active_only::Bool=true, missing::Symbol=:error,
            kwargs...) =
    _per_vertex(net -> _sna_closeness(net; missing, kwargs...),
                dnet, at; missing, active_only)

"""
    t_eigenvector(dnet::DynamicNetwork, at; active_only=true, missing=:error, kwargs...)
        -> Vector{Float64}

Eigenvector centrality of each actor at time `at` (`SNA.evcent`), computed
on the active network and scattered back with `NaN` for absent actors. Further keywords go
to `SNA.evcent`: `rescale`, `tol`, `ignore_eval`, `attr`, `diag`.

# Example
```julia
using TSNA, DynamicNetworks
dnet = DynamicNetwork(3; directed=false, observation_start=0.0, observation_end=10.0)
activate!(dnet, 0.0, 10.0; edge=(1, 2))
activate!(dnet, 0.0, 10.0; edge=(2, 3))
ev = t_eigenvector(dnet, 2.0)
ev[2] > ev[1]                       # true: the broker is most central
```
"""
t_eigenvector(dnet::DynamicNetwork, at; active_only::Bool=true, missing::Symbol=:error,
              kwargs...) =
    _per_vertex(net -> _sna_eigenvector(net; missing, kwargs...),
                dnet, at; missing, active_only)

"""
    t_pagerank(dnet::DynamicNetwork, at; damping=0.85, active_only=true, missing=:error)
        -> Vector{Float64}

PageRank of each actor at time `at` (Graphs.jl's `pagerank`; not an sna
measure), computed on the active network (absent actors receive no teleport
mass) and scattered back with `NaN` for them.

# Example
```julia
using TSNA, DynamicNetworks
dnet = DynamicNetwork(4; observation_start=0.0, observation_end=10.0)
activate!(dnet, 0.0, 10.0; edge=(1, 2))
activate!(dnet, 0.0, 10.0; edge=(3, 2))
activate!(dnet, 6.0, 10.0; vertex=4)
pr = t_pagerank(dnet, 2.0)
isnan(pr[4]), sum(pr[1:3]) ≈ 1      # (true, true)
```
"""
t_pagerank(dnet::DynamicNetwork, at; damping::Float64=0.85, active_only::Bool=true,
           missing::Symbol=:error) =
    _per_vertex(net -> _sna_pagerank(net; α=damping, missing),
                dnet, at; missing, active_only)

"""
    t_density(dnet::DynamicNetwork, at; active_only=true, missing=:error, diag=false) -> Float64

Density at time `at` (`SNA.gden`, R's `sna::gden`), over the active vertices by default. Pass
`active_only=false` to use the complete vertex universe. The missing-data
policy is checked on the original network before inactive vertices are
removed.

# Example
```julia
using TSNA, DynamicNetworks
dnet = DynamicNetwork(4; observation_start=0.0, observation_end=10.0)
activate!(dnet, 0.0, 10.0; edge=(1, 2))
activate!(dnet, 6.0, 10.0; vertex=4)
t_density(dnet, 2.0)                # 1/6: three active actors
```
"""
function t_density(dnet::DynamicNetwork, at; active_only::Bool=true,
                   missing::Symbol=:error, diag::Bool=false)
    return _sna_density(_snapshot(dnet, at; missing, active_only); missing, diag)
end

"""
    t_reciprocity(dnet::DynamicNetwork, at; method=:dyadic, active_only=true, missing=:error)
        -> Float64

Dyadic reciprocity at time `at` (`SNA.grecip`, R's `grecip(measure = "dyadic")`), over the
active vertices by default. Pass `method=:edgewise` for the fraction of edges
that are reciprocated.

# Example
```julia
using TSNA, DynamicNetworks
dnet = DynamicNetwork(3; observation_start=0.0, observation_end=10.0)
activate!(dnet, 0.0, 10.0; edge=(1, 2))
activate!(dnet, 0.0, 10.0; edge=(2, 1))
activate!(dnet, 0.0, 10.0; edge=(2, 3))
t_reciprocity(dnet, 5.0; method=:edgewise)   # 2/3
```
"""
function t_reciprocity(dnet::DynamicNetwork, at; method::Symbol=:dyadic,
                       missing::Symbol=:error, active_only::Bool=true)
    return _sna_reciprocity(_snapshot(dnet, at; missing, active_only); method, missing)
end

"""
    t_transitivity(dnet::DynamicNetwork, at; measure=:weak, type=:global, cmode=:total,
                   active_only=true, missing=:error) -> Float64

Transitivity at time `at`, over the active vertices by default.

- `type=:global` (default): `SNA.gtrans` (R's `sna::gtrans`) with its
  `measure` — `:weak` (default: the share of two-paths that are closed),
  `:strong`, or the counts `:weakcensus`/`:strongcensus`.
- `type=:average`: the mean local clustering coefficient,
  `SNA.transitivity(net; type=:average, cmode)` (`cmode=:total` is Fagiolo's
  directed coefficient, `:weak` the coefficient of the symmetrised graph).

The per-vertex `type=:local` is not a scalar and is refused; take
`SNA.transitivity(network_extract(dnet, at); type=:local)` instead.

# Example
```julia
using TSNA, DynamicNetworks
dnet = DynamicNetwork(4; directed=false, observation_start=0.0, observation_end=10.0)
activate_edges!(dnet, [(1, 2), (2, 3), (1, 3), (3, 4)], 0.0, 10.0)
t_transitivity(dnet, 5.0)                    # 0.6: 6 of 10 two-paths closed
t_transitivity(dnet, 5.0; type=:average)     # (1 + 1 + 1/3) / 3 ≈ 0.778
```
"""
function t_transitivity(dnet::DynamicNetwork, at; measure::Symbol=:weak,
                        type::Symbol=:global, cmode::Symbol=:total,
                        active_only::Bool=true, missing::Symbol=:error)
    type in (:global, :average) || throw(ArgumentError(
        "type must be :global or :average (got :$type); the per-vertex :local " *
        "coefficients are SNA.transitivity(network_extract(dnet, at); type=:local)"))
    snap = _snapshot(dnet, at; missing, active_only)
    return type == :global ? Float64(_sna_transitivity(snap; measure, missing)) :
                             Float64(_sna_average_clustering(snap; cmode, missing))
end

# =============================================================================
# Duration and Turnover Metrics
# =============================================================================

"""
    t_edge_duration(dnet::DynamicNetwork; mode=:total, aggregate=:all,
                    active_default=true, missing=:error)

Edge activity durations, after `tsna::edgeDuration`. With the defaults this
**is** `tsna::edgeDuration(nd)`: a vector with the total active time of each
edge.

Each edge's spells are merged (adjacent or overlapping spells are one
lifetime). With an observation window (`get_observation_period`) they are
clipped to it, so an open-ended spell contributes the time up to the window
end, and edges with no spell overlapping the window are left out. **Without a
window** nothing is clipped, as in tsna: an open-ended spell, and an edge with
no spell record, last `Inf`. An edge with no spell record is active
throughout when `active_default=true` (R's default).

- `mode=:total` (tsna `subject = "edges"`): one entry per edge, the sum of its
  clipped spell durations; edges in the order of `edges(dnet.network)`.
- `mode=:spell` (tsna `subject = "spells"`): one entry per merged, clipped
  spell, by edge then onset.

`aggregate` is `:all` (the vector, default), `:mean`, `:median` or `:total`.
The durations are observed, not corrected for censoring: a spell cut by the
window end is right-censored, so the mean of a sample with open spells
underestimates the mean tie lifetime (no Kaplan–Meier estimator is offered).
On `Date`/`DateTime` axes durations are in seconds; an unbounded duration is
`Inf` on every axis.

# Example
```julia
using TSNA, DynamicNetworks
dnet = DynamicNetwork(3; observation_start=0.0, observation_end=10.0)
activate!(dnet, 0.0, 5.0; edge=(1, 2))
activate!(dnet, 5.0, 8.0; edge=(1, 2))      # adjacent: one lifetime [0, 8)
activate!(dnet, 2.0, Inf; edge=(2, 3))      # open-ended: clipped at 10
t_edge_duration(dnet)                       # [8.0, 8.0]
t_edge_duration(dnet; aggregate=:mean)      # 8.0
nowin = DynamicNetwork(3)                   # no window: nothing is clipped
activate!(nowin, 2.0, 6.0; edge=(1, 2))
activate!(nowin, 6.0, Inf; edge=(2, 3))
t_edge_duration(nowin)                      # [4.0, Inf], as tsna::edgeDuration
```
"""
function t_edge_duration(dnet::DynamicNetwork{T, Time};
                         mode::Symbol=:total, aggregate::Symbol=:all,
                         active_default::Bool=true,
                         missing::Symbol=:error) where {T, Time}
    require_observed(dnet.network, missing; context="t_edge_duration")
    mode in (:spell, :total) || throw(ArgumentError("mode must be :spell or :total"))
    lo, hi = _lifetime_window(dnet)
    durations = Float64[]
    for (_, spells) in _edge_spell_table(dnet; active_default)
        _push_durations!(durations, _clip(spells, lo, hi), mode)
    end
    return _aggregate(durations, aggregate)
end

"""
    t_vertex_duration(dnet::DynamicNetwork; mode=:total, aggregate=:all,
                      active_default=true, missing=:error)

Vertex activity durations, after `tsna::vertexDuration` (see
[`t_edge_duration`](@ref) for the merging, clipping and keyword conventions):
with the defaults, one entry per vertex active in the window, the total time
it is active there (vertices in ID order; a vertex with no spell record is
active throughout). `mode=:spell` gives one entry per merged, clipped spell.
Without a window nothing is clipped, and open spells last `Inf`, as in tsna.

With a window, finite spells are clipped to it as edge spells are; tsna's
`vertexDuration` replaces only *infinite* bounds by the window and keeps
finite spells that reach outside it, so the two differ for such vertices.

# Example
```julia
using TSNA, DynamicNetworks
dnet = DynamicNetwork(3; observation_start=0.0, observation_end=10.0)
activate!(dnet, 1.0, Inf; vertex=1)
deactivate!(dnet, -Inf, Inf; vertex=3)      # never present
t_vertex_duration(dnet)                     # [9.0, 10.0]
```
"""
function t_vertex_duration(dnet::DynamicNetwork{T, Time};
                           mode::Symbol=:total, aggregate::Symbol=:all,
                           active_default::Bool=true,
                           missing::Symbol=:error) where {T, Time}
    require_observed(dnet.network, missing; context="t_vertex_duration")
    mode in (:spell, :total) || throw(ArgumentError("mode must be :spell or :total"))
    lo, hi = _lifetime_window(dnet)
    durations = Float64[]
    for (_, spells) in _vertex_spell_table(dnet; active_default)
        _push_durations!(durations, _clip(spells, lo, hi), mode)
    end
    return _aggregate(durations, aggregate)
end

function _push_durations!(durations, clipped, mode)
    isempty(clipped) && return durations
    if mode == :spell
        for s in clipped
            push!(durations, _span(s))
        end
    else
        push!(durations, sum(_span, clipped; init=0.0))
    end
    return durations
end

"""
    tied_duration(dnet::DynamicNetwork; mode=:duration, neighborhood=:out,
                  active_default=true, missing=:error) -> Vector{Float64}

For each vertex, the total time its ties are active in the observation window
(`mode=:duration`) or the number of its tie spells there (`mode=:counts`):
`tsna::tiedDuration`. `neighborhood` is `:out` (ties the vertex sends, the
default), `:in` or `:combined`; undirected networks always use `:combined`.
Edge spells are merged and clipped to the window as in
[`t_edge_duration`](@ref); vertices with no ties score 0.

**Without a window**, tsna's `tiedDuration` clips to the range
`[first, last]` of the change times (`get_change_times`), so a tie
still active after the last change time counts only up to it, and a spell
that starts at the last change time does not count. TSNA does the same. A
network with no finite change time at all raises an `ArgumentError` (tsna
would use the placeholder range (0, 1)); set a window with
`set_observation_period!` to choose the range.

# Example
```julia
using TSNA, DynamicNetworks
dnet = DynamicNetwork(3; observation_start=0.0, observation_end=10.0)
activate!(dnet, 0.0, 4.0; edge=(1, 2))
activate!(dnet, 2.0, Inf; edge=(1, 3))
tied_duration(dnet)                            # [12.0, 0.0, 0.0]
tied_duration(dnet; neighborhood=:in)          # [0.0, 4.0, 8.0]
```
"""
function tied_duration(dnet::DynamicNetwork{T, Time}; mode::Symbol=:duration,
                       neighborhood::Symbol=:out, active_default::Bool=true,
                       missing::Symbol=:error) where {T, Time}
    require_observed(dnet.network, missing; context="tied_duration")
    mode in (:duration, :counts) || throw(ArgumentError("mode must be :duration or :counts"))
    neighborhood in (:out, :in, :combined) ||
        throw(ArgumentError("neighborhood must be :out, :in or :combined"))
    is_directed(dnet) || (neighborhood = :combined)
    lo, hi = _bounds(dnet, "tied_duration")
    out = zeros(Float64, Int(nv(dnet)))
    for ((i, j), spells) in _edge_spell_table(dnet; active_default)
        clipped = _clip(spells, lo, hi)
        isempty(clipped) && continue
        amount = mode == :counts ? Float64(length(clipped)) : sum(_span, clipped; init=0.0)
        neighborhood in (:out, :combined) && (out[Int(i)] += amount)
        neighborhood in (:in, :combined) && (out[Int(j)] += amount)
    end
    return out
end

"""
    t_edge_density(dnet::DynamicNetwork; mode=:duration, agg_unit=:edge,
                   active_default=true, missing=:error) -> Float64

`tsna::tEdgeDensity` over the observation window of length `L`:

- `mode=:duration, agg_unit=:edge` (default): the fraction of the time the
  edges of the base network are active, `Σ durations / (n_edges · L)`.
- `mode=:duration, agg_unit=:dyad`: the same total divided by the number of
  dyads times `L` — the time-averaged density. The dyads are those of the base
  network, as R's `network.dyadcount` counts them: `n(n−1)` directed and
  `n(n−1)/2` undirected, plus `n` when the network allows loops; `n₁·n₂`
  on a two-mode network (`2n₁n₂` if it is directed).
- `mode=:event, agg_unit=:edge`: the number of edge spells per edge per unit
  time, `n_spells / (n_edges · L)`. (tsna 0.3.6 computes
  `n_spells / (n_edges · end − start)`, which is this value only when the
  window starts at 0; the intended formula is used here.)

**Without a window**, `L` and the clipping range are those of tsna's
`tEdgeDensity`: the range `[first, last]` of the change times. A network with
no finite change time raises an `ArgumentError` (tsna would divide by the
placeholder range (0, 1)). A range of length zero gives `NaN`, as in tsna, and
an infinite window raises an `ArgumentError`.

Returns 0 for a network with no edges, as tsna does. `mode=:event` with
`agg_unit=:dyad` is not implemented in tsna either and is refused. Under
`missing=:face` the masked dyads count as observed (R's `network.dyadcount`
would subtract the missing edges).

# Example
```julia
using TSNA, DynamicNetworks
dnet = DynamicNetwork(3; observation_start=0.0, observation_end=10.0)
activate!(dnet, 0.0, 5.0; edge=(1, 2))
activate!(dnet, 0.0, 10.0; edge=(2, 3))
t_edge_density(dnet)                           # 0.75
t_edge_density(dnet; agg_unit=:dyad)           # 0.25
```
"""
function t_edge_density(dnet::DynamicNetwork{T, Time}; mode::Symbol=:duration,
                        agg_unit::Symbol=:edge, active_default::Bool=true,
                        missing::Symbol=:error) where {T, Time}
    require_observed(dnet.network, missing; context="t_edge_density")
    mode in (:duration, :event) || throw(ArgumentError("mode must be :duration or :event"))
    agg_unit in (:edge, :dyad) || throw(ArgumentError("agg_unit must be :edge or :dyad"))
    mode == :event && agg_unit == :dyad && throw(ArgumentError(
        "t_edge_density: mode=:event with agg_unit=:dyad is not implemented " *
        "(nor is it in tsna::tEdgeDensity)"))
    n, m = Int(nv(dnet)), Int(ne(dnet))
    (n == 0 || m == 0) && return 0.0
    lo, hi = _bounds(dnet, "t_edge_density")
    len = _span(Spell(lo, hi))
    isfinite(len) || throw(ArgumentError(
        "t_edge_density needs an observation range of finite length (got $lo to " *
        "$hi); set one with set_observation_period!(dnet, start, stop)"))
    len == 0 && return NaN          # a density over no time is undefined, as in tsna
    table = _edge_spell_table(dnet; active_default)
    if mode == :event
        return sum(length(sp) for (_, sp) in table; init=0) / (m * len)
    end
    total = sum(_span(s) for (_, sp) in table for s in _clip(sp, lo, hi); init=0.0)
    units = agg_unit == :edge ? m : _dyad_count(dnet.network)
    return total / (units * len)
end

# R's network.dyadcount (without its na.omit, since masked networks are
# refused unless the caller asked for face values).
function _dyad_count(net::Network)
    n = Int(nv(net))
    if is_two_mode(net)
        n1 = Int(net.bipartite)
        return (is_directed(net) ? 2 : 1) * n1 * (n - n1)
    end
    pairs = is_directed(net) ? n * (n - 1) : n * (n - 1) ÷ 2
    return net.loops ? pairs + n : pairs
end

_dur(d) = elapsed_seconds(d)
_dur(d::Real) = Float64(d)

function _aggregate(values::Vector{Float64}, aggregate::Symbol)
    aggregate in (:all, :mean, :median, :total) ||
        throw(ArgumentError("aggregate must be :mean, :median, :total, or :all"))
    aggregate == :all && return values
    isempty(values) && return NaN
    aggregate == :mean && return mean(values)
    aggregate == :median && return median(values)
    return sum(values)
end

# Formation (onset) and dissolution (terminus) events of the merged, clipped
# edge spells. A censored bound is not an event unless include_censored=true,
# in which case it is counted at its truncated position (tsna's
# include.censored). Without a window the spells are not clipped, and only an
# infinite bound is censored (it then sits at -Inf/Inf, outside any query).
function _edge_events(dnet::DynamicNetwork{T, Time}; include_censored::Bool,
                      active_default::Bool=true) where {T, Time}
    lo, hi = _lifetime_window(dnet)
    onsets, termini = Time[], Time[]
    for (_, spells) in _edge_spell_table(dnet; active_default)
        for s in _clip(spells, lo, hi)
            (include_censored || !s.onset_censored) && push!(onsets, s.onset)
            (include_censored || !s.terminus_censored) && push!(termini, s.terminus)
        end
    end
    return sort!(onsets), sort!(termini)
end

_events_in(v, lo, hi) = searchsortedfirst(v, hi) - searchsortedfirst(v, lo)

"""
    t_edge_formation(dnet::DynamicNetwork, onset, terminus; include_censored=false,
                     missing=:error) -> Int

The number of edge **formation events** — onsets of merged edge spells — in
`[onset, terminus)`. Spells are merged per edge first (a tie active over two
adjacent waves forms once) and clipped to the observation window, if one is
set; an onset that is left-censored (stored `onset_censored`, before the
window start, or `-Inf`) is not a formation unless `include_censored=true`,
as in tsna's `include.censored = FALSE`. Without a window, as in tsna, only an
infinite onset is censored, and with `include_censored=true` it sits at
`-Inf`, outside any finite query.

`tsna::tEdgeFormation(nd, start, end, time.interval)` returns the series of
counts at each `t` in `seq(start, end, time.interval)` of onsets equal to `t`;
on integer-valued spells that is
`[t_edge_formation(dnet, t, t + 1) for t in start:end]`. Its default range,
the closed range of the change times, is what [`t_turnover`](@ref) covers on a
network without a window.

# Example
```julia
using TSNA, DynamicNetworks
dnet = DynamicNetwork(3; observation_start=0.0, observation_end=20.0)
activate!(dnet, 0.0, 5.0; edge=(1, 2))
activate!(dnet, 5.0, 9.0; edge=(1, 2))     # persists: no new formation
activate!(dnet, 12.0, 15.0; edge=(1, 2))   # re-forms
add_spell!(dnet, Spell(0.0, 3.0; onset_censored=true); edge=(2, 3))
t_edge_formation(dnet, 0.0, 20.0)                          # 2
t_edge_formation(dnet, 0.0, 20.0; include_censored=true)   # 3
```
"""
function t_edge_formation(dnet::DynamicNetwork{T, Time}, onset, terminus;
                          include_censored::Bool=false,
                          missing::Symbol=:error) where {T, Time}
    require_observed(dnet.network, missing; context="t_edge_formation")
    onset, terminus = convert(Time, onset), convert(Time, terminus)
    onset <= terminus || throw(ArgumentError("onset must be <= terminus"))
    onsets, _ = _edge_events(dnet; include_censored)
    return _events_in(onsets, onset, terminus)
end

"""
    t_edge_dissolution(dnet::DynamicNetwork, onset, terminus; include_censored=false,
                       missing=:error) -> Int

The number of edge **dissolution events** — termini of merged edge spells —
in `[onset, terminus)`, with the conventions of [`t_edge_formation`](@ref): a
right-censored terminus (stored `terminus_censored`, after the window end, or
`Inf`) is not a dissolution unless `include_censored=true`.

# Example
```julia
using TSNA, DynamicNetworks
dnet = DynamicNetwork(3; observation_start=0.0, observation_end=20.0)
activate!(dnet, 0.0, 5.0; edge=(1, 2))
activate!(dnet, 5.0, 9.0; edge=(1, 2))     # one lifetime [0, 9)
activate!(dnet, 4.0, Inf; edge=(2, 3))     # still active at the window end
t_edge_dissolution(dnet, 0.0, 20.0)        # 1
```
"""
function t_edge_dissolution(dnet::DynamicNetwork{T, Time}, onset, terminus;
                            include_censored::Bool=false,
                            missing::Symbol=:error) where {T, Time}
    require_observed(dnet.network, missing; context="t_edge_dissolution")
    onset, terminus = convert(Time, onset), convert(Time, terminus)
    onset <= terminus || throw(ArgumentError("onset must be <= terminus"))
    _, termini = _edge_events(dnet; include_censored)
    return _events_in(termini, onset, terminus)
end

"""
    t_edge_persistence(dnet::DynamicNetwork, window_size; missing=:error) -> Float64

Proportion of edges that persist (survive) across consecutive time
windows, after the temporal-correlation measure of Nicosia et al.: the
observation period is divided into windows of length `window_size`; for
each consecutive pair of windows the edges active at the start of the
first window are checked for activity at the start of the second, and the
pooled proportion `persisted / total` is returned.

Values near 1 indicate a stable network (little edge turnover); values
near 0 indicate almost complete edge replacement per window. Returns
`NaN` when there are fewer than two windows or no active edges to track.
Not a tsna function. Without an observation window the windows start at
every `window_size` step of the closed range `[first, last]` of the change
times (see [`t_turnover`](@ref)).

# Example
```julia
using TSNA, DynamicNetworks
dnet = DynamicNetwork(3; observation_start=0.0, observation_end=30.0)
activate!(dnet, 0.0, 30.0; edge=(1, 2))
activate!(dnet, 0.0, 15.0; edge=(2, 3))
t_edge_persistence(dnet, 10.0)      # 0.75: 2/2 then 1/2 survive
```
"""
function t_edge_persistence(dnet::DynamicNetwork{T, Time}, window_size;
                            missing::Symbol=:error) where {T, Time}
    require_observed(dnet.network, missing; context="t_edge_persistence")
    starts = _window_starts(dnet, window_size, "t_edge_persistence")
    length(starts) >= 2 || return NaN

    total = 0
    persisted = 0
    prev = Set{Tuple{T, T}}(active_edges(dnet, starts[1]))
    for k in 2:length(starts)
        cur = Set{Tuple{T, T}}(active_edges(dnet, starts[k]))
        total += length(prev)
        persisted += count(in(cur), prev)
        prev = cur
    end

    return total == 0 ? NaN : persisted / total
end

"""
    t_turnover(dnet::DynamicNetwork, window_size; include_censored=false, missing=:error)
        -> Vector{NamedTuple}

Formation and dissolution **event counts and rates** per window of length
`window_size` across the observation period, with the event conventions of
[`t_edge_formation`](@ref)/[`t_edge_dissolution`](@ref) (spells merged per
edge, censored bounds not counted). Every element has the same shape:
`(window_start, window_end, n_formations, n_dissolutions, formation_rate,
dissolution_rate)`, with rates per unit time (per second on calendar axes).

With an observation window the windows tile `[start, end)`, the last one cut
at the window end. **Without one** they start at `first, first + window_size,
…` up to and including the last change time — the closed range of
`get_change_times` that `tsna::tEdgeFormation(nd)` and
`tEdgeDissolution(nd)` cover by default — so events at the last change time
are counted. Without a window, `window_size = 1` on integer-valued spells
gives tsna's default series. (With a window the counts per window agree with
tsna's, but tsna's default grid is still the range of the change times, not
the window.) A network with neither a window nor a finite change time raises
an `ArgumentError`.

# Example
```julia
using TSNA, DynamicNetworks
dnet = DynamicNetwork(3; observation_start=0.0, observation_end=20.0)
activate!(dnet, 2.0, 12.0; edge=(1, 2))
activate!(dnet, 5.0, 8.0; edge=(2, 3))
w = t_turnover(dnet, 10.0)
[(x.n_formations, x.n_dissolutions) for x in w]   # [(2, 1), (0, 1)]
```
"""
function t_turnover(dnet::DynamicNetwork{T, Time}, window_size;
                    include_censored::Bool=false,
                    missing::Symbol=:error) where {T, Time}
    require_observed(dnet.network, missing; context="t_turnover")
    start_time, end_time, closed = _series_range(dnet, "t_turnover")
    _check_window(start_time, end_time, window_size)

    # Event times collected and sorted once, so each window is a pair of
    # binary searches instead of a full scan over every spell
    onsets, termini = _edge_events(dnet; include_censored)

    row(t, w_end, nf, nd, span) = (window_start=t, window_end=w_end,
                                   n_formations=nf, n_dissolutions=nd,
                                   formation_rate=span > 0 ? nf / span : NaN,
                                   dissolution_rate=span > 0 ? nd / span : NaN)
    results = typeof(row(start_time, start_time, 0, 0, 0.0))[]
    t = start_time
    while t < end_time || (closed && t == end_time)
        next = _next_window(t, window_size)
        w_end = closed ? next : min(next, end_time)
        push!(results, row(t, w_end, _events_in(onsets, t, w_end),
                           _events_in(termini, t, w_end), _dur(w_end - t)))
        t = next
    end

    return results
end

"""
    tie_decay(dnet::DynamicNetwork; method=:exponential, rate=1.0, at=nothing,
              missing=:error) -> Dict

Tie weights decayed by the time since each edge was last active at `at` (by
default the window end, or the last change time when no window is set):
`exp(-rate·Δ)` (`:exponential`) or `max(0, 1 - rate·Δ)` (`:linear`),
where Δ is the time from the end of the edge's most recent spell to `at`
(0 for currently active ties, including edges with no spell record, which
are active by default). Edges whose activity lies entirely after `at` get no
weight. Not a tsna function.

# Example
```julia
using TSNA, DynamicNetworks
dnet = DynamicNetwork(3; observation_start=0.0, observation_end=10.0)
activate!(dnet, 0.0, 10.0; edge=(1, 2))
activate!(dnet, 0.0, 5.0; edge=(2, 3))
w = tie_decay(dnet; rate=0.1)
w[(1, 2)], round(w[(2, 3)]; digits=3)   # (1.0, 0.607)
```
"""
function tie_decay(dnet::DynamicNetwork{T, Time};
                   method::Symbol=:exponential, rate::Float64=1.0,
                   at=nothing, missing::Symbol=:error) where {T, Time}
    require_observed(dnet.network, missing; context="tie_decay")
    # Default instant: the window end, else the last change time.
    at = isnothing(at) ? _bounds(dnet, "tie_decay without `at`")[2] : at
    method in (:exponential, :linear) || throw(ArgumentError("method must be :exponential or :linear"))
    isfinite(rate) && rate >= 0 || throw(ArgumentError("rate must be finite and nonnegative"))
    at = convert(Time, at)
    weights = Dict{Tuple{T, T}, Float64}()

    for (edge, spells) in _edge_spell_table(dnet)
        isempty(spells) && continue
        # Time since last activity (0 if active at `at`)
        Δ = Inf
        for s in spells
            if spell_active_at(s, at)
                Δ = 0.0
                break
            end
            s.terminus <= at && (Δ = min(Δ, _dur(at - s.terminus)))
        end
        isinf(Δ) && continue  # only spells entirely after `at`

        weights[edge] = method == :exponential ? exp(-rate * Δ) : max(0.0, 1.0 - rate * Δ)
    end

    return weights
end

# =============================================================================
# Aggregation and Time Series
# =============================================================================

const _SNA_STATS_MEASURES = (:density, :reciprocity, :transitivity, :mean_degree, :n_edges)

"""
    t_sna_stats(dnet::DynamicNetwork, times; measures=[:density, :reciprocity,
                :transitivity, :mean_degree], active_only=true, missing=:error,
                reciprocity_method=:dyadic) -> Vector{NamedTuple}

Network-level statistics of the snapshot at each time point (one row per
time, fields `time` and the requested measures, in alphabetical order).
Each snapshot is extracted once, over the active vertices by default, and
reused for all measures; `:mean_degree` counts each edge appropriately for
the network's directedness and `:n_edges` is the edge count. This is not
tsna's `tSnaStats` (which applies one sna function over a regular time grid);
see the R concordance in the documentation.

# Example
```julia
using TSNA, DynamicNetworks
dnet = DynamicNetwork(3; observation_start=0.0, observation_end=10.0)
activate!(dnet, 0.0, 10.0; edge=(1, 2))
activate!(dnet, 5.0, 10.0; edge=(2, 3))
rows = t_sna_stats(dnet, [1.0, 6.0]; measures=[:n_edges, :density])
[r.n_edges for r in rows]           # [1.0, 2.0]
```
"""
function t_sna_stats(dnet::DynamicNetwork{T, Time}, times::AbstractVector;
                     measures::Vector{Symbol}=[:density, :reciprocity,
                                               :transitivity, :mean_degree],
                     missing::Symbol=:error, active_only::Bool=true,
                     reciprocity_method::Symbol=:dyadic) where {T, Time}
    require_observed(dnet.network, missing; context="t_sna_stats")
    for m in measures
        m in _SNA_STATS_MEASURES || throw(ArgumentError(
            "unknown measure: $m (available: $(join(_SNA_STATS_MEASURES, ", ")))"))
    end
    keys_sorted = Tuple(sort(unique([:time; measures])))

    function row(t)
        snapshot = _snapshot(dnet, t; missing, active_only)
        n = Int(nv(snapshot))
        value(m) = if m == :time
            convert(Time, t)
        elseif m == :density
            Float64(_sna_density(snapshot; missing, diag=false))
        elseif m == :reciprocity
            Float64(_sna_reciprocity(snapshot; method=reciprocity_method, missing))
        elseif m == :transitivity
            Float64(_sna_transitivity(snapshot; missing))
        elseif m == :mean_degree
            n == 0 ? NaN : (is_directed(snapshot) ? ne(snapshot) / n : 2 * ne(snapshot) / n)
        else # :n_edges
            Float64(ne(snapshot))
        end
        return NamedTuple{keys_sorted}(map(value, keys_sorted))
    end
    return [row(t) for t in times]
end

"""
    window_sna_stats(dnet::DynamicNetwork, window_size; kwargs...) -> Vector{NamedTuple}

[`t_sna_stats`](@ref) sampled at the start of each window of length
`window_size` across the observation period; without a window, at every
`window_size` step of the closed range `[first, last]` of the change times
(`tsna::tSnaStats`'s default `seq(start, end)`).

# Example
```julia
using TSNA, DynamicNetworks
dnet = DynamicNetwork(3; observation_start=0.0, observation_end=20.0)
activate!(dnet, 0.0, 20.0; edge=(1, 2))
activate!(dnet, 10.0, 20.0; edge=(2, 3))
[r.n_edges for r in window_sna_stats(dnet, 10.0; measures=[:n_edges])]   # [1.0, 2.0]
```
"""
function window_sna_stats(dnet::DynamicNetwork{T, Time}, window_size;
                          missing::Symbol=:error, kwargs...) where {T, Time}
    require_observed(dnet.network, missing; context="window_sna_stats")
    times = _window_starts(dnet, window_size, "window_sna_stats")
    return t_sna_stats(dnet, times; missing, kwargs...)
end

# Window starts of a series: every `window_size` step of [start, end) of the
# observation window, or of the closed range [first, last] of the change times
# when there is no window (tsna's seq(start, end, time.interval)).
function _window_starts(dnet::DynamicNetwork{T, Time}, window_size,
                        context::AbstractString) where {T, Time}
    start_time, end_time, closed = _series_range(dnet, context)
    _check_window(start_time, end_time, window_size)
    starts = Time[]
    t = start_time
    while t < end_time || (closed && t == end_time)
        push!(starts, t)
        t = _next_window(t, window_size)
    end
    return starts
end

function _next_window(t, window_size)
    next = convert(typeof(t), t + window_size)
    next > t || throw(ArgumentError("window_size must advance the time axis"))
    return next
end

function _check_window(start_time::Time, end_time::Time, window_size) where Time
    # Unbounded ends are ±Inf, or the axis extremes on integer and calendar
    # axes (DynamicNetworks.unbounded_spell); a window reaching typemax(Int)
    # would otherwise be tiled step by step, practically forever.
    u = DynamicNetworks.unbounded_spell(Time)
    (start_time == u.onset || end_time == u.terminus ||
     (start_time isa Real && (!isfinite(start_time) || !isfinite(end_time)))) &&
        throw(ArgumentError("window statistics require a finite observation period " *
                            "(got $start_time to $end_time)"))
    _next_window(start_time, window_size)
    return nothing
end

"""
    t_aggregate(dnet::DynamicNetwork; method=:union, onset=window start,
                terminus=window end, report=false) -> Network

Collapse the dynamic network into a static one over `[onset, terminus)`, by
default the observation window. Without a window, and with neither `onset`
nor `terminus` given, the whole time axis is collapsed (R's
`network.collapse` default), and `:weighted` gives an open spell the weight
`Inf`.

- `:union` — edges active at some point in the window
- `:intersection` — edges active throughout the window
- `:weighted` — union, with each edge's total active time in the window (its
  spells merged and clipped) stored as its `:weight` attribute

This is `DynamicNetworks.network_collapse` (R's `network.collapse`, which
needs the endpoints to be active too) with an aggregation rule, and it
inherits its conversion invariants: vertex IDs are stable, so directedness,
`loops`, two-mode metadata, static attributes and the **missing-dyad mask** all
survive; spells, time-varying attributes and the observation window are dropped
by nature. Pass `report=true` for `(net, ::NetworkCore.ConversionReport)`.

# Example
```julia
using TSNA, DynamicNetworks, NetworkCore
dnet = DynamicNetwork(3; observation_start=0.0, observation_end=10.0)
activate!(dnet, 0.0, 10.0; edge=(1, 2))
activate!(dnet, 2.0, 6.0; edge=(2, 3))
ne(t_aggregate(dnet)), ne(t_aggregate(dnet; method=:intersection))   # (2, 1)
get_edge_attribute(t_aggregate(dnet; method=:weighted), :weight, 2, 3)   # 4.0
```
"""
function t_aggregate(dnet::DynamicNetwork{T, Time}; method::Symbol=:union,
                     onset=nothing, terminus=nothing,
                     report::Bool=false) where {T, Time}
    window = get_observation_period(dnet)
    if isnothing(onset) || isnothing(terminus)
        if !isnothing(window)
            onset = something(onset, window[1])
            terminus = something(terminus, window[2])
        elseif isnothing(onset) && isnothing(terminus)
            onset, terminus = _axis(dnet)
        else
            throw(ArgumentError(
                "t_aggregate: give both onset and terminus, or neither, on a network " *
                "with no observation window"))
        end
    end
    onset, terminus = convert(Time, onset), convert(Time, terminus)

    method in (:union, :intersection, :weighted) ||
        throw(ArgumentError("method must be :union, :intersection, or :weighted"))

    rule = method == :intersection ? :all : :any
    collapsed, rep = network_collapse(dnet; onset=onset, terminus=terminus,
                                      rule=rule, report=true)

    if method == :weighted
        for (edge, spells) in _edge_spell_table(dnet)
            total = sum(_span, _clip(spells, onset, terminus); init=0.0)
            if total > 0 && has_edge(collapsed, edge[1], edge[2])
                set_edge_attribute!(collapsed, :weight, edge[1], edge[2], total)
            end
        end
    end

    return report ? (collapsed, rep) : collapsed
end

# Advertise accepted policies rather than inferring them from keyword presence.
for f in (earliest_arrival, earliest_arrival!, earliest_arrival_all,
          temporal_distance, temporal_distance_matrix, reachability_matrix,
          forward_reachable_set, backward_reachable_set, temporal_path,
          t_degree, t_betweenness, t_closeness, t_eigenvector, t_pagerank,
          t_density, t_reciprocity, t_transitivity, t_edge_duration,
          t_vertex_duration, t_edge_formation, t_edge_dissolution,
          t_edge_persistence, t_turnover, tie_decay, tied_duration,
          t_edge_density, t_sna_stats,
          window_sna_stats, as_contact_sequence)
    @eval NetworkCore.missing_policies(::typeof($f)) = (:error, :face)
end

# Time-to-first-result: compile the README path (snapshot measures, paths,
# durations, events, time series) for both directedness flavours.
@setup_workload begin
    @compile_workload begin
        for directed in (true, false)
            d = DynamicNetwork(5; directed=directed, observation_start=0.0,
                               observation_end=10.0)
            activate!(d, 0.0, 5.0; edge=(1, 2))
            activate!(d, 2.0, 8.0; edge=(2, 3))
            activate!(d, 4.0, 4.0; edge=(3, 4))
            activate!(d, 6.0, 10.0; vertex=5)
            t_degree(d, 3.0); t_betweenness(d, 3.0); t_closeness(d, 3.0)
            t_eigenvector(d, 3.0); t_pagerank(d, 3.0)
            t_density(d, 3.0); t_reciprocity(d, 3.0); t_transitivity(d, 3.0)
            forward_reachable_set(d, 1, 0.0); backward_reachable_set(d, 3, 10.0)
            temporal_path(d, 1, 3, 0.0); temporal_distance_matrix(d, 0.0)
            t_edge_duration(d); t_vertex_duration(d)
            t_edge_formation(d, 0.0, 10.0); t_edge_dissolution(d, 0.0, 10.0)
            t_turnover(d, 2.0); t_edge_persistence(d, 2.0); tie_decay(d)
            tied_duration(d); t_edge_density(d)
            t_sna_stats(d, [1.0, 4.0]); window_sna_stats(d, 5.0)
            t_aggregate(d; method=:weighted); as_contact_sequence(d)
        end
    end
end

end # module
