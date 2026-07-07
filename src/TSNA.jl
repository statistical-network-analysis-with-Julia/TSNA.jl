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

Port of the R tsna package from the StatNet collection.
"""
module TSNA

using Dates
using Graphs
using Network
using NetworkDynamic
using SNA
using Statistics

# Temporal paths and reachability
export tPath, path_duration
export temporalDistance, earliestArrival
export forwardReachableSet, backwardReachableSet
export temporalPath, shortestTemporalPath

# Point-in-time measures
export tDegree, tBetweenness, tCloseness
export tEigenvector, tPagerank
export tDensity, tReciprocity, tTransitivity

# Duration and turnover
export tEdgeDuration, tVertexDuration
export tEdgeFormation, tEdgeDissolution, tTurnover
export tieDecay

# Contact sequences
export Contact, ContactSequence, as_contact_sequence

# Aggregation and time series
export tSnaStats, windowSnaStats, tAggregate

# =============================================================================
# Temporal Path Types
# =============================================================================

"""
    tPath{T, Time}

A time-respecting path through a dynamic network: `times[k]` is the
instant edge `edges[k]` is traversed, and times are non-decreasing.
"""
struct tPath{T, Time}
    vertices::Vector{T}
    times::Vector{Time}
    edges::Vector{Tuple{T, T}}

    function tPath{T, Time}(vertices::Vector{T}, times::Vector{Time},
                            edges::Vector{Tuple{T, T}}) where {T, Time}
        length(times) == length(edges) ||
            throw(ArgumentError("times and edges must have same length"))
        length(vertices) == length(edges) + 1 ||
            throw(ArgumentError("vertices must have length edges + 1"))
        new{T, Time}(vertices, times, edges)
    end
end

tPath(vertices::Vector{T}, times::Vector{Time}, edges::Vector{Tuple{T,T}}) where {T, Time} =
    tPath{T, Time}(vertices, times, edges)

Base.length(p::tPath) = length(p.edges)

function Base.show(io::IO, p::tPath)
    print(io, "tPath: ")
    for (i, v) in enumerate(p.vertices)
        print(io, v)
        if i <= length(p.times)
            print(io, " --($(p.times[i]))--> ")
        end
    end
end

"""
    path_duration(p::tPath) -> Time difference

Elapsed time between the first and last traversal of the path.
"""
function path_duration(p::tPath{T, Time}) where {T, Time}
    isempty(p.times) && return _zero_duration(Time)
    return p.times[end] - p.times[1]
end

_zero_duration(::Type{Time}) where Time<:Number = zero(Time)
_zero_duration(::Type{DateTime}) = Millisecond(0)
_zero_duration(::Type{Date}) = Day(0)

# =============================================================================
# Contact Sequence
# =============================================================================

"""
    Contact{T, Time}

A single contact (edge activation spell) in a temporal network.
"""
struct Contact{T, Time}
    source::T
    target::T
    time::Time
    duration
end

"""
    ContactSequence{T, Time}

A sequence of contacts ordered by onset time.
"""
struct ContactSequence{T, Time}
    contacts::Vector{Contact{T, Time}}
    n_vertices::Int
    directed::Bool

    function ContactSequence(contacts::Vector{Contact{T, Time}}, n::Int;
                             directed::Bool=true) where {T, Time}
        sorted = sort(contacts, by=c -> c.time)
        new{T, Time}(sorted, n, directed)
    end
end

Base.length(cs::ContactSequence) = length(cs.contacts)
Base.iterate(cs::ContactSequence, state=1) =
    state > length(cs) ? nothing : (cs.contacts[state], state + 1)

"""
    as_contact_sequence(dnet::DynamicNetwork) -> ContactSequence

Convert a dynamic network's edge spells to a contact sequence.
"""
function as_contact_sequence(dnet::DynamicNetwork{T, Time}) where {T, Time}
    contacts = Contact{T, Time}[]

    for ((i, j), spells) in dnet.edge_spells
        for spell in spells
            push!(contacts, Contact{T, Time}(i, j, spell.onset,
                                             spell.terminus - spell.onset))
        end
    end

    return ContactSequence(contacts, Int(nv(dnet)); directed=is_directed(dnet))
end

# =============================================================================
# Temporal Path Finding (interval semantics)
# =============================================================================

# Per-vertex outgoing contacts: v -> [(neighbor, spell), ...]. For
# undirected networks each spell is listed from both endpoints.
function _out_contacts(dnet::DynamicNetwork{T, Time}) where {T, Time}
    out = Dict{T, Vector{Tuple{T, Spell{Time}}}}()
    directed = is_directed(dnet)
    for ((i, j), spells) in dnet.edge_spells
        for spell in spells
            push!(get!(out, i, Tuple{T, Spell{Time}}[]), (j, spell))
            if !directed
                push!(get!(out, j, Tuple{T, Spell{Time}}[]), (i, spell))
            end
        end
    end
    return out
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
    earliestArrival(dnet, source, start_time; end_time=observation end)
        -> (arrival::Dict, parent::Dict)

Earliest-arrival times from `source` to every vertex, starting at
`start_time`, under interval semantics: an edge spell `[onset, terminus)`
is traversable at any instant in it (boarding mid-spell is allowed;
spells that began before `start_time` but are still active count).
Implemented as a Dijkstra-style label-setting search, so chains of
simultaneous (equal-onset) spells are handled correctly.

Returns the arrival-time dictionary (vertices absent = unreachable) and
the parent map `(vertex => (predecessor, boarding time))` for path
reconstruction.
"""
function earliestArrival(dnet::DynamicNetwork{T, Time}, source::T, start_time;
                         end_time=dnet.observation_period[2]) where {T, Time}
    start_time = convert(Time, start_time)
    end_time = convert(Time, end_time)
    out = _out_contacts(dnet)

    arrival = Dict{T, Time}(source => start_time)
    parent = Dict{T, Tuple{T, Time}}()
    settled = Set{T}()

    # Label-setting search (Dijkstra with linear min-scan; boarding times
    # never precede the label being settled, so labels are final)
    while true
        v = nothing
        best = nothing
        for (u, t) in arrival
            u in settled && continue
            if isnothing(best) || t < best
                v, best = u, t
            end
        end
        isnothing(v) && break
        push!(settled, v)
        t = arrival[v]

        for (w, spell) in get(out, v, Tuple{T, Spell{Time}}[])
            w in settled && continue
            depart = _board_time(spell, t, end_time)
            isnothing(depart) && continue
            if !haskey(arrival, w) || depart < arrival[w]
                arrival[w] = depart
                parent[w] = (v, depart)
            end
        end
    end

    return arrival, parent
end

"""
    temporalDistance(dnet, source, target, start_time; end_time=...)
        -> elapsed time or nothing

Elapsed time of the earliest time-respecting path from `source` to
`target` departing at `start_time` (arrival − start). Returns `nothing`
when no path exists within the window.
"""
function temporalDistance(dnet::DynamicNetwork{T, Time}, source::T, target::T,
                          start_time; end_time=dnet.observation_period[2]) where {T, Time}
    arrival, _ = earliestArrival(dnet, source, start_time; end_time=end_time)
    haskey(arrival, target) || return nothing
    return arrival[target] - convert(Time, start_time)
end

"""
    forwardReachableSet(dnet, source, start_time; end_time=...) -> Vector

Vertices reachable from `source` by a time-respecting path departing at
or after `start_time` and arriving before `end_time` (`source` included).
"""
function forwardReachableSet(dnet::DynamicNetwork{T, Time}, source::T, start_time;
                             end_time=dnet.observation_period[2]) where {T, Time}
    arrival, _ = earliestArrival(dnet, source, start_time; end_time=end_time)
    return sort(collect(keys(arrival)))
end

"""
    backwardReachableSet(dnet, target, end_time; start_time=...) -> Vector

Vertices from which `target` can be reached by a time-respecting path in
`[start_time, end_time)` — the exact dual of
[`forwardReachableSet`](@ref) (computed by forward searches, so both use
identical traversal semantics).
"""
function backwardReachableSet(dnet::DynamicNetwork{T, Time}, target::T, end_time;
                              start_time=dnet.observation_period[1]) where {T, Time}
    reachable = T[]
    for v in 1:nv(dnet)
        vT = T(v)
        if vT == target
            push!(reachable, vT)
            continue
        end
        arrival, _ = earliestArrival(dnet, vT, start_time; end_time=end_time)
        haskey(arrival, target) && push!(reachable, vT)
    end
    return reachable
end

"""
    temporalPath(dnet, source, target, start_time; end_time=...)
        -> Union{tPath, Nothing}

The **earliest-arrival** time-respecting path from `source` to `target`
departing at `start_time` (not necessarily the fewest-hops path).
Returns `nothing` when no path exists.
"""
function temporalPath(dnet::DynamicNetwork{T, Time}, source::T, target::T,
                      start_time; end_time=dnet.observation_period[2]) where {T, Time}
    arrival, parent = earliestArrival(dnet, source, start_time; end_time=end_time)
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

    return tPath(verts, times, path_edges)
end

"""
    shortestTemporalPath(dnet, source, target, start_time; kwargs...)

Alias for [`temporalPath`](@ref): the *earliest-arrival* path (kept for
API compatibility with earlier versions; note this is not the fewest-hops
path).
"""
const shortestTemporalPath = temporalPath

# =============================================================================
# Temporal Measures at a Point
# =============================================================================
#
# Snapshots retain every vertex (retain_all_vertices=true), so all
# per-vertex vectors are indexed by the dynamic network's own vertex IDs
# and always have length nv(dnet) — even when some vertices are inactive.

_snapshot(dnet, at) = network_extract(dnet, at; retain_all_vertices=true)

"""
    tDegree(dnet::DynamicNetwork, at; mode=:total) -> Vector{Float64}

Degree at time `at`, indexed by the network's vertex IDs (inactive
vertices score 0).
"""
function tDegree(dnet::DynamicNetwork{T, Time}, at; mode::Symbol=:total) where {T, Time}
    return SNA.degree_centrality(_snapshot(dnet, at); mode=mode)
end

"""
    tBetweenness(dnet::DynamicNetwork, at; normalized=false) -> Vector{Float64}

Betweenness centrality at time `at` (raw scores by default, as in SNA.jl).
"""
tBetweenness(dnet::DynamicNetwork, at; normalized::Bool=false) =
    SNA.betweenness_centrality(_snapshot(dnet, at); normalized=normalized)

"""
    tCloseness(dnet::DynamicNetwork, at) -> Vector{Float64}
"""
tCloseness(dnet::DynamicNetwork, at) = SNA.closeness_centrality(_snapshot(dnet, at))

"""
    tEigenvector(dnet::DynamicNetwork, at) -> Vector{Float64}
"""
tEigenvector(dnet::DynamicNetwork, at) = SNA.eigenvector_centrality(_snapshot(dnet, at))

"""
    tPagerank(dnet::DynamicNetwork, at; damping=0.85) -> Vector{Float64}
"""
tPagerank(dnet::DynamicNetwork, at; damping::Float64=0.85) =
    SNA.pagerank(_snapshot(dnet, at); α=damping)

"""
    tDensity(dnet::DynamicNetwork, at) -> Float64

Density at time `at`, over all `nv(dnet)` vertices (not just the active
ones).
"""
function tDensity(dnet::DynamicNetwork, at)
    return network_density(_snapshot(dnet, at))
end

"""
    tReciprocity(dnet::DynamicNetwork, at) -> Float64

Edgewise reciprocity at time `at` (fraction of edges that are
reciprocated).
"""
function tReciprocity(dnet::DynamicNetwork, at)
    return _snapshot_reciprocity(_snapshot(dnet, at))
end

function _snapshot_reciprocity(snapshot)
    !is_directed(snapshot) && return 1.0
    n_edges = ne(snapshot)
    n_edges == 0 && return 0.0

    mutual = 0
    for e in edges(snapshot)
        has_edge(snapshot, dst(e), src(e)) && (mutual += 1)
    end
    return mutual / n_edges
end

"""
    tTransitivity(dnet::DynamicNetwork, at) -> Float64

Weak transitivity at time `at` (see `SNA.transitivity`).
"""
tTransitivity(dnet::DynamicNetwork, at) = SNA.transitivity(_snapshot(dnet, at))

# =============================================================================
# Duration and Turnover Metrics
# =============================================================================

"""
    tEdgeDuration(dnet::DynamicNetwork; mode=:spell, aggregate=:mean)

Edge activity durations.

- `mode=:spell` (default, matching `tsna::edgeDuration`): each activity
  spell contributes its own duration `terminus − onset`.
- `mode=:total`: spells of the same edge are summed first (total active
  time per edge).

`aggregate` is `:mean`, `:median`, `:total`, or `:all` (raw vector).
Censoring flags are ignored (durations are observed, not corrected).
"""
function tEdgeDuration(dnet::DynamicNetwork{T, Time};
                       mode::Symbol=:spell, aggregate::Symbol=:mean) where {T, Time}
    durations = Float64[]

    if mode == :spell
        for spells in values(dnet.edge_spells), s in spells
            push!(durations, _dur(s.terminus - s.onset))
        end
    elseif mode == :total
        for spells in values(dnet.edge_spells)
            push!(durations, sum(_dur(s.terminus - s.onset) for s in spells; init=0.0))
        end
    else
        throw(ArgumentError("mode must be :spell or :total"))
    end

    return _aggregate(durations, aggregate)
end

"""
    tVertexDuration(dnet::DynamicNetwork; mode=:spell, aggregate=:mean)

Vertex activity durations (see [`tEdgeDuration`](@ref)).
"""
function tVertexDuration(dnet::DynamicNetwork{T, Time};
                         mode::Symbol=:spell, aggregate::Symbol=:mean) where {T, Time}
    durations = Float64[]

    if mode == :spell
        for spells in values(dnet.vertex_spells), s in spells
            push!(durations, _dur(s.terminus - s.onset))
        end
    elseif mode == :total
        for spells in values(dnet.vertex_spells)
            push!(durations, sum(_dur(s.terminus - s.onset) for s in spells; init=0.0))
        end
    else
        throw(ArgumentError("mode must be :spell or :total"))
    end

    return _aggregate(durations, aggregate)
end

_dur(d) = Float64(NetworkDynamic._elapsed_seconds(d))
_dur(d::Real) = Float64(d)

function _aggregate(values::Vector{Float64}, aggregate::Symbol)
    aggregate == :all && return values
    isempty(values) && return NaN
    aggregate == :mean && return mean(values)
    aggregate == :median && return median(values)
    aggregate == :total && return sum(values)
    throw(ArgumentError("aggregate must be :mean, :median, :total, or :all"))
end

"""
    tEdgeFormation(dnet::DynamicNetwork, onset, terminus) -> Int

The number of edge-spell **onset events** in `[onset, terminus)` —
formations, counted as events like `tsna::tEdgeFormationAt` (not by
comparing point samples).
"""
function tEdgeFormation(dnet::DynamicNetwork{T, Time}, onset, terminus) where {T, Time}
    onset, terminus = convert(Time, onset), convert(Time, terminus)
    count = 0
    for spells in values(dnet.edge_spells), s in spells
        onset <= s.onset < terminus && (count += 1)
    end
    return count
end

"""
    tEdgeDissolution(dnet::DynamicNetwork, onset, terminus) -> Int

The number of edge-spell **terminus events** in `[onset, terminus)` —
dissolutions counted as events. Right-censored spells (terminus flagged
censored) are excluded.
"""
function tEdgeDissolution(dnet::DynamicNetwork{T, Time}, onset, terminus) where {T, Time}
    onset, terminus = convert(Time, onset), convert(Time, terminus)
    count = 0
    for spells in values(dnet.edge_spells), s in spells
        s.terminus_censored && continue
        onset <= s.terminus < terminus && (count += 1)
    end
    return count
end

"""
    tTurnover(dnet::DynamicNetwork, window_size) -> Vector{NamedTuple}

Formation and dissolution **event counts and rates** per window of length
`window_size` across the observation period. Every element has the same
shape: `(window_start, window_end, n_formations, n_dissolutions,
formation_rate, dissolution_rate)`, with rates per unit time.
"""
function tTurnover(dnet::DynamicNetwork{T, Time}, window_size) where {T, Time}
    start_time, end_time = dnet.observation_period
    results = NamedTuple[]

    t = start_time
    while t < end_time
        w_end = min(t + window_size, end_time)
        nf = tEdgeFormation(dnet, t, w_end)
        nd = tEdgeDissolution(dnet, t, w_end)
        span = _dur(w_end - t)
        push!(results, (window_start=t, window_end=w_end,
                        n_formations=nf, n_dissolutions=nd,
                        formation_rate=span > 0 ? nf / span : NaN,
                        dissolution_rate=span > 0 ? nd / span : NaN))
        t = w_end
    end

    return results
end

"""
    tieDecay(dnet::DynamicNetwork; method=:exponential, rate=1.0, at=observation end)

Tie weights decayed by the time since each edge was last active:
`exp(-rate·Δ)` (`:exponential`) or `max(0, 1 - rate·Δ)` (`:linear`),
where Δ is the time from the end of the edge's most recent spell to `at`
(0 for currently active ties).
"""
function tieDecay(dnet::DynamicNetwork{T, Time};
                  method::Symbol=:exponential, rate::Float64=1.0,
                  at=dnet.observation_period[2]) where {T, Time}
    at = convert(Time, at)
    weights = Dict{Tuple{T, T}, Float64}()

    for (edge, spells) in dnet.edge_spells
        isempty(spells) && continue
        # Time since last activity (0 if active at `at`)
        Δ = Inf
        for s in spells
            if NetworkDynamic._spell_active_at(s, at)
                Δ = 0.0
                break
            end
            s.terminus <= at && (Δ = min(Δ, _dur(at - s.terminus)))
        end
        isinf(Δ) && continue  # only spells entirely after `at`

        w = method == :exponential ? exp(-rate * Δ) :
            method == :linear ? max(0.0, 1.0 - rate * Δ) :
            throw(ArgumentError("method must be :exponential or :linear"))
        weights[edge] = w
    end

    return weights
end

# =============================================================================
# Aggregation and Time Series
# =============================================================================

"""
    tSnaStats(dnet::DynamicNetwork, times; measures=[:density, :reciprocity,
              :transitivity, :mean_degree]) -> Vector{NamedTuple}

Snapshot statistics at each time point. Each snapshot is extracted once
and reused for all measures. `mean_degree` counts each edge appropriately
for the network's directedness.
"""
function tSnaStats(dnet::DynamicNetwork{T, Time}, times::AbstractVector;
                   measures::Vector{Symbol}=[:density, :reciprocity,
                                             :transitivity, :mean_degree]) where {T, Time}
    results = NamedTuple[]

    for t in times
        snapshot = _snapshot(dnet, t)
        n = Int(nv(snapshot))
        row = Dict{Symbol, Any}(:time => t)

        for m in measures
            row[m] = if m == :density
                network_density(snapshot)
            elseif m == :reciprocity
                _snapshot_reciprocity(snapshot)
            elseif m == :transitivity
                SNA.transitivity(snapshot)
            elseif m == :mean_degree
                n == 0 ? NaN :
                    (is_directed(snapshot) ? ne(snapshot) / n : 2 * ne(snapshot) / n)
            elseif m == :n_edges
                Float64(ne(snapshot))
            else
                throw(ArgumentError("unknown measure: $m"))
            end
        end

        push!(results, NamedTuple{Tuple(sort(collect(keys(row))))}(
            Tuple(row[k] for k in sort(collect(keys(row))))))
    end

    return results
end

"""
    windowSnaStats(dnet::DynamicNetwork, window_size; kwargs...) -> Vector{NamedTuple}

Snapshot statistics sampled at the start of each window of length
`window_size`.
"""
function windowSnaStats(dnet::DynamicNetwork{T, Time}, window_size;
                        kwargs...) where {T, Time}
    start_time, end_time = dnet.observation_period
    times = Time[]
    t = start_time
    while t < end_time
        push!(times, t)
        t += window_size
    end
    return tSnaStats(dnet, times; kwargs...)
end

"""
    tAggregate(dnet::DynamicNetwork; method=:union, onset=..., terminus=...)

Collapse the dynamic network into a static one.

- `:union` — edges ever active in the window
- `:intersection` — edges active throughout the window
- `:weighted` — union, with each edge's total active time (in the window,
  clipped) stored as its `:weight` attribute
"""
function tAggregate(dnet::DynamicNetwork{T, Time}; method::Symbol=:union,
                    onset=dnet.observation_period[1],
                    terminus=dnet.observation_period[2]) where {T, Time}
    onset, terminus = convert(Time, onset), convert(Time, terminus)

    if method == :union
        return network_collapse(dnet; onset=onset, terminus=terminus, rule=:any)
    elseif method == :intersection
        return network_collapse(dnet; onset=onset, terminus=terminus, rule=:all)
    elseif method == :weighted
        collapsed = network_collapse(dnet; onset=onset, terminus=terminus, rule=:any)
        for (edge, spells) in dnet.edge_spells
            total = 0.0
            for s in spells
                lo = max(s.onset, onset)
                hi = min(s.terminus, terminus)
                hi > lo && (total += _dur(hi - lo))
            end
            if total > 0 && has_edge(collapsed, edge[1], edge[2])
                set_edge_attribute!(collapsed, :weight, edge[1], edge[2], total)
            end
        end
        return collapsed
    else
        throw(ArgumentError("method must be :union, :intersection, or :weighted"))
    end
end

end # module
