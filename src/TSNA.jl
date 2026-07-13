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

All public functions have snake_case primary names (Julia convention)
with R-style camelCase aliases (e.g. `earliest_arrival` /
`earliestArrival`); both are exported and interchangeable.

Port of the R tsna package from the StatNet collection.
"""
module TSNA

using DataStructures: BinaryMinHeap
using Dates
using Graphs
using Networks
using NetworkDynamic
using SNA
using Statistics

# Temporal paths and reachability
export TemporalPath, tPath, path_duration
export temporal_distance, temporalDistance
export earliest_arrival, earliestArrival

# Batch / all-source temporal paths: one reused workspace instead of one set of
# scratch containers per source (TSNA.jl#1)
export TemporalPathWorkspace, earliest_arrival!
export earliest_arrival_all, earliestArrivalAll
export temporal_distance_matrix, reachability_matrix
export forward_reachable_set, forwardReachableSet
export backward_reachable_set, backwardReachableSet
export temporal_path, temporalPath
export shortest_temporal_path, shortestTemporalPath

# Point-in-time measures
export t_degree, tDegree
export t_betweenness, tBetweenness
export t_closeness, tCloseness
export t_eigenvector, tEigenvector
export t_pagerank, tPagerank
export t_density, tDensity
export t_reciprocity, tReciprocity
export t_transitivity, tTransitivity

# Duration and turnover
export t_edge_duration, tEdgeDuration
export t_vertex_duration, tVertexDuration
export t_edge_formation, tEdgeFormation
export t_edge_dissolution, tEdgeDissolution
export t_edge_persistence, tEdgePersistence
export t_turnover, tTurnover
export tie_decay, tieDecay

# Contact sequences
export Contact, ContactSequence, as_contact_sequence

# Aggregation and time series
export t_sna_stats, tSnaStats
export window_sna_stats, windowSnaStats
export t_aggregate, tAggregate

# =============================================================================
# Temporal Path Types
# =============================================================================

"""
    TemporalPath{T, Time}

A time-respecting path through a dynamic network: `times[k]` is the
instant edge `edges[k]` is traversed, and times are non-decreasing.

Also available under the R-style alias `tPath`.
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

"""
    tPath

R-style alias for [`TemporalPath`](@ref).
"""
const tPath = TemporalPath

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

Elapsed time between the first and last traversal of the path.
"""
function path_duration(p::TemporalPath{T, Time}) where {T, Time}
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
    as_contact_sequence(dnet::DynamicNetwork; missing=:error, report=false) -> ContactSequence

Convert a dynamic network's edge spells to a contact sequence: one `Contact`
per edge spell, carrying its onset and duration. Overlapping spells on the
same edge stay separate contacts; a point spell `[t,t)` becomes a contact of
zero duration.

# Conversion invariants

Preserved: the vertex count, directedness, and every edge spell's onset and
duration (so the spell set is reconstructable).

A `ContactSequence` has no slot for the rest of the dynamic network, so the
conversion is lossy by nature: **spell censoring flags**, **vertex activity
spells** (actor presence/composition), static and time-varying attributes, and
the observation window are dropped. Pass `report=true` for
`(cs, ::Networks.ConversionReport)` naming them.

A `Contact` cannot record that a dyad is *unobserved*, so a network with a
missing-dyad mask is **rejected** by default (`missing=:error`) rather than
being flattened into contacts that read as observed. Pass `missing=:face` to
convert the recorded face values anyway (the mask is then dropped, and the
report says so).
"""
function as_contact_sequence(dnet::DynamicNetwork{T, Time};
                             missing::Symbol=:error,
                             report::Bool=false) where {T, Time}
    require_observed(dnet.network, missing; context="as_contact_sequence")

    contacts = Contact{T, Time}[]

    for ((i, j), spells) in dnet.edge_spells
        for spell in spells
            push!(contacts, Contact{T, Time}(i, j, spell.onset,
                                             spell.terminus - spell.onset))
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
    record_drop!(rep, :observation_period,
                 "the observation window $(dnet.observation_period) has no " *
                 "contact-sequence counterpart")
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

# Memoized contact index: rebuilding the per-vertex contact lists on every
# path query is O(total spells), which dominated repeated-query workloads
# (e.g. backward_reachable_set runs one search per vertex). The cache is
# keyed weakly by network identity and invalidated via the network's
# `mutation_count`, which NetworkDynamic bumps on every spell mutation.
const _CONTACT_INDEX_LOCK = ReentrantLock()
const _CONTACT_INDEX_CACHE = WeakKeyDict{Any, Tuple{Int, Any}}()

function _contact_index(dnet::DynamicNetwork{T, Time}) where {T, Time}
    lock(_CONTACT_INDEX_LOCK) do
        entry = get(_CONTACT_INDEX_CACHE, dnet, nothing)
        if !isnothing(entry) && entry[1] == dnet.mutation_count
            return entry[2]::Dict{T, Vector{Tuple{T, Spell{Time}}}}
        end
        index = _out_contacts(dnet)
        _CONTACT_INDEX_CACHE[dnet] = (dnet.mutation_count, index)
        return index
    end
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
    earliest_arrival(dnet, source, start_time; end_time=observation end,
                     target=nothing) -> (arrival::Dict, parent::Dict)

Earliest-arrival times from `source` to every vertex, starting at
`start_time`, under interval semantics: an edge spell `[onset, terminus)`
is traversable at any instant in it (boarding mid-spell is allowed;
spells that began before `start_time` but are still active count).
Implemented as a heap-based Dijkstra label-setting search over a memoized
per-network contact index, so chains of simultaneous (equal-onset) spells
are handled correctly and repeated queries do not rebuild the index.

With `target` set, the search stops as soon as that vertex is settled
(its arrival time is already final); the returned dictionaries then only
cover the explored part of the network.

Returns the arrival-time dictionary (vertices absent = unreachable) and
the parent map `(vertex => (predecessor, boarding time))` for path
reconstruction.

Also available under the R-style alias `earliestArrival`.
"""
function earliest_arrival(dnet::DynamicNetwork{T, Time}, source::T, start_time;
                          end_time=dnet.observation_period[2],
                          target::Union{Nothing, T}=nothing) where {T, Time}
    ws = TemporalPathWorkspace{T, Time}()
    return earliest_arrival!(ws, dnet, source, start_time;
                             end_time=end_time, target=target)
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
                      end_time=observation end, target=nothing)
        -> (arrival::Dict, parent::Dict)

In-place [`earliest_arrival`](@ref): runs the same search, reusing `ws`'s
containers instead of allocating fresh ones.

**The returned dictionaries alias `ws`** and are overwritten by the next search
on the same workspace. Copy what you need to keep, or use
[`earliest_arrival_all`](@ref), which does that for you.
"""
function earliest_arrival!(ws::TemporalPathWorkspace{T, Time},
                           dnet::DynamicNetwork{T, Time}, source::T, start_time;
                           end_time=dnet.observation_period[2],
                           target::Union{Nothing, T}=nothing) where {T, Time}
    start_time = convert(Time, start_time)
    end_time = convert(Time, end_time)
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
                         end_time=observation end) -> Dict{T, Dict{T, Time}}

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

Also available under the R-style alias `earliestArrivalAll`.
"""
function earliest_arrival_all(dnet::DynamicNetwork{T, Time}, start_time;
                              sources=T.(1:nv(dnet)),
                              end_time=dnet.observation_period[2]) where {T, Time}
    ws = TemporalPathWorkspace{T, Time}()
    result = Dict{T, Dict{T, Time}}()
    for s in sources
        arrival, _ = earliest_arrival!(ws, dnet, T(s), start_time; end_time=end_time)
        result[T(s)] = copy(arrival)      # detach from the workspace
    end
    return result
end

"""
    earliestArrivalAll

R-style alias for [`earliest_arrival_all`](@ref).
"""
const earliestArrivalAll = earliest_arrival_all

"""
    temporal_distance_matrix(dnet, start_time; end_time=observation end)
        -> Matrix{Union{Time, Nothing}}

All-pairs temporal distances (elapsed time of the earliest time-respecting path,
`arrival − start_time`), in one batched pass. `nothing` where no path exists
within the window; zero on the diagonal.

Computed with a single reused workspace via [`earliest_arrival_all`](@ref)
rather than `nv(dnet)²` independent [`temporal_distance`](@ref) calls.
"""
function temporal_distance_matrix(dnet::DynamicNetwork{T, Time}, start_time;
                                  end_time=dnet.observation_period[2]) where {T, Time}
    n = nv(dnet)
    t0 = convert(Time, start_time)
    all_arr = earliest_arrival_all(dnet, t0; end_time=end_time)
    D = Matrix{Union{Time, Nothing}}(nothing, n, n)
    for i in 1:n
        arr = all_arr[T(i)]
        for j in 1:n
            haskey(arr, T(j)) && (D[i, j] = arr[T(j)] - t0)
        end
    end
    return D
end

"""
    reachability_matrix(dnet, start_time; end_time=observation end) -> BitMatrix

`R[i, j]` is `true` when a time-respecting path runs from `i` to `j` within the
window (the diagonal is `true`). One batched pass; see
[`earliest_arrival_all`](@ref).
"""
function reachability_matrix(dnet::DynamicNetwork{T, Time}, start_time;
                             end_time=dnet.observation_period[2]) where {T, Time}
    n = nv(dnet)
    all_arr = earliest_arrival_all(dnet, start_time; end_time=end_time)
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
    earliestArrival

R-style alias for [`earliest_arrival`](@ref).
"""
const earliestArrival = earliest_arrival

"""
    temporal_distance(dnet, source, target, start_time; end_time=...)
        -> elapsed time or nothing

Elapsed time of the earliest time-respecting path from `source` to
`target` departing at `start_time` (arrival − start). Returns `nothing`
when no path exists within the window.

Also available under the R-style alias `temporalDistance`.
"""
function temporal_distance(dnet::DynamicNetwork{T, Time}, source::T, target::T,
                           start_time; end_time=dnet.observation_period[2]) where {T, Time}
    arrival, _ = earliest_arrival(dnet, source, start_time;
                                  end_time=end_time, target=target)
    haskey(arrival, target) || return nothing
    return arrival[target] - convert(Time, start_time)
end

"""
    temporalDistance

R-style alias for [`temporal_distance`](@ref).
"""
const temporalDistance = temporal_distance

"""
    forward_reachable_set(dnet, source, start_time; end_time=...) -> Vector

Vertices reachable from `source` by a time-respecting path departing at
or after `start_time` and arriving before `end_time` (`source` included).

Also available under the R-style alias `forwardReachableSet`.
"""
function forward_reachable_set(dnet::DynamicNetwork{T, Time}, source::T, start_time;
                               end_time=dnet.observation_period[2]) where {T, Time}
    arrival, _ = earliest_arrival(dnet, source, start_time; end_time=end_time)
    return sort(collect(keys(arrival)))
end

"""
    forwardReachableSet

R-style alias for [`forward_reachable_set`](@ref).
"""
const forwardReachableSet = forward_reachable_set

"""
    backward_reachable_set(dnet, target, end_time; start_time=...) -> Vector

Vertices from which `target` can be reached by a time-respecting path in
`[start_time, end_time)` — the exact dual of
[`forward_reachable_set`](@ref) (computed by forward searches that stop
as soon as `target` is settled, so both use identical traversal
semantics).

Also available under the R-style alias `backwardReachableSet`.
"""
function backward_reachable_set(dnet::DynamicNetwork{T, Time}, target::T, end_time;
                                start_time=dnet.observation_period[1]) where {T, Time}
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
                                       end_time=end_time, target=target)
        haskey(arrival, target) && push!(reachable, vT)
    end
    return reachable
end

"""
    backwardReachableSet

R-style alias for [`backward_reachable_set`](@ref).
"""
const backwardReachableSet = backward_reachable_set

"""
    temporal_path(dnet, source, target, start_time; end_time=...)
        -> Union{TemporalPath, Nothing}

The **earliest-arrival** time-respecting path from `source` to `target`
departing at `start_time` (not necessarily the fewest-hops path).
Returns `nothing` when no path exists.

Also available under the R-style alias `temporalPath`.
"""
function temporal_path(dnet::DynamicNetwork{T, Time}, source::T, target::T,
                       start_time; end_time=dnet.observation_period[2]) where {T, Time}
    arrival, parent = earliest_arrival(dnet, source, start_time;
                                       end_time=end_time, target=target)
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

"""
    temporalPath

R-style alias for [`temporal_path`](@ref).
"""
const temporalPath = temporal_path

"""
    shortest_temporal_path(dnet, source, target, start_time; kwargs...)

Alias for [`temporal_path`](@ref): the *earliest-arrival* path (kept for
API compatibility with earlier versions; note this is not the fewest-hops
path). Also available under the R-style alias `shortestTemporalPath`.
"""
const shortest_temporal_path = temporal_path

"""
    shortestTemporalPath

R-style alias for [`shortest_temporal_path`](@ref).
"""
const shortestTemporalPath = temporal_path

# =============================================================================
# Temporal Measures at a Point
# =============================================================================
#
# Snapshots retain every vertex (retain_all_vertices=true), so all
# per-vertex vectors are indexed by the dynamic network's own vertex IDs
# and always have length nv(dnet) — even when some vertices are inactive.

_snapshot(dnet, at) = network_extract(dnet, at; retain_all_vertices=true)

"""
    t_degree(dnet::DynamicNetwork, at; mode=:total) -> Vector{Float64}

Degree at time `at`, indexed by the network's vertex IDs (inactive
vertices score 0).

Also available under the R-style alias `tDegree`.
"""
function t_degree(dnet::DynamicNetwork{T, Time}, at; mode::Symbol=:total) where {T, Time}
    return SNA.degree_centrality(_snapshot(dnet, at); mode=mode)
end

"""
    tDegree

R-style alias for [`t_degree`](@ref).
"""
const tDegree = t_degree

"""
    t_betweenness(dnet::DynamicNetwork, at; normalized=false) -> Vector{Float64}

Betweenness centrality at time `at` (raw scores by default, as in SNA.jl).

Also available under the R-style alias `tBetweenness`.
"""
t_betweenness(dnet::DynamicNetwork, at; normalized::Bool=false) =
    SNA.betweenness_centrality(_snapshot(dnet, at); normalized=normalized)

"""
    tBetweenness

R-style alias for [`t_betweenness`](@ref).
"""
const tBetweenness = t_betweenness

"""
    t_closeness(dnet::DynamicNetwork, at) -> Vector{Float64}

Closeness centrality at time `at`. Also available under the R-style alias
`tCloseness`.
"""
t_closeness(dnet::DynamicNetwork, at) = SNA.closeness_centrality(_snapshot(dnet, at))

"""
    tCloseness

R-style alias for [`t_closeness`](@ref).
"""
const tCloseness = t_closeness

"""
    t_eigenvector(dnet::DynamicNetwork, at) -> Vector{Float64}

Eigenvector centrality at time `at`. Also available under the R-style
alias `tEigenvector`.
"""
t_eigenvector(dnet::DynamicNetwork, at) = SNA.eigenvector_centrality(_snapshot(dnet, at))

"""
    tEigenvector

R-style alias for [`t_eigenvector`](@ref).
"""
const tEigenvector = t_eigenvector

"""
    t_pagerank(dnet::DynamicNetwork, at; damping=0.85) -> Vector{Float64}

PageRank at time `at`. Also available under the R-style alias `tPagerank`.
"""
t_pagerank(dnet::DynamicNetwork, at; damping::Float64=0.85) =
    SNA.pagerank(_snapshot(dnet, at); α=damping)

"""
    tPagerank

R-style alias for [`t_pagerank`](@ref).
"""
const tPagerank = t_pagerank

"""
    t_density(dnet::DynamicNetwork, at) -> Float64

Density at time `at`, over all `nv(dnet)` vertices (not just the active
ones).

Also available under the R-style alias `tDensity`.
"""
function t_density(dnet::DynamicNetwork, at)
    return network_density(_snapshot(dnet, at))
end

"""
    tDensity

R-style alias for [`t_density`](@ref).
"""
const tDensity = t_density

"""
    t_reciprocity(dnet::DynamicNetwork, at) -> Float64

Edgewise reciprocity at time `at` (fraction of edges that are
reciprocated).

Also available under the R-style alias `tReciprocity`.
"""
function t_reciprocity(dnet::DynamicNetwork, at)
    return _snapshot_reciprocity(_snapshot(dnet, at))
end

"""
    tReciprocity

R-style alias for [`t_reciprocity`](@ref).
"""
const tReciprocity = t_reciprocity

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
    t_transitivity(dnet::DynamicNetwork, at) -> Float64

Weak transitivity at time `at` (see `SNA.transitivity`).

Also available under the R-style alias `tTransitivity`.
"""
t_transitivity(dnet::DynamicNetwork, at) = SNA.transitivity(_snapshot(dnet, at))

"""
    tTransitivity

R-style alias for [`t_transitivity`](@ref).
"""
const tTransitivity = t_transitivity

# =============================================================================
# Duration and Turnover Metrics
# =============================================================================

"""
    t_edge_duration(dnet::DynamicNetwork; mode=:spell, aggregate=:mean)

Edge activity durations.

- `mode=:spell` (default, matching `tsna::edgeDuration`): each activity
  spell contributes its own duration `terminus − onset`.
- `mode=:total`: spells of the same edge are summed first (total active
  time per edge).

`aggregate` is `:mean`, `:median`, `:total`, or `:all` (raw vector).
Censoring flags are ignored (durations are observed, not corrected).

Also available under the R-style alias `tEdgeDuration`.
"""
function t_edge_duration(dnet::DynamicNetwork{T, Time};
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
    tEdgeDuration

R-style alias for [`t_edge_duration`](@ref).
"""
const tEdgeDuration = t_edge_duration

"""
    t_vertex_duration(dnet::DynamicNetwork; mode=:spell, aggregate=:mean)

Vertex activity durations (see [`t_edge_duration`](@ref)).

Also available under the R-style alias `tVertexDuration`.
"""
function t_vertex_duration(dnet::DynamicNetwork{T, Time};
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

"""
    tVertexDuration

R-style alias for [`t_vertex_duration`](@ref).
"""
const tVertexDuration = t_vertex_duration

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
    t_edge_formation(dnet::DynamicNetwork, onset, terminus) -> Int

The number of edge-spell **onset events** in `[onset, terminus)` —
formations, counted as events like `tsna::tEdgeFormationAt` (not by
comparing point samples).

Also available under the R-style alias `tEdgeFormation`.
"""
function t_edge_formation(dnet::DynamicNetwork{T, Time}, onset, terminus) where {T, Time}
    onset, terminus = convert(Time, onset), convert(Time, terminus)
    count = 0
    for spells in values(dnet.edge_spells), s in spells
        onset <= s.onset < terminus && (count += 1)
    end
    return count
end

"""
    tEdgeFormation

R-style alias for [`t_edge_formation`](@ref).
"""
const tEdgeFormation = t_edge_formation

"""
    t_edge_dissolution(dnet::DynamicNetwork, onset, terminus) -> Int

The number of edge-spell **terminus events** in `[onset, terminus)` —
dissolutions counted as events. Right-censored spells (terminus flagged
censored) are excluded.

Also available under the R-style alias `tEdgeDissolution`.
"""
function t_edge_dissolution(dnet::DynamicNetwork{T, Time}, onset, terminus) where {T, Time}
    onset, terminus = convert(Time, onset), convert(Time, terminus)
    count = 0
    for spells in values(dnet.edge_spells), s in spells
        s.terminus_censored && continue
        onset <= s.terminus < terminus && (count += 1)
    end
    return count
end

"""
    tEdgeDissolution

R-style alias for [`t_edge_dissolution`](@ref).
"""
const tEdgeDissolution = t_edge_dissolution

"""
    t_edge_persistence(dnet::DynamicNetwork, window_size) -> Float64

Proportion of edges that persist (survive) across consecutive time
windows, after the temporal-correlation measure of Nicosia et al.: the
observation period is divided into windows of length `window_size`; for
each consecutive pair of windows the edges active at the start of the
first window are checked for activity at the start of the second, and the
pooled proportion `persisted / total` is returned.

Values near 1 indicate a stable network (little edge turnover); values
near 0 indicate almost complete edge replacement per window. Returns
`NaN` when there are fewer than two windows or no active edges to track.

Also available under the R-style alias `tEdgePersistence`.
"""
function t_edge_persistence(dnet::DynamicNetwork{T, Time}, window_size) where {T, Time}
    start_time, end_time = dnet.observation_period

    starts = Time[]
    t = start_time
    while t < end_time
        push!(starts, t)
        t += window_size
    end
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
    tEdgePersistence

R-style alias for [`t_edge_persistence`](@ref).
"""
const tEdgePersistence = t_edge_persistence

"""
    t_turnover(dnet::DynamicNetwork, window_size) -> Vector{NamedTuple}

Formation and dissolution **event counts and rates** per window of length
`window_size` across the observation period. Every element has the same
shape: `(window_start, window_end, n_formations, n_dissolutions,
formation_rate, dissolution_rate)`, with rates per unit time.

Also available under the R-style alias `tTurnover`.
"""
function t_turnover(dnet::DynamicNetwork{T, Time}, window_size) where {T, Time}
    start_time, end_time = dnet.observation_period

    # Event times collected and sorted once, so each window is a pair of
    # binary searches instead of a full scan over every spell
    onsets = Time[]
    termini = Time[]
    for spells in values(dnet.edge_spells), s in spells
        push!(onsets, s.onset)
        s.terminus_censored || push!(termini, s.terminus)
    end
    sort!(onsets)
    sort!(termini)
    events_in(v, lo, hi) = searchsortedfirst(v, hi) - searchsortedfirst(v, lo)

    results = NamedTuple[]
    t = start_time
    while t < end_time
        w_end = min(t + window_size, end_time)
        nf = events_in(onsets, t, w_end)
        nd = events_in(termini, t, w_end)
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
    tTurnover

R-style alias for [`t_turnover`](@ref).
"""
const tTurnover = t_turnover

"""
    tie_decay(dnet::DynamicNetwork; method=:exponential, rate=1.0, at=observation end)

Tie weights decayed by the time since each edge was last active:
`exp(-rate·Δ)` (`:exponential`) or `max(0, 1 - rate·Δ)` (`:linear`),
where Δ is the time from the end of the edge's most recent spell to `at`
(0 for currently active ties).

Also available under the R-style alias `tieDecay`.
"""
function tie_decay(dnet::DynamicNetwork{T, Time};
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

"""
    tieDecay

R-style alias for [`tie_decay`](@ref).
"""
const tieDecay = tie_decay

# =============================================================================
# Aggregation and Time Series
# =============================================================================

"""
    t_sna_stats(dnet::DynamicNetwork, times; measures=[:density, :reciprocity,
                :transitivity, :mean_degree]) -> Vector{NamedTuple}

Snapshot statistics at each time point. Each snapshot is extracted once
and reused for all measures. `mean_degree` counts each edge appropriately
for the network's directedness.

Also available under the R-style alias `tSnaStats`.
"""
function t_sna_stats(dnet::DynamicNetwork{T, Time}, times::AbstractVector;
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
    tSnaStats

R-style alias for [`t_sna_stats`](@ref).
"""
const tSnaStats = t_sna_stats

"""
    window_sna_stats(dnet::DynamicNetwork, window_size; kwargs...) -> Vector{NamedTuple}

Snapshot statistics sampled at the start of each window of length
`window_size`.

Also available under the R-style alias `windowSnaStats`.
"""
function window_sna_stats(dnet::DynamicNetwork{T, Time}, window_size;
                          kwargs...) where {T, Time}
    start_time, end_time = dnet.observation_period
    times = Time[]
    t = start_time
    while t < end_time
        push!(times, t)
        t += window_size
    end
    return t_sna_stats(dnet, times; kwargs...)
end

"""
    windowSnaStats

R-style alias for [`window_sna_stats`](@ref).
"""
const windowSnaStats = window_sna_stats

"""
    t_aggregate(dnet::DynamicNetwork; method=:union, onset=..., terminus=..., report=false)

Collapse the dynamic network into a static one.

- `:union` — edges ever active in the window
- `:intersection` — edges active throughout the window
- `:weighted` — union, with each edge's total active time (in the window,
  clipped) stored as its `:weight` attribute

This is `NetworkDynamic.network_collapse` with an aggregation rule, and it
inherits its conversion invariants: vertex IDs are stable, so directedness,
`loops`, two-mode metadata, static attributes and the **missing-dyad mask** all
survive; spells, time-varying attributes and the observation window are dropped
by nature. Pass `report=true` for `(net, ::Networks.ConversionReport)`.

Also available under the R-style alias `tAggregate`.
"""
function t_aggregate(dnet::DynamicNetwork{T, Time}; method::Symbol=:union,
                     onset=dnet.observation_period[1],
                     terminus=dnet.observation_period[2],
                     report::Bool=false) where {T, Time}
    onset, terminus = convert(Time, onset), convert(Time, terminus)

    method in (:union, :intersection, :weighted) ||
        throw(ArgumentError("method must be :union, :intersection, or :weighted"))

    rule = method == :intersection ? :all : :any
    collapsed, rep = network_collapse(dnet; onset=onset, terminus=terminus,
                                      rule=rule, report=true)

    if method == :weighted
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
    end

    return report ? (collapsed, rep) : collapsed
end

"""
    tAggregate

R-style alias for [`t_aggregate`](@ref).
"""
const tAggregate = t_aggregate

end # module
