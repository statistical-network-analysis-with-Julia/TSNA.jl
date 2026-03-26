"""
    TSNA.jl - Temporal Social Network Analysis

Provides descriptive analysis tools for dynamic networks including
temporal centrality, reachability, path analysis, and duration metrics.

Port of the R tsna package from the StatNet collection.
"""
module TSNA

using Graphs
using LinearAlgebra
using Network
using NetworkDynamic
using SNA
using Statistics
using StatsBase

# Temporal paths
export tPath, tReachable
export temporalDistance, forwardReachableSet, backwardReachableSet
export temporalPath, shortestTemporalPath

# Temporal centrality
export tDegree, tBetweenness, tCloseness
export tEigenvector, tPagerank

# Temporal network measures
export tDensity, tReciprocity, tTransitivity
export tEdgeDuration, tVertexDuration
export tEdgePersistence, tTurnover

# Contact sequences
export Contact, ContactSequence, as_contact_sequence
export tieDecay

# Aggregation
export tSnaStats, tAggregate
export windowSnaStats

# =============================================================================
# Temporal Path Types
# =============================================================================

"""
    tPath{T, Time}

A temporal path through a dynamic network.

A valid temporal path must have non-decreasing times along the path.

# Fields
- `vertices::Vector{T}`: Sequence of vertices
- `times::Vector{Time}`: Times at which each transition occurs
- `edges::Vector{Tuple{T, T}}`: Edges traversed
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
    path_duration(p::tPath) -> Time

Get the total duration of a temporal path.
"""
function path_duration(p::tPath{T, Time}) where {T, Time}
    isempty(p.times) && return zero(Time)
    return p.times[end] - p.times[1]
end

# =============================================================================
# Contact Sequence
# =============================================================================

"""
    Contact{T, Time}

A single contact (edge activation) in a temporal network.
"""
struct Contact{T, Time}
    source::T
    target::T
    time::Time
    duration::Time
end

"""
    ContactSequence{T, Time}

A sequence of contacts ordered by time.
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

Convert a dynamic network to a contact sequence.
"""
function as_contact_sequence(dnet::DynamicNetwork{T, Time}) where {T, Time}
    contacts = Contact{T, Time}[]

    for ((i, j), spells) in dnet.edge_spells
        for spell in spells
            push!(contacts, Contact{T, Time}(i, j, spell.onset, spell.terminus - spell.onset))
        end
    end

    return ContactSequence(contacts, nv(dnet); directed=is_directed(dnet))
end

# =============================================================================
# Temporal Measures at a Point
# =============================================================================

"""
    tDegree(dnet::DynamicNetwork, at::Time; mode=:total) -> Vector{Int}

Compute degree centrality at a specific time point.

# Arguments
- `dnet`: Dynamic network
- `at`: Time point
- `mode`: :in, :out, or :total for directed networks
"""
function tDegree(dnet::DynamicNetwork{T, Time}, at::Time;
                 mode::Symbol=:total) where {T, Time}
    snapshot = network_extract(dnet, at)

    if mode == :total
        return [Graphs.degree(snapshot, v) for v in vertices(snapshot)]
    elseif mode == :in
        return [Graphs.indegree(snapshot, v) for v in vertices(snapshot)]
    elseif mode == :out
        return [Graphs.outdegree(snapshot, v) for v in vertices(snapshot)]
    else
        throw(ArgumentError("mode must be :in, :out, or :total"))
    end
end

"""
    tDensity(dnet::DynamicNetwork, at::Time) -> Float64

Compute network density at a specific time point.
"""
function tDensity(dnet::DynamicNetwork{T, Time}, at::Time) where {T, Time}
    snapshot = network_extract(dnet, at)
    n = nv(snapshot)
    n <= 1 && return 0.0

    m = ne(snapshot)
    max_edges = is_directed(snapshot) ? n * (n - 1) : n * (n - 1) ÷ 2
    return m / max_edges
end

"""
    tReciprocity(dnet::DynamicNetwork, at::Time) -> Float64

Compute reciprocity at a specific time point (directed networks).
"""
function tReciprocity(dnet::DynamicNetwork{T, Time}, at::Time) where {T, Time}
    snapshot = network_extract(dnet, at)
    !is_directed(snapshot) && return 1.0

    n_edges = ne(snapshot)
    n_edges == 0 && return 0.0

    mutual = 0
    for e in edges(snapshot)
        if has_edge(snapshot, dst(e), src(e))
            mutual += 1
        end
    end

    return mutual / n_edges
end

"""
    tTransitivity(dnet::DynamicNetwork, at::Time) -> Float64

Compute transitivity (clustering coefficient) at a specific time point.
"""
function tTransitivity(dnet::DynamicNetwork{T, Time}, at::Time) where {T, Time}
    snapshot = network_extract(dnet, at)
    return Graphs.global_clustering_coefficient(snapshot)
end

"""
    tBetweenness(dnet::DynamicNetwork, at::Time; normalized=true) -> Vector{Float64}

Compute betweenness centrality at a specific time point.
"""
function tBetweenness(dnet::DynamicNetwork{T, Time}, at::Time;
                      normalized::Bool=true) where {T, Time}
    snapshot = network_extract(dnet, at)
    bc = Graphs.betweenness_centrality(snapshot; normalize=normalized)
    return bc
end

"""
    tCloseness(dnet::DynamicNetwork, at::Time) -> Vector{Float64}

Compute closeness centrality at a specific time point.
"""
function tCloseness(dnet::DynamicNetwork{T, Time}, at::Time) where {T, Time}
    snapshot = network_extract(dnet, at)
    return Graphs.closeness_centrality(snapshot)
end

"""
    tEigenvector(dnet::DynamicNetwork, at::Time) -> Vector{Float64}

Compute eigenvector centrality at a specific time point.
"""
function tEigenvector(dnet::DynamicNetwork{T, Time}, at::Time) where {T, Time}
    snapshot = network_extract(dnet, at)
    return Graphs.eigenvector_centrality(snapshot)
end

"""
    tPagerank(dnet::DynamicNetwork, at::Time; damping=0.85) -> Vector{Float64}

Compute PageRank at a specific time point.
"""
function tPagerank(dnet::DynamicNetwork{T, Time}, at::Time;
                   damping::Float64=0.85) where {T, Time}
    snapshot = network_extract(dnet, at)
    return Graphs.pagerank(snapshot, damping)
end

# =============================================================================
# Temporal Path Finding
# =============================================================================

"""
    temporalDistance(dnet::DynamicNetwork, source, target, start_time; kwargs...) -> Time

Find the earliest arrival time from source to target starting at start_time.
Returns Inf if no path exists.
"""
function temporalDistance(dnet::DynamicNetwork{T, Time}, source::T, target::T,
                          start_time::Time) where {T, Time}
    source == target && return start_time

    # BFS-style search through time
    n = nv(dnet)
    earliest_arrival = fill(typemax(Time), n)
    earliest_arrival[source] = start_time

    # Sort edge spells by onset time
    edge_events = Tuple{Time, T, T, Time}[]  # (onset, src, dst, terminus)
    for ((i, j), spells) in dnet.edge_spells
        for spell in spells
            if spell.onset >= start_time
                push!(edge_events, (spell.onset, i, j, spell.terminus))
            end
        end
    end
    sort!(edge_events, by=x -> x[1])

    # Process events in time order
    for (onset, i, j, terminus) in edge_events
        if earliest_arrival[i] <= onset
            # Can traverse this edge
            if onset < earliest_arrival[j]
                earliest_arrival[j] = onset
            end
        end
    end

    return earliest_arrival[target] == typemax(Time) ? convert(Time, Inf) : earliest_arrival[target]
end

"""
    forwardReachableSet(dnet::DynamicNetwork, source, start_time) -> Set{T}

Find all vertices reachable from source starting at start_time.
"""
function forwardReachableSet(dnet::DynamicNetwork{T, Time}, source::T,
                             start_time::Time) where {T, Time}
    n = nv(dnet)
    reachable = Set{T}([source])
    earliest_arrival = fill(typemax(Time), n)
    earliest_arrival[source] = start_time

    # Sort edges by onset
    edge_events = Tuple{Time, T, T}[]
    for ((i, j), spells) in dnet.edge_spells
        for spell in spells
            if spell.onset >= start_time
                push!(edge_events, (spell.onset, i, j))
            end
        end
    end
    sort!(edge_events, by=x -> x[1])

    for (onset, i, j) in edge_events
        if earliest_arrival[i] <= onset
            if onset < earliest_arrival[j]
                earliest_arrival[j] = onset
                push!(reachable, j)
            end
        end
    end

    return reachable
end

"""
    backwardReachableSet(dnet::DynamicNetwork, target, end_time) -> Set{T}

Find all vertices that can reach target by end_time.
"""
function backwardReachableSet(dnet::DynamicNetwork{T, Time}, target::T,
                              end_time::Time) where {T, Time}
    n = nv(dnet)
    reachable = Set{T}([target])
    latest_departure = fill(typemin(Time), n)
    latest_departure[target] = end_time

    # Sort edges by terminus (decreasing)
    edge_events = Tuple{Time, T, T}[]
    for ((i, j), spells) in dnet.edge_spells
        for spell in spells
            if spell.terminus <= end_time
                push!(edge_events, (spell.terminus, i, j))
            end
        end
    end
    sort!(edge_events, by=x -> -x[1])  # Decreasing order

    for (terminus, i, j) in edge_events
        if latest_departure[j] >= terminus
            if terminus > latest_departure[i]
                latest_departure[i] = terminus
                push!(reachable, i)
            end
        end
    end

    return reachable
end

"""
    shortestTemporalPath(dnet::DynamicNetwork, source, target, start_time) -> Union{tPath, Nothing}

Find a shortest temporal path from source to target.
"""
function shortestTemporalPath(dnet::DynamicNetwork{T, Time}, source::T, target::T,
                              start_time::Time) where {T, Time}
    source == target && return tPath([source], Time[], Tuple{T,T}[])

    n = nv(dnet)
    earliest_arrival = fill(typemax(Time), n)
    predecessor = fill((T(0), T(0), typemin(Time)), n)  # (prev_vertex, via_edge_src, time)
    earliest_arrival[source] = start_time

    # Sort edges by onset
    edge_events = Tuple{Time, T, T, Time}[]
    for ((i, j), spells) in dnet.edge_spells
        for spell in spells
            if spell.onset >= start_time
                push!(edge_events, (spell.onset, i, j, spell.terminus))
            end
        end
    end
    sort!(edge_events, by=x -> x[1])

    for (onset, i, j, terminus) in edge_events
        if earliest_arrival[i] <= onset && onset < earliest_arrival[j]
            earliest_arrival[j] = onset
            predecessor[j] = (i, i, onset)
        end
    end

    earliest_arrival[target] == typemax(Time) && return nothing

    # Reconstruct path
    vertices = T[target]
    times = Time[]
    edges = Tuple{T,T}[]

    current = target
    while current != source
        prev, _, time = predecessor[current]
        prev == 0 && return nothing
        pushfirst!(vertices, prev)
        pushfirst!(times, time)
        pushfirst!(edges, (prev, current))
        current = prev
    end

    return tPath(vertices, times, edges)
end

# =============================================================================
# Duration and Persistence Metrics
# =============================================================================

"""
    tEdgeDuration(dnet::DynamicNetwork; aggregate=:mean) -> Union{Float64, Dict}

Compute edge duration statistics.

# Arguments
- `aggregate`: :mean, :median, :total, or :all for per-edge values
"""
function tEdgeDuration(dnet::DynamicNetwork{T, Time};
                       aggregate::Symbol=:mean) where {T, Time}
    durations = Dict{Tuple{T,T}, Time}()

    for (edge, spells) in dnet.edge_spells
        total_dur = sum(spell.terminus - spell.onset for spell in spells; init=zero(Time))
        durations[edge] = total_dur
    end

    isempty(durations) && return 0.0

    if aggregate == :all
        return durations
    elseif aggregate == :mean
        return mean(values(durations))
    elseif aggregate == :median
        return median(collect(values(durations)))
    elseif aggregate == :total
        return sum(values(durations))
    else
        throw(ArgumentError("aggregate must be :mean, :median, :total, or :all"))
    end
end

"""
    tVertexDuration(dnet::DynamicNetwork; aggregate=:mean) -> Union{Float64, Vector}

Compute vertex activity duration statistics.
"""
function tVertexDuration(dnet::DynamicNetwork{T, Time};
                         aggregate::Symbol=:mean) where {T, Time}
    durations = zeros(Time, nv(dnet))

    for (v, spells) in dnet.vertex_spells
        durations[v] = sum(spell.terminus - spell.onset for spell in spells; init=zero(Time))
    end

    if aggregate == :all
        return durations
    elseif aggregate == :mean
        return mean(durations)
    elseif aggregate == :median
        return median(durations)
    elseif aggregate == :total
        return sum(durations)
    else
        throw(ArgumentError("aggregate must be :mean, :median, :total, or :all"))
    end
end

"""
    tEdgePersistence(dnet::DynamicNetwork, window_size::Time) -> Float64

Compute the proportion of edges that persist across time windows.
"""
function tEdgePersistence(dnet::DynamicNetwork{T, Time}, window_size::Time) where {T, Time}
    start_time, end_time = dnet.observation_period
    n_windows = ceil(Int, (end_time - start_time) / window_size)
    n_windows < 2 && return 1.0

    persisted = 0
    total = 0

    for w in 1:(n_windows - 1)
        t1 = start_time + (w - 1) * window_size
        t2 = start_time + w * window_size

        edges_t1 = active_edges(dnet, t1)
        edges_t2 = active_edges(dnet, t2)

        for e in edges_t1
            total += 1
            if e in edges_t2
                persisted += 1
            end
        end
    end

    return total == 0 ? 1.0 : persisted / total
end

"""
    tTurnover(dnet::DynamicNetwork, window_size::Time) -> NamedTuple

Compute edge turnover (formation and dissolution rates).
"""
function tTurnover(dnet::DynamicNetwork{T, Time}, window_size::Time) where {T, Time}
    start_time, end_time = dnet.observation_period
    n_windows = ceil(Int, (end_time - start_time) / window_size)
    n_windows < 2 && return (formation_rate=0.0, dissolution_rate=0.0)

    formations = 0
    dissolutions = 0
    at_risk_form = 0
    at_risk_diss = 0

    for w in 1:(n_windows - 1)
        t1 = start_time + (w - 1) * window_size
        t2 = start_time + w * window_size

        edges_t1 = Set(active_edges(dnet, t1))
        edges_t2 = Set(active_edges(dnet, t2))

        # New edges at t2
        new_edges = setdiff(edges_t2, edges_t1)
        formations += length(new_edges)

        # Lost edges
        lost_edges = setdiff(edges_t1, edges_t2)
        dissolutions += length(lost_edges)

        # At-risk counts
        n = nv(dnet)
        max_possible = is_directed(dnet) ? n * (n - 1) : n * (n - 1) ÷ 2
        at_risk_form += max_possible - length(edges_t1)
        at_risk_diss += length(edges_t1)
    end

    formation_rate = at_risk_form == 0 ? 0.0 : formations / at_risk_form
    dissolution_rate = at_risk_diss == 0 ? 0.0 : dissolutions / at_risk_diss

    return (
        formation_rate=formation_rate,
        dissolution_rate=dissolution_rate,
        n_formations=formations,
        n_dissolutions=dissolutions
    )
end

"""
    tieDecay(dnet::DynamicNetwork; method=:exponential) -> Float64

Estimate tie decay rate from observed data.
"""
function tieDecay(dnet::DynamicNetwork{T, Time}; method::Symbol=:exponential) where {T, Time}
    durations = Float64[]

    for (edge, spells) in dnet.edge_spells
        for spell in spells
            push!(durations, spell.terminus - spell.onset)
        end
    end

    isempty(durations) && return 0.0

    if method == :exponential
        # MLE for exponential distribution rate parameter
        return 1.0 / mean(durations)
    elseif method == :halflife
        return log(2) / mean(durations)
    else
        throw(ArgumentError("method must be :exponential or :halflife"))
    end
end

# =============================================================================
# Aggregation and Time Series
# =============================================================================

"""
    tSnaStats(dnet::DynamicNetwork, times::AbstractVector{Time};
              stats::Vector{Symbol}=[:density, :reciprocity]) -> Dict{Symbol, Vector}

Compute SNA statistics at multiple time points.
"""
function tSnaStats(dnet::DynamicNetwork{T, Time}, times::AbstractVector{Time};
                   stats::Vector{Symbol}=[:density, :reciprocity]) where {T, Time}
    result = Dict{Symbol, Vector{Float64}}()

    for stat in stats
        result[stat] = Float64[]
    end

    for t in times
        snapshot = network_extract(dnet, t)
        n = nv(snapshot)

        for stat in stats
            val = if stat == :density
                n <= 1 ? 0.0 : ne(snapshot) / (is_directed(snapshot) ? n*(n-1) : n*(n-1)÷2)
            elseif stat == :reciprocity
                tReciprocity(dnet, t)
            elseif stat == :transitivity
                n >= 3 ? Graphs.global_clustering_coefficient(snapshot) : 0.0
            elseif stat == :n_edges
                Float64(ne(snapshot))
            elseif stat == :n_vertices
                Float64(n)
            elseif stat == :mean_degree
                n == 0 ? 0.0 : 2.0 * ne(snapshot) / n
            else
                @warn "Unknown statistic: $stat"
                NaN
            end
            push!(result[stat], val)
        end
    end

    return result
end

"""
    windowSnaStats(dnet::DynamicNetwork, window_size::Time;
                   stats::Vector{Symbol}=[:density]) -> Dict{Symbol, Vector}

Compute SNA statistics in sliding windows.
"""
function windowSnaStats(dnet::DynamicNetwork{T, Time}, window_size::Time;
                        stats::Vector{Symbol}=[:density],
                        step::Time=window_size) where {T, Time}
    start_time, end_time = dnet.observation_period
    times = collect(start_time:step:(end_time - window_size))

    result = Dict{Symbol, Vector{Float64}}()
    for stat in stats
        result[stat] = Float64[]
    end

    for t in times
        window_net = network_extract(dnet, t, t + window_size; rule=:any)
        n = nv(window_net)

        for stat in stats
            val = if stat == :density
                n <= 1 ? 0.0 : ne(window_net) / (is_directed(window_net) ? n*(n-1) : n*(n-1)÷2)
            elseif stat == :n_edges
                Float64(ne(window_net))
            elseif stat == :mean_degree
                n == 0 ? 0.0 : 2.0 * ne(window_net) / n
            else
                NaN
            end
            push!(result[stat], val)
        end
    end

    return result
end

"""
    tAggregate(dnet::DynamicNetwork; method=:union) -> Network

Aggregate dynamic network to static using specified method.

# Methods
- `:union`: Include edge if ever active
- `:intersection`: Include edge if always active
- `:weighted`: Create weighted network by total activation time
"""
function tAggregate(dnet::DynamicNetwork{T, Time}; method::Symbol=:union) where {T, Time}
    if method == :union
        return network_collapse(dnet)
    elseif method == :intersection
        # Edge must be active throughout observation period
        collapsed = Network{T}(; n=nv(dnet.network), directed=is_directed(dnet.network))
        obs_start, obs_end = dnet.observation_period

        for (edge, spells) in dnet.edge_spells
            # Check if any spell covers entire observation period
            if any(s.onset <= obs_start && s.terminus >= obs_end for s in spells)
                add_edge!(collapsed, edge[1], edge[2])
            end
        end

        return collapsed
    elseif method == :weighted
        collapsed = Network{T}(; n=nv(dnet.network), directed=is_directed(dnet.network))

        for (edge, spells) in dnet.edge_spells
            if !isempty(spells)
                add_edge!(collapsed, edge[1], edge[2])
                total_time = sum(s.terminus - s.onset for s in spells)
                set_edge_attribute!(collapsed, edge[1], edge[2], :weight, total_time)
            end
        end

        return collapsed
    else
        throw(ArgumentError("method must be :union, :intersection, or :weighted"))
    end
end

end # module
