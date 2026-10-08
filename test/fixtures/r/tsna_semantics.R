# Golden fixture for TSNA.jl: durations, formation/dissolution events, paths
# and per-vertex snapshot measures against R tsna / sna / networkDynamic.
#
#   Rscript test/fixtures/r/tsna_semantics.R > test/fixtures/tsna_semantics.toml
#
# Four groups of random dynamic networks, all with R's defaults:
#
# 1. [values.case_*], cases 1-60: one-mode networks on 6 vertices (directed and
#    undirected) observed over net.obs.period = [0, 12]. Edge activations have
#    integer onsets in -3..13 and durations in {0 (point spell), 1, 2, 3, 5,
#    Inf}, so spells are adjacent, overlapping, open-ended or start before the
#    window; some edges and vertices never receive a spell (R's active.default
#    = TRUE makes them active throughout). Recorded:
#    - tsna::edgeDuration (subject "edges" and "spells"), labelled by dyad;
#    - tsna::vertexDuration, and, computed here from networkDynamic's stored
#      vertex spells, the durations clipped to [0, 12] that TSNA documents;
#    - tsna::tEdgeFormation / tEdgeDissolution over t = 0..12, with
#      include.censored FALSE and TRUE;
#    - sna::degree (cmode "freeman"), betweenness, closeness and gden of
#      network.extract(nd, at = t), labelled by the original vertex ids;
#    - tsna::tiedDuration (duration and counts; out/in/combined) and
#      tsna::tEdgeDensity (duration per edge and per dyad, event per edge).
# 2. cases 61-80: the same over [0, 12] on two-mode networks (undirected and
#    directed) and on directed networks with self-loops, where tEdgeDensity's
#    dyad count is network.dyadcount's n1*n2 (2*n1*n2 directed) and n^2.
#    No snapshot measures (sna's measures do not take two-mode networks).
# 3. [values.nowin_*]: networks with NO net.obs.period -- one-mode, two-mode
#    and looped -- plus four fixed networks with an open spell starting at the
#    last change time, a dissolution at the first change time, and no finite
#    spell bound at all. Recorded: edgeDuration (Inf for open spells),
#    vertexDuration, the default tEdgeFormation/tEdgeDissolution series over
#    the closed range of the change times, tiedDuration, tEdgeDensity, and
#    tPath(direction = "fwd") from every vertex at the default start. tsna's
#    tPath is unreliable when spell times are negative (it warns so, and its
#    answers change when every time is shifted by a constant), so the paths
#    are computed on a copy shifted to non-negative times; the elapsed times
#    tdist do not depend on the shift.
# 4. [values.panel_*]: samplk-style panels (three waves on 18 vertices) built
#    with networkDynamic(network.list = ...): edgeDuration, the default event
#    series, tiedDuration, tEdgeDensity, tSnaStats(gden, grecip) at the waves
#    and tPath from vertex 1.
suppressMessages({library(tsna); library(sna)})
seed <- 20261003L
set.seed(seed)

fmt <- function(x) {
  if (is.na(x)) return("NaN")
  if (is.infinite(x)) return(if (x > 0) "Inf" else "-Inf")
  trimws(formatC(x, digits = 17, format = "g"))
}
q <- function(s) paste0('"', s, '"')
toml_strings <- function(v) if (length(v) == 0) "[]" else paste0("[", paste(q(v), collapse = ", "), "]")
toml_ints <- function(v) paste0("[", paste(v, collapse = ", "), "]")
# A bare TOML number: TOML spells the special values nan, inf and -inf.
tnum <- function(x) {
  if (is.na(x)) return("nan")
  if (is.infinite(x)) return(if (x > 0) "inf" else "-inf")
  fmt(x)
}
toml_nums <- function(v) paste0("[", paste(sapply(v, tnum), collapse = ", "), "]")
quiet <- function(expr) { invisible(capture.output(x <- suppressWarnings(suppressMessages(expr)))); x }

n <- 6L
n_onemode <- 60L
n_special <- 20L
n_nowin <- 40L
n_fixed <- 4L        # the fixed networks of group 3
n_panels <- 10L
times_vertex <- c(1, 4.5, 7, 10)

cat('name = "tsna_semantics"\n\n[provenance]\n')
cat(sprintf('r_version = "%s"\n', paste(R.version$major, R.version$minor, sep = ".")))
cat(sprintf('tsna_version = "%s"\n', as.character(packageVersion("tsna"))))
cat(sprintf('sna_version = "%s"\n', as.character(packageVersion("sna"))))
cat(sprintf('networkDynamic_version = "%s"\n', as.character(packageVersion("networkDynamic"))))
cat(sprintf("seed = %d\n", seed))
cat('script = "test/fixtures/r/tsna_semantics.R"\n')
cat(sprintf('date = "%s"\n', format(Sys.Date())))
cat(sprintf('dataset = "%d one-mode and %d two-mode/looped random dynamic networks on %d vertices with net.obs.period [0, 12]; %d networks without net.obs.period (one-mode, two-mode, looped, and 4 fixed edge cases); %d samplk-style 3-wave panels on 18 vertices"\n',
            n_onemode, n_special, n, n_nowin, n_panels))
cat("\n[tolerance]\n# Durations and counts are exact; sna's floating-point centralities are compared at 1e-10.\ndurations = 0.0\ncentrality = 1e-10\n\n[values]\n")
cat(sprintf("n_cases = %d\n", n_onemode + n_special))
cat(sprintf("n_onemode = %d\n", n_onemode))
cat(sprintf("n_nowin = %d\n", n_nowin + n_fixed))
cat(sprintf("n_panels = %d\n", n_panels))
cat(sprintf("times = [%s]\n", paste(times_vertex, collapse = ", ")))

## ---- shared helpers --------------------------------------------------------

# Candidate dyads of a network: within-mode pairs excluded on two-mode
# networks, i == j only with loops.
random_pairs <- function(n, directed, bip, loops, k) {
  pairs <- t(replicate(k, sample(1:n, 2, replace = loops)))
  if (bip > 0) {
    pairs <- cbind(sample(1:bip, k, replace = TRUE), sample((bip + 1):n, k, replace = TRUE))
    if (directed) pairs <- t(apply(pairs, 1, function(p) if (runif(1) < 0.5) rev(p) else p))
  }
  if (!directed) pairs <- t(apply(pairs, 1, sort))
  unique(pairs)
}

random_ops <- function(nw, n, pairs) {
  ops <- character(0)
  for (v in 1:n) {
    for (i in seq_len(sample(0:2, 1, prob = c(0.4, 0.4, 0.2)))) {
      on <- sample(-2:11, 1); te <- on + sample(c(1, 2, 4, 8, Inf), 1)
      nw <- activate.vertices(nw, onset = on, terminus = te, v = v)
      ops <- c(ops, sprintf("av %d 0 %s %s", v, fmt(on), fmt(te)))
    }
  }
  for (eid in seq_len(nrow(pairs))) {
    for (i in seq_len(sample(0:3, 1, prob = c(0.15, 0.35, 0.3, 0.2)))) {
      on <- sample(-3:13, 1); te <- on + sample(c(0, 1, 2, 3, 5, Inf), 1)
      nw <- activate.edges(nw, onset = on, terminus = te, e = eid)
      ops <- c(ops, sprintf("ae %d %d %s %s", pairs[eid, 1], pairs[eid, 2], fmt(on), fmt(te)))
    }
  }
  list(nw = nw, ops = ops)
}

lab_fun <- function(directed) function(a, b) if (directed) paste(a, b) else paste(min(a, b), max(a, b))

# edgeDuration (edges and spells), labelled by dyad, from tsna itself.
edge_durations <- function(nw, pairs, directed) {
  lab <- lab_fun(directed)
  del <- as.data.frame(nw)
  ed <- edgeDuration(nw)
  ids <- sort(unique(del$edge.id))
  ed_lab <- sapply(seq_along(ids), function(k) sprintf("%s %s", lab(pairs[ids[k], 1], pairs[ids[k], 2]), fmt(ed[k])))
  es <- edgeDuration(nw, subject = "spells")
  es_lab <- sapply(seq_len(nrow(del)), function(k)
    sprintf("%s %s %s", lab(del$tail[k], del$head[k]), fmt(del$onset[k]), fmt(es[k])))
  list(ed = ed_lab, es = es_lab)
}

vertex_durations <- function(nw) {
  vd <- vertexDuration(nw)
  vsl <- get.vertex.activity(nw, as.spellList = TRUE)
  vids <- sort(unique(vsl$vertex.id))
  sapply(seq_along(vids), function(k) sprintf("%d %s", vids[k], fmt(vd[k])))
}

tied_records <- function(nw, directed) {
  nb <- if (directed) c("out", "in", "combined") else "combined"
  list(nb = nb,
       dur = sapply(nb, function(k) paste(sapply(tiedDuration(nw, neighborhood = k), fmt), collapse = ",")),
       cnt = sapply(nb, function(k) paste(tiedDuration(nw, mode = "counts", neighborhood = k), collapse = ",")))
}

emit_common <- function(directed, bip, loops, pairs, ops) {
  cat(sprintf("directed = %s\n", if (directed) "true" else "false"))
  cat(sprintf("bipartite = %d\n", bip))
  cat(sprintf("loops = %s\n", if (loops) "true" else "false"))
  cat(sprintf("edges = %s\n", toml_strings(paste(pairs[, 1], pairs[, 2]))))
  cat(sprintf("ops = %s\n", toml_strings(ops)))
}

## ---- groups 1 and 2: observed over [0, 12] ---------------------------------

for (g in seq_len(n_onemode + n_special)) {
  special <- g > n_onemode
  if (!special) {
    directed <- g %% 3 != 0; bip <- 0L; loops <- FALSE
    nw <- network.initialize(n, directed = directed)
    pairs <- unique(t(replicate(8, sample(1:n, 2))))
    if (!directed) pairs <- unique(t(apply(pairs, 1, sort)))
  } else {
    kind <- (g - n_onemode) %% 4
    directed <- kind != 0; bip <- if (kind %in% c(0, 1)) sample(2:4, 1) else 0L
    loops <- kind %in% c(2, 3)
    nw <- network.initialize(n, directed = directed, bipartite = if (bip > 0) bip else FALSE, loops = loops)
    pairs <- random_pairs(n, directed, bip, loops, 9)
  }
  for (r in seq_len(nrow(pairs))) add.edge(nw, pairs[r, 1], pairs[r, 2])
  ro <- random_ops(nw, n, pairs); nw <- ro$nw; ops <- ro$ops
  nw %n% "net.obs.period" <- list(observations = list(c(0, 12)), mode = "discrete",
                                  time.increment = 1, time.unit = "step")
  network.vertex.names(nw) <- 1:n

  durs <- edge_durations(nw, pairs, directed)
  vd_lab <- vertex_durations(nw)
  # TSNA clips vertex spells to the window as edge spells are clipped; tsna's
  # vertexDuration only replaces infinite bounds. The clipped durations,
  # computed from networkDynamic's stored spells with its spells.overlap:
  vclip <- character(0); vinside <- character(0)
  for (v in 1:n) {
    sp <- get.vertex.activity(nw, v = v)[[1]]
    if (is.null(sp)) next
    fin <- c(sp[is.finite(sp)])      # every finite bound, kept or not
    keep <- sapply(seq_len(nrow(sp)), function(r) spells.overlap(c(0, 12), sp[r, ]))
    if (!any(keep)) next
    kept <- sp[keep, , drop = FALSE]
    vclip <- c(vclip, sprintf("%d %s", v, fmt(sum(pmin(kept[, 2], 12) - pmax(kept[, 1], 0)))))
    vinside <- c(vinside, sprintf("%d %s", v, if (all(fin >= 0 & fin <= 12)) "true" else "false"))
  }

  form <- as.integer(tEdgeFormation(nw, start = 0, end = 12))
  diss <- as.integer(tEdgeDissolution(nw, start = 0, end = 12))
  form_c <- as.integer(tEdgeFormation(nw, start = 0, end = 12, include.censored = TRUE))
  diss_c <- as.integer(tEdgeDissolution(nw, start = 0, end = 12, include.censored = TRUE))

  snap <- character(0)
  if (!special) {
    gm <- if (directed) "digraph" else "graph"
    for (t in times_vertex) {
      x <- network.extract(nw, at = t)
      k <- network.size(x)
      if (k == 0) { snap <- c(snap, sprintf("%s|||||", fmt(t))); next }
      ids_t <- as.integer(x %v% "vertex.names")
      dg <- degree(x, gmode = gm, cmode = if (directed) "freeman" else "degree")
      bt <- betweenness(x, gmode = gm)
      cl <- closeness(x, gmode = gm)
      gd <- gden(x, mode = gm)
      snap <- c(snap, sprintf("%s|%s|%s|%s|%s|%s", fmt(t), paste(ids_t, collapse = ","),
                              paste(sapply(dg, fmt), collapse = ","),
                              paste(sapply(bt, fmt), collapse = ","),
                              paste(sapply(cl, fmt), collapse = ","), fmt(gd)))
    }
  }

  tied <- tied_records(nw, directed)
  dens <- c(tEdgeDensity(nw), tEdgeDensity(nw, agg.unit = "dyad"), tEdgeDensity(nw, mode = "event"))

  cat(sprintf("\n[values.case_%d]\n", g))
  emit_common(directed, bip, loops, pairs, ops)
  cat(sprintf("edge_duration = %s\n", toml_strings(durs$ed)))
  cat(sprintf("spell_duration = %s\n", toml_strings(durs$es)))
  cat(sprintf("vertex_duration = %s\n", toml_strings(vd_lab)))
  cat(sprintf("vertex_duration_clipped = %s\n", toml_strings(vclip)))
  cat(sprintf("vertex_inside = %s\n", toml_strings(vinside)))
  cat(sprintf("formation = %s\n", toml_ints(form)))
  cat(sprintf("dissolution = %s\n", toml_ints(diss)))
  cat(sprintf("formation_censored = %s\n", toml_ints(form_c)))
  cat(sprintf("dissolution_censored = %s\n", toml_ints(diss_c)))
  cat(sprintf("snapshots = %s\n", toml_strings(snap)))
  cat(sprintf("tied_neighborhoods = %s\n", toml_strings(tied$nb)))
  cat(sprintf("tied_duration = %s\n", toml_strings(tied$dur)))
  cat(sprintf("tied_counts = %s\n", toml_strings(tied$cnt)))
  cat(sprintf("dyad_count = %s\n", tnum(network.dyadcount(nw))))
  cat(sprintf("edge_density = %s\n", toml_nums(dens)))
}

## ---- group 3: no net.obs.period --------------------------------------------

# The fixed networks: (a) a tie that forms at the last change time and stays
# open; (b) a left-censored tie that dissolves at the first change time; (c) a
# base edge with no spell record next to one active over (-Inf, Inf), so that
# there is no finite spell bound at all; (d) an open tie forming at the last
# change time, as an animation would see it.
fixed <- list(
  list(n = 3, directed = TRUE, pairs = rbind(c(1, 2), c(2, 3)),
       ops = c("ae 1 2 2 6", "ae 2 3 6 Inf")),
  list(n = 2, directed = TRUE, pairs = rbind(c(1, 2)),
       ops = c("ae 1 2 -Inf 5", "ae 1 2 8 10")),
  list(n = 3, directed = TRUE, pairs = rbind(c(1, 2), c(2, 3)),
       ops = c("ae 2 3 -Inf Inf")),
  list(n = 3, directed = FALSE, pairs = rbind(c(1, 2), c(2, 3)),
       ops = c("ae 1 2 0 5", "ae 2 3 5 Inf"))
)
stopifnot(length(fixed) == n_fixed)

for (g in seq_len(n_nowin + n_fixed)) {
  if (g <= n_nowin) {
    kind <- g %% 5
    directed <- kind != 0 && kind != 3
    bip <- if (kind == 3) sample(2:4, 1) else 0L
    loops <- kind == 4
    nn <- n
    nw <- network.initialize(nn, directed = directed, bipartite = if (bip > 0) bip else FALSE, loops = loops)
    pairs <- random_pairs(nn, directed, bip, loops, 8)
    for (r in seq_len(nrow(pairs))) add.edge(nw, pairs[r, 1], pairs[r, 2])
    ro <- random_ops(nw, nn, pairs); nw <- ro$nw; ops <- ro$ops
  } else {
    f <- fixed[[g - n_nowin]]
    directed <- f$directed; bip <- 0L; loops <- FALSE; nn <- f$n; pairs <- f$pairs; ops <- f$ops
    nw <- network.initialize(nn, directed = directed)
    for (r in seq_len(nrow(pairs))) add.edge(nw, pairs[r, 1], pairs[r, 2])
    for (op in ops) {
      p <- strsplit(op, " ")[[1]]
      eid <- get.edgeIDs(nw, v = as.integer(p[2]), alter = as.integer(p[3]))[1]
      nw <- activate.edges(nw, onset = as.numeric(p[4]), terminus = as.numeric(p[5]), e = eid)
    }
  }
  network.vertex.names(nw) <- 1:nn
  changes <- get.change.times(nw)

  durs <- edge_durations(nw, pairs, directed)
  vd_lab <- vertex_durations(nw)
  if (length(changes) > 0) {
    fs <- tEdgeFormation(nw)
    series_start <- start(fs)[1]
    form <- as.integer(fs)
    diss <- as.integer(tEdgeDissolution(nw))
    form_c <- as.integer(tEdgeFormation(nw, include.censored = TRUE))
    diss_c <- as.integer(tEdgeDissolution(nw, include.censored = TRUE))
  } else {
    series_start <- NA; form <- diss <- form_c <- diss_c <- integer(0)
  }
  tied <- quiet(tied_records(nw, directed))
  dens <- quiet(c(tEdgeDensity(nw), tEdgeDensity(nw, agg.unit = "dyad"), tEdgeDensity(nw, mode = "event")))
  start_path <- if (length(changes) > 0) min(changes) else 0
  shift <- max(0, -start_path)
  nw_paths <- network.initialize(nn, directed = directed, bipartite = if (bip > 0) bip else FALSE,
                                 loops = loops)
  for (r in seq_len(nrow(pairs))) add.edge(nw_paths, pairs[r, 1], pairs[r, 2])
  for (op in ops) {
    p <- strsplit(op, " ")[[1]]
    on <- as.numeric(p[4]) + shift; te <- as.numeric(p[5]) + shift
    if (p[1] == "av") {
      nw_paths <- activate.vertices(nw_paths, onset = on, terminus = te, v = as.integer(p[2]))
    } else {
      eid <- get.edgeIDs(nw_paths, v = as.integer(p[2]), alter = as.integer(p[3]))[1]
      nw_paths <- activate.edges(nw_paths, onset = on, terminus = te, e = eid)
    }
  }
  paths <- character(0)
  for (v in 1:nn) {
    tp <- try(quiet(tPath(nw_paths, v = v, direction = "fwd", start = start_path + shift)),
              silent = TRUE)
    if (inherits(tp, "try-error")) next
    paths <- c(paths, sprintf("%d|%s", v, paste(sapply(tp$tdist, fmt), collapse = ",")))
  }

  cat(sprintf("\n[values.nowin_%d]\n", g))
  cat(sprintf("n = %d\n", nn))
  emit_common(directed, bip, loops, pairs, ops)
  cat(sprintf("change_times = %s\n", toml_nums(changes)))
  cat(sprintf("edge_duration = %s\n", toml_strings(durs$ed)))
  cat(sprintf("spell_duration = %s\n", toml_strings(durs$es)))
  cat(sprintf("vertex_duration = %s\n", toml_strings(vd_lab)))
  cat(sprintf("series_start = %s\n", tnum(series_start)))
  cat(sprintf("formation = %s\n", toml_ints(form)))
  cat(sprintf("dissolution = %s\n", toml_ints(diss)))
  cat(sprintf("formation_censored = %s\n", toml_ints(form_c)))
  cat(sprintf("dissolution_censored = %s\n", toml_ints(diss_c)))
  cat(sprintf("tied_neighborhoods = %s\n", toml_strings(tied$nb)))
  cat(sprintf("tied_duration = %s\n", toml_strings(tied$dur)))
  cat(sprintf("tied_counts = %s\n", toml_strings(tied$cnt)))
  cat(sprintf("edge_density = %s\n", toml_nums(dens)))
  cat(sprintf("path_start = %s\n", tnum(start_path)))
  cat(sprintf("path_shift = %s\n", tnum(shift)))
  cat(sprintf("paths = %s\n", toml_strings(paths)))
}

## ---- group 4: samplk-style panels --------------------------------------------

next_wave <- function(prev, nn, directed, p_new) {
  m <- prev * (matrix(runif(nn * nn), nn) < 0.7)
  m <- pmax(m, matrix(runif(nn * nn), nn) < p_new)
  diag(m) <- 0
  if (!directed) { m[lower.tri(m)] <- 0; m <- pmax(m, t(m)) }
  m
}
edge_strings <- function(m, directed) {
  idx <- which(m == 1, arr.ind = TRUE)
  if (!directed) idx <- idx[idx[, 1] < idx[, 2], , drop = FALSE]
  if (nrow(idx) == 0) return(character(0))
  idx <- idx[order(idx[, 1], idx[, 2]), , drop = FALSE]
  paste(idx[, 1], idx[, 2])
}
for (g in seq_len(n_panels)) {
  nn <- 18L; directed <- g %% 4 != 0
  waves <- list(next_wave(matrix(0, nn, nn), nn, directed, 0.15))
  for (k in 2:3) waves[[k]] <- next_wave(waves[[k - 1]], nn, directed, 0.05)
  nets <- lapply(waves, function(m) network(m, directed = directed))
  nd <- quiet(networkDynamic(network.list = nets))
  el <- as.matrix.network.edgelist(nd)
  lab <- lab_fun(directed)
  del <- as.data.frame(nd)
  ed <- edgeDuration(nd)
  ids <- sort(unique(del$edge.id))
  ed_lab <- sapply(seq_along(ids), function(k) {
    e <- nd$mel[[ids[k]]]
    sprintf("%s %s", lab(e$outl, e$inl), fmt(ed[k]))
  })
  form <- as.integer(tEdgeFormation(nd)); diss <- as.integer(tEdgeDissolution(nd))
  tied <- tied_records(nd, directed)
  dens <- c(tEdgeDensity(nd), tEdgeDensity(nd, agg.unit = "dyad"))
  gd <- as.numeric(tSnaStats(nd, "gden"))
  gr <- as.numeric(tSnaStats(nd, "grecip"))
  tp <- tPath(nd, v = 1, start = 0)$tdist

  cat(sprintf("\n[values.panel_%d]\n", g))
  cat(sprintf("directed = %s\n", if (directed) "true" else "false"))
  cat(sprintf("n = %d\n", nn))
  for (k in 1:3) cat(sprintf("wave_%d = %s\n", k, toml_strings(edge_strings(waves[[k]], directed))))
  cat(sprintf("edge_duration = %s\n", toml_strings(ed_lab)))
  cat(sprintf("formation = %s\n", toml_ints(form)))
  cat(sprintf("dissolution = %s\n", toml_ints(diss)))
  cat(sprintf("tied_neighborhoods = %s\n", toml_strings(tied$nb)))
  cat(sprintf("tied_duration = %s\n", toml_strings(tied$dur)))
  cat(sprintf("tied_counts = %s\n", toml_strings(tied$cnt)))
  cat(sprintf("edge_density = %s\n", toml_nums(dens)))
  cat(sprintf("gden = %s\n", toml_nums(gd)))
  cat(sprintf("grecip = %s\n", toml_nums(gr)))
  cat(sprintf("tpath_1 = %s\n", toml_nums(tp)))
}
