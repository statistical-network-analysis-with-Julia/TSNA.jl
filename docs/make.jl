using Documenter
using TSNA

DocMeta.setdocmeta!(TSNA, :DocTestSetup, :(using TSNA); recursive=true)

makedocs(
    sitename = "TSNA.jl",
    modules = [TSNA],
    authors = "Statistical Network Analysis with Julia",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical = "https://Statistical-network-analysis-with-Julia.github.io/TSNA.jl",
        edit_link = "main",
    ),
    repo = Documenter.Remotes.GitHub("Statistical-network-analysis-with-Julia", "TSNA.jl"),
    pages = [
        "Home" => "index.md",
        "Getting Started" => "getting_started.md",
        "User Guide" => [
            "Temporal Centrality" => "guide/centrality.md",
            "Temporal Paths" => "guide/paths.md",
            "Duration Metrics" => "guide/metrics.md",
        ],
        "API Reference" => [
            "Centrality" => "api/centrality.md",
            "Paths" => "api/paths.md",
            "Metrics" => "api/metrics.md",
        ],
    ],
    # STRICT. Undefined bindings, bad cross-references, duplicate docs and
    # malformed markdown are build ERRORS, so they cannot silently accumulate
    # again (a docs build that passes while warning is one that will rot).
    #
    # `checkdocs = :exports` is the one deliberate exclusion: every *exported*
    # name must be documented, but internal machinery (materialized/private
    # types, `Base`/`Graphs` method extensions, inner constructors) need not be
    # -- filler docstrings for names a user never types are worse than none.
    warnonly = false,
    checkdocs = :exports,
)

deploydocs(
    repo = "github.com/statistical-network-analysis-with-Julia/TSNA.jl.git",
    devbranch = "main",
    versions = [
        "stable" => "dev", # serve dev docs at /stable until a release is tagged
        "dev" => "dev",
    ],
    push_preview = true,
)
