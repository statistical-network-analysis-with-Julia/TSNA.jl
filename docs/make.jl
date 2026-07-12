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
    warnonly = [:missing_docs, :docs_block],
)

deploydocs(
    repo = "github.com/Statistical-network-analysis-with-Julia/TSNA.jl.git",
    devbranch = "main",
    versions = [
        "stable" => "dev", # serve dev docs at /stable until a release is tagged
        "dev" => "dev",
    ],
    push_preview = true,
)
