using TOML

# Reconstruct local source paths for tests, documentation and benchmarks.
const package_dir = dirname(@__DIR__)
const workspace_dir = dirname(package_dir)
const visited = Set{String}()
function checkout_sources(path)
    path = dirname(abspath(joinpath(path, "Project.toml")))
    path in visited && return
    push!(visited, path)
    for environment in (path, joinpath(path, "docs"), joinpath(path, "benchmark"), joinpath(path, "test"))
        project_file = joinpath(environment, "Project.toml")
        isfile(project_file) || continue
        project = TOML.parsefile(project_file)
        for source in values(get(project, "sources", Dict()))
            haskey(source, "path") || continue
            target = dirname(abspath(joinpath(environment, source["path"], "Project.toml")))
            target == path && continue
            dirname(target) == workspace_dir ||
                error("Local dependency must be a sibling checkout: $target")
            if !isdir(target)
                repo = "https://github.com/statistical-network-analysis-with-Julia/$(basename(target)).git"
                run(`git clone --depth 1 --quiet $repo $target`)
            end
            checkout_sources(target)
        end
    end
end
checkout_sources(package_dir)
