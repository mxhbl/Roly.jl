using Documenter
using Roly

makedocs(
    sitename = "Roly.jl",
    modules = [Roly],
    pages = [
        "Home" => "index.md",
        "Workflow" => "workflow.md",
        "Custom particle species" => "custom_species.md",
        "Orientation and twists" => "orientation.md",
        "Bounded and unbounded rules" => "growth.md",
        "API reference" => "api.md",
    ],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
    ),
    checkdocs = :none,
)

deploydocs(
    repo = "github.com/goodrichgroup/Roly.jl.git",
    devbranch = "main",
    push_preview = true,
)
