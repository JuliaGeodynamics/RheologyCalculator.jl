push!(LOAD_PATH, dirname(@__DIR__))

using Documenter
using DocumenterVitepress
using RheologyCalculator
using RheologyCalculator.RheologyModels

@info "Making documentation..."
makedocs(;
    sitename = "RheologyCalculator.jl",
    authors  = "Albert de Montserrat and Boris Kaus",
    modules  = [RheologyCalculator, RheologyCalculator.RheologyModels],
    format   = DocumenterVitepress.MarkdownVitepress(;
        repo = "github.com/JuliaGeodynamics/RheologyCalculator.jl",
        devbranch = "main",
        devurl = "dev",
        sidebar_drawer = true,
    ),
    warnonly = Documenter.except(:footnote),
    checkdocs = :exports,
    draft    = false,
    source   = "src",
    build    = "build",
    pages    = [
        "Home" => "index.md",
        "Composites" => "composites.md",
        "Rheology" => "rheology.md",
        "Solving" => "stress.md",
        "Elastic correction" => "strain_rate_correction.md",
        "Sparsity detection" => "sparsity_detection.md",
        "API" => "api.md",
    ],
)

DocumenterVitepress.deploydocs(;
    repo       = "github.com/JuliaGeodynamics/RheologyCalculator.jl",
    target     = joinpath(@__DIR__, "build"),
    branch     = "gh-pages",
    devbranch  = "main",
    push_preview = true,
)
