using ImageStacking
using Documenter

using_directives = :(using ImageStacking)
DocMeta.setdocmeta!(ImageStacking, :DocTestSetup, using_directives; recursive=true)

is_ci_env = (get(ENV, "CI", nothing) == true)
@info "is_ci_env == $is_ci_env"

makedocs(;
    sitename="ImageStacking.jl",
    modules=[ImageStacking],
    doctest=false,
    checkdocs = :exports,
    format=Documenter.HTML(;
        prettyurls = is_ci_env,
        canonical = "https://esitohi.github.io/ImageStacking.jl",
        edit_link = "main",
        assets = String[],
    ),
    pages=[
        "Home" => "index.md"
        "Stacking methods" => "stacking.md"
        "Normalization" => "normalization.md"
        #= TODO: create more pages
        "API" => Any[
            "Stacking" => "api/stacking.md"
        ]
        =#
    ],
)

deploydocs(;
    repo="github.com/esitohi/ImageStacking.jl.git",
    devbranch="main",
)
