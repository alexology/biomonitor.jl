using Documenter
using Biomonitor

DocMeta.setdocmeta!(Biomonitor, :DocTestSetup, :(using Biomonitor); recursive=true)

makedocs(;
    modules=[Biomonitor],
    sitename="Biomonitor.jl",
    authors="Alex Laini",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://alexology.github.io/Biomonitor.jl",
        edit_link="main",
        size_threshold_warn=150 * 2^10,
    ),
    pages=[
        "Home" => "index.md",
        "From biomonitoR to Biomonitor.jl" => "r_to_julia.md",
        "Validation" => "validation.md",
        "API" => "api.md",
    ],
    checkdocs=:exports,
)

deploydocs(; repo="github.com/alexology/Biomonitor.jl", devbranch="main")
