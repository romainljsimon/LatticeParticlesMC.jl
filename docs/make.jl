using LatticeParticlesMC
using Documenter

DocMeta.setdocmeta!(LatticeParticlesMC, :DocTestSetup, :(using LatticeParticlesMC); recursive=true)

makedocs(;
    modules=[LatticeParticlesMC],
    authors="Romain Simon, Leonardo Galliano",
    sitename="LatticeParticlesMC.jl",
    format=Documenter.HTML(;
        canonical="https://romainljsimon.github.io/LatticeParticlesMC.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/romainljsimon/LatticeParticlesMC.jl",
    devbranch="main",
)
