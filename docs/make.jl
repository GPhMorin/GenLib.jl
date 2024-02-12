using Documenter, DataFrames, DataStructures, GenLib

makedocs(sitename="GenLib.jl",
         pages = [
            "Home" => "index.md",
            "tutorials.md",
            "reference.md",
            "bibliography.md"
            ],
         format = Documenter.HTML(
            prettyurls = get(ENV, "CI", nothing) == "true"
            )
        )
deploydocs(
   repo = "github.com/GPhMorin/GenLib.jl.git"
)