using Documenter, GenLib, DataFrames

makedocs(sitename="GenLib.jl",
         pages = [
            "Home" => "index.md",
            "tutorials.md",
            "reference.md",
            "bibliography.md",
            "release.md"
         ],
         format = Documenter.HTML(
            prettyurls = get(ENV, "CI", nothing) == "true"
            )
        )