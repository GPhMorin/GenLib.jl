var documenterSearchIndex = {"docs":
[{"location":"release/#Release-Notes","page":"Release Notes","title":"Release Notes","text":"","category":"section"},{"location":"release/#Version-[v0.1.0](https://github.com/GPhMorin/GenLib.jl/releases/tag/v0.1.0)-2024-02-11","page":"Release Notes","title":"Version v0.1.0 - 2024-02-11","text":"","category":"section"},{"location":"release/","page":"Release Notes","title":"Release Notes","text":"Initial release.","category":"page"},{"location":"tutorials/#Tutorials","page":"Tutorials","title":"Tutorials","text":"","category":"section"},{"location":"tutorials/#Installing-this-Package","page":"Tutorials","title":"Installing this Package","text":"","category":"section"},{"location":"tutorials/","page":"Tutorials","title":"Tutorials","text":"To install, using REPL:","category":"page"},{"location":"tutorials/","page":"Tutorials","title":"Tutorials","text":"julia> using Pkg; Pkg.add(url=\"https://github.com/GPhMorin/GenLib.jl\")","category":"page"},{"location":"tutorials/","page":"Tutorials","title":"Tutorials","text":"Or in the Pkg REPL mode (]):","category":"page"},{"location":"tutorials/","page":"Tutorials","title":"Tutorials","text":"add https://github.com/GPhMorin/GenLib.jl","category":"page"},{"location":"tutorials/","page":"Tutorials","title":"Tutorials","text":"You may then use the library with using GenLib or, to mimic the behavior of R's GENLIB, import GenLib as gen.","category":"page"},{"location":"tutorials/#Loading-a-Pedigree","page":"Tutorials","title":"Loading a Pedigree","text":"","category":"section"},{"location":"tutorials/#From-a-DataFrame","page":"Tutorials","title":"From a DataFrame","text":"","category":"section"},{"location":"tutorials/","page":"Tutorials","title":"Tutorials","text":"Here is a pedigree with full first-degree cousins.","category":"page"},{"location":"tutorials/","page":"Tutorials","title":"Tutorials","text":"import GenLib as gen\nusing DataFrames\ninds = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]\nfathers = [0, 0, 0, 1, 1, 0, 3, 3, 6, 6]\nmothers = [0, 0, 0, 2, 2, 0, 4, 4, 5, 5]\nsexes = [1, 2, 1, 2, 2, 1, 2, 1, 1, 2]\ndf = DataFrame([inds, fathers, mothers, sexes], [:ind, :father, :mother, :sex])\nped = gen.genealogy(df)","category":"page"},{"location":"tutorials/#From-a-CSV-File","page":"Tutorials","title":"From a CSV File","text":"","category":"section"},{"location":"tutorials/","page":"Tutorials","title":"Tutorials","text":"The original GENLIB package for R contains two sample pedigrees. They are also available in GenLib.jl as genea140 and geneaJi.","category":"page"},{"location":"tutorials/","page":"Tutorials","title":"Tutorials","text":"import GenLib as gen\ngenea140 = gen.genea140\nped = gen.genealogy(genea140)","category":"page"},{"location":"tutorials/","page":"Tutorials","title":"Tutorials","text":"import GenLib as gen\ngeneaJi = gen.geneaJi\nped = gen.genealogy(geneaJi)","category":"page"},{"location":"tutorials/#Accessing-an-Individual","page":"Tutorials","title":"Accessing an Individual","text":"","category":"section"},{"location":"tutorials/","page":"Tutorials","title":"Tutorials","text":"A pedigree is an ordered dictionary where the key is the ID and the value corresponds to the Individual.","category":"page"},{"location":"tutorials/","page":"Tutorials","title":"Tutorials","text":"The individual's parents and children are accessed by reference.","category":"page"},{"location":"tutorials/","page":"Tutorials","title":"Tutorials","text":"import GenLib as gen\ngenea140 = gen.genea140;\nped = gen.genealogy(genea140);\nped[33724]\nped[33724].mother\nped[33724].father\nped[33724].children\nped[33724].children[2].father","category":"page"},{"location":"tutorials/#Getting-Founders-and-Probands","page":"Tutorials","title":"Getting Founders and Probands","text":"","category":"section"},{"location":"tutorials/","page":"Tutorials","title":"Tutorials","text":"This is done using the founder and pro functions, respectively.","category":"page"},{"location":"tutorials/","page":"Tutorials","title":"Tutorials","text":"import GenLib as gen\ngenea140 = gen.genea140;\nped = gen.genealogy(genea140);\nfounder = gen.founder(ped)\npro = gen.pro(ped)","category":"page"},{"location":"tutorials/#Finding-Most-Recent-Common-Ancestors","page":"Tutorials","title":"Finding Most Recent Common Ancestors","text":"","category":"section"},{"location":"tutorials/","page":"Tutorials","title":"Tutorials","text":"import GenLib as gen\ngenea140 = gen.genea140;\nped = gen.genealogy(genea140);\npro = gen.pro(ped);\npro1 = pro[1]\npro2 = pro[2]\nmrcas, meioses = gen.findMRCA(ped, [pro1, pro2]);\nmrcas\nmeioses","category":"page"},{"location":"tutorials/#Computing-Genetic-Contributions","page":"Tutorials","title":"Computing Genetic Contributions","text":"","category":"section"},{"location":"tutorials/","page":"Tutorials","title":"Tutorials","text":"This is done with the gc function.","category":"page"},{"location":"tutorials/","page":"Tutorials","title":"Tutorials","text":"import GenLib as gen\ngenea140 = gen.genea140;\nped = gen.genealogy(genea140);\ncontributions = gen.gc(ped)\nsum(contributions, dims=2)","category":"page"},{"location":"tutorials/#Computing-Kinship-Coefficients","page":"Tutorials","title":"Computing Kinship Coefficients","text":"","category":"section"},{"location":"tutorials/","page":"Tutorials","title":"Tutorials","text":"This is done using one of the phi functions.","category":"page"},{"location":"tutorials/#Pairwise-Coefficient","page":"Tutorials","title":"Pairwise Coefficient","text":"","category":"section"},{"location":"tutorials/","page":"Tutorials","title":"Tutorials","text":"Let's take the two siblings above for example.","category":"page"},{"location":"tutorials/","page":"Tutorials","title":"Tutorials","text":"import GenLib as gen\ngenea140 = gen.genea140;\nped = gen.genealogy(genea140);\npro1 = 10033\npro2 = 113470\npro1 = ped[pro1]\npro2 = ped[pro2]\ngen.phi(pro1, pro2)","category":"page"},{"location":"tutorials/#Square-Matrix","page":"Tutorials","title":"Square Matrix","text":"","category":"section"},{"location":"tutorials/","page":"Tutorials","title":"Tutorials","text":"import GenLib as gen\ngeneaJi = gen.geneaJi;\nped = gen.genealogy(geneaJi);\ngen.phi(ped)","category":"page"},{"location":"reference/#Reference","page":"Reference","title":"Reference","text":"","category":"section"},{"location":"reference/","page":"Reference","title":"Reference","text":"genea140\ngeneaJi\npop140\nIndividual\nSTATE\nbranching\nchildren\nfindDistance(::OrderedDict{Int64, Individual}, ::Vector{Int64}, ::Int64)\nfindFounders\nfindMRCA\nfounder\ngc\ngenealogy\ngenout\nocc\nphi\npro\nrec","category":"page"},{"location":"reference/#GenLib.genea140","page":"Reference","title":"GenLib.genea140","text":"Genealogical information for 140 individuals from the Quebec Reference Sample.\n\nAccording to the R GENLIB documentation, genea140 corresponds to \"a genealogical corpus made of 41523 individuals from the province of Quebec, Canada. A total of 140 individuals have been sampled in seven sub-populations, listed in pop140, and their genealogies were reconstructed as far back as possible using the BALSAC population register and the Early Quebec Population Register.\"\n\n\n\n\n\n","category":"constant"},{"location":"reference/#GenLib.geneaJi","page":"Reference","title":"GenLib.geneaJi","text":"A highly inbred pedigree.\n\nAccording to the R GENLIB documentation, geneaJi corresponds to \"a modified version of a pedigree of two Jicaque Indians studied by Chapman & Jacquard (1971).\"\n\n\n\n\n\n","category":"constant"},{"location":"reference/#GenLib.pop140","page":"Reference","title":"GenLib.pop140","text":"Population of origin of the 140 Quebec samples.\n\nAccording to the R GENLIB documentation, pop140 corresponds to \"140 individuals from the genealogical corpus from Quebec (…) sampled from 7 different opulations from 5 regions: Quebec City, Montreal, Saguenay, North Shore, Gaspesia. In Gaspesia we find 3 different populations: French-Canadians, Acadians and Loyalists.\"\n\n\n\n\n\n","category":"constant"},{"location":"reference/#GenLib.Individual","page":"Reference","title":"GenLib.Individual","text":"mutable struct Individual\n    ID::Int64\n    father::Union{Nothing, Individual}\n    mother::Union{Nothing, Individual}\n    index::Int64\n    children::Vector{Individual}\n    sex::Int64\n    state::STATE\n    probability::Float64\n    sort::Int64\n    ancestor::Bool\n    descendant::Bool\n    occurrence::Int64\nend\n\nThe unit structure of a pedigree.\n\n\n\n\n\n","category":"type"},{"location":"reference/#GenLib.STATE","page":"Reference","title":"GenLib.STATE","text":"@enum STATE begin\n    PROBAND\n    FOUNDER\n    UNEXPLORED\nend\n\nAn enumeration of the Individual's state.\n\n\n\n\n\n","category":"type"},{"location":"reference/#GenLib.branching","page":"Reference","title":"GenLib.branching","text":"branching(pedigree::OrderedDict{Int64, Individual}; pro::Vector{Int64} = pro(pedigree), ancestors::Vector{Int64} = founder(pedigree))\n\nReturn a pedigree that filters individuals who are in the paths between select probands and ancestors.\n\n\n\n\n\n","category":"function"},{"location":"reference/#GenLib.children","page":"Reference","title":"GenLib.children","text":"children(pedigree::OrderedDict{Int64, Individual}, ID::Int64)\n\nReturn the children of an individual.\n\n\n\n\n\n","category":"function"},{"location":"reference/#GenLib.findDistance-Tuple{OrderedDict{Int64, Individual}, Vector{Int64}, Int64}","page":"Reference","title":"GenLib.findDistance","text":"findDistance(pedigree::OrderedDict{Int64, Individual}, IDs::Vector{Int64}, ancestorID::Int64)\n\nReturn the distance between two individuals and their ancestor.\n\n\n\n\n\n","category":"method"},{"location":"reference/#GenLib.findFounders","page":"Reference","title":"GenLib.findFounders","text":"findFounders(pedigree::OrderedDict{Int64}, IDs::Vector{Int64})\n\nReturn a vector of founders from whom the IDs descend.\n\n\n\n\n\n","category":"function"},{"location":"reference/#GenLib.findMRCA","page":"Reference","title":"GenLib.findMRCA","text":"findMRCA(pedigree::OrderedDict{Int64, Individual}, IDs::Vector{Int64})\n\nReturn a tuple of type ::Tuple{::Vector{Int64}, ::Matrix{Int64}} consisting in the vector of individuals' most recent common ancestors (MRCAs) and the matrix of meioses between each proband and each MRCA.\n\n\n\n\n\n","category":"function"},{"location":"reference/#GenLib.founder","page":"Reference","title":"GenLib.founder","text":"founder(pedigree::OrderedDict{Int64, Individual})\n\nReturn a vector of founder IDs in alphabetical order.\n\n\n\n\n\n","category":"function"},{"location":"reference/#GenLib.gc","page":"Reference","title":"GenLib.gc","text":"gc(pedigree::OrderedDict{Int64, Individual}; pro::Vector{Int64} = pro(pedigree), ancestors::Vector{Int64} = founder(pedigree))\n\nReturn a matrix of the genetic contribution of each ancestor (columns) to each proband (rows).\n\n\n\n\n\n","category":"function"},{"location":"reference/#GenLib.genealogy","page":"Reference","title":"GenLib.genealogy","text":"genealogy(dataframe::DataFrame)\n\nReturn an ordered pedigree of individuals from a DataFrame.\n\n\n\n\n\ngenealogy(filename::String)\n\nReturn an ordered pedigree of individuals from a CSV file.\n\n\n\n\n\n","category":"function"},{"location":"reference/#GenLib.genout","page":"Reference","title":"GenLib.genout","text":"genout(pedigree::OrderedDict{Int64, Individual}, sorted::Bool = false)\n\nReturn a pedigree as a DataFrame.\n\nIf sorted is false (the default), then the individuals will appear in the same order as in the pedigree.\n\nIf sorted is true, then the individuals will appear in alphabetical ID order.\n\n\n\n\n\n","category":"function"},{"location":"reference/#GenLib.occ","page":"Reference","title":"GenLib.occ","text":"occ(pedigree::OrderedDict{Int64, Individual}; pro::Vector{Int64} = pro(genealogy), ancestors::Vector{Int64} = founder(genealogy), typeOcc::String = \"IND\")\n\nReturn a matrix of ancestors' occurrences.\n\nIf typeOcc is \"IND\" (default), then the matrix corresponds to the occurrence per individual. If typeOcc is \"TOTAL\", then the matrix corresponds to the total occurrence.\n\n\n\n\n\n","category":"function"},{"location":"reference/#GenLib.phi","page":"Reference","title":"GenLib.phi","text":"phi(individualᵢ::Individual, individualⱼ::Individual, Ψ::Union{Nothing, Matrix{Float64}} = nothing)\n\nReturn the kinship coefficient between two individuals.\n\nA matrix of the founders' kinships may optionally be provided.\n\n\n\n\n\nphi(pedigree::OrderedDict{Int64, Individual}, rowIDs::Vector{Int64}, columnIDs::Vector{Int64})\n\nReturn a rectangle matrix of kinship coefficients between row IDs and column IDs. The kinship of someone with themself is replaced with their inbreeding.\n\n\n\n\n\nphi(pedigree::OrderedDict{Int64, Individual}, Ψ::Matrix{Float64})\n\nReturn a square matrix of pairwise kinship coefficients between all probands given the founders' kinships.\n\nAn implementation of the recursive-cut algorithm presented in Kirkpatrick et al., 2019.\n\n\n\n\n\nphi(pedigree::OrderedDict{Int64, Individual}, probandIDs::Vector{Int64} = pro(pedigree); verbose::Bool = false)\n\nReturn the square matrix of the pairwise kinship coefficients of a set of probands.\n\nIf no probands are given, return the square matrix for all probands in the pedigree.\n\nAn implementation of the recursive-cut algorithm presented in Kirkpatrick et al., 2019.\n\n\n\n\n\n","category":"function"},{"location":"reference/#GenLib.pro","page":"Reference","title":"GenLib.pro","text":"pro(pedigree::OrderedDict{Int64, Individual})\n\nReturn a vector of proband IDs in alphabetical order.\n\n\n\n\n\n","category":"function"},{"location":"reference/#GenLib.rec","page":"Reference","title":"GenLib.rec","text":"rec(pedigree::OrderedDict{Int64, Individual}, probandIDs::Vector{Int64} = pro(genealogy), ancestorIDs::Vector{Int64} = founder(genealogy))\n\nReturn the number of descendants of each ancestor.\n\n\n\n\n\n","category":"function"},{"location":"bibliography/#Bibliography","page":"Bibliography","title":"Bibliography","text":"","category":"section"},{"location":"bibliography/","page":"Bibliography","title":"Bibliography","text":"Gauvin, H., Lefebvre, JF., Moreau, C. et al. GENLIB: an R package for the analysis of genealogical data. BMC Bioinformatics 16, 160 (2015). https://doi.org/10.1186/s12859-015-0581-5","category":"page"},{"location":"bibliography/","page":"Bibliography","title":"Bibliography","text":"KARIGL, G. (1981), A recursive algorithm for the calculation of identity coefficients. Annals of Human Genetics, 45: 299-305. https://doi.org/10.1111/j.1469-1809.1981.tb00341.x","category":"page"},{"location":"bibliography/","page":"Bibliography","title":"Bibliography","text":"Brent Kirkpatrick, Shufei Ge, Liangliang Wang, Efficient computation of the kinship coefficients, Bioinformatics, Volume 35, Issue 6, March 2019, Pages 1002–1008, https://doi.org/10.1093/bioinformatics/bty725","category":"page"},{"location":"#GenLib.jl","page":"Home","title":"GenLib.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Tools for pedigree analysis.","category":"page"},{"location":"","page":"Home","title":"Home","text":"A pure Julia port of R's GENLIB genetics and genealogical library.","category":"page"},{"location":"#Package-Features","page":"Home","title":"Package Features","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Basic functions that use the same syntax as R's GENLIB, sometimes faster than the original implementations;\nA function that computes the most recent common ancestors several magnitudes faster than gen.findMRCA in R's GENLIB;\nThe fastest available implementation to compute pairwise kinship coefficients, based on the algorithm by Karigl, 1981;\nThe fastest available implementation to compute a square matrix of kinship coefficients, based on the algorithm by Kirkpatrick et al., 2019.","category":"page"},{"location":"","page":"Home","title":"Home","text":"The Tutorials explain how to get started using GenLib.","category":"page"},{"location":"","page":"Home","title":"Home","text":"See the Index for the complete list of documented functions and types.","category":"page"},{"location":"","page":"Home","title":"Home","text":"The Bibliography lists the sources used for implementing the algorithms.","category":"page"},{"location":"#main-index","page":"Home","title":"Index","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Pages = [\"reference.md\"]","category":"page"}]
}
