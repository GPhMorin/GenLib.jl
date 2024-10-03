var documenterSearchIndex = {"docs":
[{"location":"tutorials/#Tutorials","page":"Tutorials","title":"Tutorials","text":"","category":"section"},{"location":"tutorials/#Installing-this-Package","page":"Tutorials","title":"Installing this Package","text":"","category":"section"},{"location":"tutorials/","page":"Tutorials","title":"Tutorials","text":"To install, using REPL:","category":"page"},{"location":"tutorials/","page":"Tutorials","title":"Tutorials","text":"julia> using Pkg; Pkg.add(\"GenLib\")","category":"page"},{"location":"tutorials/","page":"Tutorials","title":"Tutorials","text":"Or in the Pkg REPL mode (]):","category":"page"},{"location":"tutorials/","page":"Tutorials","title":"Tutorials","text":"add GenLib","category":"page"},{"location":"tutorials/","page":"Tutorials","title":"Tutorials","text":"You may then use the library with using GenLib or, to mimic the behavior of R's GENLIB, import GenLib as gen.","category":"page"},{"location":"tutorials/#Loading-a-Pedigree","page":"Tutorials","title":"Loading a Pedigree","text":"","category":"section"},{"location":"tutorials/#From-a-DataFrame","page":"Tutorials","title":"From a DataFrame","text":"","category":"section"},{"location":"tutorials/","page":"Tutorials","title":"Tutorials","text":"Here is a pedigree with full first-degree cousins.","category":"page"},{"location":"tutorials/","page":"Tutorials","title":"Tutorials","text":"import GenLib as gen\nusing DataFrames\ninds = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]\nfathers = [0, 0, 0, 1, 1, 0, 3, 3, 6, 6]\nmothers = [0, 0, 0, 2, 2, 0, 4, 4, 5, 5]\nsexes = [1, 2, 1, 2, 2, 1, 2, 1, 1, 2]\ndf = DataFrame([inds, fathers, mothers, sexes], [:ind, :father, :mother, :sex])\nped = gen.genealogy(df)","category":"page"},{"location":"tutorials/#From-a-CSV-File","page":"Tutorials","title":"From a CSV File","text":"","category":"section"},{"location":"tutorials/","page":"Tutorials","title":"Tutorials","text":"The original GENLIB package for R contains two sample pedigrees. They are also available in GenLib.jl as GenLib.genea140 and GenLib.geneaJi.","category":"page"},{"location":"tutorials/","page":"Tutorials","title":"Tutorials","text":"import GenLib as gen\ngenea140 = gen.genea140\nped = gen.genealogy(genea140)","category":"page"},{"location":"tutorials/","page":"Tutorials","title":"Tutorials","text":"import GenLib as gen\ngeneaJi = gen.geneaJi\nped = gen.genealogy(geneaJi)","category":"page"},{"location":"tutorials/#Accessing-an-Individual","page":"Tutorials","title":"Accessing an Individual","text":"","category":"section"},{"location":"tutorials/","page":"Tutorials","title":"Tutorials","text":"A pedigree is an ordered dictionary where the key is the ID and the value corresponds to the GenLib.Individual.","category":"page"},{"location":"tutorials/","page":"Tutorials","title":"Tutorials","text":"The individual's parents and children are accessed by reference.","category":"page"},{"location":"tutorials/","page":"Tutorials","title":"Tutorials","text":"import GenLib as gen\ngenea140 = gen.genea140;\nped = gen.genealogy(genea140);\nped[33724]\nped[33724].mother\nped[33724].father\nped[33724].children\nped[33724].children[2].father","category":"page"},{"location":"tutorials/#Getting-Founders-and-Probands","page":"Tutorials","title":"Getting Founders and Probands","text":"","category":"section"},{"location":"tutorials/","page":"Tutorials","title":"Tutorials","text":"This is done using the GenLib.founder and GenLib.pro functions, respectively.","category":"page"},{"location":"tutorials/","page":"Tutorials","title":"Tutorials","text":"import GenLib as gen\ngenea140 = gen.genea140;\nped = gen.genealogy(genea140);\nfounder = gen.founder(ped)\npro = gen.pro(ped)","category":"page"},{"location":"tutorials/#Finding-Most-Recent-Common-Ancestors","page":"Tutorials","title":"Finding Most Recent Common Ancestors","text":"","category":"section"},{"location":"tutorials/","page":"Tutorials","title":"Tutorials","text":"import GenLib as gen\ngenea140 = gen.genea140;\nped = gen.genealogy(genea140);\npro = gen.pro(ped);\npro1 = pro[1]\npro2 = pro[2]\ngenMatrix = gen.findMRCA(ped, [pro1, pro2]);\ngenMatrix.individuals\ngenMatrix.ancestors\ngenMatrix.meioses","category":"page"},{"location":"tutorials/#Computing-Genetic-Contributions","page":"Tutorials","title":"Computing Genetic Contributions","text":"","category":"section"},{"location":"tutorials/","page":"Tutorials","title":"Tutorials","text":"This is done with the GenLib.gc function.","category":"page"},{"location":"tutorials/","page":"Tutorials","title":"Tutorials","text":"import GenLib as gen\ngenea140 = gen.genea140;\nped = gen.genealogy(genea140);\ncontributions = gen.gc(ped)\nsum(contributions, dims=2)","category":"page"},{"location":"tutorials/#Computing-Kinship-Coefficients","page":"Tutorials","title":"Computing Kinship Coefficients","text":"","category":"section"},{"location":"tutorials/","page":"Tutorials","title":"Tutorials","text":"This is done using one of the GenLib.phi functions.","category":"page"},{"location":"tutorials/#Pairwise-Coefficient","page":"Tutorials","title":"Pairwise Coefficient","text":"","category":"section"},{"location":"tutorials/","page":"Tutorials","title":"Tutorials","text":"Let's take the two siblings above for example.","category":"page"},{"location":"tutorials/","page":"Tutorials","title":"Tutorials","text":"import GenLib as gen\ngenea140 = gen.genea140;\nped = gen.genealogy(genea140);\npro1 = ped[10033]\npro2 = ped[113470]\ngen.phi(pro1, pro2)","category":"page"},{"location":"tutorials/#Square-Matrix","page":"Tutorials","title":"Square Matrix","text":"","category":"section"},{"location":"tutorials/","page":"Tutorials","title":"Tutorials","text":"import GenLib as gen\ngeneaJi = gen.geneaJi\nped = gen.genealogy(geneaJi)\ngen.phi(ped)","category":"page"},{"location":"reference/#Reference","page":"Reference","title":"Reference","text":"","category":"section"},{"location":"reference/#Structures","page":"Reference","title":"Structures","text":"","category":"section"},{"location":"reference/","page":"Reference","title":"Reference","text":"GenLib.Individual\nGenLib.Pedigree\nGenLib.GenMatrix","category":"page"},{"location":"reference/#GenLib.Individual","page":"Reference","title":"GenLib.Individual","text":"struct Individual <: AbstractIndividual\n    ID::Int64\n    father::Union{Nothing, Individual}\n    mother::Union{Nothing, Individual}\n    children::Vector{Individual}\n    sex::Int64\n    rank::Int64\nend\n\nThe unit structure of a GenLib.Pedigree.\n\n\n\n\n\n","category":"type"},{"location":"reference/#GenLib.Pedigree","page":"Reference","title":"GenLib.Pedigree","text":"struct Pedigree{T}\n\nA minimal structure wrapping an OrderedDict with individuals accessed by ID.\n\n\n\n\n\n","category":"type"},{"location":"reference/#GenLib.GenMatrix","page":"Reference","title":"GenLib.GenMatrix","text":"struct GenMatrix\n    individuals::Vector{Int64}\n    ancestors::Vector{Int64}\n    meioses::Matrix{Int64}\nend\n\nA matrix that goes with individuals as rows and ancestors as columns.\n\n\n\n\n\n","category":"type"},{"location":"reference/#Available-data","page":"Reference","title":"Available data","text":"","category":"section"},{"location":"reference/","page":"Reference","title":"Reference","text":"GenLib.genea140\nGenLib.geneaJi\nGenLib.pop140","category":"page"},{"location":"reference/#GenLib.genea140","page":"Reference","title":"GenLib.genea140","text":"Genealogical information for 140 individuals from the Quebec Reference Sample.\n\nAccording to the R GENLIB documentation, genea140 corresponds to \"a genealogical corpus made of 41523 individuals from the province of Quebec, Canada. A total of 140 individuals have been sampled in seven sub-populations, listed in pop140, and their genealogies were reconstructed as far back as possible using the BALSAC population register and the Early Quebec Population Register.\"\n\nLoading the Dataset\n\nimport GenLib as gen\ngenea140 = gen.genea140\nped = gen.genealogy(genea140)\n\n\n\n\n\n","category":"constant"},{"location":"reference/#GenLib.geneaJi","page":"Reference","title":"GenLib.geneaJi","text":"A highly inbred pedigree.\n\nAccording to the R GENLIB documentation, geneaJi corresponds to \"a modified version of a pedigree of two Jicaque Indians studied by Chapman & Jacquard (1971).\"\n\nLoading the Dataset\n\nimport GenLib as gen\ngeneaJi = gen.geneaJi\nped = gen.genealogy(geneaJi)\n\n\n\n\n\n","category":"constant"},{"location":"reference/#GenLib.pop140","page":"Reference","title":"GenLib.pop140","text":"Population of origin of the 140 Quebec samples.\n\nAccording to the R GENLIB documentation, pop140 corresponds to \"140 individuals from the genealogical corpus from Quebec (…) sampled from 7 different populations from 5 regions: Quebec City, Montreal, Saguenay, North Shore, Gaspesia. In Gaspesia we find 3 different populations: French-Canadians, Acadians and Loyalists.\"\n\nLoading the Dataset\n\nimport GenLib as gen\npop140 = gen.pop140\npop = gen.population(pop140)\n\n\n\n\n\n","category":"constant"},{"location":"reference/#Create-a-pedigree","page":"Reference","title":"Create a pedigree","text":"","category":"section"},{"location":"reference/","page":"Reference","title":"Reference","text":"GenLib.genealogy","category":"page"},{"location":"reference/#GenLib.genealogy","page":"Reference","title":"GenLib.genealogy","text":"genealogy(dataframe::DataFrame; sort::Bool = true)\n\nReturn an ordered pedigree of individuals from a DataFrame.\n\nExample\n\nimport GenLib as gen\nusing DataFrames\ninds = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]\nfathers = [0, 0, 0, 1, 1, 0, 3, 3, 6, 6]\nmothers = [0, 0, 0, 2, 2, 0, 4, 4, 5, 5]\nsexes = [1, 2, 1, 2, 2, 1, 2, 1, 1, 2]\ndf = DataFrame([inds, fathers, mothers, sexes], [:ind, :father, :mother, :sex])\nped = gen.genealogy(df)\n\n\n\n\n\ngenealogy(filename::String)\n\nReturn an ordered pedigree of individuals from a CSV file.\n\nExample\n\nimport GenLib as gen\ngenea140 = gen.genea140\nped = gen.genealogy(genea140)\n\n\n\n\n\n","category":"function"},{"location":"reference/#Extract-a-pedigree","page":"Reference","title":"Extract a pedigree","text":"","category":"section"},{"location":"reference/","page":"Reference","title":"Reference","text":"GenLib.branching","category":"page"},{"location":"reference/#GenLib.branching","page":"Reference","title":"GenLib.branching","text":"branching(pedigree::Pedigree; pro::Union{Vector{Int64}, Nothing} = nothing,\nancestors::Union{Vector{Int64}, Nothing} = nothing)\n\nReturn a pedigree that filters individuals who are in the paths between select probands and ancestors.\n\n\n\n\n\n","category":"function"},{"location":"reference/#Output-a-pedigree","page":"Reference","title":"Output a pedigree","text":"","category":"section"},{"location":"reference/","page":"Reference","title":"Reference","text":"GenLib.genout","category":"page"},{"location":"reference/#GenLib.genout","page":"Reference","title":"GenLib.genout","text":"genout(pedigree::Pedigree, sorted::Bool = false)\n\nReturn a pedigree as a DataFrame.\n\nIf sorted is false (the default), then the individuals will appear in the same order as in the pedigree.\n\nIf sorted is true, then the individuals will appear in alphabetical ID order.\n\nimport GenLib as gen\ngeneaJi = gen.geneaJi\nped = gen.genealogy(geneaJi)\ngen.genout(ped)\n\n\n\n\n\n","category":"function"},{"location":"reference/#Identify-individuals","page":"Reference","title":"Identify individuals","text":"","category":"section"},{"location":"reference/","page":"Reference","title":"Reference","text":"GenLib.founder\nGenLib.pro\nGenLib.sibship\nGenLib.children\nGenLib.findFounders\nGenLib.findMRCA\nGenLib.ancestor","category":"page"},{"location":"reference/#GenLib.founder","page":"Reference","title":"GenLib.founder","text":"founder(pedigree::Pedigree)\n\nReturn a vector of founder IDs in alphabetical order.\n\nExample\n\nimport GenLib as gen\ngenea140 = gen.genea140\nped = gen.genealogy(genea140)\nfounders = gen.founder(ped)\n\n\n\n\n\n","category":"function"},{"location":"reference/#GenLib.pro","page":"Reference","title":"GenLib.pro","text":"pro(pedigree::Pedigree)\n\nReturn a vector of proband IDs in alphabetical order.\n\nExample\n\nimport GenLib as gen\ngenea140 = gen.genea140\nped = gen.genealogy(genea140)\nprobands = gen.pro(ped)\n\n\n\n\n\n","category":"function"},{"location":"reference/#GenLib.sibship","page":"Reference","title":"GenLib.sibship","text":"sibship(pedigree::Pedigree, IDs::Vector{Int64}; halfSibling = true)\n\nReturn the siblings of specified individuals.\n\nIf halfSibling is true, half-siblings will be included.\n\n\n\n\n\n","category":"function"},{"location":"reference/#GenLib.children","page":"Reference","title":"GenLib.children","text":"children(pedigree::Pedigree, ID::Int64)\n\nReturn the children of an individual.\n\n\n\n\n\n","category":"function"},{"location":"reference/#GenLib.findFounders","page":"Reference","title":"GenLib.findFounders","text":"findFounders(pedigree::Pedigree, IDs::Vector{Int64})\n\nReturn a vector of founders from whom the IDs descend.\n\n\n\n\n\n","category":"function"},{"location":"reference/#GenLib.findMRCA","page":"Reference","title":"GenLib.findMRCA","text":"findMRCA(pedigree::Pedigree, IDs::Vector{Int64})\n\nReturn a GenLib.GenMatrix of meioses between individuals and their most recent common ancestors (MRCAs).\n\nExample\n\nimport GenLib as gen\ngenea140 = gen.genea140\nped = gen.genealogy(genea140)\npro = gen.pro(ped)\npro1 = pro[1]\npro2 = pro[2]\ngenMatrix = gen.findMRCA(ped, [pro1, pro2])\n\n\n\n\n\n","category":"function"},{"location":"reference/#GenLib.ancestor","page":"Reference","title":"GenLib.ancestor","text":"ancestor(pedigree::Pedigree, ID::Int64)\n\nReturn a vector of an individual's ancestors.\n\n\n\n\n\nancestor(pedigree::Pedigree, IDs::Vector{Int64})\n\nReturn a vector of several individual's ancestors.\n\n\n\n\n\n","category":"function"},{"location":"reference/#Describe","page":"Reference","title":"Describe","text":"","category":"section"},{"location":"reference/","page":"Reference","title":"Reference","text":"GenLib.nomen\nGenLib.nowomen\nGenLib.noind\nGenLib.depth\nGenLib.completeness\nGenLib.rec\nGenLib.occ\nGenLib.findDistance","category":"page"},{"location":"reference/#GenLib.nomen","page":"Reference","title":"GenLib.nomen","text":"nomen(pedigree::Pedigree)\n\nReturn the number of men in the pedigree.\n\n\n\n\n\n","category":"function"},{"location":"reference/#GenLib.nowomen","page":"Reference","title":"GenLib.nowomen","text":"nowomen(pedigree::Pedigree)\n\nReturn the number of women in the pedigree.\n\n\n\n\n\n","category":"function"},{"location":"reference/#GenLib.noind","page":"Reference","title":"GenLib.noind","text":"noind(pedigree::Pedigree)\n\nReturn the number of individuals in the pedigree.\n\n\n\n\n\n","category":"function"},{"location":"reference/#GenLib.depth","page":"Reference","title":"GenLib.depth","text":"depth(pedigree::Pedigree)\n\nReturn the maximum depth of a pedigree.\n\n\n\n\n\n","category":"function"},{"location":"reference/#GenLib.completeness","page":"Reference","title":"GenLib.completeness","text":"completeness(pedigree::Pedigree; pro::Vector{Int64} = pro(pedigree),\ngenNo::Vector{Int64} = Int64[], type::String = \"MEAN\")\n\nReturn a dataframe with the completeness at each generation (one row per generation).\n\ngenNo: A vector of the generations to output. The probands are at generation 0.\n\ntype: If \"MEAN\", the mean completeness for each generation. If \"IND\", the completeness for each generation for each proband.\n\n\n\n\n\n","category":"function"},{"location":"reference/#GenLib.rec","page":"Reference","title":"GenLib.rec","text":"rec(pedigree::Pedigree, probandIDs::Vector{Int64} = pro(genealogy),\nancestorIDs::Vector{Int64} = founder(genealogy))\n\nReturn the number of descendants of each ancestor.\n\n\n\n\n\n","category":"function"},{"location":"reference/#GenLib.occ","page":"Reference","title":"GenLib.occ","text":"occ(pedigree::Pedigree; pro::Vector{Int64} = pro(genealogy),\nancestors::Vector{Int64} = founder(genealogy), typeOcc::String = \"IND\")\n\nReturn a matrix of ancestors' occurrences.\n\nIf typeOcc is \"IND\" (default), then the matrix corresponds to the occurrence per individual. If typeOcc is \"TOTAL\", then the matrix corresponds to the total occurrence.\n\nExample\n\nimport GenLib as gen\ngeneaJi = gen.geneaJi\nped = gen.genealogy(geneaJi)\nocc = gen.occ(ped, typeOcc = \"TOTAL\")\n\n\n\n\n\n","category":"function"},{"location":"reference/#GenLib.findDistance","page":"Reference","title":"GenLib.findDistance","text":"findDistance(pedigree::Pedigree, IDs::Vector{Int64}, ancestorID::Int64)\n\nReturn the distance between two individuals and their ancestor.\n\n\n\n\n\n","category":"function"},{"location":"reference/#Compute","page":"Reference","title":"Compute","text":"","category":"section"},{"location":"reference/","page":"Reference","title":"Reference","text":"GenLib.gc\nGenLib.phi\nGenLib.phiMean\nGenLib.f","category":"page"},{"location":"reference/#GenLib.gc","page":"Reference","title":"GenLib.gc","text":"gc(pedigree::Pedigree; pro::Vector{Int64} = pro(pedigree),\n    ancestors::Vector{Int64} = founder(pedigree))\n\nReturn a matrix of the genetic contribution of each ancestor (columns) to each proband (rows).\n\nExample\n\nimport GenLib as gen\ngenea140 = gen.genea140\nped = gen.genealogy(genea140)\ncontributions = gen.gc(ped)\n\n\n\n\n\n","category":"function"},{"location":"reference/#GenLib.phi","page":"Reference","title":"GenLib.phi","text":"phi(individualᵢ::Individual, individualⱼ::Individual)\n\nReturn the kinship coefficient between two individuals.\n\nAdapted from Karigl, 1981.\n\nExample\n\nimport GenLib as gen\ngenea140 = gen.genea140\nped = gen.genealogy(genea140)\npro1 = ped[10033]\npro2 = ped[113470]\ngen.phi(pro1, pro2)\n\n\n\n\n\nphi(individualᵢ::IndexedIndividual, individualⱼ::IndexedIndividual, Ψ::Matrix{Float64})\n\nReturn the kinship coefficient between two individuals given a matrix of the founders' kinships.\n\nAdapted from Karigl, 1981, and Kirkpatrick et al., 2019.\n\n\n\n\n\nphi(pedigree::Pedigree, probandIDs::Vector{Int64} = pro(pedigree);\n    verbose::Bool = false)\n\nReturn a square matrix of pairwise kinship coefficients between probands.\n\nIf no probands are given, return the square matrix for all probands in the pedigree.\n\nThe algorithm computes pairwise kinships in parallel, a hybrid between the algorithms of Karigl, 1981, and Kirkpatrick et al., 2019.\n\nIf verbose is true, print the information about the cut vertices. If compute is true (the default), compute the kinship matrix. If it is false, only print the information about the cut vertices.\n\nExample\n\nimport GenLib as gen\ngeneaJi = gen.geneaJi\nped = gen.genealogy(geneaJi)\ngen.phi(ped)\n\n\n\n\n\nphi(individualᵢ::IndexedIndividual, individualⱼ::IndexedIndividual,\nϕ::T) where T <: KinshipMatrix\n\nReturn the kinship coefficient between two individuals given a dictionary of the individuals' kinships.\n\nAdapted from Karigl, 1981, and Kirkpatrick et al., 2019.\n\n\n\n\n\nphi(individualᵢ::IndexedIndividual, individualⱼ::IndexedIndividual,\nϕ::T) where T <: KinshipMatrix\n\nReturn the kinship coefficient between two individuals given a dictionary of the individuals' kinships.\n\nAdapted from Karigl, 1981, and Kirkpatrick et al., 2019.\n\n\n\n\n\n","category":"function"},{"location":"reference/#GenLib.phiMean","page":"Reference","title":"GenLib.phiMean","text":"function phiMean(ϕ::Matrix{Float64})\n\nReturn the mean kinship from a given kinship matrix.\n\n\n\n\n\nfunction phiMean(ϕ::VectorKinshipMatrix)\n\nReturn the mean kinship from a given sparse kinship matrix.\n\n\n\n\n\nfunction phiMean(ϕ::DictKinshipMatrix)\n\nReturn the mean kinship from a given sparse kinship matrix.\n\n\n\n\n\n","category":"function"},{"location":"reference/#GenLib.f","page":"Reference","title":"GenLib.f","text":"f(pedigree::Pedigree, IDs::Vector{Int64})\n\nReturn the coefficients of inbreeding of a vector of individuals.\n\n\n\n\n\n","category":"function"},{"location":"bibliography/#Bibliography","page":"Bibliography","title":"Bibliography","text":"","category":"section"},{"location":"bibliography/#Gauvin-et-al.,-2015","page":"Bibliography","title":"Gauvin et al., 2015","text":"","category":"section"},{"location":"bibliography/","page":"Bibliography","title":"Bibliography","text":"Gauvin, H., Lefebvre, JF., Moreau, C. et al. GENLIB: an R package for the analysis of genealogical data. BMC Bioinformatics 16, 160 (2015). https://doi.org/10.1186/s12859-015-0581-5","category":"page"},{"location":"bibliography/#Karigl,-1981","page":"Bibliography","title":"Karigl, 1981","text":"","category":"section"},{"location":"bibliography/","page":"Bibliography","title":"Bibliography","text":"KARIGL, G. (1981), A recursive algorithm for the calculation of identity coefficients. Annals of Human Genetics, 45: 299-305. https://doi.org/10.1111/j.1469-1809.1981.tb00341.x","category":"page"},{"location":"bibliography/#Kirkpatrick-et-al.,-2019","page":"Bibliography","title":"Kirkpatrick et al., 2019","text":"","category":"section"},{"location":"bibliography/","page":"Bibliography","title":"Bibliography","text":"Brent Kirkpatrick, Shufei Ge, Liangliang Wang, Efficient computation of the kinship coefficients, Bioinformatics, Volume 35, Issue 6, March 2019, Pages 1002–1008, https://doi.org/10.1093/bioinformatics/bty725","category":"page"},{"location":"#GenLib.jl","page":"Home","title":"GenLib.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Tools for pedigree analysis.","category":"page"},{"location":"","page":"Home","title":"Home","text":"A pure Julia port of R's GENLIB genetics and genealogical library.","category":"page"},{"location":"#Package-Features","page":"Home","title":"Package Features","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Basic functions that use the same syntax as R's GENLIB, sometimes faster than the original implementations;\nA function that computes the most recent common ancestors several magnitudes faster than gen.findMRCA in R's GENLIB;\nThe fastest available implementation to compute pairwise kinship coefficients, based on the algorithm by Karigl, 1981;\nThe fastest available implementation to compute a square matrix of kinship coefficients, based on the algorithms by Karigl, 1981 and Kirkpatrick et al., 2019.","category":"page"},{"location":"","page":"Home","title":"Home","text":"The Tutorials explain how to get started using GenLib.","category":"page"},{"location":"","page":"Home","title":"Home","text":"See the Index for the complete list of documented functions and types.","category":"page"},{"location":"","page":"Home","title":"Home","text":"The Bibliography lists the sources used for implementing the algorithms.","category":"page"},{"location":"#main-index","page":"Home","title":"Index","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Pages = [\"reference.md\"]","category":"page"}]
}
