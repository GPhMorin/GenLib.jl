# Tutorials

## Installation

To install, using REPL:

```julia-repl
julia> using Pkg; Pkg.add(url="https://github.com/GPhMorin/GenLib.jl")
```

Or in the Pkg REPL mode (`]`):

```pkg
add https://github.com/GPhMorin/GenLib.jl
```

You may then use the library with `using GenLib` or, to mimic
the behavior of R's GENLIB, `import GenLib as gen`.

## Loading a Pedigree

### From a DataFrame

Here is a pedigree with full first-degree cousins.

```@example
import GenLib as gen
using DataFrames
inds = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
fathers = [0, 0, 0, 1, 1, 0, 3, 3, 6, 6]
mothers = [0, 0, 0, 2, 2, 0, 4, 4, 5, 5]
sexes = [1, 2, 1, 2, 2, 1, 2, 1, 1, 2]
df = DataFrame([inds, fathers, mothers, sexes], [:ind, :father, :mother, :sex])
ped = gen.genealogy(df)
```

### From a CSV File

The original GENLIB package for R contains two sample pedigrees.
They are also available in GenLib.jl as [`genea140`](@ref) and [`geneaJi`](@ref).

```@example
import GenLib as gen
genea140 = gen.genea140
ped = gen.genealogy(genea140)
```

```@example
import GenLib as gen
geneaJi = gen.geneaJi
ped = gen.genealogy(geneaJi)
```