# Tutorials

## Installing this Package

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

## Accessing an Individual

A pedigree is an ordered dictionary where the key is the ID
and the value corresponds to the [`Individual`](@ref).

The individual's parents and children are accessed by reference.

```@repl
import GenLib as gen
genea140 = gen.genea140;
ped = gen.genealogy(genea140);
ped[33724]
ped[33724].mother
ped[33724].father
ped[33724].children
ped[33724].children[2].father
```

## Getting Founders and Probands

This is done using the [`founder`](@ref) and [`pro`](@ref) functions, respectively.

```@repl
import GenLib as gen
genea140 = gen.genea140;
ped = gen.genealogy(genea140);
founder = gen.founder(ped)
pro = gen.pro(ped)
```

## Finding Most Recent Common Ancestors

```@repl
import GenLib as gen
genea140 = gen.genea140;
ped = gen.genealogy(genea140);
pro = gen.pro(ped);
pro1 = pro[1]
pro2 = pro[2]
mrcas, meioses = gen.findMRCA(ped, [pro1, pro2]);
mrcas
meioses
```

## Computing Genetic Contributions

This is done with the [`gc`](@ref) function.

```@repl
import GenLib as gen
genea140 = gen.genea140;
ped = gen.genealogy(genea140);
contributions = gen.gc(ped)
sum(contributions, dims=2)
```

## Computing Kinship Coefficients

This is done using one of the [`phi`](@ref) functions.

### Pairwise Coefficient

Let's take the two siblings above for example.

```@example
import GenLib as gen
genea140 = gen.genea140;
ped = gen.genealogy(genea140);
pro1 = 10033
pro2 = 113470
pro1 = ped[pro1]
pro2 = ped[pro2]
gen.phi(pro1, pro2)
```

### Square Matrix

```@example
import GenLib as gen
geneaJi = gen.geneaJi;
ped = gen.genealogy(geneaJi);
gen.phi(ped)
```
