# GenLib.jl

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://gphmorin.github.io/GenLib.jl/dev)
[![Build Status](https://github.com/GPhMorin/GenLib.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/GPhMorin/GenLib.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![codecov](https://codecov.io/gh/GPhMorin/GenLib.jl/graph/badge.svg?token=3A5C7F4H87)](https://codecov.io/gh/GPhMorin/GenLib.jl)

*Tools for pedigree analysis.*

A pure Julia port of [R's GENLIB](https://cran.r-project.org/web/packages/GENLIB/index.html) genetics and genealogical library.

## Package Features

* Basic functions that use the same syntax as R's GENLIB, sometimes faster than the original implementations;
* A function that computes the most recent common ancestors several magnitudes faster than `gen.findMRCA` in R's GENLIB;
* The fastest available implementation to compute pairwise kinship coefficients, based on the algorithm by [Karigl, 1981](https://gphmorin.github.io/GenLib.jl/dev/bibliography/#Karigl,-1981);
* The fastest available implementation to compute a square matrix of kinship coefficients, based on the algorithms by [Karigl, 1981](https://gphmorin.github.io/GenLib.jl/dev/bibliography/#Karigl,-1981) and [Kirkpatrick et al., 2019](https://gphmorin.github.io/GenLib.jl/dev/bibliography/#Kirkpatrick-et-al.,-2019).

## Documentation

The [Tutorials](https://gphmorin.github.io/GenLib.jl/dev/tutorials/) explain how to get started using GenLib.

See the [Index](https://gphmorin.github.io/GenLib.jl/dev/#main-index) for the complete list of documented functions and types.

The [Bibliography](https://gphmorin.github.io/GenLib.jl/dev/bibliography/) lists the sources used for implementing the algorithms.


## Installation

```julia
julia> ]
pkg> add GenLib
```

## Examples

```julia
julia> import GenLib as gen

julia> ped = gen.genealogy(gen.geneaJi)
A pedigree with:
29 individuals;
44 parent-child relations;
15 men;
14 women;
3 subjects;
8 generations.

julia> gen.phi(ped)
3×3 Matrix{Float64}:
 0.591797   0.371094   0.0722656
 0.371094   0.591797   0.0722656
 0.0722656  0.0722656  0.535156

julia> gen.pro(ped)
3-element Vector{Int64}:
  1
  2
 29

julia> ped[1]
ind: 1
father: 4
mother: 3
sex: 1

julia> ped[1].father
ind: 4
father: 8
mother: 6
sex: 1

julia> ped[1].father.ID
4

julia> ped[1].father.children
2-element Vector{GenLib.Individual}:
 ind: 1
father: 4
mother: 3
sex: 1
 ind: 2
father: 4
mother: 3
sex: 1
```
