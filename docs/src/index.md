# GenLib.jl
*Tools for pedigree analysis.*

A pure Julia port of [R's GENLIB](https://cran.r-project.org/web/packages/GENLIB/index.html) genetics and genealogical library.

## Package Features
- Basic functions that use the same syntax as R's GENLIB, sometimes faster than the original implementations;
- A function that computes the most recent common ancestors several magnitudes faster than `gen.findMRCA` in R's GENLIB;
- The fastest available implementation to compute pairwise kinship coefficients, based on the algorithm by [Karigl, 1981](@ref);
- The fastest available implementation to compute a square matrix of kinship coefficients, based on the algorithm by [Kirkpatrick et al., 2019](@ref)

The [Tutorials](@ref) explain how to get started using GenLib.

See the [Index](@ref main-index) for the complete list of documented functions and types.

The [Bibliography](@ref) lists the sources used for implementing the algorithms.

### [Index](@id main-index)

```@index
Pages = ["reference.md"]
```