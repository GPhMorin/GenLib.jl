# GenLib.jl

[![Build Status](https://github.com/GPhMorin/GenLib.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/GPhMorin/GenLib.jl/actions/workflows/CI.yml?query=branch%3Amain)

An unofficial, pure Julia port of R's GENLIB genetics and genealogical library.

To install, using REPL: `using Pkg; Pkg.add(url="https://github.com/GPhMorin/GenLib.jl")` or `add https://github.com/GPhMorin/GenLib.jl` in the Pkg REPL mode. You may then use the library with `using GenLib` or `import GenLib`.

If you use this project in a publication, please cite the original R package:

> Gauvin, H., Lefebvre, JF., Moreau, C. *et al.* GENLIB: an R package for the analysis of genealogical data. *BMC Bioinformatics* **16**, 160 (2015). https://doi.org/10.1186/s12859-015-0581-5

```bibtex
@article{gauvin2015genlib,
  title={GENLIB: an R package for the analysis of genealogical data},
  author={Gauvin, H{\'e}lo{\"\i}se and Lefebvre, Jean-Fran{\c{c}}ois and Moreau, Claudia and Lavoie, Eve-Marie and Labuda, Damian and V{\'e}zina, H{\'e}l{\`e}ne and Roy-Gagnon, Marie-H{\'e}l{\`e}ne},
  journal={BMC bioinformatics},
  volume={16},
  number={1},
  pages={1--10},
  year={2015},
  publisher={BioMed Central}
}
```

If you are a developer, you may clone this repository as `path/to/home/.julia/dev/GenLib/` and using: `import Pkg; Pkg.activate("path/to/home/.julia/dev/GenLib/")`. You may then use your local installation with the normal `using GenLib` or `import GenLib`.
