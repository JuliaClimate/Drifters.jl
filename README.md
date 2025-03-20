# Drifters.jl

[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaClimate.github.io/Drifters.jl/dev)
[![codecov](https://codecov.io/gh/JuliaClimate/Drifters.jl/branch/master/graph/badge.svg?token=qgMClXyzVB)](https://codecov.io/gh/JuliaClimate/Drifters.jl)
[![CI](https://github.com/JuliaClimate/Drifters.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/JuliaClimate/Drifters.jl/actions/workflows/ci.yml)

[![DOI](https://joss.theoj.org/papers/10.21105/joss.02813/status.svg)](https://doi.org/10.21105/joss.02813)
[![DOI](https://zenodo.org/badge/208676176.svg)](https://zenodo.org/badge/latestdoi/208676176)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/JuliaClimate/Drifters.jl/master)

**Drifters.jl** computes point displacements over a gridded domain. It is geared towards the analysis of Climate, Ocean, etc models (`Arakawa C-grids` are natively supported) and the simulation of material transports within the Earth System (e.g. plastics or planktons in the Ocean; dusts or chemicals in the Atmosphere). 

[<img src="https://user-images.githubusercontent.com/20276764/90925145-ca90cc80-e3be-11ea-8eed-559307dcb925.png" width="40%">](https://youtu.be/tsdf4fmYt1k) [<img src="https://user-images.githubusercontent.com/20276764/90924860-41799580-e3be-11ea-96bd-9a5784d00ecc.png" width="40%">](https://youtu.be/82HPnYBtoVo)

Inter-operability with popular climate model grids via [MeshArrays.jl](https://github.com/JuliaClimate/MeshArrays.jl) is an important aspect. The package can read and write individual displacement collection files, including those generated by the [MIT general circulation model](https://mitgcm.readthedocs.io/en/latest/?badge=latest). `Drifters`'s initial test suite is based on global ocean model simulations called [ECCO (v4r2)](https://eccov4.readthedocs.io/en/latest/) and [CBIOMES (alpha)](https://cbiomes.readthedocs.io/en/latest/) (see [Forget et al. 2015](https://doi.org/10.5194/gmd-8-3071-2015)).

### Installation

```
using Pkg
Pkg.add("Drifters")
Pkg.test("Drifters")
```

### Examples

The movies highlighted above are based on the [Global Climatology](https://juliaclimate.github.io/Drifters.jl/dev/examples/global_ocean_circulation.html) and [Three Dimensional](https://juliaclimate.github.io/Drifters.jl/dev/examples/three_dimensional_ocean.html) examples. 

The basic examples shown below are from the [Simple 3-D Flow](https://juliaclimate.github.io/Drifters.jl/dev/examples/solid_body_rotation.html) and [Simple 2-D Flow](https://juliaclimate.github.io/Drifters.jl/dev/examples/random_flow_field.html) notebooks. For more, please report to [the docs](https://JuliaClimate.github.io/Drifters.jl/dev).

[<img src="https://user-images.githubusercontent.com/20276764/84766999-b801ad80-af9f-11ea-922a-610ad8a257dc.png"  width="45%">](https://www.youtube.com/watch?v=W5DNqJG9jt0) <img src="https://user-images.githubusercontent.com/20276764/94491485-595ee900-01b6-11eb-95e6-c2cacb812f46.png" width="30%"> 

<img src="https://github.com/JuliaClimate/Drifters.jl/raw/master/examples/figs/RandomFlow.gif" width="40%"> <img src="https://user-images.githubusercontent.com/20276764/94755085-4b010080-0361-11eb-94fd-ff1b68b99f94.png" width="40%"> 

### Would like to contribute?

Thank you for your interest in contributing to this package! 

Please start by filing a ticket in the GitHub issue tracker if you think you’ve found a bug, found a potential fix for a bug, would like to suggest additional features, would like guidance on how to add features, or could help improve the documentation for example. 

You can also send a GitHub pull request directly but filing a ticket first is often preferable e.g. to avoid duplicating efforts or confusion. For an example of how to work with pull requests, please refer to the [MITgcm documentation](https://mitgcm.readthedocs.io/en/latest/contributing/contributing.html).
 

