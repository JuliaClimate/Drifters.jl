# Select testsets via env var, e.g.: JULIA_TESTSETS=time_axes,sde julia runtests_modular.jl
# Omit the env var (or set to "all") to run everything.
const _RUN = split(get(ENV, "JULIA_TESTSETS", "all"), ",")
include_maybe(name) = ("all" in _RUN || name in _RUN) && include("testsets/$name.jl")

#use this environment variable to bypass downloads / API calls that require server access
_SKIP_DOWNLOADS = parse(Bool,get(ENV, "SKIP_DOWNLOADS", "false"))

using Test, Documenter, Drifters, Suppressor, CairoMakie
import StochasticDiffEq, MITgcm, Climatology

import Drifters: MeshArrays, NetCDF, CSV, DataFrames, JLD2

skip_docs_testing=true

if !_SKIP_DOWNLOADS
    MITgcm.getdata("mitgcmsmall")
    Climatology.get_ecco_velocity_if_needed()
    Climatology.get_occa_velocity_if_needed()
    Climatology.get_ecco_variable_if_needed("THETA")
    Climatology.get_ecco_variable_if_needed("SALT")
    MeshArrays.GridLoad(MeshArrays.GridSpec(ID=:LLC90))
    MeshArrays.GridLoad(MeshArrays.GridSpec(ID=:onedegree))
end

include_maybe("time_axes")
include_maybe("sde")
include_maybe("oscar")
include_maybe("ecco")
include_maybe("occa")
include_maybe("simple")
include_maybe("downloads")
include_maybe("global")
include_maybe("various")
_SKIP_DOWNLOADS ? nothing : include_maybe("doctests")
