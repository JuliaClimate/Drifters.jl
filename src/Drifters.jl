module Drifters

using MeshArrays, CyclicArrays, OrdinaryDiffEq, DataFrames, Random
import NetCDF, CSV

function data_path end
function read_data_ECCO end

"""
    _SDEProblem(a,b,c,d)
    
Placeholder for the `StochasticDiffEq` extension, which
calls `StochasticDiffEq.SDEProblem(f,g,u₀,tspan)`.

```
using Drifters, StochasticDiffEq
?Drifters._SDEProblem
```
"""
function _SDEProblem(a,b,c,d)
    @warn "requires using StochasticDiffEq"
    error("Placeholder function for _SDEProblem should not be used")
end

include("API.jl")
include("compute.jl")
include("data_wrangling.jl")
include("toy_models.jl")
include("initial_positions.jl")
include("Downloads.jl")

include("example_ECCO.jl")
include("example_OCCA.jl")
include("example_GOM.jl")
include("example_Oscar.jl")

export Individuals, ∫!, solve!, DataFrame, groupby
export FlowFields, convert_to_FlowFields, to_C_grid!
export uvwArrays, uvArrays, uvwMeshArrays, uvMeshArrays
export dxdt!, dxy_dt_CyclicArray, dxy_dt_replay
export postprocess_MeshArray, add_lonlat!, postprocess_xy, interp_to_xy
export nearest_to_xy, randn_lonlat, interp_to_lonlat
export gcdist, stproj, stproj_inv
export random_flow_field, vortex_flow_field
export DriftersDataset
InDiPlot=DriftersDataset; export InDiPlot #backward compatibility

__init__() = begin
    datadeps.__init__datadeps()
end

end # module
