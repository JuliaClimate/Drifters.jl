module IndividualDisplacements

using MeshArrays, CyclicArrays, OrdinaryDiffEq, DataFrames
using Random, Artifacts, LazyArtifacts

include("API.jl")
include("compute.jl")
include("data_wrangling.jl")
include("toy_models.jl")
include("various.jl")

DiffEqBase.solve!(I::Individuals,args...)=∫!(I::Individuals,args...)
DataFrames.groupby(I::Individuals,args...) = groupby(I.🔴,args...)
DataFrames.DataFrame(I::Individuals) = I.🔴

export Individuals, ∫!, solve!, DataFrame, groupby
export FlowFields, convert_to_FlowFields
export F_Array3D, F_Array2D, F_MeshArray3D, F_MeshArray2D
export dxdt!, dxy_dt_CyclicArray, dxy_dt_replay
export postprocess_MeshArray, add_lonlat!, postprocess_xy, interp_to_xy
export nearest_to_xy, randn_lonlat, interp_to_lonlat
export gcdist, stproj, stproj_inv

export random_flow_field, vortex_flow_field

end # module
