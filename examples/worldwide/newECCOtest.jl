# For run the script


using Pkg; Pkg.activate("/Users/yysong/Documents/julia/Julia_Project_Lagrangian/")
Pkg.develop(path="/Users/yysong/Desktop/Githubs/Drifters.jl")
using Drifters, MeshArrays, Climatology, MITgcm
ECCOmodule = Drifters.ECCO
CSV = Drifters.CSV
DataFrames = Drifters.DataFrames
#include("examples/worldwide/newECCOtest.jl")
out_path=joinpath("/Users/yysong/Desktop/data/ENSO/Lagrangian/Tropical_1998Jan_backward/")
!ispath(out_path) ? mkdir(out_path) : nothing
kk=1

#I=main_comp_SO(kk)



#function main_comp_SO(kk)

println(kk)

## Parameters

file_base = "initial_1_1"
backward_time = true
backward_time ? file_base=file_base*"_◀◀" : file_base=file_base*"_▶▶"

k=0
np=30000
ny=2 #number of years
nm=12 #number of months
dpth="/Users/yysong/Desktop/data/OCCA2/OCCA2HR1_diags"

## Set up FlowFields

#P,D=ECCOmodule.init_FlowFields(k=k, backward_time=backward_time,
 #       time_unit=:DateTime, dpth=dpth, datasets=:OCCA2, climatology=false, year0=1980)
P,D=ECCOmodule.init_FlowFields(k=k, backward_time=backward_time,time_unit=:DateTime, dpth=dpth, datasets=:OCCA2)

## Initialize individuals

#df=Drifters.init.init_regional_3d(np, D, lons=(-180.0,180.0), lats=(-30.0,30.0), zs=1:25)
df=Drifters.init.initial_positions(D.Γ, np, (-180.0,180.0), (-30.0,30.0), 1:25)
#df = Drifters.init.init_regional_3d(np , D, zs=1:25)  # or Gulf Stream defaults: omit kwargs
S = ECCOmodule.init_storage(np, 100, length(D.Γ.RC), 50)
I = Individuals(P, df.x, df.y, df.z, df.fid,
       (D=merge(D,S), ∫=ECCOmodule.custom∫,
        🔧=ECCOmodule.custom🔧,
        🔴=deepcopy(ECCOmodule.custom🔴)))

my∫! = ECCOmodule.custom∫!

## Initial integration (first half-month from 1998-01-15)

T=Drifters.DateTime(1998,1,1)
D.🔄(P, D, T)
## Monthly simulation loop (2 years = 24 months)

function step!(I::Individuals)
        I.D.🔄(I)
        ECCOmodule.custom∫!(I)
end

    
[step!(I) for y=1:ny, m=1:nm]

## Save output

#file_output=joinpath(out_path, file_base*".csv")
#CSV.write(file_output, I.🔴)

#return I
#end
