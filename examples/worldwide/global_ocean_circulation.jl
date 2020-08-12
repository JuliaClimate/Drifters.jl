
# # Global Ocean Simulation
#
#md # [![](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/JuliaClimate/IndividualDisplacements.jl/web1?filepath=docs/src/notebooks/global_ocean_circulation.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/notebooks/global_ocean_circulation.ipynb)
#
# Particles moving with climatological monthly mean flow at selected depth level
# (e.g. `k=10` for 95 m) based on an ocean state estimate (ECCO v4 r2 from https://ecco-group.org).
# For additional documentation e.g. see :
# [1](https://JuliaClimate.github.io/MeshArrays.jl/dev/),
# [2](https://JuliaClimate.github.io/IndividualDisplacements.jl/dev/),
# [3](https://docs.juliadiffeq.org/latest/solvers/ode_solve.html),
# [4](https://en.wikipedia.org/wiki/Displacement_(vector))
#
# [![simulated particle movie (5m)](https://user-images.githubusercontent.com/20276764/84766999-b801ad80-af9f-11ea-922a-610ad8a257dc.png)](https://youtu.be/W5DNqJG9jt0)

#nb # %% {"slideshow": {"slide_type": "subslide"}, "cell_type": "markdown"}
# ## 1. Get Software & Iput Files
#
# - packages + helper functions
# - grid and velocity files

using IndividualDisplacements, MeshArrays, OrdinaryDiffEq
using Statistics, DataFrames, MITgcmTools, OceanStateEstimation

p=dirname(pathof(IndividualDisplacements))
include(joinpath(p,"../examples/recipes_plots.jl"))
include(joinpath(p,"../examples/helper_functions.jl"))
get_grid_if_needed()

#velocity files:
pp="$p/../examples/nctiles_climatology"
q=dirname(pathof(OceanStateEstimation))
qq="$q/../examples/nctiles_climatology"
!isfile(pp*".csv") ? run(`cp $qq.csv $pp.csv`) : nothing
!isdir(pp) ? run(`mkdir $pp`) : nothing
!isdir(pp*"/UVELMASS") ? get_from_dataverse("UVELMASS",pp) : nothing
!isdir(pp*"/VVELMASS") ? get_from_dataverse("VVELMASS",pp) : nothing

#nb # %% {"slideshow": {"slide_type": "subslide"}, "cell_type": "markdown"}
# ## 2. Set Up Parameters & Inputs
#
# - depth level and duration
# - read grid variables
# - read & normalize velocities

k=10 #choice of vertical level
ny=10 #number of simulated years (20 for k>20)
r_reset = 0.01 #fraction of the particles reset per month (0.05 for k<=10)

#read grid and set up connections between subdomains
γ=GridSpec("LatLonCap",joinpath(p,"../examples/GRID_LLC90/"))
Γ=GridLoad(γ)
Γ=merge(Γ,IndividualDisplacements.NeighborTileIndices_cs(Γ))

#initialize u0,u1 etc
𝑃=read_uvetc(k,0.0,Γ,joinpath(p,"../examples/nctiles_climatology/"));

#nb # %% {"slideshow": {"slide_type": "subslide"}, "cell_type": "markdown"}
# ### Single interpolation & trajectory Test
#
# Interpolate Velocity (normalized)
(u0,du)=initialize_lonlat(Γ,-160.1,35.1; msk=Γ["hFacC"][:,k]);
⬡!(du,u0,𝑃,0.0);
#u,v,f=du[:]

# Solve for trajectory

prob = ODEProblem(⬡!,u0,𝑃.𝑇,𝑃)
sol = solve(prob,Tsit5(),reltol=1e-4,abstol=1e-4);
#sol = solve(prob,Euler(),dt=1e6)
#size(sol)

#nb # %% {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ## 3. Main Computation Loop
#
# - initial particle positions randomly over Global Ocean
# - initial integration from time 0 to 0.5 month
# - update velocity fields & repeat for ny years

# ## 3.1 Initialization & Initial Solution
#
# Initial Positions:

#(u0,du)=initialize_gridded(𝑃,10)

#(lon, lat) = randn_lonlat(20000)
#(u0,du)=initialize_lonlat(Γ,lon,lat; msk=Γ["hFacC"][:,k])

#Or

if false
    lo0, lo1 = (-160.0, -150.0)
    la0, la1 = (35.0, 45.0)
    n = 100
    lon = lo0 .+ (lo1 - lo0) .* rand(n)
    lat = la0 .+ (la1 - la0) .* rand(n)
    (u0, du) = initialize_lonlat(Γ, lon, lat; msk = Γ["hFacC"][:, k])
else
    (lon, lat) = randn_lonlat(100000)
    (u0, du) = initialize_lonlat(Γ, lon, lat; msk = Γ["hFacC"][:, k])
end

# Arrays For Storing Results
u0_store = deepcopy(u0)
n_store = size(u0_store,2);

# Fraction of the particles reset per month

#r_reset = 0.05
n_reset = Int(round(r_reset*n_store));
#k_reset = rand(1:size(u0_store,2), n_reset)

#nb # %% {"slideshow": {"slide_type": "subslide"}, "cell_type": "markdown"}
# Solve for all trajectories for first 1/2 month

prob = ODEProblem(⬡!,u0,𝑃.𝑇,𝑃)
sol = solve(prob,Euler(),dt=5*86400.0);
#size(sol)

#nb # %% {"slideshow": {"slide_type": "subslide"}, "cell_type": "markdown"}
# Map `i,j` to `lon,lat` coordinates and convert to `DataFrames`

df=postprocess_lonlat(sol,𝑃);
println(size(df))

#nb # %% {"slideshow": {"slide_type": "subslide"}, "cell_type": "markdown"}
# Define iteration function

function iter!(df,𝑃,u0)
    #need an inplace version to update 𝑃 contents?
    update_uvetc!(k,𝑃.𝑇[2],𝑃)
    prob = ODEProblem(⬡!,u0,𝑃.𝑇,𝑃)
    sol = solve(prob,Euler(),dt=5*86400.0)
    tmp = postprocess_lonlat(sol,𝑃)

    k_reset = rand(1:size(u0_store,2), n_reset)
    k_new = rand(1:size(u0_store,2), n_reset)
    t_reset = Int(size(tmp,1)/n_store)-1

    tmp[k_reset.+t_reset*n_store,2:end].=NaN #reset a random subset of particles
    nt=size(sol,3)
    nf=size(sol,2)
    append!(df,tmp[nf+1:end,:])

    u0[:,:] = deepcopy(sol[:,:,end])
    u0[:,k_reset].=deepcopy(u0_store[:,k_new]) #reset a random subset of particles
end

#nb # %% {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ## 3.2 Repeat For `ny*12` Months
#
# _A fraction of the particles, randomly selected, is reset every month to maintain a relatively homogeneous coverage of the Global Ocean by the fleet of particles._

u0 = deepcopy(sol[:,:,end])
dt=86400.0*365.0/12.0
𝑃.𝑇[1]=-dt/2.0
𝑃.𝑇[2]=dt/2.0
for y=1:ny
    for m=1:12
        iter!(df,𝑃,u0)
    end
    println(𝑃.𝑇)
    println(size(df))
end

#nb # %% {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ## 4. Plot trajectories
#

#nb # %% {"slideshow": {"slide_type": "subslide"}, "cell_type": "markdown"}
# Examples Using Various Plotting Packages are provided below.

#p=dirname(pathof(IndividualDisplacements))
#nn=1000

#include(joinpath(p,"../examples/recipes_plots.jl"))
#plt=PlotBasic(df,nn,180.)

#include(joinpath(p,"../examples/recipes_pyplot.jl"))
#PyPlot.figure(); PlotMapProj(df,nn)

#include(joinpath(p,"../examples/recipes_makie.jl"))
#AbstractPlotting.inline!(true) #for Juno, set to false
#scene=PlotMakie(df,nn,180.0)
##Makie.save("LatLonCap300mDepth.png", scene)

#nb # %% {"slideshow": {"slide_type": "skip"}, "cell_type": "markdown"}
# Here we create `lon`, `lat`, and `DL` (log10 of bottom depth) to use in plot background:

lon=[i for i=-179.5:1.0:179.5, j=-89.5:1.0:89.5]
lat=[j for i=-179.5:1.0:179.5, j=-89.5:1.0:89.5]
(f,i,j,w,_,_,_)=InterpolationFactors(Γ,vec(lon),vec(lat))

DL=log10.(Interpolate(Γ["Depth"],f,i,j,w))
DL[findall((!isfinite).(DL))].=NaN
DL=reshape(DL,size(lon));

#nb # %% {"slideshow": {"slide_type": "subslide"}, "cell_type": "markdown"}
# Then we generate plot or movie using `GeoMakie.jl` ...

if false
#using ArgoData
#p = include(joinpath(dirname(pathof(ArgoData)),"movies.jl"));
df.y=2000 .+ df.t ./86400/365
yy=collect(2000:0.05:2000+ny)
scene = ProjMap(DL,colorrange=(2.,4.))
ProjScatterMovie(scene,df,yy,"GlobalDomain_fleet_k"*"$k"*"_v1.mp4",dt=1.0,mrksz=5e3)
end

#nb # %% {"slideshow": {"slide_type": "subslide"}, "cell_type": "markdown"}
# ... or using `Plots.jl`:

if false
contourf(lon[:,1],lat[1,:],transpose(DL),clims=(1.5,5),c = :ice, colorbar=false)

dt=0.0001
t0=0.0
t1=𝑃.𝑇[2]

t=t1
df_t = df[ (df.t.>t-dt).&(df.t.<=t) , :]
scatter!(df_t.lon,df_t.lat,markersize=1.5,c=:red,leg=:none,
    xlims=(-180.0,180.0),ylims=(-90.0,90.0),marker = (:circle, stroke(0)))

t=t0
df_t = df[ (df.t.>t-dt).&(df.t.<=t) , :]
scatter!(df_t.lon,df_t.lat,markersize=1.5,c=:yellow,leg=:none,
    xlims=(-180.0,180.0),ylims=(-90.0,90.0),marker = (:dot, stroke(0)))
end
