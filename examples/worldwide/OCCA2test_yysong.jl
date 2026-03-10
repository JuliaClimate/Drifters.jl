 # For run the script
"""
using Pkg; Pkg.activate("/Users/yysong/Documents/julia/Julia_Project_Lagrangian/")
using Drifters, MeshArrays, Climatology, MAT, CSV
include("worldwide/OCCA2_yysong.jl")
using Main.OCCA2_yysong
include("worldwide/OCCA2test_yysong.jl")
out_path=joinpath("/Users/yysong/Desktop/Ciara/global_ocean_tmp_SE_2014Jan_backward_tropics")
!ispath(out_path) ? mkdir(out_path) : nothing
kk=1

main_comp_SO(kk)

"""
# For plotting

"""
using Drifters, GLMakie
#include("worldwide/global_ocean_circulation.jl") remove 

pth=out_path
fil=joinpath(pth,"initial_1_1_◀◀.csv")
df=CSV.read(fil,DataFrame)

? how to get I to plot when I is not saved in the csv file?
x=DriftersDataset( data=(I=I,df=I.🔴,), options=(plot_type=:global_plot1,) )

file_output_mp4="/Users/yysong/Desktop/figures/Julia_practice/Ciara_2014Jan_backward_tropics.mp4"
record(fig, file_output_mp4, -50:nt, framerate = 25) do t
    tt[]=max(t,0)
end

"""  

function main_comp_SO(kk)

println(kk)

## Initial Position Files


file_base = "initial_1_1"
backward_time = true
backward_time ? file_base=file_base*"_◀◀" : file_base=file_base*"_▶▶"

k=0
t=1.076922e9
np_0=20000
ny=1 #number of years
nm=6 #number of months

P,D=init_FlowFields(k=k,t=t,backward_time=backward_time)

#----------------- mask df -----------------------------------------------#

γ=MeshArrays.GridSpec(ID=:LLC90)
Γ=MeshArrays.GridLoad(γ,option="full")
lon0=-180.0; lon1=180.0; lat0=-30.0; lat1=30.0; 
masks_h = (Γ.XC.>=lon0)*(Γ.XC.<=lon1)*(Γ.YC.>=lat0)*(Γ.YC.<=lat1)
lon_0 = Dict()
lat_0 = Dict()
rac_0 = Dict()
for i in 1:5
    
    pos = findall(x -> x != 0.0 && !isnan(x), masks_h[i])
    if !isempty(pos)
        lon_0[i] = D.Γ.XC[i][pos]
        lat_0[i] = D.Γ.YC[i][pos]
        rac_0[i] = D.Γ.RAC[i][pos]
    end
end


lon_1 = round.(vcat(values(lon_0)...),digits=10)
lat_1 = round.(vcat(values(lat_0)...),digits=10)
rac_1 = round.(vcat(values(rac_0)...),digits=10)
lat_rg = sort(unique(lat_1))
lon_rg = sort(unique(lon_1))
lonint = 1.0 # interval between longitudes
latint = 0.6 # interval between latitudes
nlon1 = length(lon_1);
# create random lons, lats by each grid
all_lon_2 = Vector{Vector{Float64}}()
all_lat_2 = Vector{Vector{Float64}}()
for i = 1:nlon1
    np_egrid = Int(div(2*np_0*rac_1[i],sum(rac_1)))
    lon_2 = rand(np_egrid) .* lonint .+ lon_1[i] .- lonint/2    
    lat_2 = rand(np_egrid) .* latint .+ lat_1[i] .- latint/2
    push!(all_lon_2, lon_2)
    push!(all_lat_2, lat_2)
end
lon = vcat(all_lon_2...)
lat = vcat(all_lat_2...)

"""
fig = Figure(size = (1000, 800), backgroundcolor = :grey80)
ax = Axis(fig[1,1])
Makie.scatter!(lon,lat,color=:red, markersize = 10, label = "particles")
Makie.scatter!(lon_1,lat_1,color=:blue, markersize = 13, label = "original")
axislegend(ax, position = :rt)
fig
"""

(_,_,_,_,f,x,y)=Drifters.InterpolationFactors(D.Γ,lon,lat)

zs=11:25  # Note: z index, not depth values
z=zs[1] .+rand(length(lon))*(zs[end]-zs[1])

xx=Int.(round.(x))
yy=Int.(round.(y))
zz=Int.(round.(z))

batch_T=zeros(size(x,1),50)
dx,dy=(x - floor.(x) .+ 0.5,y - floor.(y) .+ 0.5);
i_c = Int32.(floor.(x)) .+ 1;
j_c = Int32.(floor.(y)) .+ 1;
nr=50
# m=findall( (f.!==0).*((!isnan).(x)) )   #CartesianIndex
m=findall( ((f.!==0).*((!isnan).(x)))[:] )  # Vector
# horizontal bilinear interpolation 
for k in 1:nr, jj in m
    tmp0=(1.0-dx[jj])*(1.0-dy[jj])*D.exmsk[f[jj],k][i_c[jj],j_c[jj]]+
    (dx[jj])*(1.0-dy[jj])*D.exmsk[f[jj],k][i_c[jj]+1,j_c[jj]]+
    (1.0-dx[jj])*(dy[jj])*D.exmsk[f[jj],k][i_c[jj],j_c[jj]+1]+
    (dx[jj])*(dy[jj])*D.exmsk[f[jj],k][i_c[jj]+1,j_c[jj]+1]
    #
    tmp1=(1.0-dx[jj])*(1.0-dy[jj])*D.θ1[f[jj],k][i_c[jj],j_c[jj]]+
    (dx[jj])*(1.0-dy[jj])*D.θ1[f[jj],k][i_c[jj]+1,j_c[jj]]+
    (1.0-dx[jj])*(dy[jj])*D.θ1[f[jj],k][i_c[jj],j_c[jj]+1]+
    (dx[jj])*(dy[jj])*D.θ1[f[jj],k][i_c[jj]+1,j_c[jj]+1]
    batch_T[jj,k]=tmp1/tmp0
end
# vertical linear interpolation
local_T=0*x;
for p in m
    k1=floor(z[p]+0.5)
        a2=(z[p]+0.5)-k1
        k2=Int(min(max(k1+1,1),nr))
        k1=Int(min(max(k1,1),nr))
        local_T[p]=(1-a2)*batch_T[p,k1]+a2*batch_T[p,k2]
end

df = DataFrame(fid=f[:], x=x[:], y=y[:], z=z[:])
XC=exchange(D.Γ.XC) #add 1 lon point at each edge
YC=exchange(D.Γ.YC) #add 1 lat point at each edge
df.fid[findall(df.fid.==0)].=1 # replace 0 indices
add_lonlat!(df,XC,YC)  # add lon and lat in df

# remove rows which are zeros
df = filter(row -> all(!=(0), row), eachrow(df))

# count new np and make it as multiple of 100
np = floor(Int, size(df, 1) / 100) * 100
println(np)
# reserve the first np rows
df = first(df, np)

#------------------------ df done -----------------------------------------#

S = init_storage(np,100,length(D.Γ.RC),50)
I = Individuals(P,df.x,df.y,df.z,df.fid,
    (D=merge(D,S),∫=custom∫,🔧=custom🔧,🔴=deepcopy(custom🔴)))

T=(0.0,I.P.T[2])
custom∫!(I,T)

[step!(I) for y=1:ny, m=1:nm]

file_output=joinpath(out_path,file_base*".csv")
CSV.write(file_output, Float32.(I.🔴))

I

end

function step!(I::Individuals)
    I.D.🔄(I)
    custom∫!(I)
end
