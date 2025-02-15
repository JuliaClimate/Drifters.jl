module Oscar

import Drifters: uvArrays, Individuals, FlowFields
import Drifters: dxdt!, ∫!, update_FlowFields!, postprocess_xy

using Glob, DataFrames, CSV, NetCDF

n_part=1000
reset_rate = 0.05 #per day
dT=86400.0
nt=30

list_files(path="data",year=2021)=glob("oscar_currents_final_$(year)*.nc",path)

"""
The Oscar data product
 
https://podaac.jpl.nasa.gov/dataset/OSCAR_L4_OC_FINAL_V2.0
 
averaged surface currents are provided on a global 0.25 x 0.25 degree grid
 
An example netcdf granule is here:
https://archive.podaac.earthdata.nasa.gov/podaac-ops-cumulus-protected/OSCAR_L4_OC_FINAL_V2.0/oscar_currents_final_20220504.nc

To download files, in shell : 

```
podaac-data-downloader -c OSCAR_L4_OC_FINAL_V2.0 -d ./data --start-date 2021-01-01T00:00:00Z --end-date 2021-02-01T00:00:00Z -e ""
```

Source :

- https://github.com/podaac/data-subscriber
- https://podaac.jpl.nasa.gov/dataset/OSCAR_L4_OC_FINAL_V2.0#capability-modal-download
"""

## tools

preproc(x) = begin
  y=zeros(size(x))
  [(x[a]>-999 ? y[a]=x[a] : nothing) for a in eachindex(x)]
  permutedims(y[:,:])
end

function read_oscar(fil,G,dT)
  u=preproc(ncread(fil,"u")); v=preproc(ncread(fil,"v"))
  #ug=preproc(ncread(fil,"ug")); vg=preproc(ncread(fil,"vg"))
  lon=ncread(fil,"lon"); lat=ncread(fil,"lat"); time=ncread(fil,"time")
  O=(u=u,v=v,lon=lon,lat=lat,time=time)

  uC=O.u./G.dxF; vC=O.v./G.dyF
  uW=0.5*(circshift(uC, (1,0))+uC) #staggered u converted to grid point units (m/s -> 1/s)
  vS=0.5*(circshift(vC, (0,1))+vC) #staggered v converted to grid point units (m/s -> 1/s)

  P=(u=uW,v=vS,dT=dT,T=[0,dT])
  (O,P)
end

function update_FlowFields!(I::Individuals,G,dT,list)
  dT=I.P.T[2]-I.P.T[1]
  m0=Int(I.P.T[2]/dT)
  println(m0)
  _,P=read_oscar(list[m0],G,dT)
  I.P.u0.=P.u; I.P.v0.=P.v
  _,P=read_oscar(list[m0+1],G,dT)
  I.P.u1.=P.u; I.P.v1.=P.v
  I.P.T.=dT*[m0;m0+1]
end

function reset_📌!(I::Individuals,frac::Number,📌::Array)
  np=length(I.🆔)
  n_reset = Int(round(frac*np))
  k_reset = rand(1:np, n_reset)
  l_reset = rand(1:np, n_reset)
  I.📌[k_reset]=deepcopy(📌[l_reset])
  #isempty(I.🔴.ID) ? m=maximum(I.🆔) : m=max(maximum(I.🔴.ID),maximum(I.🆔))
  #I.🆔[k_reset]=collect(1:n_reset) .+ m
end

lonlat(x,y)=(x*0.25,y*0.25-90.0)

function custom🔧(sol,F::uvArrays,D::NamedTuple;id=missing,T=missing)
  df=postprocess_xy(sol,F,D,id=id,T=T)
  
  ll=lonlat.(df.x,df.y)
  df.lon=Float32.([l[1] for l in ll])
  df.lat=Float32.([l[2] for l in ll])

#  vel=0*df.x; du=zeros(2); [vel[k]=sqrt(sum(dxdt!(du,[df.x[k],df.y[k]],F,T[2]).^2)) for k in 1:length(vel)];
  dxdt=0*df.x; du=zeros(2); [dxdt[k]=dxdt!(du,[df.x[k],df.y[k]],F,T[2])[1] for k in 1:length(dxdt)]
  dydt=0*df.y; du=zeros(2); [dydt[k]=dxdt!(du,[df.x[k],df.y[k]],F,T[2])[2] for k in 1:length(dydt)]
  df.dxdt=Float32.(dxdt); df.dydt=Float32.(dydt)

  df.x.=Float32.(df.x); df.y.=Float32.(df.y); df.t.=Float32.(df.t)

  return df
end

custom🔴 = DataFrame(ID=Int[], x=Float32[], y=Float32[], t=Float32[], lon=Float32[], lat=Float32[], dxdt=Float32[], dydt=Float32[])

## grid factors and flow fields normalization

grid()=begin
  xc=0.0:0.25:359.75
  yc=-89.75:0.25:89.75
  yg=-89.875:0.25:89.875
  rSphere = 6370.0*1000
  dxF = rSphere*deg2rad.(cosd.(yc)*0.25)
  dyF = rSphere*deg2rad.(diff(yg))
  x=sind.(yg[2:end])-sind.(yg[1:end-1])
  RAC = rSphere*rSphere*1.0*deg2rad.(abs.(x))

  (
    XC=xc*ones(1,length(yc)),
    YC=permutedims(yc*ones(1,length(xc))),
    dxF=permutedims(dxF*ones(1,length(xc))),
    dyF=permutedims(dyF*ones(1,length(xc))),
    RAC=permutedims(RAC*ones(1,length(xc))),
    )
end

init(n_part,x1)=0.0.+x1*rand(n_part)

"""
    main_loop(;  input_files=list_files(),
      output_file=joinpath("movies","oscar_v06.csv"), do_save=false)

```
using Drifters, GLMakie

I=Drifters.Oscar.main_loop()

using Proj, MeshArrays, GeoJSON
lon0=-160.0; proj=Proj.Transformation(MA_preset=2,lon0=lon0)
options=(plot_type=:Oscar_plot,proj=proj,lon0=-160,add_background=true,add_polygons=true)

J=DriftersDataset( data=(df=I.🔴,), options=options)
plot(J)
"""
function main_loop(;  input_files=list_files(),
                      output_file=joinpath("movies","oscar_v06.csv"), 
                      do_save=false)

## initialize grid and flow fields

G=grid()
O,P=read_oscar(input_files[1],G,dT)
F=FlowFields(u=P.u,v=P.v,period=P.T)
siz=size(P.u)
  
## initial particle positions

x0=init(2*n_part,siz[1]); y0=init(2*n_part,siz[2])
v0=0*x0; du=zeros(2); [v0[k]=sqrt(sum(dxdt!(du,[x0[k],y0[k]],F,1.0).^2)) for k in 1:length(x0)];
x0=x0[findall(v0.>0)[1:n_part]]
y0=y0[findall(v0.>0)[1:n_part]]

📌=[[x0[i],y0[i]] for i in 1:length(x0)]

## solve

FF=uvArrays{Float64}(zeros(siz),zeros(siz),zeros(siz),zeros(siz),P.T)
II=Individuals(FF,x0,y0,(🔧=custom🔧,🔴=custom🔴))
update_FlowFields!(II,G,dT,input_files)
for tt=1:nt
  ∫!(II)
  update_FlowFields!(II,G,dT,input_files)
  reset_📌!(II,reset_rate,📌)
end
#II_t=groupby(II.🔴,:t)

do_save ? CSV.write(output_file,II.🔴) : nothing

II
end

end
