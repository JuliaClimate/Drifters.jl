module OCCA

import Drifters: data_path, uvwMeshArrays, FlowFields
import Drifters: postprocess_MeshArray, add_lonlat!, interp_to_xy

##

import NetCDF
function rd(filename, varname,n)
   fil = NetCDF.open(filename, varname)
   siz = size(fil)
   tmp = zeros(siz[1:2]...,n)
   [tmp .+= fil[:,:,1:n,t] for t=1:12]
   tmp ./= 12.0
   tmp[findall(tmp.<-1e22)] .= 0.0
   return tmp
end

read_data(v0,data_path,n)=begin
   file_name=
   if v0=="u"
      "DDuvel.0406clim.nc"
   elseif v0=="v"
      "DDvvel.0406clim.nc"
   elseif v0=="w"
      "DDwvel.0406clim.nc"
   elseif v0=="theta"
      "DDtheta.0406clim.nc"
   elseif v0=="salt"
      "DDsalt.0406clim.nc"
   else
      "unknown"
   end
   fileIn=joinpath(data_path,file_name)
   rd(fileIn,v0,n)
end

##

import Drifters.DataFrames: DataFrame
import Drifters.MeshArrays as MeshArrays
import Drifters.MeshArrays: MeshArray, exchange
import Drifters.OrdinaryDiffEq: solve, Tsit5, remake
import Drifters.OrdinaryDiffEq: ODEProblem, EnsembleProblem

function setup(;backward_in_time::Bool=false,nmax=Inf)

   γ=MeshArrays.GridSpec(ID=:onedegree)
   Γ=MeshArrays.GridLoad(γ;option="full")
   n=length(Γ.RC)
   isfinite(nmax) ? n=min(n,Int(nmax)) : nothing

   g=Γ.XC.grid
   func=(u -> MeshArrays.update_location_PeriodicDomain!(u,g))

   jj=[:hFacC, :hFacW, :hFacS, :DXG, :DYG, :RAC, :RAZ, :RAS]
   ii=findall([!in(i,jj) for i in keys(Γ)])
   Γ=(; zip(Symbol.(keys(Γ)[ii]), values(Γ)[ii])...)

   backward_in_time ? s=-1.0 : s=1.0
   s=Float32(s)
   pth=data_path(:OCCA)
  
   u=s*read(read_data("u",pth,n),MeshArray(γ,Float32,n))
   v=s*read(read_data("v",pth,n),MeshArray(γ,Float32,n))

   w=s*read_data("w",pth,n)
   w=-cat(w,zeros(360, 160),dims=3)
   w[:,:,1] .=0.0
   w=read(w,MeshArray(γ,Float32,n+1))

   θ=read(read_data("theta",pth,n),MeshArray(γ,Float32,n))
   #𝑆=read(read_data("salt",pth,n),MeshArray(γ,Float32,n))

   for i in eachindex(u)
      u[i]=u[i]./Γ.DXC[1]
      v[i]=v[i]./Γ.DYC[1]
   end

   for i in eachindex(u)
      u[i]=circshift(u[i],[-180 0])
      v[i]=circshift(v[i],[-180 0])
      θ[i]=circshift(θ[i],[-180 0])
#      𝑆[i]=circshift(𝑆[i],[-180 0])
   end

   for i in eachindex(w)
      w[i]=w[i]./Γ.DRC[min(i[2]+1,n)]
      w[i]=circshift(w[i],[-180 0])
   end

   tmpx=circshift(Γ.XC[1],[-180 0])
   tmpx[1:180,:]=tmpx[1:180,:] .- 360.0
   Γ.XC[1]=tmpx

   tmpx=circshift(Γ.XG[1],[-180 0])
   tmpx[1:180,:]=tmpx[1:180,:] .- 360.0
   Γ.XG[1]=tmpx
   Γ.Depth[1]=circshift(Γ.Depth[1],[-180 0])

   t0=0.0; t1=86400*366*2.0;

   (u,v)=exchange(u,v)
   w=exchange(w)

   P=FlowFields(u,u,v,v,w,w,[t0,t1],func)

   XC=exchange(Γ.XC)
   YC=exchange(Γ.YC)

   iso=MeshArrays.isosurface(θ,15,Γ)
   iso[findall(isnan.(iso))].=0.
   iso=exchange(iso)

   D = (iso=iso, XC=XC, YC=YC, RF=Γ.RF, RC=Γ.RC,ioSize=(360,160,n), Γ=Γ)

   return P,D

end

custom🔴 = DataFrame(ID=Int[], fid=Int[], x=Float64[], y=Float64[],
   k=Float64[], z=Float64[], iso=Float64[], t=Float64[],
   lon=Float64[], lat=Float64[], dlon=Float64[], dlat=Float64[], 
   year=Float64[], col=Symbol[])

function custom🔧(sol,P::uvwMeshArrays,D::NamedTuple;id=missing,T=missing)
   df=postprocess_MeshArray(sol,P,D,id=id,T=T)
   add_lonlat!(df,D.XC,D.YC)
   df.dlon=0*df.lon
   df.dlat=0*df.lat

   #add year (convenience time axis for plotting)
   df.year=df.t ./86400/365

   #add depth (i.e. the 3rd, vertical, coordinate)
   k=[[sol.u[i].u[1][3] for i in 1:size(sol,3)];[sol.u[i].u[end][3] for i in 1:size(sol,3)]]
   nz=length(D.RC)
   df.k=min.(max.(k[:],Ref(0.0)),Ref(nz)) #level
   k=Int.(floor.(df.k)); w=(df.k-k);
   df.z=D.RF[1 .+ k].*(1 .- w)+D.RF[2 .+ k].*w #depth

   #add selected isotherm depth
   df.iso=interp_to_xy(df,D.iso)

   #add color = f(iso-z)
   c=fill(:gold,length(df.iso))
   c[findall(df.iso.<df.z)].=:violet
   df.col=c

   #to plot e.g. Pacific Ocean transports, shift longitude convention?
   df.lon[findall(df.lon .< 0.0 )] = df.lon[findall(df.lon .< 0.0 )] .+360.0
   return df
end

function custom∫(prob)
	u0 = prob.u0
	prob_func(prob,i,repeat) = remake(prob,u0=u0[i])
	indiv_prob = ODEProblem(prob.f,u0[1],prob.tspan,prob.p)
	ensemble_prob = EnsembleProblem(indiv_prob,prob_func=prob_func,safetycopy=false)
	solve(ensemble_prob, Tsit5(), reltol=1e-8, abstol=1e-8, 
         trajectories=length(u0),saveat=365/12*86400.0)
end

end
