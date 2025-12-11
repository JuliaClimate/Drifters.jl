"""
    convert_to_FlowFields(U::Array{T,2},V::Array{T,2},t1::T) where T

Convert a pair of U,V arrays (staggered C-grid velocity field in 2D) to
a `uvMeshArrays` struct ready for integration of individual displacements
from time `t0=0` to time `t1`.
"""
function convert_to_FlowFields(U::Array{T,2},V::Array{T,2},t1::T) where T
    np,nq=size(U)
    Γ=MeshArrays.Grids_simple.periodic_domain(np,nq)

    g=Γ.XC.grid
    u=MeshArray(g,[U])
    v=MeshArray(g,[V])
    (u,v)==MeshArrays.exchange_main(u,v,1)
    func=(u -> MeshArrays.update_location_dpdo!(u,g))

    uvMeshArrays{eltype(u.MA)}(u.MA,u.MA,v.MA,v.MA,[0,t1],func)
end

"""
    postprocess_MeshArray(sol,P::FlowFields,D::NamedTuple; id=missing, T=missing)

Copy `sol` to a `DataFrame` & map position to lon,lat coordinates
using "exchanged" D.XC, D.YC via `add_lonlat!`
"""
function postprocess_MeshArray(sol,P::FlowFields, D::NamedTuple; id=missing, T=missing)
    ismissing(id) ? id=collect(1:size(sol,2)) : nothing
    ismissing(T) ? T=P.T : nothing
    
    if isa(sol,EnsembleSolution)
        nd=length(sol.u[1][:,1])
        np=length(sol)
        x=[[sol.u[i][:,1][1] for i in 1:np];[sol.u[i][:,end][1] for i in 1:np]]
        y=[[sol.u[i][:,1][2] for i in 1:np];[sol.u[i][:,end][2] for i in 1:np]]
        fIndex=[[sol.u[i][:,1][nd] for i in 1:np];[sol.u[i][:,end][nd] for i in 1:np]]
        t=[fill(T[1],np);fill(T[2],np)]
        id=[id[:,1];id[:,1]]
    else
        nt=length(sol.u)
        x=sol[1,:]
        y=sol[2,:]
        fIndex=sol[end,:]
        t=T[1] .+ (T[2]-T[1]) * collect(0:nt-1) / (nt-1)
        id=fill(id[1],nt)
    end

    P.u0.grid.nFaces==1 ? fIndex=ones(size(x)) : nothing
    
    df = DataFrame(ID=id[:], x=x[:], y=y[:], fid=fIndex[:], t=t[:])

    return df
#    return id,x,y,fIndex,t
end

"""
    update_location!(pos,P)

Update `pos` via `P.update_location!` if needed.
"""
function update_location!(pos,P)
    g=P.u0.grid
    while MeshArrays.location_is_out(pos,g)
        P.update_location!(pos)
    end
end

"""
    add_lonlat!(df::DataFrame,XC,YC)

Add lon & lat to dataframe using "exchanged" XC, YC
"""
function add_lonlat!(df::DataFrame,XC::MeshArray_wh,YC::MeshArray_wh)
    x = cosd.(YC.MA) * cosd.(XC.MA)
    y = cosd.(YC.MA) * sind.(XC.MA)
    z = sind.(YC.MA)

    df.lon=0*df.x
    df.lat=0*df.x
    for i in eachindex(df.lon)
        xx=interp_to_xy(df.x[i],df.y[i],df.fid[i],x)
        yy=interp_to_xy(df.x[i],df.y[i],df.fid[i],y)
        zz=interp_to_xy(df.x[i],df.y[i],df.fid[i],z)
        df.lat[i] = asind(zz/sqrt(xx^2+yy^2+zz^2))
        df.lon[i] = atand(yy, xx)
    end

    return df
end

"""
    add_lonlat!(df::DataFrame,XC,YC,func::Function)

Add lon & lat to dataframe using "exchanged" XC, YC after updating 
subdomain indices (via func) if needed (MeshArrays.location_is_out)
"""
function add_lonlat!(df::DataFrame,XC::MeshArray_wh,YC::MeshArray_wh,func::Function)
    g=XC.MA.grid
    u=zeros(3)

    for i in eachindex(df.x)
        u[:]=[df.x[i];df.y[i];df.fid[i]]
        while MeshArrays.location_is_out(u,g)
            func(u)
            df.x[i]=u[1]
            df.y[i]=u[2]
            df.fid[i]=u[3]
        end
    end

    add_lonlat!(df,XC,YC)

    return df
end

"""
    postprocess_xy(sol,P::FlowFields,D::NamedTuple; id=missing, T=missing)

Copy `sol` to a `DataFrame` & map position to x,y coordinates,
and define time axis for a simple doubly periodic domain
"""
function postprocess_xy(sol,P::FlowFields,D::NamedTuple; id=missing, T=missing)
    ismissing(id) ? id=collect(1:size(sol,2)) : nothing
    ismissing(T) ? T=P.T : nothing

    isa(P.u0,MeshArray) ? (nx,ny)=P.u0.grid.ioSize[1:2] : (nx,ny)=size(P.u0)[1:2]

    if isa(sol,EnsembleSolution)
        nd=length(sol.u[1][1])
        np=length(sol)
        x=[mod.([sol[i][1,1] for i in 1:np],Ref(nx));
            mod.([sol[i][1,end] for i in 1:np],Ref(nx))];
        y=[mod.([sol[i][2,1] for i in 1:np],Ref(ny));
            mod.([sol[i][2,end] for i in 1:np],Ref(ny))]
        t=[fill(T[1],np);fill(T[2],np)]
        id=[id[:,1];id[:,1]]
    else
        nt=length(sol.u)
        x=mod.(sol[1,:],Ref(nx))
        y=mod.(sol[2,:],Ref(ny))
        t=T[1] .+ (T[2]-T[1]) * collect(0:nt-1) / (nt-1)
        id=fill(id[1],nt)
    end

    return DataFrame(ID=id[:], t=t[:], x=x[:], y=y[:])
end

"""
    randn_lonlat(nn=1,seed=missing)

Randomly distributed longitude, latitude positions on the sphere.
"""
function randn_lonlat(nn=1;seed=missing)
    !ismissing(seed) ? rng = MersenneTwister(1234) : rng = MersenneTwister()
    tmp = randn(rng, Float64, (nn, 3))
    tmpn = tmp ./ sqrt.(sum(tmp.^2, dims=2))
    lon = rad2deg.(atan.(tmpn[:,2], tmpn[:,1]))
    lat = 90.0 .- rad2deg.(acos.(tmpn[:,3]))
    return lon, lat
end

"""
    nearest_to_xy(α::MeshArray,x,y,f)

Value of α at eachindex of the grid cell center nearest to `x,y` on subdomain array / facet `f`
"""
nearest_to_xy(α::MeshArray,x,y,f) = [α[Int(f[i]),1][ Int(round(x[i] .+ 0.5)), Int(round(y[i] .+ 0.5)) ] for i in eachindex(x)]

"""
    nearest_to_xy(α::Array,x,y)

Value of α at eachindex of the grid cell center nearest to `x,y`
"""
nearest_to_xy(α::Array,x,y) = [α[ Int(round(x[i] .+ 0.5)), Int(round(y[i] .+ 0.5)) ] for i in eachindex(x)]

"""
    interp_to_lonlat(X::MeshArray,Γ::NamedTuple,lon,lat)

Use MeshArrays.Interpolate() to interpolate to arbitrary positions (e.g., a regular grid for plotting).

# Extended help

```jldoctest
using Drifters
import Drifters: MeshArrays
γ=MeshArrays.GridSpec(ID=:LLC90)
Γ=MeshArrays.GridLoad(γ,option="full")

lon=[i for i=20.:20.0:380., j=-70.:10.0:70.]
lat=[j for i=20.:20.0:380., j=-70.:10.0:70.]
tmp1=interp_to_lonlat(Γ.Depth,Γ,lon,lat)

(f,i,j,w,_,_,_)=MeshArrays.InterpolationFactors(Γ,vec(lon),vec(lat))
IntFac=(lon=lon,lat=lat,f=f,i=i,j=j,w=w)
tmp1=interp_to_lonlat(Γ.Depth,IntFac)
    
prod(isapprox(maximum(tmp1),5896.,atol=1.0))

# output

true
```
"""
function interp_to_lonlat(X::MeshArray,Γ::NamedTuple,lon,lat)
    (f,i,j,w,_,_,_)=InterpolationFactors(Γ,vec(lon),vec(lat))
    return reshape(Interpolate(X,f,i,j,w),size(lon))
end

"""
    interp_to_lonlat(X::MeshArray,IntFac::NamedTuple)

Use MeshArrays.Interpolate() to interpolate to arbitrary positions (e.g., a regular grid for plotting).
"""
function interp_to_lonlat(X::MeshArray,IntFac::NamedTuple)
    (; f,i,j,w,lon,lat) = IntFac
    return reshape(Interpolate(X,f,i,j,w),size(lon))
end

"""
    interp_to_xy(df::DataFrame,Zin::MeshArray)

Interpolate "exchanged" / "hallo-included" Zin to df[!,:x], df[!,:y] on df[!,:fid]
"""
function interp_to_xy(df::DataFrame,Zin::MeshArray_wh)
    x=df[!,:x];
    y=df[!,:y];
    f=Int.(df[!,:fid]);
    dx,dy=(x - floor.(x) .+ 0.5,y - floor.(y) .+ 0.5);
    i_c = Int32.(floor.(x)) .+ 1;
    j_c = Int32.(floor.(y)) .+ 1;

    Z=zeros(length(x),4)
    [Z[k,:]=Zin.MA[f[k]][i_c[k]:i_c[k]+1,j_c[k]:j_c[k]+1][:] for k in 1:length(i_c)]

    return (1.0 .-dx).*(1.0 .-dy).*Z[:,1]+dx.*(1.0 .-dy).*Z[:,2] +
           (1.0 .-dx).*dy.*Z[:,3]+dx.*dy.*Z[:,4]
end

function interp_to_xy(x,y,f,Zin::MeshArray)
    dx,dy=(x - floor(x) + 0.5,y - floor(y) + 0.5)
    i_c = Int(floor(x)) + 1
    j_c = Int(floor(y)) + 1
    ff=Int(f)

    return (1.0-dx)*(1.0-dy)*Zin[ff][i_c,j_c] +
    dx*(1.0-dy)*Zin[ff][i_c+1,j_c] +
    (1.0-dx)*dy*Zin[ff][i_c,j_c+1] +
    dx*dy*Zin[ff][i_c+1,j_c+1]
end

"""
    stproj(XC,YC;XC0=0.0,YC0=90.0)

Apply to XC,YC (longitude, latitude) the stereographic projection
which would put XC0,YC0 (longitude, latitude) at x,y=0,0
"""
function stproj(XC,YC,XC0=0.0,YC0=90.0)

# compute spherical coordinates:
phi=pi/180*XC; theta=pi/180*(90-YC);
phi0=pi/180*XC0; theta0=pi/180*(90-YC0);

# compute cartesian coordinates:
x=sin(theta)*cos(phi)
y=sin(theta)*sin(phi)
z=cos(theta)

# bring chosen longitude to x>0,y=0:
xx=cos(phi0)*x+sin(phi0)*y
yy=-sin(phi0)*x+cos(phi0)*y
zz=z
# bring chosen point to South Pole:
x=cos(theta0)*xx-sin(theta0)*zz
y=yy
z=sin(theta0)*xx+cos(theta0)*zz

# stereographic projection from South Pole:
xx=x/(1+z)
yy=y/(1+z)

return xx,yy
end

"""
    stproj_inv(xx,yy;XC0=0.0,YC0=90.0)

Apply to xx,yy (stereographic projection coordinates) the reverse 
of the stereographic projection which puts XC0,YC0 (longitude, 
latitude) at x,y=0,0
"""
function stproj_inv(xx,yy,XC0=0.0,YC0=90.0)
    phi0=pi/180*XC0; theta0=pi/180*(90-YC0);

    # Reverse stereographic projection from North Pole:
    (x,y,z)=2*[xx yy (1-xx^2-yy^2)/2]./(1+xx^2+yy^2)

    # bring chosen point back from South Pole:
    xx=cos(theta0)*x+sin(theta0)*z
    yy=y
    zz=-sin(theta0)*x+cos(theta0)*z

    # bring chosen longitude back from x>0,y=0:
    x=cos(phi0)*xx-sin(phi0)*yy
    y=sin(phi0)*xx+cos(phi0)*yy
    z=zz

    # compute spherical coordinates:
    theta=atan(sqrt(x^2+y^2)/z)
    phi=atan(y, x)
    
    XC=180/pi*phi
    theta>=0 ? YC=90-180/pi*theta : YC=-90-180/pi*theta

    return XC,YC
end
