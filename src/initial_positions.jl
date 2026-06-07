
module init

using Drifters, MeshArrays, DataFrames, CSV
import Drifters: randn_lonlat

#deprecated function name:
init_positions(np ::Int; filename="global_ocean_circulation.csv") =
    read_initial_positions(np; filename=filename)

"""
    read_initial_positions(np ::Int; filename="global_ocean_circulation.csv")

"""
function read_initial_positions(np ::Int; filename="global_ocean_circulation.csv")
    if filename=="global_ocean_circulation.csv"
        p=dirname(pathof(Drifters))
        fil=joinpath(p,"../examples/worldwide/global_ocean_circulation.csv")
    else
        fil=filename
    end
    return DataFrame(CSV.File(fil))[1:np,:]
end

"""
    init_global_randn(np ::Int , D::NamedTuple)

Randomly distribute `np` points over the Earth, within `D.msk` 
region, and return position in grid index space (`i,j,subdomain`).
"""
function init_global_randn(np ::Int , D::NamedTuple)
    (lon, lat) = randn_lonlat(maximum([2*np 10]))
    (_,_,_,_,f,x,y)=InterpolationFactors(D.Γ,lon,lat)
    m=findall( (f.!==0).*((!isnan).(x)) )
    n=findall(nearest_to_xy(D.msk,x[m],y[m],f[m]).==1.0)[1:np]
    xyf=permutedims([x[m[n]] y[m[n]] f[m[n]]])
    return DataFrame(x=xyf[1,:],y=xyf[2,:],fid=xyf[3,:])
end

"""
    init_regional_3d(np::Int, D::NamedTuple;
        lons=(-81.0, -79.0), lats=(26.0, 28.0), zs=0:27, n_try=3)

- Randomly distribute `np` ocean particles within a lon/lat/depth box. 
- Use `InterpolationFactors` to assign face indices (fid) and grid coordinates (x,y)
- Call `msk_and_z` to :
    - Filter out points such that `D.msk.==0`. 
    - Distribute vertical position based on `zs`.
- Return a DataFrame with columns `x, y, z, fid, lon, lat`.

_note : Default `lons`, `lats`, and `zs` match the former `init gulf stream` / Florida Strait setup. For a wider domain (e.g. tropical band), specify something like `lons=(-180.0,180.0), lats=(-30.0,30.0), zs=1:25`._

_note : initially oversampling by a factor of `n try>1` is needed to find `np` valid points when the mask has zero valued points. Raises an error if fewer than `np` valid are found._
"""
function init_regional_3d(np::Int, D::NamedTuple;
        lons=(-81.0, -79.0), lats=(26.0, 28.0), zs=0:27, n_try=3)
    n_tot = n_try * np
    lon = lons[1] .+ (lons[2] - lons[1]) .* rand(n_tot)
    lat = lats[1] .+ (lats[2] - lats[1]) .* rand(n_tot)
    (_,_,_,_,f,x,y) = Drifters.InterpolationFactors(D.Γ, lon, lat)
    xy=(lon=lon,lat=lat,x=x,y=y,f=f)
    msk_and_z(xy,np; msk=D.msk, zs=zs)
end

function msk_at_xy(xy,msk)
    (;lon,lat,x,y,f) = xy
    np=length(lon)
    M = zeros(np)
    m = findall((f .!== 0) .* ((.!isnan).(x)))
    M[m] = Drifters.nearest_to_xy(msk, x[m], y[m], f[m])
    M
end

"""
    msk_and_z(xy::NamedTuple, np::Int; msk::MeshArray, zs=0:27)

- Filter out `np` points such that `msk.==0` out of `xy`. 
- Distribute vertical position based on `zs`.

_note : initially oversampling by a factor of `n try>1` is needed to find `np` valid points when the mask has zero valued points. Raises an error if fewer than `np` valid are found._
"""
function msk_and_z(xy::NamedTuple,np::Int; msk::MeshArray, zs=0:27)
    (;lon,lat,x,y,f) = xy

    M = msk_at_xy(xy,msk)
    ocean = findall(M .== 1.0)

    length(ocean) < np && error(
        "Only $(length(ocean)) ocean points found (need $np). " *
        "Widen the lon/lat range or reduce np.")

    sel = ocean[1:np]
    xyf = permutedims([x[sel] y[sel] f[sel]])
    z = zs[1] .+ rand(np) .* (zs[end] - zs[1])

    return DataFrame(x=xyf[1,:], y=xyf[2,:], z=z, fid=xyf[3,:], lon=lon[sel], lat=lat[sel])
end

"""
    initial_positions(Γ; nf=10000, lon_rng=(-160.0,-159.0), lat_rng=(30.0,31.0))

Randomly assign initial positions in longitude,latitude ranges. Positions are 
expressed in, normalized, grid point units (x,y in the 0,nx and 0,ny range). 
To convert from longitude,latitude here we take advantage of the regularity 
of the 1 degree grid being used -- for a more general alternative, see the 
global ocean example.

**Warning**: this function assigns all particles to face 1 and uses a simple
linear lon/lat-to-grid mapping. It only works for small regions on face 1.
For regional or global initialization on the LLC90 grid, use `init_regional_3d`.
"""
function initial_positions(Γ::NamedTuple, nf=10000, lon_rng=(-160.0,-159.0), lat_rng=(30.0,31.0), level=1)
   lon=lon_rng[1] .+(lon_rng[2]-lon_rng[1]).*rand(nf)
   lat=lat_rng[1] .+(lat_rng[2]-lat_rng[1]).*rand(nf)
   x=lon .+ (21. - Γ.XC[1][21,1])
   y=lat .+ (111. - Γ.YC[1][1,111])

   return DataFrame(:x => x, :y => y, :z => fill(level,nf),:fid => fill(1,nf))
end

end
