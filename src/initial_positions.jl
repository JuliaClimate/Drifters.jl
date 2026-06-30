
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
    initial_positions_2d(nchunk::Int, Γ; mask=nothing)

Randomly distribute `np` points over the Earth, within `mask` 
region, and return position in grid index space (`i,j,subdomain`).
"""
function initial_positions_2d(nchunk::Int, Γ; mask=nothing)
    (lon, lat) = randn_lonlat(maximum([2*nchunk 10]))
    (_,_,_,_,f,x,y) = InterpolationFactors(Γ, lon, lat)

    M = Bool.((f.!==0).*((!isnan).(x)))
    for i in 1:length(M)
        if M[i]
            Bool((x[i]<=0).+(x[i]>=Γ.XC.fSize[f[i]][1])) ? (M[i]=false) : nothing
            Bool((y[i]<=0).+(y[i]>=Γ.XC.fSize[f[i]][2])) ? (M[i]=false) : nothing
        end
    end
    
    m = findall(M)
    mask_vals = mask !== nothing ? msk_at_xy(x[m], y[m], mask) : ones(length(m))
    
    return DataFrame(
        x=x[m], 
        y=y[m], 
        fid=f[m], 
        lon=lon[m], 
        lat=lat[m], 
        M=mask_vals
    )
end

function initial_positions_2d(np::Int, Γ, nchunk; mask=nothing)
	result = DataFrame(
		x=fill(NaN, np), 
		y=fill(NaN, np), 
		fid=fill(0, np), 
		lon=fill(NaN, np), 
		lat=fill(NaN, np), 
		M=fill(NaN, np)
	)
	
	idx = 1  # pointer to current position
	
	while idx <= np
		xy0 = initial_positions_2d(nchunk, Γ; mask=mask)
		xy_valid = xy0[xy0.M .== 1.0, :]
		
		# How many rows can we fit before reaching np?
		n_to_add = min(nrow(xy_valid), np - idx + 1)
		
		# Write in place
		result[idx:(idx + n_to_add - 1), :] = xy_valid[1:n_to_add, :]
		idx += n_to_add
	end
	
	return result
end

#for backwards compatibility :
init_global_randn(np ::Int , D::NamedTuple)=initial_positions_2d(np, D.Γ)

"""
    init_regional_3d(np::Int, D::NamedTuple;
        lons=(-81.0, -79.0), lats=(26.0, 28.0), zs=0:27, n_try=3)

- Randomly distribute `np` ocean particles within a lon/lat/depth box. 
- Use `InterpolationFactors` to assign face indices (fid) and grid coordinates (x,y)
- Call `msk_filter` to filter out points such that `D.msk.==0`. 
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
    xy=DataFrame(lon=lon,lat=lat,x=x[:],y=y[:],fid=f[:])
    xy=msk_filter(xy,np,D.msk)
    xy.z = zs[1] .+ rand(np) .* (zs[end] - zs[1])
    xy
end

function msk_at_xy(xy,msk)
    (;lon,lat,x,y,fid) = xy
    np=length(lon)
    M = zeros(np)
    m = findall((fid .!== 0) .* ((.!isnan).(x)))
    M[m] = Drifters.nearest_to_xy(msk, x[m], y[m], fid[m])
    M
end

"""
    msk_filter(xy::NamedTuple, np::Int, msk::MeshArray)

- Filter out `np` points such that `msk.==0` out of `xy`. 

_note : initially oversampling by a factor of `n try>1` is needed to find `np` valid points when the mask has zero valued points. Raises an error if fewer than `np` valid are found._
"""
function msk_filter(xy::DataFrame,np::Int,msk::MeshArray)
    (;lon,lat,x,y,fid) = xy

    M = msk_at_xy(xy,msk)
    ocean = findall(M .== 1.0)

    length(ocean) < np && error(
        "Only $(length(ocean)) ocean points found (need $np). " *
        "Widen the lon/lat range or reduce np.")

    sel = ocean[1:np]
    return DataFrame(x=x[sel], y=y[sel], fid=fid[sel], lon=lon[sel], lat=lat[sel])
end

"""
    deprecated_initial_positions(Γ; nf=10000, lon_rng=(-160.0,-159.0), lat_rng=(30.0,31.0))

Randomly assign initial positions in longitude,latitude ranges. Positions are 
expressed in, normalized, grid point units (x,y in the 0,nx and 0,ny range). 
To convert from longitude,latitude here we take advantage of the regularity 
of the 1 degree grid being used -- for a more general alternative, see the 
global ocean example.

**Warning**: this function assigns all particles to face 1 and uses a simple
linear lon/lat-to-grid mapping. It only works for small regions on face 1.
For regional or global initialization on the LLC90 grid, use `init_regional_3d`.
"""
function deprecated_initial_positions(Γ::NamedTuple, nf=10000, lon_rng=(-160.0,-159.0), lat_rng=(30.0,31.0), level=1)
   lon=lon_rng[1] .+(lon_rng[2]-lon_rng[1]).*rand(nf)
   lat=lat_rng[1] .+(lat_rng[2]-lat_rng[1]).*rand(nf)
   x=lon .+ (21. - Γ.XC[1][21,1])
   y=lat .+ (111. - Γ.YC[1][1,111])

   return DataFrame(:x => x, :y => y, :z => fill(level,nf),:fid => fill(1,nf))
end

end
