
lon180(x)=Float64(x>180.0 ? x-360.0 : x)
lon360(x)=Float64(x<0.0 ? x+360.0 : x)

function background()
    dx=0.1
    lon,lat,basemap=demo.get_basemap()
    fig = Figure(size = (1200, 800), backgroundcolor = :grey80)
    ax = Axis(fig[1, 1])
    im=image!(ax,lon[1,1]..lon[end,1],lat[1,1]..lat[1,end],basemap)
    #hidedecorations!(ax)
    fig,ax
end

"""
    global_plot1_yys(🔴::DataFrame)

Plot initial and final positions, superimposed on a globalmap of ocean depth log.

```
using Drifters, GLMakie
include("worldwide/global_ocean_circulation.jl")

x=DriftersDataset( data=(I=𝐼,df=tmp_🔴,), options=(plot_type=:global_plot1_yys,) )
fig,tt=plot(x)
fig

file_output_mp4=tempname()*".mp4"
record(fig, file_output_mp4, -50:nt, framerate = 25) do t
    tt[]=max(t,0)
end
```
"""
function global_plot1(🔴::DataFrame;
	time=0,xlims=(0.0,360.0),ylims=(-80.0,90.0),
    colormap=:linear_wcmr_100_45_c42_n256,
	colorrange=(-1300,00),add_colorbar=false)

    fig,ax=background()

    sort!(🔴, [:t, :ID])
    🔴 = unique(🔴, [:t, :ID])

    nt=length(unique(🔴.t))
    🔴_by_t=groupby(🔴, :t)
    all_ids = Int.(🔴_by_t[1].ID)
    np=length(all_ids)
	ii= np <= 10000 ? collect(1:np) : sort(randperm(np)[1:10000])
    rows_ii = ii
    tmp1=🔴_by_t[1][rows_ii,:lon].!==🔴_by_t[nt][rows_ii,:lon]
    tmp2=🔴_by_t[1][rows_ii,:lat].!==🔴_by_t[nt][rows_ii,:lat]
    jj_rows=rows_ii[findall(tmp1.*tmp2)]

    time==0 ? tt=Observable(nt) : tt=Observable(time)    
    for tx in -12:0 
        ttt=@lift(max(1,$tt+tx))
        lon_tt=@lift(lon360.(🔴_by_t[$ttt][jj_rows,:lon]))
        lat_tt=@lift(🔴_by_t[$ttt][jj_rows,:lat])
        d_tt=@lift(max.(🔴_by_t[$ttt][jj_rows,:d],Ref(-1200)))
        scatter!(ax,lon_tt,lat_tt,markersize=4.0,
        color=d_tt,colorrange=colorrange,colormap=colormap)
    end

    lon_t1=lon360.(🔴_by_t[1][jj_rows,:lon])
    lat_t1=🔴_by_t[1][jj_rows,:lat]
    scatter!(ax,lon_t1,lat_t1,markersize=1.0,color=:lightblue)

    limits!(ax,xlims...,ylims...)

    add_colorbar ? Colorbar(fig[1,2],colorrange=colorrange,colormap=colormap) : nothing

    return fig,tt
end

## Oscar example

function Oscar_plot(df=[]; plot_type="Oscar_plot", lon=[], lat=[], 
			color=:red, colormap=:thermal, colorrange=(0,10), markersize=0.5, 
			add_background=false,add_polygons=false,
			proj=[],lon0=0.0)
    
	f = Figure(size=(1200,800))
	ttl="Drifters.jl + Oscar model"
	ax = f[1, 1] = Axis(f, aspect = DataAspect(), title = ttl)
	pr_ax = (proj==[] ? ax : MeshArrays.ProjAxis(ax; proj=proj,lon0=lon0))
  
	if add_background
		BG=get_background()
		surf = surface!(pr_ax,BG.lon,BG.lat,0*BG.lat; color=BG.DD,
			  colorrange=BG.colorrange, colormap=BG.colormap, shading = NoShading)
	end
	if add_polygons
		pol=read_pol()
		lines!(pr_ax; polygons=pol,color=:black,linewidth=0.5)
	end
	proj==[] ? nothing : MeshArrays.grid_lines!(pr_ax;color=:grey10,linewidth=0.5)
		
	if isa(lon,Observable)
		xy=@lift(proj.($lon,$lat))
		sc=scatter!(ax,xy,markersize=markersize,color=color,colormap=colormap,colorrange=colorrange)
		isa(color,Symbol) ? nothing : Colorbar(f[1,2], sc, height = Relative(0.5))
	elseif !isempty(lon)
		xy=proj.(lon,lat)
		sc=scatter!(ax,xy,markersize=markersize,color=color,colormap=colormap,colorrange=colorrange)
		isa(color,Symbol) ? nothing : Colorbar(f[1,2], sc, height = Relative(0.5))
	else
		I_t=groupby(df,:t)
		times=collect(1:length(I_t))
		[scatter!(pr_ax,I_t[t].lon,I_t[t].lat,color=:blue,markersize=1) for t in times]
	  	scatter!(pr_ax,I_t[times[1]].lon,I_t[times[1]].lat,markersize=2,color=:red)
	end
  
	f
  end

function get_background()
	γ=MeshArrays.GridSpec(ID=:LLC90)
	λ=MeshArrays.interpolation_setup()
	hFacC=MeshArrays.GridLoadVar("hFacC",γ)
	μ=MeshArrays.land_mask(hFacC[:,1])

	Depth=MeshArrays.GridLoadVar("Depth",γ)
	DD=MeshArrays.Interpolate(μ*Depth,λ.f,λ.i,λ.j,λ.w)
	DD=reshape(DD,size(λ.lon))
	(lon=λ.lon,lat=λ.lat,DD=DD,colorrange=[-6000.0,6000.0],colormap=:grays)
end

function read_pol()
	fil=MeshArrays.demo.download_polygons("countries.geojson")
	pol=MeshArrays.read_polygons(fil)
end
