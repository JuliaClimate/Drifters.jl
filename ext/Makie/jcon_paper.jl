
## Drifters as plotted in JuliaCon Proceedings paper

EarthRadius=6371e3 #in meters
res=1/2 #resolution if test case
lola(x,y)=(-100+x*res,17+y*res) #convert x/y to lon/lat

"""
    plot_drifters_jcon(gdf ; prefix="",pol=[],xlims=(-180.0,180.0),ylims=(-90.0,90.0),vmax=10.0)

```
include("LoopCurrent_replay.jl")
LoopC=DriftersDataset( data=(gdf=gdf,), options=(plot_type="jcon_drifters",
				prefix=prefix,xlims=(-98,-78),ylims=(18,31),pol=pol) )
plot(LoopC)
```
"""
function plot_drifters_jcon(gdf ; 	plot_type="jcon_drifters", prefix="",pol=[],
									xlims=(-180.0,180.0),ylims=(-90.0,90.0),vmax=10.0)
    fi00()=Figure(size=(900,600),fontsize=24, textcolor=:grey90)
	fi0=with_theme(fi00,theme_dark()) 
    ax0=Axis(fi0[1,1],xlabel="longitude",ylabel="latitude",
        title=prefix*"surface drifter trajectories and speed (m/s)")
    cr=(0,2.5); cm=:speed

    m=ones(maximum([size(D,1) for D in gdf]))
    for D in gdf
		if in("longitude",names(D))
			x=D.longitude
			y=D.latitude
		else
			tmp=lola.(D.x,D.y)
			x=[x[1] for x in tmp]
			y=[x[2] for x in tmp]
		end

        if length(x) > 10
            dt=(in("t",names(D)) ? diff(D.t)[1] : diff(D.time)[1].value/1000)
            v=EarthRadius/dt*[gcdist(x[i],x[i+1],y[i],y[i+1]) for i in [1:length(x)-1 ; length(x)-1]]

            m.=1.0
            sum(v.>vmax)>0 ? m[findall(v.>vmax)].=NaN : nothing
            n=1:length(v)
            lines!(x.*m[n],y.*m[n],color=v.*m[n],colorrange=cr, colormap=cm)
        end
    end

    xlims!(xlims...); ylims!(ylims...)
	isa(pol,Array) ? nothing : lines!(pol,color=:mediumpurple,linewidth=4)
	Colorbar(fi0[1,2], colorrange=cr, colormap=cm, height = Relative(0.65))
	
	fi0
end
