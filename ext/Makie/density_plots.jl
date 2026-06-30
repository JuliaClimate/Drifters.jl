
function Drifters_density2d_plot(df; time=missing)
	gdf=groupby(df,:year)
    tt=Int(floor(length(gdf)/2))
    t=(!ismissing(time) ? time : tt)

	lon,lat,n=Drifters.histogram2d(DataFrame(gdf[1]),lon=-179:2:179,lat=-89:2:89)
	_,_,n_mid=Drifters.histogram2d(DataFrame(gdf[t]),lon=-179:2:179,lat=-89:2:89)
	_,_,n_end=Drifters.histogram2d(DataFrame(gdf[end]),lon=-179:2:179,lat=-89:2:89)
	fi=Figure(size=(900,400)); 
	ax1=Axis(fi[1,1],title="start"); heatmap!(ax1,lon,lat,log10.(n))
	ax2=Axis(fi[1,2],title="time $t"); heatmap!(ax2,lon,lat,log10.(n_mid))
	ax3=Axis(fi[1,3],title="end"); heatmap!(ax3,lon,lat,log10.(n_end))
	#lo=[-100,0]; la=[10,50]; [(ax.limits = (lo,la)) for ax in (ax1,ax2,ax3)]
	fi
end

#plot(DriftersDataset(I,options=(plot_type=:density_z)))
function Drifters_density1d_plot(df; time=missing)
    np=length(unique(df.ID))
	gdf=groupby(df,:year)
    tt=Int(floor(length(gdf)/2))
    t=(!ismissing(time) ? time : tt)

	z,n=Drifters.histogram1d(-gdf[1].d,-10:10:5000)
	_,n_mid=Drifters.histogram1d(-gdf[t].d,-10:10:5000)
	_,n_end=Drifters.histogram1d(-gdf[end].d,-10:10:5000)
	fi=Figure(size=(900,400)); 
	ax1=Axis(fi[1,1],yreversed=true); barplot!(ax1,z,n,direction=:x)
	ax2=Axis(fi[1,2],yreversed=true); barplot!(ax2,z,n_mid,direction=:x)
	ax3=Axis(fi[1,3],yreversed=true); barplot!(ax3,z,n_end,direction=:x)
	de=[0,500]; nu=[0,2*np]; [(ax.limits = (nu,de)) for ax in (ax1,ax2,ax3)]
	fi
end


