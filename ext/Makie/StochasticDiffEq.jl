
begin
	## Plotting

	to_depth(z; scale=1.0)=-scale*z
	
	function plot_paths(;za=missing,zb=missing,dz=100,lw=0.5,scale=1.0,anom=false)
		np=size(za,1)
		pp=1:dz:np
		fig=Figure(size=(300,500)); ax=Axis(fig[1,1])
		if anom==true
			[lines!(to_depth(za[p,:].-za[p,1],scale=scale),color=:blue,linewidth=1) for p in pp]
		else
			[lines!(to_depth(za[p,:],scale=scale),color=:blue,linewidth=lw) for p in pp]
			ylims!(-scale,0.0)
		end
		if !ismissing(zb)
			ax=Axis(fig[1,2])
			if anom==true
				[lines!(to_depth(za[p,:].-za[p,1],scale=scale),color=:blue,linewidth=1) for p in pp]
			else
				[lines!(to_depth(zb[p,:],scale=scale),color=:red,linewidth=lw) for p in pp]		
				ylims!(-scale,0.0)
			end
		end
		fig
	end

	function plot_stats(st;T=missing)
		fig=Figure(size=(300,500)); ax=Axis(fig[1,1])
		dz=0.1; k=[0+dz*(i0-1) for i0 in 1:10].+dz/2
		z=to_depth(k)
		lines!(st[:,1],z)#,ylim=(-0.1,1.1))
		lines!(st[:,1].+st[:,2],z)
		lines!(st[:,1].-st[:,2],z)
		if !ismissing(T)
			kk=((1:length(T)).-0.5)./length(T)
			zz=to_depth(kk)
			lines!(T,zz,color=:black,linestyle=:dash)
		end
		fig
	end	

    function plot_EulerianModel(T,T0)
		kk=((1:length(T)).-0.5)./length(T)
		zz=to_depth(kk)
		fig=Figure(size=(300,500)); ax=Axis(fig[1,1])
	    lines!(T0,zz); lines!(T,zz)
	    return fig
	end
	
	"Plotting functions"
end
