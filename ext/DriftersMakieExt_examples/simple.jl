
"""
    plot_start_end(I::Individuals)

Plot the initial and final positions as scatter plot in `lon,lat` or `x,y` plane.
"""
function plot_start_end(I::Individuals)
🔴_by_t = Drifters.DataFrames.groupby(I.🔴, :t)
set_theme!(theme_light())
fig=Figure(size = (900, 600))
ms=4
if hasproperty(🔴_by_t[1],:lon)
	a = Axis(fig[1, 1],xlabel="longitude",ylabel="latitude")		
	scatter!(a,🔴_by_t[1].lon,🔴_by_t[1].lat,color=:green2,markersize=ms)
	scatter!(a,🔴_by_t[end].lon,🔴_by_t[end].lat,color=:red,markersize=ms)
else
	a = Axis(fig[1, 1],xlabel="longitude",ylabel="latitude")		
	scatter!(a,🔴_by_t[1].x,🔴_by_t[1].y,color=:green2,markersize=ms)
	scatter!(a,🔴_by_t[end].x,🔴_by_t[end].y,color=:red,markersize=ms)
end
return fig
end

"""
    simple_plot1(I,ϕ)

```
using Drifters, CairoMakie
include("basics/random_flow_field.jl")
x=DriftersDataset( data=(I=I,ϕ=ϕ), options=(plot_type=:simple_plot1,) )
plot(x)
```
"""
function simple_plot1(I,ϕ)
	I_t = groupby(I, :t)
	nt=length(I_t)

	time = Observable(nt)
	xp=@lift( I_t[$time].x )
	yp=@lift( I_t[$time].y )

	siz=size(ϕ)
	xx=-0.5 .+ collect(1:siz[1])
	yy=-0.5 .+ collect(1:siz[2])
	ll=collect(-1.0:0.2:1.0)*maximum(ϕ)
	
	fig=Figure(size = (900, 600)); set_theme!(theme_light())
	a = Axis(fig[1, 1],xlabel="x",ylabel="y", title="Positions and Streamfunction")		
	contourf!(a, xx,yy, ϕ, levels=ll, colormap=:grayC)
	scatter!(a,I_t[1].x,I_t[1].y,color=:green2)
	scatter!(a,xp,yp,color=:red)

	fig
end

##

"""
    simple_plot2(I)

```
using Drifters, CairoMakie
include("basics/solid_body_rotation.jl")
x=DriftersDataset( data=(I=I,), options=(plot_type=:simple_plot2,) )
plot(x)
```
"""
function simple_plot2(I)
	nt=length(I.🔴.x)
	
	time = Observable(nt)
	xx=@lift( [I.🔴.x[1:$time];fill(NaN,nt-$time)] )
	yy=@lift( [I.🔴.y[1:$time];fill(NaN,nt-$time)] )
	zz=@lift( [I.🔴.z[1:$time];fill(NaN,nt-$time)] )
	
	set_theme!(theme_light())
	f=Figure(size = (900, 600))
	a = Axis3(f[1, 1],xlabel="x",ylabel="y",zlabel="z",
		title="Solid body rotation / Spiral example")		
	lines!(a,xx,yy,zz,linewidth=1.0,color=:black)
	scatter!(a,[I.🔴.x[1]],[I.🔴.y[1]],[I.🔴.z[1]],color=:red)
	scatter!(a,[I.🔴.x[nt]],[I.🔴.y[nt]],[I.🔴.z[nt]],color=:green)

	f
end
