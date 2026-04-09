module ex_SDE

using Statistics
import Drifters.SciMLBase: DiscreteCallback
import Distributions: Histogram, fit, Normal, cdf, pdf

## helper functions for the example

DEFAULT_KAPPA = [0.0, 0.019719386, 0.030045405, 0.03708679, 0.04244723, 0.036685575, 0.018587558, 0.0039903168, 
	0.001206897, 0.0006229491, 0.00027001562, 0.000102254, 3.624377f-5, 1.6010854f-5, 1.0576804f-5, 
	1.0f-5, 1.0f-5, 1.0f-5, 1.0f-5, 1.0f-5, 1.0f-5, 1.0f-5, 1.0f-5, 1.0f-5, 1.0f-5, 1.0f-5, 1.0f-5, 1.0f-5, 
	1.0f-5, 1.0f-5, 1.0f-5, 1.0f-5, 1.0f-5, 1.0f-5, 1.0f-5, 1.0f-5, 1.0f-5, 1.0f-5, 1.0f-5, 1.0f-5, 1.0f-5, 
	1.0f-5, 1.0f-5, 1.0f-5, 1.0f-5, 1.0f-5, 1.0f-5, 1.0f-5, 1.0f-5, 1.0f-5]
# GGL90Kappa in early winter, from a grid cell. 

"""
    fold_tails(z)

"""
function fold_tails(z)
    while !isempty(findall( xor.(z.>1.0,z.<0.0) ))
        z[findall(z.<0.0)].=-z[findall(z.<0.0)]
        z[findall(z.>1.0)].=(2.0 .- z[findall(z.>1.0)])
    end
end

## 

_at(x, t) = x(t)                 # for functions/functors
_at(x::Number, t) = x            # for constants

"""
function kappa_erf(u,p,t)

Use 1-erf as diffusivity, and the factor of dW is
        g = sqrt(2kappa)
A mirror image above the sea surface is created 
to handle the surface boundary correctly. 
"""
function kappa_erf(u,p,t)
    μp, σp = p
    μ = _at(μp, t)
    σ = _at(σp, t)
    d = Normal(μ, σ)
    return (1.0 .-cdf(d, u))
end

g_erf(u,p,t) = sqrt.(kappa_erf(u,p,t).*2.)

"""
function f_gauss(u,p,t)

Use 1-erf as diffusivity, and the drift velocity
        f = dkappa/du
A mirror image above the sea surface is created 
to handle the surface boundary correctly. 
"""
function f_gauss(u,p,t)
    μp, σp = p
    μ = _at(μp, t)
    σ = _at(σp, t)
    d = Normal(μ, σ)
    return (pdf(d, -u).-pdf(d, u))
end

"""
function kappa_piecewise(u,p,t)

Diffusivity is one above depth - thickness,
is zero below depth,
and is a linear transition between the two in 
the middle. 
"""
function kappa_piecewise(u,p,t)    
    depth_p,thickness_p,kappa0_p = p
    depth = _at(depth_p,t)
    thickness = _at(thickness_p,t)
    kappa0 = _at(kappa0_p,t)
    kappa = kappa0*ifelse.(u.>depth,
        0,
        ifelse.(
            u.<depth-thickness,
            1,
            (depth .- u)./thickness
            ))
    return kappa
end

#this represent 2K=g*g Equation 
g_piecewise(u,p,t) = sqrt.(kappa_piecewise(u,p,t).*2.)

g_piecewise_3d(u,p,t) = [0.0,0.0,g_piecewise(u[3],p,t)]

"""
function f_piecewise(u,p,t)

Drift corresponding to  kappa_piecewise,
just a constant drift backinto the mixed layer
of 1/thickness. 1 has unit of diffusivity. 
"""
function f_piecewise(u,p,t)
    depth_p,thickness_p,kappa0_p = p
    depth = _at(depth_p,t)
    thickness = _at(thickness_p,t)
    kappa0 = _at(kappa0_p,t)
    #this represents f=dK/dz
    drift = -kappa0*Float64.((u.>depth-thickness).&(u.<depth))./thickness
end

f_piecewise_3d(u,p,t) = [0.0,0.0,f_piecewise(u[3],p,t)]

"""
    particle_density(sol; nbins=20, xmin=nothing, xmax=nothing, normalize=true)

Compute time-depth histogram of particle distribution. 

```
SDE = Base.get_extension(Drifters, :DriftersStochasticDiffEqExt)

fig=Figure()
for a in 1:2
    zi=(a==1 ? IC.u₀a : IC.u₀b)
    z,sol=SDE.solve_paths(zi)
    ρ,x_centers=ex_SDE.particle_density(sol)
    t_centers=1:size(ρ,2)

    ax=Axis(fig[1,a]);
    heatmap!(t_centers,x_centers,permutedims(ρ))
    ylims!(0,1); ax.yreversed[]=true
end
fig
```
"""
function particle_density(sol; nbins=20, xmin=nothing, xmax=nothing, normalize=true)
    ts = sol.t
    us = sol.u 

    allmins = minimum(minimum, us)
    allmaxs = maximum(maximum, us)
    xmin = isnothing(xmin) ? allmins : xmin
    xmax = isnothing(xmax) ? allmaxs : xmax

    x_edges = collect(range(xmin, xmax; length=nbins+1))
    x_centers = (x_edges[1:end-1] .+ x_edges[2:end]) ./ 2
    dx = x_edges[2] - x_edges[1]

    nt = length(ts)
    ρ = zeros(nbins, nt)

    for (k, uk) in enumerate(us)
        h = fit(Histogram, uk, x_edges)
        counts = Float64.(h.weights)
        if normalize
            ρ[:, k] .= counts ./ (sum(counts) * dx)
        else
            ρ[:, k] .= counts
        end
    end

    return ρ,x_centers
end

"""
function surface_reflect

A function for creating discrete callback function to 
prevent particles
leaving the surface
"""
function surface_reflect()
    condition(u,t,integrator) = true
    function affect!(integrator)
        integrator.u .= abs.(integrator.u)   
    end
    cb = DiscreteCallback(condition, affect!)
    return cb
end

"""
function surface_reflect

A function for creating discrete callback
functions to prevent particles leaving
the surface and the bottom of the mixed layer,
similar to fold_tails. 
h(t) is the mixed layer depth. 
"""
function surface_and_bottom_reflect(h::Function)
    condition(u,t,integrator) = true
    function affect!(integrator)
        bottom_depth = h(integrator.t)
        integrator.u .= abs.(integrator.u)
        integrator.u .= bottom_depth.-abs.(bottom_depth.-integrator.u)
    end
    cb = DiscreteCallback(condition, affect!)
    return cb
end

"""
    function mix_neighbors(za,ca,zb,cb,p)

Mix properties between cells that are within the same grid cell (`dz`).
"""
function mix_neighbors(za,ca,zb,cb,p)
    dz=0.1
    for i0 in 1:10
        z0=dz*(i0-1)
        z1=dz*i0
        ia=findall((za.>z0).*(za.<=z1))
        ib=findall((zb.>z0).*(zb.<=z1))
        tmp=mean([ca[ia];cb[ib]])
        ca[ia].=(1-p)*ca[ia] .+ p*tmp
        cb[ib].=(1-p)*cb[ib] .+ p*tmp
    end
end

function gridded_stats(IC::NamedTuple)
    (; u₀a,u₀b,ca,cb,np) = IC
    gridded_stats(IC.u₀a,IC.ca,IC.u₀b,IC.cb)
end
    
function gridded_stats(za,ca,zb,cb)
    out=zeros(10,2)
    dz=0.1
    t=size(za,2)
    for i0=1:10
        z0=0+dz*(i0-1)
        ia=findall( (za[:,t].>z0).*(za[:,t].<=z0+dz) );
        ib=findall( (zb[:,t].>z0).*(zb[:,t].<=z0+dz) );
        tmp=[ca[ia,t];cb[ib,t]]
        out[i0,1]=mean(tmp)
        out[i0,2]=std(tmp)
    end
    out
end

# initial conditions
function initial_conditions(np=10000)	
	u₀a=0.5*rand(np)
	ca=zeros(np)
	u₀b=0.5 .+ 0.5*rand(np)
	cb=ones(np)
    (u₀a=u₀a,u₀b=u₀b,ca=ca,cb=cb,np=np)
end

function ocean_reflect()
    condition(u,t,integrator) = true
    function affect!(integrator)
        integrator.u .= -abs.(integrator.u)   # reflect all components (vector/array-safe)
    end
    cb = DiscreteCallback(condition, affect!)
    return cb
end

function build_interpolator(rf,kappa)
	kap = [kappa[2];kappa[2:end]; 0.0]
	rff = [1000;rf[2:end]]# it is all the same diffusivity above surface
	itp = linear_interpolation(reverse(rff), reverse(kap))
	return itp
end

function build_kappa_profile(itp)
	kappa_func(u,p,t) = itp.(u)
	f_func(u,p,t) = Float64[Interpolations.gradient(itp, z)[1] for z in u]
	g_func(u,p,t) = sqrt.(max.(itp.(u), 0.0) .* 2.0)
	return kappa_func,f_func,g_func
end

function create_sde_updator(kappa,rf)
	itp = build_interpolator(rf,kappa)
	κ,𝒻,ℊ = build_kappa_profile(itp)
	function updator(u_0,tf;dt=10)
		tspan = (0.0,tf)
		p = nothing
		prob = SDEProblem(𝒻,ℊ,u_0,tspan,p)
		sol=solve(prob,EulerHeun(),dt=dt,callback = ocean_reflect());
		return sol.u[end]
	end
	return updator
end

# parameter configuration

function default_parameters()
    tspan = (0.0,1.0) #time span in s
    dt=1e-3 #output frequency in s
    configuration=:piecewise #configration for K profile
    mldepth(t)=0.5 #MLD depth in normalized units (0 to 1)
    thickness(t)=0.1 #transition layer thickness in normalized units (0 to 1)
    seafloor(t)=1.0 #sea floor depth in normalized units (0 to 1)
    mlkappa(t)=0.001 #diffusivity in m2/s
    depthscale(t)=2000 #depth scale in m
    params=(tspan=tspan,dt=dt,mlkappa=mlkappa,depthscale=depthscale,
        mldepth=mldepth,thickness=thickness,seafloor=seafloor,
        configuration=configuration)
end

## Eulerian Model for comparison

function EulerianModel(nt=1)
    N=20
    dt=1e-4
    dx=1.0/2/N
    T=[zeros(N);ones(N)]
    T0=deepcopy(T)
    for tt in 1:nt
        dTr=(circshift(T,-1)-T); dTr[end]=0;
        dTl=(T-circshift(T,+1)); dTl[1]=0;
        T.+=(dTr-dTl)*dt/dx/dx
    end
    return T,T0
end

end
