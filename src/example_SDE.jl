module ex_SDE

using Statistics
import Drifters.SciMLBase: DiscreteCallback
import Distributions: Histogram, fit, Normal, cdf, pdf

## helper functions for the example

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
    kappa = similar(u, Float64)
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

"""
function particle_density

get a t-x heatmap for the parcel distribution. 
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
