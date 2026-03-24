module ex_SDE

using Statistics
using Distributions
if Base.find_package("DifferentialEquations") !== nothing
    @eval using DifferentialEquations
end

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

_at(x, t) = x(t)                 # for functions/functors
_at(x::Number, t) = x            # for constants
"""
function g_erf(u,p,t)

Use 1-erf as diffusivity, and g = \sqrt(2\kappa)
"""
function g_erf(u,p,t)
    μp, σp = p
    μ = _at(μp, t)
    σ = _at(σp, t)
    d = Normal(μ, σ)
    return sqrt.((1.0 .-cdf(d, u)).*2)
end

"""
function f_gauss(u,p,t)

Use 1-erf as diffusivity, and f = d\kappa/du
"""
function f_gauss(u,p,t)
    μp, σp = p
    μ = _at(μp, t)
    σ = _at(σp, t)
    d = Normal(μ, σ)
    return -pdf(d, u)
end

"""
function hovmoller_density

get a t-x heatmap for the parcel distribution. 
"""
function hovmoller_density(sol; nbins=20, xmin=nothing, xmax=nothing, normalize=true)
    ts = sol.t
    us = sol.u                      # each entry: vector of particle positions at time ts[k]

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
