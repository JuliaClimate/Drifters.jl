# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,jl:light
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.4'
#       jupytext_version: 1.2.4
#   kernelspec:
#     display_name: Julia 1.3.1
#     language: julia
#     name: julia-1.3
# ---

# # Single Particle Simulation Example
#
# In this example we simulate the trajectory of a particle in an idealized representation of an ocean eddy. To start, we use solid body rotation around a central point. As an exercise left to the user, directions are provided e.g. to add a convergence / divergence term. 
#
# - 1. setup the software and initialize example
# - 2. simulate trajectories & plot results
# - 3. experiment with parameters (user)

# + {"slideshow": {"slide_type": "subslide"}, "cell_type": "markdown"}
# ### For More Documentation 
#
# - <https://docs.juliadiffeq.org/latest> 
# - <https://en.wikipedia.org/wiki/Displacement_(vector)>
# - <https://juliaclimate.github.io/IndividualDisplacements.jl/dev>
# - <https://juliaclimate.github.io/MeshArrays.jl/dev>
#

# + {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ## 1.1 Import Software, i.e. `Julia` Packages
# -

using OrdinaryDiffEq, Plots
using IndividualDisplacements, MeshArrays

# + {"slideshow": {"slide_type": "subslide"}, "cell_type": "markdown"}
# ## 1.2  Gridded Domain Example
#
# - define `SetPeriodicDomain` function, which uses `MeshArrays.jl`
# - call `SetPeriodicDomain` function with a chosen grid size; e.g. `np=16`

# +
np=16

Γ=setup_periodic_domain(np);

# + {"slideshow": {"slide_type": "subslide"}, "cell_type": "markdown"}
# ## 1.3 Define Time Period & Velocity Field
#
# - define time range
# - define velocity field(s)
# - store in `𝑃` (dictionary) with grid variables

# + {"slideshow": {"slide_type": "skip"}}
#time range
t0=0.0
t1=0.95*2*pi
#t1=2.95*2*pi

#solid-body rotation around central location
i=Int(np/2+1)
u=-(Γ["YG"].-Γ["YG"][1][i,i])
v=(Γ["XG"].-Γ["XG"][1][i,i])

#add some convergence to / divergence from central location
d=0.0 
#d=-0.10
u=u+d*(Γ["XG"].-Γ["XG"][1][i,i])
v=v+d*(Γ["YG"].-Γ["YG"][1][i,i])

#store everything in a dictionnary
𝑃=Dict("u0" => u, "u1" => u, "v0" => v, "v1" => v, "t0" => t0, "t1" => t1)
𝑃=merge(𝑃,Γ);

# + {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ## 1.4 Define initial position and time period
# -

u0=np*[1/3,1/3]
du=fill(0.0,2)
𝑇 = (𝑃["t0"],𝑃["t1"]);

# + {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ## 2.1 solve for particle trajectory
#
# - `ODEProblem` formulates the differential equation along with the time period `𝑇`, parameters `𝑃`
# - `solve` then performs the integration over `𝑇`, starting from `u0`
#
# _Try `?ODEProblem` or `?solve` for additional documentation_

# +
prob = ODEProblem(⬡,u0,𝑇,𝑃)
sol = solve(prob,Tsit5(),reltol=1e-8)

x,y=sol[1,:],sol[2,:]
nt=length(x)

# + {"slideshow": {"slide_type": "slide"}, "cell_type": "markdown"}
# ## 2.2 visualize particle trajectory
#
# - define `myplot` convenience function
# - generate animation using `myplot`
#
# _Want a single plot? try uncommenting one line at a time in the next cell below_

# + {"slideshow": {"slide_type": "-"}}
myplot(i)=plot(x[1:i],y[1:i],linewidth=2,arrow = 2,
    title="Solid body rotation / Spiral example",leg=false,
    xaxis="x",yaxis="y",xlims=(0,np),ylims=(0,np))


# + {"slideshow": {"slide_type": "subslide"}}
p=Int(ceil(nt/100))
anim = @animate for i ∈ 1:p:nt
    myplot(i)
end
pth=tempdir()*"/"
gif(anim, pth*"SolidBodyRotation.gif", fps = 15)
# + {"slideshow": {"slide_type": "subslide"}}
plt=myplot(nt)
scatter!(plt,[u0[1]],[u0[2]])
savefig(plt,pth*"SolidBodyRotation.png")
# -

