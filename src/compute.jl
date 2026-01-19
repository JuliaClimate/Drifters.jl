
#needed to avoid allocations:
Flo=Union{Float32,Float64}
mydt(tim::Flo,T::Array{Float32,1})=(tim-T[1])/(T[2]-T[1])
mydt(tim::Flo,T::Array{Float64,1})=(tim-T[1])/(T[2]-T[1])

#needed to avoid allocations:
fSize(f::Array{Tuple{Int,Int}},i::Int) = (f[i][1],f[i][2])
fview(f::Array{Array{Float32,2},1},i::Int) = view(f[i],:,:)
fview(f::Array{Array{Float64,2},1},i::Int) = view(f[i],:,:)
fview(f::Array{Array{Float32,2},2},i::Int,j::Int) = view(f[i,j],:,:)
fview(f::Array{Array{Float64,2},2},i::Int,j::Int) = view(f[i,j],:,:)

"""
    dxdt!(du,u,p::uvwMeshArrays,tim)

Interpolate velocity from gridded fields (3D; with halos) to position `u`
(`x,y,z,fIndex`) to compute the derivative of position v time  `du_dt`.

# Extended help

```jldoctest; output = false
using Drifters
u,v,w,pos,func=vortex_flow_field(format=:MeshArray)
F=FlowFields(u,u,v,v,0*w,1*w,[0,3*pi],func)
I=Individuals(F,pos...)
∫!(I)

using CairoMakie
J=DriftersDataset( data=(I=I,), options=(plot_type=:simple_plot2,) )
f=plot(J)

ref=[9.35, 7.94, 1.29, 1.0] #hide
prod(isapprox.(I.📌,ref,atol=1.0)) # hide

# output

true
```
"""
function dxdt!(du::Array{T,1},u::Array{T,1},P::uvwMeshArrays,tim) where T
    dt=mydt(tim,P.T)
    g=P.u0.grid
    #
    while MeshArrays.location_is_out(u,g)
        P.update_location!(u)
    end

    x=u[1]
    y=u[2]
    z=u[3]
    fIndex = Int(u[4])
    (nx,ny)=fSize(g.fSize,fIndex)
    nz=size(P.u0,2)
    #
    (dx,dy,dz)=(x - floor(x),y - floor(y),z - floor(z))
    i_c = Int32(floor(x)) + 2
    j_c = Int32(floor(y)) + 2
    k_c = Int32(floor(z)) + 1
    #
    (i_w,i_e)=(i_c,i_c+1)
    (j_s,j_n)=(j_c,j_c+1)
    (k_l,k_r)=(k_c,k_c+1)
    #
    k_c=min(k_c,nz)
    k_l=min(k_l,nz)
    k_r=min(k_r,nz)
    k_c=max(k_c,1)
    k_l=max(k_l,1)
    k_r=max(k_r,1)
    #
    onemdt=(1.0-dt)
    onemdx=(1.0-dx)
    onemdy=(1.0-dy)
    onemdz=(1.0-dz)
    #interpolate u to position and time
    u0=fview(P.u0.f,fIndex,k_c)
    u1=fview(P.u1.f,fIndex,k_c)
    du[1]=onemdx*onemdt*u0[i_w,j_c]+
    dx*onemdt*u0[i_e,j_c]+
    onemdx*dt*u1[i_w,j_c]+
    dx*dt*u1[i_e,j_c]
    #interpolate v to position and time
    v0=fview(P.v0.f,fIndex,k_c)
    v1=fview(P.v1.f,fIndex,k_c)
    du[2]=onemdy*onemdt*v0[i_c,j_s]+
    dy*onemdt*v0[i_c,j_n]+
    onemdy*dt*v1[i_c,j_s]+
    dy*dt*v1[i_c,j_n]
    #interpolate w to position and time
    du[3]=onemdz*onemdt*P.w0.f[fIndex,k_l][i_c,j_c]+
    dz*onemdt*P.w0.f[fIndex,k_r][i_c,j_c]+
    onemdz*dt*P.w1.f[fIndex,k_l][i_c,j_c]+
    dz*dt*P.w1.f[fIndex,k_r][i_c,j_c]
    #leave face index unchanged
    du[4]=0.0
    #
    return du
end

"""
    dxdt!(du,u,p::uvMeshArrays,tim)

Interpolate velocity from gridded fields (2D; with halos) to position `u`
(`x,y,fIndex`) to compute the derivative of position v time  `du_dt`.

# Extended help

```jldoctest; output = false
using Drifters
u,v,w,pos,func=random_flow_field(format=:MeshArray)
F=FlowFields(u,u,v,v,[0,1.0],func)
I=Individuals(F,pos...)
∫!(I)

isa(I.📌,Vector)

# output

true
```
"""
function dxdt!(du::Array{T,1},u::Array{T,1},P::uvMeshArrays,tim) where T
    #compute positions in index units
    dt=mydt(tim,P.T)
    g=P.u0.grid
    #
    while MeshArrays.location_is_out(u,g)
        P.update_location!(u)
    end

    x=u[1]
    y=u[2]
    fIndex = Int(u[3])
    (nx,ny)=fSize(g.fSize,fIndex)
    #
    (dx,dy)=(x - floor(x),y - floor(y))
    i_c = Int32(floor(x)) + 2
    j_c = Int32(floor(y)) + 2
    #
    (i_w,i_e)=(i_c,i_c+1)
    (j_s,j_n)=(j_c,j_c+1)
    #interpolate u to position and time
    u0=fview(P.u0.f,fIndex)
    u1=fview(P.u1.f,fIndex)
    du[1]=(1.0-dx)*(1.0-dt)*u0[i_w,j_c]+
    dx*(1.0-dt)*u0[i_e,j_c]+
    (1.0-dx)*dt*u1[i_w,j_c]+
    dx*dt*u1[i_e,j_c]
    #interpolate v to position and time
    v0=fview(P.v0.f,fIndex)
    v1=fview(P.v1.f,fIndex)
    du[2]=(1.0-dy)*(1.0-dt)*v0[i_c,j_s]+
    dy*(1.0-dt)*v0[i_c,j_n]+
    (1.0-dy)*dt*v1[i_c,j_s]+
    dy*dt*v1[i_c,j_n]
    #leave face index unchanged
    du[3]=0.0
    #
    return du
end

"""
    location_is_out(u,nx,ny)

Test if location `u` is outside the gridded domain (`0,nx`, `0,ny`).
"""
location_is_out(u,nx,ny) = (u[1]>nx)&&(u[1]<0)&&(u[2]>ny)&&(u[2]<0)

"""
    dxdt!(du,u,P::uvwArrays,tim)

Interpolate velocity from gridded fields (3D; NO halos) to position `u`
(`x,y,z`) to compute the derivative of position v time  `du_dt`.

# Extended help

```jldoctest; output = false
using Drifters
u,v,w,pos=vortex_flow_field(format=:Array)
F=FlowFields(u,u,v,v,0*w,1*w,[0,3*pi])
I=Individuals(F,pos...)
∫!(I)

ref=[9.35, 7.93, 1.28] # hide
prod(isapprox.(I.📌,ref,atol=1.0)) # hide

# output

true
```
"""
function dxdt!(du::Array{T,1},u::Array{T,1},P::uvwArrays,tim) where T
    #compute positions in index units
    dt=mydt(tim,P.T)
    (nx,ny,nz) = size(P.u0)
    #
    while location_is_out(u,nx,ny)
        u[1]=( u[1]>nx ? u[1]-nx : u[1] )
        u[1]=( u[1]<0 ? u[1]+nx : u[1] )
        u[2]=( u[2]>ny ? u[2]-ny : u[2] )
        u[2]=( u[2]<0 ? u[2]+ny : u[2] )
    end
    #
    x=mod(u[1],nx)
    y=mod(u[2],ny)
    z=mod(u[3],nz)
    #
    dx=x - floor(x)
    dy=y - floor(y)
    dz=z - floor(z)
    #
    i_c = Int32(floor(x)) + 1
    j_c = Int32(floor(y)) + 1
    k_c = Int32(floor(z)) + 1
    #
    i_c==0 ? i_c=nx : nothing
    j_c==0 ? j_c=ny : nothing
    #k_c==0 ? k_c=nz : nothing
    #
    (i_w,i_e)=(i_c,i_c+1)
    x>=nx-1 ? (i_w,i_e)=(nx,1) : nothing
    (j_s,j_n)=(j_c,j_c+1)
    y>=ny-1 ? (j_s,j_n)=(ny,1) : nothing
    (k_l,k_r)=(k_c,k_c+1)
    #z>=nz-1 ? (k_l,k_r)=(nz,1) : nothing

    #interpolate u to position and time
    du[1]=(1.0-dx)*(1.0-dt)*P.u0[i_w,j_c,k_c]+
    dx*(1.0-dt)*P.u0[i_e,j_c,k_c]+
    (1.0-dx)*dt*P.u1[i_w,j_c,k_c]+
    dx*dt*P.u1[i_e,j_c,k_c]
    #interpolate v to position and time
    du[2]=(1.0-dy)*(1.0-dt)*P.v0[i_c,j_s,k_c]+
    dy*(1.0-dt)*P.v0[i_c,j_n,k_c]+
    (1.0-dy)*dt*P.v1[i_c,j_s,k_c]+
    dy*dt*P.v1[i_c,j_n,k_c]
    #interpolate w to position and time
    du[3]=(1.0-dz)*(1.0-dt)*P.w0[i_c,j_c,k_l]+
    dz*(1.0-dt)*P.w0[i_c,j_c,k_r]+
    (1.0-dz)*dt*P.w1[i_c,j_c,k_l]+
    dz*dt*P.w1[i_c,j_c,k_r]
    #
    return du
end

function dxdt!(du::Array{T,2},u::Array{T,2},P::uvwArrays,tim) where T
    [dxdt!(du[i],u[i],P,tim) for i=1:size(u,2)]
end

"""
    dxdt!(du,u,P::uvArrays,tim)

Interpolate velocity from gridded fields (2D; NO halos) to position `u`
(`x,y`) to compute the derivative of position v time  `du_dt`.

# Extended help

```jldoctest; output = false
using Drifters
u,v,w,pos=random_flow_field(format=:Array)
F=FlowFields(u,u,v,v,[0,1.0])
I=Individuals(F,pos...)
∫!(I)

isa(I.📌,Vector)

# output

true
```
"""
function dxdt!(du::Array{T,1},u::Array{T,1},P::uvArrays,tim) where T
    dt=mydt(tim,P.T)
    (nx,ny) = size(P.u0)
    #
    while location_is_out(u,nx,ny)
        u[1]=( u[1]>nx ? u[1]-nx : u[1] )
        u[1]=( u[1]<0 ? u[1]+nx : u[1] )
        u[2]=( u[2]>ny ? u[2]-ny : u[2] )
        u[2]=( u[2]<0 ? u[2]+ny : u[2] )
    end
    #
    x=mod(u[1],nx)
    y=mod(u[2],ny)
    #
    dx=x - floor(x)
    dy=y - floor(y)
    #
    i_c = Int32(floor(x)) + 1
    j_c = Int32(floor(y)) + 1
    #
    i_c==0 ? i_c=nx : nothing
    j_c==0 ? j_c=ny : nothing
    #
    (i_w,i_e)=(i_c,i_c+1)
    x>=nx-1 ? (i_w,i_e)=(nx,1) : nothing
    (j_s,j_n)=(j_c,j_c+1)
    y>=ny-1 ? (j_s,j_n)=(ny,1) : nothing
    #interpolate u to position and time
    du[1]=(1.0-dx)*(1.0-dt)*P.u0[i_w,j_c]+
    dx*(1.0-dt)*P.u0[i_e,j_c]+
    (1.0-dx)*dt*P.u1[i_w,j_c]+
    dx*dt*P.u1[i_e,j_c]
    #interpolate v to position and time
    du[2]=(1.0-dy)*(1.0-dt)*P.v0[i_c,j_s]+
    dy*(1.0-dt)*P.v0[i_c,j_n]+
    (1.0-dy)*dt*P.v1[i_c,j_s]+
    dy*dt*P.v1[i_c,j_n]
    #
    return du
end

"""
    dxy_dt_CyclicArray(du,u,P::NamedTuple,tim)

Nearest neighbor (?) velocity from gridded fields (2D; NO halos but
not needed when CyclicArrays is used to extend valid indice ranges).

# Extended help

_notes:_ spatial interpolation & temporal interpolation are lacking

```jldoctest; output = false
using Drifters
p=dirname(pathof(Drifters))
include(joinpath(p,"../examples/more/example_CyclicArray.jl"))
(x,y)=cyclicarray_example()

ref=[330.5 290.5]
prod(isapprox.([x[end] y[end]],ref,atol=1.0))

# output

true
```
"""
function dxy_dt_CyclicArray(du::Array{T,2},u::Array{T,2},P::NamedTuple,tim) where T
    np=size(du,2)
    xi,yi=(u[1,:],u[2,:])

    #not needed:
    #xi=floor.(P.xg[1,Int.(sign.(xi).*floor.(abs.(xi)))])+rem.(xi,1)
    #yi=floor.(P.yg[Int.(sign.(yi).*floor.(abs.(yi))),1])+rem.(yi,1)

    i=Int.(floor.(P.xg[1,Int.(floor.(xi))]))
    j=Int.(floor.(P.yg[Int.(floor.(yi)),1]))
    du[1,:]=[P.u[i[ii],j[ii]] for ii in 1:np]
    du[2,:]=[P.v[i[ii],j[ii]] for ii in 1:np]
    return du
end

"""
    dxy_dt_replay(du,u,p::DataFrame,t)

Interpolate velocity from MITgcm float_trajectories output and return
position increment `du`.

# Extended help

```jldoctest; output = false
using Drifters
p=dirname(pathof(Drifters))
include(joinpath(p,"../examples/more/detailed_look.jl"))
prod(isapprox.(sol[:,end],ref[:,end],atol=1.0))

# output

true
```
"""
function dxy_dt_replay(du,u,p::DataFrame,t)
    tt=t/3600.0
    tt0=Int32(floor(tt))
    w=tt-tt0
    du[1]=(1.0-w)*p[tt0+1,:uVel]+w*p[tt0+2,:uVel]
    du[2]=(1.0-w)*p[tt0+1,:vVel]+w*p[tt0+2,:vVel]
end
