
## Flow field parameters

using Dates

"""
    abstract type FlowFields

Data structure that provide access to flow fields (gridded as arrays) which will be 
used to interpolate velocities to individual locations later on (once embedded in
an `Individuals` struct). 

Following the C-grid convention also used in `MITgcm` (https://mitgcm.readthedocs.io) 
flow fields are expected to be staggered as follows: grid cell i,j has its center 
located at i-1/2,j-1/2 while the corresponding `u[i,j]` (resp. `v[i,j]) is 
located at i-1,j-1/2 (resp. i-1/2,j-1). 

By convention :
- time is in seconds but `T` can also be provided as a vector of `DateTime`
- position is in grid units. Local coordinates cover `(0 to nx, 0 to ny)` if the grid is of size `(nx, ny)`.
- velocity fields must be normalized accordingly to grid units (in 1/s rather than m/s)
  being used in a `FlowFields` constructors (whether `Array` or `MeshArray`)

```
uvArrays(u0,u1,v0,v1,T)
uvwArrays(u0,u1,v0,v1,w0,w1,T)
uvMeshArrays(u0,u1,v0,v1,T,update_location!)
uvwMeshArrays(u0,u1,v0,v1,w0,w1,T,update_location!)
```

Using the `FlowFields` constructor which gets selected by the type of `u0` etc.

For exmaple, with Array type in 3D :

```
F=FlowFields(u,u,v,v,0*w,1*w,[0.0,10.0])
```

Or, with `MeshArray`` type in 2D :

```
F=FlowFields(u,u,v,v,[0.0,10.0],func)
```

as shown in the online documentation examples.

"""
abstract type FlowFields end

struct uvArrays{Ty} <: FlowFields
    u0::Array{Ty,2}
    u1::Array{Ty,2}
    v0::Array{Ty,2}
    v1::Array{Ty,2}
    T::Union{Array{Ty},Array{DateTime}}
end

function FlowFields(u0::Array{Ty,2},u1::Array{Ty,2},
    v0::Array{Ty,2},v1::Array{Ty,2},T::Union{Array,Tuple}) where Ty
    #ensure T is a vector and enforce type
    T=time_set_type.(collect(T),Ty)
    #check array size concistency
    tst=prod([(size(u0)==size(tmp)) for tmp in (u1,v0,v1)])
    !tst ? error("inconsistent array sizes") : nothing
    #call constructor
    uvArrays(u0,u1,v0,v1,T)
end

struct uvwArrays{Ty} <: FlowFields
    u0::Array{Ty,3}
    u1::Array{Ty,3}
    v0::Array{Ty,3}
    v1::Array{Ty,3}
    w0::Array{Ty,3}
    w1::Array{Ty,3}
    T::Union{Array{Ty},Array{DateTime}}
end

"""
    FlowFields(;    u::Union{Array,Tuple}=[], v::Union{Array,Tuple}=[], w::Union{Array,Tuple}=[], 
                    period::Union{Array,Tuple}=[])

Simplified FlowFields constructor : using keywords, with time invariant flow fields.

```
uC, vC, _ = SimpleFlowFields(16)
to_C_grid!(uC,dims=1)
to_C_grid!(uC,dims=2)
F=FlowFields(u=uC,v=vC,period=(0,10.))
```
"""
function FlowFields(; u::Union{Array,Tuple}=[], v::Union{Array,Tuple}=[], w::Union{Array,Tuple}=[], 
    period::Union{Array,Tuple}=[])
    (isa(u,Tuple)||length(u[:])==2) ? (u0=u[1]; u1=u[2]) : (u0=u; u1=u)
    (isa(v,Tuple)||length(v[:])==2) ? (v0=v[1]; v1=v[2]) : (v0=v; v1=v)
    (isa(w,Tuple)||length(w[:])==2) ? (w0=w[1]; w1=w[2]) : (w0=w; w1=w)
    if isempty(period)
        @warn "period needs to be defined"
    end
    if !isempty(u0) && !isempty(v0)
        if !isempty(w0)
            FlowFields(u0,u1,v0,v1,w0,w1,period)
        else
            FlowFields(u0,u1,v0,v1,period)
        end
    else
        []
    end
end

"""
    to_C_grid!(x;dims=0)

```
if gridtype==:centered
    to_C_grid!(u0,dims=1)
    to_C_grid!(u1,dims=1)
    to_C_grid!(v0,dims=2)
    to_C_grid!(v1,dims=2)
    if !isempty(w0)
        to_C_grid!(w0,dims=3)
        to_C_grid!(w1,dims=3)
    end
end
```
"""
to_C_grid!(x;dims=0) = begin
    if (dims==1)&&(ndims(x)==2)
        x.=0.5*(circshift(x, (1,0))+x)
    elseif (dims==2)&&(ndims(x)==2)
        x.=0.5*(circshift(x, (0,1))+x)
    elseif dims==1
        x.=0.5*(circshift(x, (1,0,0))+x)
    elseif dims==2
        x.=0.5*(circshift(x, (0,1,0))+x)
    elseif dims==3
        x.=0.5*(circshift(x, (0,0,1))+x)
    end
end

function FlowFields(u0::Array{Ty,3},u1::Array{Ty,3},v0::Array{Ty,3},v1::Array{Ty,3},
    w0::Array{Ty,3},w1::Array{Ty,3},T::Union{Array,Tuple}) where Ty
    #ensure T is a vector and enforce type
    T=time_set_type.(collect(T),Ty)
    #check array size concistency
    tst=prod([(size(u0)==size(tmp)) for tmp in (u1,v0,v1)])
    tst=tst*prod([(size(u0)==size(tmp).-(0,0,1)) for tmp in (w0,w1)])
    !tst ? error("inconsistent array sizes") : nothing
    #call constructor
    uvwArrays(u0,u1,v0,v1,w0,w1,T)
end

struct uvMeshArrays{Ty} <: FlowFields
    u0::AbstractMeshArray{Ty,1}
    u1::AbstractMeshArray{Ty,1}
    v0::AbstractMeshArray{Ty,1}
    v1::AbstractMeshArray{Ty,1}
    T::Union{Array{Ty},Array{DateTime}}
    update_location!::Function
end

function FlowFields(u0::MeshArrays.MeshArray_wh,u1::MeshArrays.MeshArray_wh,
    v0::MeshArrays.MeshArray_wh,v1::MeshArrays.MeshArray_wh,
    T::Union{Array,Tuple},update_location!::Function)

    FlowFields(u0.MA,u1.MA,v0.MA,v1.MA,T,update_location!)
end

function FlowFields(u0::AbstractMeshArray{Ty,1},u1::AbstractMeshArray{Ty,1},
    v0::AbstractMeshArray{Ty,1},v1::AbstractMeshArray{Ty,1},
    T::Union{Array,Tuple},update_location!::Function) where Ty
    #ensure T is a vector and enforce type
    T=time_set_type.(collect(T),Ty)
    #check array size concistency
    tst=prod([(size(u0)==size(tmp))*(u0.fSize==tmp.fSize) for tmp in (u1,v0,v1)])
    !tst ? error("inconsistent array sizes") : nothing
    #call constructor
    uvMeshArrays(u0,u1,v0,v1,T,update_location!)
end

struct uvwMeshArrays{Ty} <: FlowFields
    u0::AbstractMeshArray{Ty,2}
    u1::AbstractMeshArray{Ty,2}
    v0::AbstractMeshArray{Ty,2}
    v1::AbstractMeshArray{Ty,2}
    w0::AbstractMeshArray{Ty,2}
    w1::AbstractMeshArray{Ty,2}
    T::Union{Array{Ty},Array{DateTime}}
    update_location!::Function
end

function FlowFields(u0::MeshArrays.MeshArray_wh,u1::MeshArrays.MeshArray_wh,
    v0::MeshArrays.MeshArray_wh,v1::MeshArrays.MeshArray_wh,
    w0::MeshArrays.MeshArray_wh,w1::MeshArrays.MeshArray_wh,
    T::Union{Array,Tuple},update_location!::Function)

    FlowFields(u0.MA,u1.MA,v0.MA,v1.MA,w0.MA,w1.MA,T,update_location!)
end

function FlowFields(u0::AbstractMeshArray{Ty,2},u1::AbstractMeshArray{Ty,2},
    v0::AbstractMeshArray{Ty,2},v1::AbstractMeshArray{Ty,2},
    w0::AbstractMeshArray{Ty,2},w1::AbstractMeshArray{Ty,2},
    T::Union{Array,Tuple},update_location!::Function) where Ty
    #ensure T is a vector and enforce type
    T=time_set_type.(collect(T),Ty)
    #check array size consistency
    tst=prod([(size(u0)==size(tmp))*(u0.fSize==tmp.fSize) for tmp in (u1,v0,v1)])
    tst=tst*prod([(size(u0)==size(tmp).-(0,1))*(u0.fSize==tmp.fSize) for tmp in (w0,w1)])
    !tst ? error("inconsistent array sizes") : nothing

    #call constructor
    uvwMeshArrays(u0,u1,v0,v1,w0,w1,T,update_location!)
end

"""
    defaults for Individuals constructor
"""

default_solver(prob::ODEProblem) = solve(prob,Tsit5(),reltol=1e-8,abstol=1e-8)

function ensemble_solver(prob::ODEProblem;solver=Tsit5(),reltol=1e-8,abstol=1e-8,safetycopy=false)
	u0 = prob.u0
	prob_func(prob,i,repeat) = remake(prob,u0=u0[i])
	indiv_prob = ODEProblem(prob.f,u0[1],prob.tspan,prob.p)
#	indiv_prob = _SDEProblem(prob.f,u0[1],prob.tspan,prob.p)
	ensemble_prob = EnsembleProblem(indiv_prob,prob_func=prob_func,safetycopy=safetycopy)
	solve(ensemble_prob, solver, reltol=reltol, abstol=abstol, trajectories=length(u0))
end

a=fill(0.0,1,1)
default_flowfields = uvArrays{Float64}(a,a,a,a,[0. 1.])
default_recorder = DataFrame(ID=Int[], x=Float64[], y=Float64[], t=Float64[])
default_postproc = (x->x)

"""
    struct Individuals{T,N}

- Data:           📌 (position),   🔴(record), 🆔 (ID), P (`FlowFields`)
- Functions:      🚄 (velocity),   ∫ (integration), 🔧(postprocessing)
- NamedTuples:    D (diagnostics),      M (metadata)

The velocity function 🚄 typically computes velocity at individual positions (📌 to start) within the 
specified space-time domain by interpolating gridded variables (provided via P). Individual trajectories 
are computed by integrating (∫) interpolated velocities through time. Normally, integration is done by 
calling ∫! which updates 📌 at the end and records results in 🔴 via 🔧. Ancillary data, for use in 
🔧 for example, can be provided in D and metadata stored in M.

Unicode cheatsheet:

- 📌=`\\:pushpin:<tab>`,          🔴=`\\:red_circle:<tab>`, 🆔=`\\:id:<tab>`
- 🚄=`\\:bullettrain_side:<tab>`, ∫=`\\int<tab>`,          🔧=`\\:wrench:<tab>`
- P=`\\itP<tab>`,                 D=`\\itD<tab>`,           M=`\\itM<tab>`

Simple constructors that use `FlowFields` to choose adequate defaults:

- Individuals(F::uvArrays,x,y)
- Individuals(F::uvwArrays,x,y,z)
- Individuals(F::uvMeshArrays,x,y,fid)
- Individuals(F::uvwMeshArrays,x,y,z,fid)

Further customization is achievable via keyword constructors:

```
df=DataFrame( ID=[], x=[], y=[], z=[], t = [])
I=Individuals{Float64,2}(📌=zeros(3,10),🆔=1:10,🔴=deepcopy(df))
I=Individuals(📌=zeros(3,2),🆔=collect(1:2),🔴=deepcopy(df))
```

Or via the plain text (or no-unicode) constructors:

```
df=DataFrame( ID=[], x=[], y=[], z=[], t = [])
I=(position=zeros(3,2),ID=1:2,record=deepcopy(df))
I=Individuals(I)
```
"""
Base.@kwdef struct Individuals{Ty,N}
   📌  ::Array{Ty,N} = Array{Ty,N}(undef, Tuple(Int.(zeros(1,N)))) #\:pushpin:<tab>
   🔴  ::DataFrame = similar(default_recorder) #\:red_circle:<tab>
   🆔   ::Array{Int,1} = Array{Int,1}(undef, 0) #\:id:<tab>
   🚄  ::Function = dxdt! #\:bullettrain_side:<tab>
   ∫   ::Function = default_solver #\int<tab>
   🔧  ::Function = default_postproc #\:wrench:<tab>
   P   ::FlowFields = default_flowfields #\itP<tab>
   D   ::NamedTuple = NamedTuple() #\itD<tab>
   M   ::NamedTuple = NamedTuple() #\itM<tab>
end

function Individuals(NT::NamedTuple)

    haskey(NT,:position) ? 📌=NT.position : 📌=Array{Float64,2}(undef, Tuple(Int.(zeros(1,2))))
    haskey(NT,:record) ? 🔴=NT.record : 🔴=similar(default_recorder)
    haskey(NT,:ID) ? 🆔=NT.ID : 🆔=collect(1:size(📌,2))    
    haskey(NT,:velocity) ? 🚄=NT.velocity : 🚄=dxdt!
    haskey(NT,:integration) ? ∫=NT.integration : ∫=default_solver
    haskey(NT,:postprocessing) ? 🔧=NT.postprocessing : 🔧=default_postproc
    haskey(NT,:parameters) ? P=NT.parameters : P=default_flowfields
    haskey(NT,:diagnostics) ? D=NT.diagnostics : D=NamedTuple()
    haskey(NT,:metadata) ? M=NT.metadata : M=NamedTuple()
    isa(📌,UnitRange) ? 📌=collect(📌) : nothing
    haskey(NT,:type) ? T=NT.type : T=eltype(📌)

    Individuals{T,ndims(📌)}(📌=📌,🔴=🔴,🆔=🆔,🚄=🚄,∫=∫,🔧=🔧,P=P,D=D,M=M)    
end

function Individuals(F::uvArrays,x,y, NT::NamedTuple = NamedTuple())
    📌=permutedims([[x[i];y[i]] for i in eachindex(x)])
    if length(📌)==1
        📌=📌[1]
        ∫=default_solver 
    else
        ∫=ensemble_solver
    end
    T=eltype(📌)

    🔴 = DataFrame(ID=Int[], x=Float64[], y=Float64[], t=Float64[])
    haskey(NT,:🔴) ? 🔴=NT.🔴 : nothing

    🔧 = postprocess_xy
    haskey(NT,:🔧) ? 🔧=NT.🔧 : nothing

    🆔=collect(1:size(📌,2))
    haskey(NT,:🆔) ? 🆔=NT.🆔 : nothing

    haskey(NT,:∫) ? ∫=NT.∫ : nothing

    D=NamedTuple()
    haskey(NT,:D) ? D=NT.D : nothing
    
    Individuals{T,ndims(📌)}(P=F,📌=📌,🔴=🔴,🆔=🆔,🚄=dxdt!,∫=∫,🔧=🔧,D=D)
end

function Individuals(F::uvwArrays,x,y,z, NT::NamedTuple = NamedTuple())
    📌=permutedims([[x[i];y[i];z[i]] for i in eachindex(x)])
    if length(📌)==1
        📌=📌[1]
        ∫=default_solver 
    else
        ∫=ensemble_solver
    end
    T=eltype(📌)

    🔴 = DataFrame(ID=Int[], x=Float64[], y=Float64[], z=Float64[], t=Float64[])
    haskey(NT,:🔴) ? 🔴=NT.🔴 : nothing

    function 🔧(sol,F::uvwArrays,D::NamedTuple;id=missing,T=missing)
        df=postprocess_xy(sol,F,D,id=id,T=T)
        if isa(sol,EnsembleSolution)
            np=length(sol)
            z=[[sol[i][1,3] for i in 1:np];[sol[3][1,end] for i in 1:np]]
        else
            z=sol[3,:]
        end
        df.z=z[:]
        return df
    end
    haskey(NT,:🔧) ? 🔧=NT.🔧 : nothing

    🆔=collect(1:size(📌,2))
    haskey(NT,:🆔) ? 🆔=NT.🆔 : nothing

    haskey(NT,:∫) ? ∫=NT.∫ : nothing

    D=NamedTuple()
    haskey(NT,:D) ? D=NT.D : nothing
    
    Individuals{T,ndims(📌)}(P=F,📌=📌,🔴=🔴,🆔=🆔,🚄=dxdt!,∫=∫,🔧=🔧,D=D)
end

function Individuals(F::uvMeshArrays,x,y,fid, NT::NamedTuple = NamedTuple())
    📌=permutedims([[x[i];y[i];fid[i]] for i in eachindex(x)])
    if length(📌)==1
        📌=📌[1]
        ∫=default_solver 
    else
        ∫=ensemble_solver
    end
    T=eltype(📌)

    🔴 = DataFrame(ID=Int[], x=Float64[], y=Float64[], fid=Int64[], t=Float64[])
    haskey(NT,:🔴) ? 🔴=NT.🔴 : nothing

    🔧 = postprocess_MeshArray
    haskey(NT,:🔧) ? 🔧=NT.🔧 : nothing

    🆔=collect(1:size(📌,2))
    haskey(NT,:🆔) ? 🆔=NT.🆔 : nothing

    haskey(NT,:∫) ? ∫=NT.∫ : nothing

    D=NamedTuple()
    haskey(NT,:D) ? D=NT.D : nothing

    Individuals{T,ndims(📌)}(P=F,📌=📌,🔴=🔴,🆔=🆔,🚄=dxdt!,∫=∫,🔧=🔧,D=D)
end

function Individuals(F::uvwMeshArrays,x,y,z,fid, NT::NamedTuple = NamedTuple())
    📌=permutedims([[x[i];y[i];z[i];fid[i]] for i in eachindex(x)])
    if length(📌)==1
        📌=📌[1]
        ∫=default_solver 
    else
        ∫=ensemble_solver
    end
    T=eltype(📌)

    🔴 = DataFrame(ID=Int[], x=Float64[], y=Float64[], z=Float64[], fid=Int64[], t=Union{Float64,DateTime}[])
    haskey(NT,:🔴) ? 🔴=NT.🔴 : nothing

    function 🔧(sol,F::uvwMeshArrays,D::NamedTuple;id=missing,T=missing)
        df=postprocess_MeshArray(sol,F,D,id=id,T=T)
        if isa(sol,EnsembleSolution)
            np=length(sol)
            z=[[sol.u[i][1][3] for i in 1:np];[sol.u[i][end][3] for i in 1:np]]
        else
            z=sol[3,:]
        end
        df.z=z[:]
        return df
    end
    haskey(NT,:🔧) ? 🔧=NT.🔧 : nothing

    🆔=collect(1:size(📌,2))
    haskey(NT,:🆔) ? 🆔=NT.🆔 : nothing

    haskey(NT,:∫) ? ∫=NT.∫ : nothing

    D=NamedTuple()
    haskey(NT,:D) ? D=NT.D : nothing

    Individuals{T,ndims(📌)}(P=F,📌=📌,🔴=🔴,🆔=🆔,🚄=dxdt!,∫=∫,🔧=🔧,D=D)
end

function time_in_seconds(T)
    t=if isa(T,DateTime)
        86400*datetime2julian.(T)
    else
        T
    end
end

function time_in_DateTime(T)
    t=if isa(T,DateTime)
        T
    else
        julian2datetime.(T./86400)
    end
end

function time_set_type(T,Ty)
    t=if isa(T,DateTime)
        T
    else
        convert(Ty,T)
    end
end

"""
    ∫!(I::Individuals,T::Tuple)

Displace simulated individuals continuously through space over time period T starting from position 📌. 

- This is typically achieved by computing the cumulative integral of velocity experienced by each individual along its trajectory (∫ 🚄 dt).
- The current default is `solve(prob,Tsit5(),reltol=1e-8,abstol=1e-8)` but all solver options from the [OrdinaryDiffEq.jl](https://github.com/SciML/OrdinaryDiffEq.jl) package are available.
- After this, `∫!` is also equipped to postprocess results recorded into 🔴 via the 🔧 workflow, and the last step in `∫!` consists in updating 📌 to be ready for continuing in a subsequent call to `∫!`.
"""
function ∫!(I::Individuals,T::Tuple)
    (; 🚄,📌,P, D, 🔧, 🆔, 🔴, ∫) = I

    TT=time_in_seconds.(T)
    prob = ODEProblem(🚄,📌, TT ,P)
#    prob = _SDEProblem(🚄,📌, TT ,P)
    sol = ∫(prob)

    tmp = 🔧(sol,P,D, id=🆔, T=T)
    isempty(🔴) ? np =0 : np=length(🆔)
    append!(🔴,tmp[np+1:end,:],promote=true)

    if isa(sol,EnsembleSolution)
        np=length(sol)
        📌[:] = deepcopy([sol[:,i].u[end] for i in 1:np])
        if isa(P,uvwMeshArrays)||isa(P,uvMeshArrays)
            [update_location!(i,P) for i in I.📌]
        end
    else
        nd=length(size(sol))
        nd==3 ? 📌[:,:] = deepcopy(sol[:,:,end]) : 📌[:] = deepcopy(sol[:,end])
    end

    sol
end

∫!(I::Individuals,T::Array) = ∫!(I::Individuals,(T[1],T[2]))

"""
    ∫!(I::Individuals)

Call ∫!(I::Individuals,I.P.T)
"""
∫!(I::Individuals) = ∫!(I::Individuals,I.P.T)

## Convenience Methods (size,show,similar)

Base.size(A::Individuals) = size(A.📌)

function Base.show(io::IO, I::Individuals)
    (; 🚄,📌,P, D, M, 🔧, 🆔, 🔴, ∫) = I
    printstyled(io, "  📌 details     = ",color=:normal)
    printstyled(io, "$(size(📌)) $(typeof(I).parameters[1])\n",color=:blue)
    printstyled(io, "  🔴 details     = ",color=:normal)
    printstyled(io, "$(size(🔴)) $(names(🔴))\n",color=:blue)
    printstyled(io, "  🆔 range       = ",color=:normal)
    printstyled(io, "$(extrema(🆔))\n",color=:blue)
    printstyled(io, "  🚄 function    = ",color=:normal)
    printstyled(io, "$(🚄)\n",color=:blue)
    printstyled(io, "  ∫  function    = ",color=:normal)
    printstyled(io, "$(∫)\n",color=:blue)
    printstyled(io, "  🔧 function    = ",color=:normal)
    printstyled(io, "$(🔧)\n",color=:blue)
    printstyled(io, "  P  details     = ",color=:normal)
    printstyled(io, "$(fieldnames(typeof(P)))\n",color=:blue)
  return
end

function Base.similar(I::Individuals)
    (; 🚄,📌,P, D, M, 🔧, 🆔, 🔴, ∫) = I
    T = typeof(I).parameters[1]
    N = ndims(I.📌)
    return Individuals{T,N}(📌=similar(📌),🔴=similar(🔴),🆔=similar(🆔),
                          🚄=🚄, ∫=∫, 🔧=🔧, P=P, D=D, M=M)
end

"""
    Base.diff(I::Individuals)

Difference in grid unit coordinates (dx,dy) between final and initial positions.
"""
function Base.diff(I::Individuals)
    f(x)=last(x).-first(x)
    🔴_by_ID = groupby(I.🔴, :ID)
    return combine(🔴_by_ID,nrow,:x => f => :dx,:y => f => :dy)
end

"""
    gcdist(I::Individuals)

Great circle distance (gcd in radians) between final and initial positions.
"""
function gcdist(I::Individuals)
    🔴_by_ID = groupby(I.🔴, :ID)
    tmp = combine(🔴_by_ID, 
    :lon => first => :lo1,:lon => last => :lo2,
    :lat => first => :la1,:lat => last => :la2)

    tmp.gcd=[gcdist(tmp.lo1[i],tmp.lo2[i],tmp.la1[i],tmp.la2[i]) for i in 1:size(tmp,1)]
    return tmp
end

"""
    gcdist(lo1,lo2,la1,la2)

Great circle distance (gcd in radians) between two positions.
"""
gcdist(lo1,lo2,la1,la2) = acos(min(1.0,max(-1.0,sind(la1)*sind(la2)+cosd(la1)*cosd(la2)*cosd(lo1-lo2))))

EarthRadius=6371e3 #in meters

"""
    abstract type AbstractDriftersDataset

Data structure used for plotting. See the docs for examples.
"""
abstract type AbstractDriftersDataset <: Any end

Base.@kwdef struct DriftersDataset <: AbstractDriftersDataset
    path :: String = tempdir()
    name :: String = "unknown"
    options :: NamedTuple = NamedTuple()
    data :: NamedTuple = NamedTuple()
end

#to avoid unicode try these:
OrdinaryDiffEq.OrdinaryDiffEqCore.solve!(I::Individuals,args...)=∫!(I::Individuals,args...)
DataFrames.groupby(I::Individuals,args...) = groupby(I.🔴,args...)
DataFrames.DataFrame(I::Individuals) = I.🔴
