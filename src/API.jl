
"""
    defaults for Individuals constructor
"""

day=86400.0
mon=365/12*day
OneMonth=[-0.5*mon,0.5*mon]

solver_default(prob) = solve(prob,Euler(),dt=day)
param_default = ( 𝑇=OneMonth , 🔄=(x->x), u0=[], u1=[], v0=[], v1=[])
rec_default = DataFrame(fill(Float64, 7),[:ID, :x, :y, :t, :lon, :lat, :fid])
postprocess_default = (x->x)

"""
    struct Individuals{T}

- Data:           📌 (position),   🔴(record),           🆔 (ID)
- Functions:      🚄 (velocity),   ∫ (integration), 🔧(postprocessing)
- NamedTuples:    𝑃  (parameters), 𝐷 (diagnostics),      𝑀 (metadata)

Default keyword constructor example:

```
df=DataFrame( ID=[], x=[], y=[], z=[], t = [])
𝐼=Individuals{Float64}(📌=zeros(3,10),🆔=1:10,🔴=deepcopy(df))
𝐼=Individuals(📌=zeros(3,2),🆔=collect(1:2),🔴=deepcopy(df))
```

Plain text (or no-unicode) constructor example:

```
df=DataFrame( ID=[], x=[], y=[], z=[], t = [])
I=(position=zeros(3,2),ID=1:2,record=deepcopy(df))
I=Individuals(I)
```

Keyword cheatsheet:

- 📌=`\\:pushpin:<tab>`,        🔴=`\\:red_circle:<tab>`, 🆔=`\\:id:<tab>`
- 🚄=`\\bullettrain_side<tab>`, ∫=`\\int<tab>`,           🔧=`\\wrench<tab>`
- 𝑃=`\\itP<tab>`,               𝐷=`\\itD<tab>`,            𝑀 =`\\itM<tab>`
"""
Base.@kwdef struct Individuals{T}
   📌  ::Array{T,2} = Array{T,2}(undef, Tuple(Int.(zeros(1,2)))) #\:pushpin:<tab>
   🔴  ::DataFrame = rec_default #\:red_circle:<tab>
   🆔   ::Array{Int,1} = Array{Int,1}(undef, 0) #\:id:<tab>
   🚄  ::Function = dxy_dt #\bullettrain_side<tab>
   ∫   ::Function = solver_default #\int<tab>
   🔧  ::Function = postprocess_default #\wrench<tab>
   𝑃   ::NamedTuple = param_default #\itP<tab>
   𝐷   ::NamedTuple = NamedTuple() #\itD<tab>
   𝑀   ::NamedTuple = NamedTuple() #\itM<tab>vec
end

"""
    Individuals(NT::NamedTuple)

Constructor that uses a NamedTuple with only plain text keywords (i.e. no-unicode needed).

```
df=DataFrame( ID=[], x=[], y=[], z=[], t = [])
I=(position=zeros(3,2),ID=1:2,record=deepcopy(df))
I=Individuals(I)
```
"""
function Individuals(NT::NamedTuple)

    haskey(NT,:position) ? 📌=NT.position : 📌=Array{Float64,2}(undef, Tuple(Int.(zeros(1,2))))
    haskey(NT,:record) ? 🔴=NT.record : 🔴=rec_default
    haskey(NT,:ID) ? 🆔=NT.ID : 🆔=Array{Int,1}(undef, 0)
    haskey(NT,:velocity) ? 🚄=NT.velocity : 🚄=dxy_dt
    haskey(NT,:integration) ? ∫=NT.integration : ∫=solver_default
    haskey(NT,:postprocessing) ? 🔧=NT.postprocessing : 🔧=postprocess_default
    haskey(NT,:parameters) ? 𝑃=NT.parameters : 𝑃=param_default
    haskey(NT,:diagnostics) ? 𝐷=NT.diagnostics : 𝐷=NamedTuple()
    haskey(NT,:metadata) ? 𝑀=NT.metadata : 𝑀=NamedTuple()
    isa(📌,UnitRange) ? 📌=collect(📌) : nothing
    haskey(NT,:type) ? T=NT.type : T=eltype(📌)

    Individuals{T}(📌=📌,🔴=🔴,🆔=🆔,🚄=🚄,∫=∫,🔧=🔧,𝑃=𝑃,𝐷=𝐷,𝑀=𝑀)    
end

"""
    ∫!(𝐼::Individuals,𝑇::Tuple)

Displace simulated individuals continuously through space over time period 𝑇 starting from position 📌. 

- This is typically achieved by computing the cumulative integral of velocity experienced by each individual along its trajectory (∫ 🚄 dt).
- The current default is `solve(prob,Euler(),dt=day)` but all solver options from the [OrdinaryDiffEq.jl](https://github.com/SciML/OrdinaryDiffEq.jl) package are available.
- After this, `∫!` is also equiped to postprocess results recorded into 🔴 via the 🔧 workflow, and the last step in `∫!` consiste in updating 📌 to be ready for continuing in a subsequent call to `∫!`.
"""
function ∫!(𝐼::Individuals,𝑇::Tuple)
    @unpack 🚄,📌,𝑃, 🔧, 🆔, 🔴, ∫ = 𝐼

    prob = ODEProblem(🚄,📌, 𝑇 ,𝑃)
    sol = ∫(prob)

    tmp = 🔧(sol,𝑃, id=🆔, 𝑇=𝑇)

    isempty(🔴) ? np =0 : np=length(🆔)
    append!(🔴,tmp[np+1:end,:])

    📌[:,:] = deepcopy(sol[:,:,end])
end

## Convenience Methods (size,show,similar)

Base.size(A::Individuals) = size(A.📌)

function Base.show(io::IO, 𝐼::Individuals) where {T}
    @unpack 🚄,📌,𝑃, 𝐷, 𝑀, 🔧, 🆔, 🔴, ∫ = 𝐼
    printstyled(io, "  📌 details     = ",color=:normal)
    printstyled(io, "$(size(📌)) $(typeof(𝐼).parameters[1])\n",color=:blue)
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
    printstyled(io, "  𝑃  details     = ",color=:normal)
    printstyled(io, "$(keys(𝑃))\n",color=:blue)
  return
end

function Base.similar(𝐼::Individuals)
    @unpack 🚄,📌,𝑃, 𝐷, 𝑀, 🔧, 🆔, 🔴, ∫ = 𝐼
    T = typeof(𝐼).parameters[1]
    return Individuals{T}(📌=similar(📌),🔴=similar(🔴),🆔=similar(🆔),
                          🚄=🚄, ∫=∫, 🔧=🔧, 𝑃=𝑃, 𝐷=𝐷, 𝑀=𝑀)
end

function Base.diff(𝐼::Individuals)
    f(x)=last(x).-first(x)
    🔴_by_ID = groupby(𝐼.🔴, :ID)
    return combine(🔴_by_ID,nrow,:lat => f => :dlat,:lon => f => :dlon)
end

