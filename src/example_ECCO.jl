module ECCO

import MeshArrays, DataDeps, CSV, JLD2

import Drifters: postprocess_MeshArray, add_lonlat!, OrdinaryDiffEq
import OrdinaryDiffEq: solve, Tsit5, ODEProblem
import Drifters: update_location!, Individuals, uvMeshArrays, uvwMeshArrays
import Drifters: FlowFields, ensemble_solver, data_path, read_data_ECCO, order
import Drifters.DataFrames: DataFrame
import Drifters.MeshArrays as MeshArrays
import Drifters.MeshArrays: gcmgrid, MeshArray, exchange

export init_FlowFields, init_storage
export custom∫, custom🔧, custom🔴, custom∫!
#export reset_📌!, init_z_if_needed

"""
    reset_📌!(I::Individuals,frac::Number,📌::Array)

Randomly select a fraction (frac) of the particles and reset 
their positions (I.📌) to a random subset of the specified 📌.
"""
function reset_📌!(I::Individuals,frac::Number,📌::Array)
    np=length(I.🆔)
    n_reset = Int(round(frac*np))
    k_reset = rand(1:np, n_reset)
    l_reset = rand(1:np, n_reset)
    I.📌[k_reset]=deepcopy(📌[l_reset])
    #isempty(I.🔴.ID) ? m=maximum(I.🆔) : m=max(maximum(I.🔴.ID),maximum(I.🆔))
    #I.🆔[k_reset]=collect(1:n_reset) .+ m
end

"""
    init_storage(np)

- np=number of individuals to store
- nn=number of individuals per chunk
- nl=number of vertical levels to store
- nt=number of time steps to store
"""
function init_storage(np,nn,nl,nt=2)
    (   prof_T=NaN*zeros(np,nl,nt),batch_T=zeros(2*nn,nl),local_T=zeros(2*nn),
        prof_S=NaN*zeros(np,nl,nt),batch_S=zeros(2*nn,nl),local_S=zeros(2*nn)
        ) 
end

"""
    setup_FlowFields(k::Int,Γ::NamedTuple,func::Function,pth::String)

Define `FlowFields` data structure along with specified grid (`Γ` NamedTuple), 
function `func` (e.g., `(u -> MeshArrays.update_location_llc!(u,Γ)))`, 
and file location (`pth`).
    
_Note: the initial implementation approximates month durations to 
365 days / 12 months for simplicity and sets P.T to [-mon/2,mon/2]_
"""
function setup_FlowFields(k::Int,Γ::NamedTuple,func::Function,pth::String,backward_time=false)
    XC=exchange(Γ.XC) #add 1 lon point at each edge
    YC=exchange(Γ.YC) #add 1 lat point at each edge
    iDXC=1. ./Γ.DXC
    iDYC=1. ./Γ.DYC
    γ=Γ.XC.grid
    mon=86400.0*365.0/12.0
    
    if k==0
        msk=Γ.hFacC
        msk=1.0*(msk .> 0.0)
        (_,nr)=size(msk)
        exmsk=similar(msk)
        for k=1:nr
            exmsk[:,k]=exchange(msk[:,k]).MA
        end
        P=FlowFields(MeshArray(γ,Float32,nr),MeshArray(γ,Float32,nr),
        MeshArray(γ,Float32,nr),MeshArray(γ,Float32,nr),
        MeshArray(γ,Float32,nr+1),MeshArray(γ,Float32,nr+1),
        [-mon/2,mon/2],func)
    else
        msk=Γ.hFacC[:, k]
        msk=1.0*(msk .> 0.0)
        exmsk=exchange(msk).MA
        P=FlowFields(MeshArray(γ,Float32),MeshArray(γ,Float32),
        MeshArray(γ,Float32),MeshArray(γ,Float32),[-mon/2,mon/2],func)    
    end
    
    D = (🔄 = update_FlowFields!, pth=pth,
         XC=XC, YC=YC, iDXC=iDXC, iDYC=iDYC,
         k=k, msk=msk, exmsk=exmsk, 
         θ0=similar(msk), θ1=similar(msk),
         S0=similar(msk), S1=similar(msk))

    #add parameters related to gridded domain decomposition
    D = merge(D , MeshArrays.NeighborTileIndices_cs(Γ))

    tmp=(Γ=Γ, backward_time=backward_time)
    D=merge(D,tmp)

    return P,D
end

"""
    update_FlowFields!(P::uvMeshArrays,D::NamedTuple,t::Float64)

Update flow field arrays (in P), P.T, and ancillary variables (in D) 
according to the chosen time `t` (in `seconds`). 

_Note: for now, it is assumed that (1) the time interval `dt` between 
consecutive records is diff(P.T), (2) monthly climatologies are used 
with a periodicity of 12 months, (3) vertical P.k is selected_
"""
function update_FlowFields!(P::uvMeshArrays,D::NamedTuple,t::AbstractFloat)
    dt=P.T[2]-P.T[1]

    m0=Int(floor((t+dt/2.0)/dt))
    m1=m0+1
    t0=m0*dt-dt/2.0
    t1=m1*dt-dt/2.0

    m0=mod(m0,12)
    m0==0 ? m0=12 : nothing
    m1=mod(m1,12)
    m1==0 ? m1=12 : nothing

    velocity_factor=1.0
    if D.backward_time
        velocity_factor=-1.0
        m0=13-m0
        m1=13-m1
    end

    (U,V)=read_velocities(P.u0.grid,m0,D.pth)
    u0=velocity_factor*U[:,D.k]; v0=velocity_factor*V[:,D.k]
    u0[findall(isnan.(u0))]=0.0; v0[findall(isnan.(v0))]=0.0 #mask with 0s rather than NaNs
    u0=u0.*D.iDXC; v0=v0.*D.iDYC; #normalize to grid units
    (u0,v0)=MeshArrays.exchange_main(u0,v0,1) #add 1 point at each edge for u and v

    (U,V)=read_velocities(P.u0.grid,m1,D.pth)
    u1=velocity_factor*U[:,D.k]; v1=velocity_factor*V[:,D.k]
    u1[findall(isnan.(u1))]=0.0; v1[findall(isnan.(v1))]=0.0 #mask with 0s rather than NaNs
    u1=u1.*D.iDXC; v1=v1.*D.iDYC; #normalize to grid units
    (u1,v1)=MeshArrays.exchange_main(u1,v1,1) #add 1 point at each edge for u and v

    P.u0[:]=Float32.(u0.MA[:])
    P.u1[:]=Float32.(u1.MA[:])
    P.v0[:]=Float32.(v0.MA[:])
    P.v1[:]=Float32.(v1.MA[:])

    θ0=read_data_ECCO(m0,"THETA",P,D)
    θ0[findall(isnan.(θ0))]=0.0 #mask with 0s rather than NaNs
    D.θ0[:]=Float32.(θ0[:,1])

    θ1=read_data_ECCO(m1,"THETA",P,D)
    θ1[findall(isnan.(θ1))]=0.0 #mask with 0s rather than NaNs
    D.θ1[:]=Float32.(θ1[:,1])

    S0=read_data_ECCO(m0,"SALT",P,D)
    S0[findall(isnan.(S0))]=0.0 #mask with 0s rather than NaNs
    D.S0[:]=Float32.(S0[:,1])

    S1=read_data_ECCO(m1,"SALT",P,D)
    S1[findall(isnan.(S1))]=0.0 #mask with 0s rather than NaNs
    D.S1[:]=Float32.(S1[:,1])

    D.θ0[:]=exchange(D.θ0).MA
    D.θ1[:]=exchange(D.θ1).MA
    D.S0[:]=exchange(D.S0).MA
    D.S1[:]=exchange(D.S1).MA

    P.T[:]=[t0,t1]

end

"""
    update_FlowFields!(P::uvwMeshArrays,D::NamedTuple,t::Float64)

Update flow field arrays (in P), P.T, and ancillary variables (in D) 
according to the chosen time `t` (in `seconds`). 

_Note: for now, it is assumed that (1) the time interval `dt` between 
consecutive records is diff(P.T), (2) monthly climatologies are used 
with a periodicity of 12 months, (3) vertical P.k is selected_
"""
function update_FlowFields!(P::uvwMeshArrays,D::NamedTuple,t::AbstractFloat)
    dt=P.T[2]-P.T[1]

    m0=Int(floor((t+dt/2.0)/dt))
    m1=m0+1
    t0=m0*dt-dt/2.0
    t1=m1*dt-dt/2.0

    m0=mod(m0,12)
    m0==0 ? m0=12 : nothing
    m1=mod(m1,12)
    m1==0 ? m1=12 : nothing

    velocity_factor=1.0
    if D.backward_time
        velocity_factor=-1.0
        m0=13-m0
        m1=13-m1
    end

    (_,nr)=size(D.Γ.hFacC)

    (U,V)=read_velocities(P.u0.grid,m0,D.pth)
    u0=velocity_factor*U; v0=velocity_factor*V
    u0[findall(isnan.(u0))]=0.0; v0[findall(isnan.(v0))]=0.0 #mask with 0s rather than NaNs
    for k=1:nr
        u0[:,k]=u0[:,k].*D.iDXC; v0[:,k]=v0[:,k].*D.iDYC; #normalize to grid units
        (tmpu,tmpv)=MeshArrays.exchange_main(u0[:,k],v0[:,k],1) #add 1 point at each edge for u and v
        u0[:,k]=tmpu.MA
        v0[:,k]=tmpv.MA
    end
    w0=velocity_factor*read_data_ECCO(m0,"WVELMASS",joinpath(D.pth,"WVELMASS"),P.u0.grid,:)
    w0[findall(isnan.(w0))]=0.0 #mask with 0s rather than NaNs

    (U,V)=read_velocities(P.u0.grid,m1,D.pth)
    u1=velocity_factor*U; v1=velocity_factor*V
    u1[findall(isnan.(u1))]=0.0; v1[findall(isnan.(v1))]=0.0 #mask with 0s rather than NaNs
    for k=1:nr
        u1[:,k]=u1[:,k].*D.iDXC; v1[:,k]=v1[:,k].*D.iDYC; #normalize to grid units
        (tmpu,tmpv)=MeshArrays.exchange_main(u1[:,k],v1[:,k],1) #add 1 point at each edge for u and v
        u1[:,k]=tmpu.MA
        v1[:,k]=tmpv.MA
    end
    w1=velocity_factor*read_data_ECCO(m1,"WVELMASS",joinpath(D.pth,"WVELMASS"),P.u0.grid,:)
    w1[findall(isnan.(w1))]=0.0 #mask with 0s rather than NaNs

    P.u0[:,:]=Float32.(u0[:,:])
    P.u1[:,:]=Float32.(u1[:,:])
    P.v0[:,:]=Float32.(v0[:,:])
    P.v1[:,:]=Float32.(v1[:,:])
    for k=1:nr
        tmpw=MeshArrays.exchange_main(-w0[:,k]).MA
        P.w0[:,k]=Float32.(tmpw./D.Γ.DRC[k])
        tmpw=MeshArrays.exchange_main(-w1[:,k]).MA
        P.w1[:,k]=Float32.(tmpw./D.Γ.DRC[k])
    end
    P.w0[:,1]=0*MeshArrays.exchange_main(-w0[:,1]).MA
    P.w1[:,1]=0*MeshArrays.exchange_main(-w1[:,1]).MA
    P.w0[:,nr+1]=0*MeshArrays.exchange_main(-w0[:,1]).MA
    P.w1[:,nr+1]=0*MeshArrays.exchange_main(-w1[:,1]).MA

    θ0=read_data_ECCO(m0,"THETA",joinpath(D.pth,"THETA"),P.u0.grid,:)
    θ0[findall(isnan.(θ0))]=0.0 #mask with 0s rather than NaNs
    D.θ0[:,:]=Float32.(θ0[:,:])

    θ1=read_data_ECCO(m1,"THETA",joinpath(D.pth,"THETA"),P.u0.grid,:)
    θ1[findall(isnan.(θ1))]=0.0 #mask with 0s rather than NaNs
    D.θ1[:,:]=Float32.(θ1[:,:])

    S0=read_data_ECCO(m0,"SALT",joinpath(D.pth,"SALT"),P.u0.grid,:)
    S0[findall(isnan.(S0))]=0.0 #mask with 0s rather than NaNs
    D.S0[:,:]=Float32.(S0[:,:])

    S1=read_data_ECCO(m1,"SALT",joinpath(D.pth,"SALT"),P.u0.grid,:)
    S1[findall(isnan.(S1))]=0.0 #mask with 0s rather than NaNs
    D.S1[:,:]=Float32.(S1[:,:])

    for k=1:nr
        D.θ0[:,k]=exchange(D.θ0[:,k]).MA
        D.θ1[:,k]=exchange(D.θ1[:,k]).MA
        D.S0[:,k]=exchange(D.S0[:,k]).MA
        D.S1[:,k]=exchange(D.S1[:,k]).MA
    end

    P.T[:]=[t0,t1]
end

"""
    update_FlowFields!(I::Individuals)

Update flow field arrays in I.P, update I.P.T, and ancillary variables in I.D
such that I.P.T includes time `t_ϵ=I.P.T[2]+eps(I.P.T[2])`.

_Note: for now, it is assumed that (1) the time interval `dt` between 
consecutive records is diff(P.T), (2) monthly climatologies are used 
with a periodicity of 12 months, (3) vertical P.k is selected_
"""
function update_FlowFields!(I::Individuals)
    t_ϵ=I.P.T[2]+eps(I.P.T[2])
    I.D.🔄(I.P,I.D,t_ϵ)
end

"""
    read_velocities(γ::gcmgrid,t::Int,pth::String)

Read velocity components `u,v` from files in `pth`for time `t`
"""
function read_velocities(γ::gcmgrid,t::Int,pth::String)
    u=read_data_ECCO(t,"UVELMASS",joinpath(pth,"UVELMASS"),γ,:)
    v=read_data_ECCO(t,"VVELMASS",joinpath(pth,"VVELMASS"),γ,:)
    return u,v
end

"""
    init_FlowFields(;k=1)

Set up Global Ocean particle simulation in 2D with seasonally varying flow field.
"""
function init_FlowFields(; k=1, backward_time=false)
  
  #read grid and set up connections between subdomains
  γ=MeshArrays.GridSpec(ID=:LLC90)
  Γ=MeshArrays.GridLoad(γ,option="full")
  f(x,y)=Float32.(MeshArrays.GridLoadVar(x,y))
  tmp=( DXC=f("DXC",γ),DYC=f("DYC",γ),hFacC=f("hFacC",γ),
        Depth=f("Depth",γ),RC=f("RC",γ),DRC=f("DRC",γ))
  Γ=merge(Γ,tmp)
  Γ=merge(Γ,MeshArrays.NeighborTileIndices_cs(Γ))
  func=(u -> MeshArrays.update_location_llc!(u,Γ))

  #initialize u0,u1 etc
  P,D=setup_FlowFields(k,Γ,func,data_path(:ECCO),backward_time)
  D.🔄(P,D,0.0)

  #add background map for plotting
  λ=get_interp_coefficients(Γ)
  ODL=OceanDepthLog(λ,Γ)
  
  #(optional) fraction of the particles reset per month (e.g., 0.05 for k<=10)
  r_reset = 1/12/4

  #add parameters for use in reset!
  tmp=(frac=r_reset, ODL=ODL)
  D=merge(D,tmp)

  return P,D
end

function get_interp_coefficients(Γ)
    fil=joinpath(MeshArrays.mydatadep("interp_halfdeg"),"interp_coeffs_halfdeg.jld2")
    λ=JLD2.load(fil)
    λ=MeshArrays.Dict_to_NamedTuple(λ)
end

function OceanDepthLog(λ,Γ)
    DL=MeshArrays.Interpolate(λ.μ*Γ.Depth,λ.f,λ.i,λ.j,λ.w)
    DL=reshape(DL,size(λ.lon))
    DL[findall(DL.<0)].=0
    DL=transpose(log10.(DL))
    DL[findall((!isfinite).(DL))].=NaN
    (lon=λ.lon[:,1],lat=λ.lat[1,:],fld=DL,rng=(1.5,5))
end

custom∫(prob) = ensemble_solver(prob,solver=Tsit5(),reltol=1e-5,abstol=1e-5)

custom🔴 = DataFrame(ID=Int[], fid=Int[], x=Float64[], y=Float64[], 
lon=Float64[], lat=Float64[], z=Float64[], d=Float64[], 
θ=Float64[], SSθ=Float64[], S=Float64[], SSS=Float64[], year=Float64[], t=Float64[])

function custom🔧(sol,𝐹::uvwMeshArrays,D::NamedTuple;id=missing,T=missing)

    df=postprocess_MeshArray(sol,𝐹,D,id=id,T=T)
    np=length(sol.u)
    z=[[sol.u[i][:,1][3] for i in 1:np];[sol.u[i][:,end][3] for i in 1:np]]
    df.z=z[:]
    df.year=df.t ./86400/365
    add_lonlat!(df,D.XC,D.YC)
    k=Int.(floor.(z)); w=(z-k);
	df.d=D.Γ.RF[1 .+ k].*(1 .- w)+D.Γ.RF[2 .+ k].*w

    #for k in 1:nr
    # D.batch_T[:,k]=interp_to_xy(df,D.θ1[:,k])./interp_to_xy(df,D.exmsk[:,k])
    #end

    x=df[!,:x];
    y=df[!,:y];
    f=Int.(df[!,:fid]);
    dx,dy=(x - floor.(x) .+ 0.5,y - floor.(y) .+ 0.5);
    i_c = Int32.(floor.(x)) .+ 1;
    j_c = Int32.(floor.(y)) .+ 1;
    
    nr=size(D.exmsk,2)

    #need time interpolation (df.t)
    for k in 1:nr, jj in 1:length(i_c)
        tmp0=(1.0-dx[jj])*(1.0-dy[jj])*D.exmsk[f[jj],k][i_c[jj],j_c[jj]]+
        (dx[jj])*(1.0-dy[jj])*D.exmsk[f[jj],k][i_c[jj]+1,j_c[jj]]+
        (1.0-dx[jj])*(dy[jj])*D.exmsk[f[jj],k][i_c[jj],j_c[jj]+1]+
        (dx[jj])*(dy[jj])*D.exmsk[f[jj],k][i_c[jj]+1,j_c[jj]+1]
        #
        tmp1=(1.0-dx[jj])*(1.0-dy[jj])*D.θ1[f[jj],k][i_c[jj],j_c[jj]]+
        (dx[jj])*(1.0-dy[jj])*D.θ1[f[jj],k][i_c[jj]+1,j_c[jj]]+
        (1.0-dx[jj])*(dy[jj])*D.θ1[f[jj],k][i_c[jj],j_c[jj]+1]+
        (dx[jj])*(dy[jj])*D.θ1[f[jj],k][i_c[jj]+1,j_c[jj]+1]
        D.batch_T[jj,k]=tmp1/tmp0
        #
        tmp1=(1.0-dx[jj])*(1.0-dy[jj])*D.S1[f[jj],k][i_c[jj],j_c[jj]]+
        (dx[jj])*(1.0-dy[jj])*D.S1[f[jj],k][i_c[jj]+1,j_c[jj]]+
        (1.0-dx[jj])*(dy[jj])*D.S1[f[jj],k][i_c[jj],j_c[jj]+1]+
        (dx[jj])*(dy[jj])*D.S1[f[jj],k][i_c[jj]+1,j_c[jj]+1]
        D.batch_S[jj,k]=tmp1/tmp0
    end

    #need time interpolation (df.t)
    for p=1:size(df,1)
        #k=max(Int(floor(z[p])),1)
        #local_T[p]=batch_T[p,k]
        k1=floor(z[p]+0.5)
        a2=(z[p]+0.5)-k1
        k2=Int(min(max(k1+1,1),nr))
        k1=Int(min(max(k1,1),nr))
        D.local_T[p]=(1-a2)*D.batch_T[p,k1]+a2*D.batch_T[p,k2]
        D.local_S[p]=(1-a2)*D.batch_S[p,k1]+a2*D.batch_S[p,k2]
    end

    df.SSθ=D.batch_T[:,1]
    df.θ=D.local_T[:]
    df.SSS=D.batch_S[:,1]
    df.S=D.local_S[:]

    if D.frac==0.0
        t=Int(round(0.5+df.t[end]/(T[2]-T[1])))
        nn=Int(size(D.batch_T,1)/2)
        #df.ID[1]==1 ? println(t) : nothing        
        D.prof_T[df.ID[nn+1:end],:,t].=D.batch_T[nn+1:end,:]
        D.prof_S[df.ID[nn+1:end],:,t].=D.batch_S[nn+1:end,:]
    end

    return df
end

function custom∫!(I::Individuals,T)
    (; 🚄,📌,P,D,🔧,🆔,🔴,∫) = I

    vel=0*vec(📌)
    [🚄(vel[i],📌[i],P,T[1]) for i in 1:length(vel)]
    nd=ndims(vel[1])-1
    vel=[sqrt(sum(vel[ii][1:nd].^2)) for ii in eachindex(vel)]
    vel=[(ii,vel[ii]) for ii=1:length(vel)]
    sort!(vel, by = x -> x[2])
    ii=[vel[ii][1] for ii=1:length(vel)]
    nn=Int(size(D.batch_T,1)/2)
    ni=Int(ceil(length(ii)/nn))

    nt=6
    dt=(I.P.T[2]-I.P.T[1])/nt
    nj=Int(round(ni*min(T[2]/86400/365,1)))
    
    tmp=deepcopy(custom🔴)
    for i=1:ni
        mm=min(nn,length(ii)-nn*(i-1))
#        println("i="*string(i))
        jj=ii[nn*(i-1) .+ collect(1:mm)]
      for tt in 1:nt
        TT=[I.P.T[1]+(tt-1)*dt, I.P.T[1]+tt*dt]
        prob = ODEProblem(🚄,permutedims(📌[jj]), TT ,P)
        sol = ∫(prob)
        append!(tmp, 🔧(sol,P,D,id=🆔[jj], T=TT))
        📌[jj] = deepcopy([sol[i].u[end] for i in 1:mm])
        if i<=nj
         if isa(P,uvwMeshArrays)||isa(P,uvMeshArrays)
             [update_location!(i,P) for i in 📌[jj]]
         end
        end
      end
    end
    sort!(tmp, order(:t))

    isempty(🔴) ? np =0 : np=length(🆔)
    append!(🔴,tmp[np+1:end,:])
#    isempty(🔴) ? append!(🔴,tmp) : 🔴[:,:]=tmp[:,:]
    return true
end

custom∫!(x) = custom∫!(x,(x.P.T[1],x.P.T[2]))

end #module ECCO_FlowFields
