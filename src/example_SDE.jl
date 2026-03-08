module ex_SDE

using Statistics

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
