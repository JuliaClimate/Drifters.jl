module SDE_example

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

end
