
module DriftersMITgcmExt

using MITgcm, Drifters
import Drifters: read_data_ECCO, read_data_mdsio, NetCDF

function read_data_ECCO(m0,v0,path,grid,k)
    MITgcm.read_nctiles(joinpath(path,v0),v0,grid,I=(:,:,k,m0))
end
read_data_ECCO(m0,v0,P,D)=read_data_ECCO(m0,v0, joinpath(D.pth,v0) , P.u0.grid , D.k )

read_data_mdsio(path, varname) = MITgcm.read_mdsio(path, varname)

end
