module DriftersClimatologyExt

using Climatology, Drifters
import Drifters: data_path

function data_path(name::Symbol)
    #use this environment variable to bypass downloads / API calls that require server access
    _SKIP_DOWNLOADS = parse(Bool,get(ENV, "SKIP_DOWNLOADS", "false"))
    if name==:OCCA
        _SKIP_DOWNLOADS ? nothing : Climatology.get_occa_velocity_if_needed()
        Climatology.ScratchSpaces.OCCA
    elseif name==:ECCO
        _SKIP_DOWNLOADS ? nothing : Climatology.get_ecco_velocity_if_needed()
        _SKIP_DOWNLOADS ? nothing : Climatology.get_ecco_variable_if_needed("THETA")
        _SKIP_DOWNLOADS ? nothing : Climatology.get_ecco_variable_if_needed("SALT")
        Climatology.ScratchSpaces.ECCO        
    else
        tempdir()
    end
end

end
