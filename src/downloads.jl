module datadeps

using DataDeps, Dataverse, Glob
import DataDeps: @datadep_str

"""
    unpackDV(filepath)

Like DataDeps's `:unpack` but using `Dataverse.untargz` and remove the `.tar.gz` file.
"""    
function unpackDV(filepath; do_warn=false)
    tmp_path=Dataverse.untargz(filepath)
    tmp_path2=joinpath(tmp_path,basename(filepath)[1:end-7])
    tmp_path=(ispath(tmp_path2) ? tmp_path2 : tmp_path)
    if isdir(tmp_path)
        [mv(p,joinpath(dirname(filepath),basename(p))) for p in glob("*",tmp_path)]
        [println(joinpath(dirname(filepath),basename(p))) for p in glob("*",tmp_path)]
        rm(filepath)
    else
        rm(filepath)
        mv(tmp_path,joinpath(dirname(filepath),basename(tmp_path)))
    end
    do_warn ? println("done with unpackDV for "*filepath) : nothing
end

"""
    __init__source_code()

Register data dependency with DataDep.
"""
function __init__datadeps()
    register(DataDep("global_ocean_circulation_inputs","inputs to global_ocean_circulation",
        ["https://zenodo.org/records/14029202/files/global_ocean_circulation_inputs.tar.gz"],
        ["6a53a695299811acad35b62fe1bd58cd3b0d94e361d557a1f268ad5e0c220e52"],      
        post_fetch_method=unpackDV))
    register(DataDep("global_ocean_circulation_outputs","outputs of global_ocean_circulation",
        ["https://zenodo.org/records/14029202/files/global_ocean_circulation_outputs.tar.gz"],
        ["2c308e84d598afba9bb9f6d20736dfb76d62d8e10befea6ac6fb1127da8ed609"],      
        post_fetch_method=unpackDV))
    register(DataDep("flt_example","files for MITgcm/pkg/flt example",
        ["https://github.com/gaelforget/flt_example/archive/v1.0.tar.gz"],
        ["3306d37f5e70d0fa253a5b87b527f727891b84179d8411971bbc10f3164d5fde"],      
        post_fetch_method=unpackDV))
    register(DataDep("Gulf_of_Mexico","Sample Drifter Dataset for example in JuliaConProceedings2023",
        ["https://zenodo.org/records/14229203/files/Drifters_example.jld2"],
        ["8c6ae8865821ac261316bfea3f31d1d362f210f806a85bfe5bb9775f3f91c2fa"],
        ))
    register(DataDep("Oscar_2021_small","Sample output from Drifters.jl applied to Oscar data set",
        ["https://zenodo.org/records/14876667/files/Drifters_Oscar_small.csv"],
        ["d0eb67039ec77cacf88313c0827d9a591f2675aa485dbc97bee72123e4bba30a"],
        ))
end

"""
    getdata(nam::String)

Add data to the scratch space folder. Known options for `nam` include 
"global_ocean_circulation_inputs", 
"""
function getdata(nam::String)
    withenv("DATADEPS_ALWAYS_ACCEPT"=>true) do
        if nam=="global_ocean_circulation_inputs"
            datadep"global_ocean_circulation_inputs"
        elseif nam=="global_ocean_circulation_outputs"
            datadep"global_ocean_circulation_outputs"
        elseif nam=="flt_example"
            datadep"flt_example"
        elseif nam=="Gulf_of_Mexico"
            datadep"Gulf_of_Mexico"
        elseif nam=="Oscar_2021_small"
            datadep"Oscar_2021_small"
        else
            println("unknown dataset")
        end
    end
end

end #module datadeps
