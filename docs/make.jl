using Documenter, Literate, PlutoSliderServer, Drifters

import MITgcm; MITgcm.getdata("mitgcmsmall")

using Climatology

#download data dependencies if needed
Drifters.datadeps.getdata("flt_example")
Climatology.get_ecco_velocity_if_needed();
Climatology.get_occa_velocity_if_needed();
Climatology.get_ecco_variable_if_needed("THETA")
Climatology.get_ecco_variable_if_needed("SALT")

# generate tutorials and how-to guides using Literate
src = joinpath(@__DIR__, "src/")
lit = joinpath(@__DIR__, "../examples/more/")
notebooks = joinpath(src, "notebooks")

execute = false # Set to true for executing notebooks and documenter!
nb = true      # Set to true to generate the notebooks

lst1 = ["detailed_look","particle_cloud"]
lst2 = ["detailed_look","particle_cloud"]
tst1(x) = !isempty(lst1) && Bool(sum(isequal.(x, lst1)))
tst2(x) = !isempty(lst2) && Bool(sum(isequal.(x, lst2)))

for (root, _, files) in walkdir(lit), file in files
    splitext(file)[2] == ".jl" || continue
	tst1(splitext(file)[1]) || continue
    ipath = joinpath(root, file)
    opath = splitdir(replace(ipath, lit=>src))[1]
    Literate.markdown(ipath, opath, documenter = execute)
    nb && Literate.notebook(ipath, notebooks, execute = execute*tst2(splitext(file)[1]))
end

# Documentation structure
ismd(f) = splitext(f)[2] == ".md"
pages(folder) = [joinpath(folder, f) for f in readdir(joinpath(src, folder)) if ismd(f)]

makedocs(;
    sitename = "Drifters.jl",
    repo = Remotes.GitHub("JuliaClimate", "Drifters.jl"),
    authors="JuliaClimate <gforget@mit.edu>",
    format = Documenter.HTML(),
    pages = [
	"Introduction" => "index.md",
        "User Guide" => "workflow.md",
        "Examples" => "examples.md", 
        "Tool Box" => "API.md"],
    doctest = false,
    warnonly = [:cross_references,:missing_docs],
    modules = [Drifters]
)

pth_in = joinpath(@__DIR__, "..","examples")
pth_out = joinpath(@__DIR__, "build","examples")
lst=("global_ocean_circulation.jl","three_dimensional_ocean.jl",
     "solid_body_rotation.jl","random_flow_field.jl",
     "interactive_UI.jl","Oscar_model.jl","StochasticDiffEq.jl")
subpth=("worldwide","worldwide","basics","basics","worldwide","worldwide","basics")
for ii in 1:length(lst)
    fil_in=joinpath(pth_in,subpth[ii],lst[ii])
    fil_out=joinpath(pth_out,lst[ii][1:end-2]*"html")
    PlutoSliderServer.export_notebook(fil_in)
    mv(fil_in[1:end-2]*"html",fil_out)
    #cp(fil_in[1:end-2]*"html",fil_out)
end

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/JuliaClimate/Drifters.jl.git",
)
