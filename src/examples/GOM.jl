
import JLD2

function Gulf_of_Mexico_setup()
    file=joinpath(datadeps.getdata("Gulf_of_Mexico"),"Drifters_example.jld2")
    jld=JLD2.jldopen(file)
    pol=jld["polygons"]
    obs=jld["drifters_real"]
    u=jld["u"]; v=jld["v"]
    x0=jld["x0"]; y0=jld["y0"]
    dT=jld["dT"]; nt=jld["nt"]
    T=(0.0,dT)
    (x0,y0,u=u,v=v,T=T,dT=dT,nt=nt,pol=pol,obs=obs)
end
