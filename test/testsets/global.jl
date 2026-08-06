@testset "global" begin
    p0=Drifters.datadeps.getdata("global_ocean_circulation_inputs")
    ECCOmodule=Drifters.ECCO
    P,D=ECCOmodule.init_FlowFields(k=1)

    file_input=joinpath(p0,"initial_10_1.csv")
    df = Drifters.read_initial_positions(10,filename=file_input)
    I=Individuals(P,df.x,df.y,df.f,(D=D,))

    zer=(eltype(I.P.T)==Drifters.DateTime ? Drifters.DateTime(2000,1,1) : 0.0)
    T=(zer,I.P.T[2])
    D.🔄(P,D,T[1])

    ∫!(I,T)

    add_lonlat!(I.🔴,D.XC,D.YC)
    add_lonlat!(I.🔴,D.XC,D.YC,P.update_location!)
    tmp=interp_to_xy(I.🔴,D.YC)
    gcdist(I)

    @test prod(abs.(tmp).<90.0)

    tmp1=randn_lonlat(10)
    tmp2=stproj_inv(stproj(30.0,30.0)...)
    @test prod(isapprox.(tmp2,30.0,atol=1.0))
end
