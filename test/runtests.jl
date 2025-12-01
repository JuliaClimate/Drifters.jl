using Test, Documenter, Drifters, Suppressor, CairoMakie

import Drifters: MeshArrays, NetCDF, CSV, DataFrames, JLD2

import MITgcm; MITgcm.getdata("mitgcmsmall")

import Climatology
Climatology.get_ecco_velocity_if_needed()
Climatology.get_occa_velocity_if_needed()
Climatology.get_ecco_variable_if_needed("THETA")
Climatology.get_ecco_variable_if_needed("SALT")

MeshArrays.GridLoad(MeshArrays.GridSpec(ID=:LLC90))
MeshArrays.GridLoad(MeshArrays.GridSpec(ID=:onedegree))

@testset "Oscar" begin
    fil=joinpath(Drifters.datadeps.getdata("Oscar_2021_small"),"Drifters_Oscar_small.csv")
    df=CSV.read(fil,DataFrame)
    J=DriftersDataset( data=(df=df,), options=(plot_type=:Oscar_plot,))
    fig=plot(J)
    @test isa(fig,Figure)
    grid=Drifters.Oscar.grid()
    @test isa(grid,NamedTuple)
end

@testset "ECCO" begin
    ECCOmodule = Drifters.ECCO
    CSV = Drifters.CSV
    DataFrames = Drifters.DataFrames
    Individuals = Drifters.Individuals

    k=0
    P,D=ECCOmodule.init_FlowFields(k=k); np=100
    df0 = Drifters.init.init_global_randn(np , D)
    df = Drifters.init.init_gulf_stream(np , D)
    S = ECCOmodule.init_storage(np,100,length(D.Γ.RC),50)
    I = Individuals(P,df.x,df.y,df.z,df.fid,
        (D=merge(D,S),∫=ECCOmodule.custom∫,🔧=ECCOmodule.custom🔧,🔴=deepcopy(ECCOmodule.custom🔴)))
    my∫! = ECCOmodule.custom∫!
    T=(0.0,I.P.T[2])
    my∫!(I,T)
    @test isa(I,Individuals)

    tmp_🔴=I.🔴
    nt=length(unique(tmp_🔴.t))	
    xlims=(-85.0,5.0)
    ylims=(20.0,67.0)

    x=Drifters.DriftersDataset( data=(I=I,df=tmp_🔴,), options=(plot_type=:global_plot1,) )
    fig,tt=CairoMakie.plot(x)
    @test isa(fig,CairoMakie.Figure)
end

@testset "OCCA" begin
    OCCAmodule=Drifters.OCCA
	initial_positions=Drifters.init.initial_positions
	P,D=OCCAmodule.setup(nmax=5)
	nf=100; lo=(-160.0,-150.0); la=(30.0,40.0); level=2.5;
	df=initial_positions(D.Γ, nf, lo, la, level)
	I=Individuals(P,df.x,df.y,df.z,df.fid,(🔴=OCCAmodule.custom🔴,🔧=OCCAmodule.custom🔧, D=D))
	T=(0.0,10*86400.0)
	∫!(I,T)

    fig=CairoMakie.plot( DriftersDataset( data=(I=I,), options=(plot_type=:plot_start_end,) ) )
    @test isa(fig,CairoMakie.Figure)
end

@testset "simple" begin
    function SimpleFlowFields(nx,dx)
        XC = dx*(collect(1:2*nx) .- 0.5)
        YC = dx*(collect(1:nx) .- 0.5)        
        fac=0.1
        f(x, y) = sin(x) + cos(y) #streamfunction
        ϕ = fac*[f(x, y) for x in XC,y in YC] #streamfunction
        uC = -fac*[sin(y) for x in XC,y in YC] #dphi/dy at cell center
        vC = -fac*[cos(x) for x in XC,y in YC] #-dphi/dx at cell center
        return uC, vC, ϕ
    end
    
    nx=16; dx= π/nx; T=(0.,10.)
    uC, vC, ϕ = SimpleFlowFields(nx,dx)
    F=FlowFields(u=uC/dx,v=vC/dx,period=T)

    np,nq=size(F.u0)
    x=np*(0.4 .+ 0.2*rand(100))
    y=nq*(0.4 .+ 0.2*rand(100))
    I=Individuals(F,x,y)
    solve!(I,T)

    fig=CairoMakie.plot( DriftersDataset( data=(I=I,ϕ=ϕ), options=(plot_type=:simple_plot1,) ) )

    @test isa(fig,CairoMakie.Figure)
end

@testset "downloads" begin
    p0=Drifters.datadeps.getdata("global_ocean_circulation_inputs")
    Drifters.datadeps.getdata("flt_example")
    @test ispath(p0)
end

@testset "global" begin
    p0=Drifters.datadeps.getdata("global_ocean_circulation_inputs")
    ECCOmodule=Drifters.ECCO
    P,D=ECCOmodule.init_FlowFields()
    file_input=joinpath(p0,"initial_10_1.csv")
    df = Drifters.init.init_positions(10,filename=file_input)
    I=Individuals(P,df.x,df.y,df.f,(D=D,))
    T=(0.0,I.P.T[2])
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

@testset "various" begin
    u,v,w,pos=random_flow_field(format=:Array)
    F=FlowFields(u,u,v,v,[0,1.0])
    I=Individuals(F,pos...)
    ∫!(I)
    
    @suppress show(I)
    diff(I)
    size(I)
    J=similar(I)
    @test isa(J,Individuals)

    𝐺=convert_to_FlowFields(u,v,10.0)
    tmp2=nearest_to_xy(𝐺.u0,3.,3.,1.)
    @test isa(tmp2,Array)
    tmp3=nearest_to_xy(F.u0,3.,3.)
    @test isa(tmp3,Array)

    uC, vC, _ = random_flow_field(np=16)
    F=FlowFields(u=uC,v=vC,period=(0,10.))
    @test isa(F,uvArrays)

    df=DataFrame( ID=[], x=[], y=[], z=[], t = [])
    I=(position=zeros(3,2),ID=1:2,record=deepcopy(df))
    I=Individuals(I)
    @test isa(I,Individuals)

    GM=Drifters.Gulf_of_Mexico_setup()
    @test isa(GM.T,Tuple)
end

@testset "doctests" begin
    doctest(Drifters; manual = false)
end
