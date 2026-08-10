@testset "ECCO" begin
    ECCOmodule = Drifters.ECCO
    Individuals = Drifters.Individuals

    k=0
    P,D=ECCOmodule.init_FlowFields(k=k)
    if true
        println("generating initial condition")
        np=100
        df0 = Drifters.init.initial_positions_2d(np , D.Γ)
        df = Drifters.init.init_regional_3d(np , D)
        fil=tempname()*".csv"; CSV.write(fil,df)
    else
        println("rereading initial condition")
        np=1000
        fil="ini_pos_ecco.csv"
        df=CSV.read(fil,DataFrame)
    end

    df0 = Drifters.init.initial_positions_2d(np , D.Γ)
    df = Drifters.init.init_regional_3d(np , D)
    S = ECCOmodule.init_storage(np,100,length(D.Γ.RC),50)
    I = Individuals(P,df.x,df.y,df.z,df.fid,
        (D=merge(D,S),∫=ECCOmodule.custom∫,🔧=ECCOmodule.custom🔧,🔴=deepcopy(ECCOmodule.custom🔴)))

    zer=(eltype(I.P.T)==Drifters.DateTime ? Drifters.DateTime(2000,1,1) : 0.0)
    T=(zer,I.P.T[2])
    D.🔄(P,D,T[1])
    ECCOmodule.custom∫!(I,T)
    @test isa(I,Individuals)

    tmp_🔴=I.🔴
    nt=length(unique(tmp_🔴.t))
    xlims=(-85.0,5.0)
    ylims=(20.0,67.0)

    x=Drifters.DriftersDataset( data=(I=I,df=tmp_🔴,), options=(plot_type=:global_plot1,) )
    fig,tt=CairoMakie.plot(x)
    @test isa(fig,CairoMakie.Figure)
end
