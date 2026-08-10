@testset "OCCA" begin
    OCCAmodule=Drifters.OCCA
	P,D=OCCAmodule.setup(nmax=5)
	nf=100; lo=(-160.0,-150.0); la=(30.0,40.0); level=2.5;
	df=Drifters.init.deprecated_initial_positions(D.Γ, nf, lo, la, level)
	I=Individuals(P,df.x,df.y,df.z,df.fid,(🔴=OCCAmodule.custom🔴,🔧=OCCAmodule.custom🔧, D=D))
	T=(0.0,10*86400.0)
	∫!(I,T)

    fig=CairoMakie.plot( DriftersDataset( data=(I=I,), options=(plot_type=:plot_start_end,) ) )
    @test isa(fig,CairoMakie.Figure)
end
