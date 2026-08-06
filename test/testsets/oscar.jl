@testset "Oscar" begin
    fil=joinpath(Drifters.datadeps.getdata("Oscar_2021_small"),"Drifters_Oscar_small.csv")
    df=CSV.read(fil,DataFrame)
    J=DriftersDataset( data=(df=df,), options=(plot_type=:Oscar_plot,))
    fig=plot(J)
    @test isa(fig,Figure)
    grid=Drifters.Oscar.grid()
    @test isa(grid,NamedTuple)
end
