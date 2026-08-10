@testset "simple" begin
    function SimpleFlowFields(nx,dx)
        XC = dx*(collect(1:2*nx) .- 0.5)
        YC = dx*(collect(1:nx) .- 0.5)
        fac=0.1
        f(x, y) = sin(x) + cos(y)
        ϕ = fac*[f(x, y) for x in XC,y in YC]
        uC = -fac*[sin(y) for x in XC,y in YC]
        vC = -fac*[cos(x) for x in XC,y in YC]
        return uC, vC, ϕ
    end

    nx=16; dx= π/nx;
    uC, vC, ϕ = SimpleFlowFields(nx,dx)

    D0=Drifters.DateTime(2000,1,1)
    D1=Drifters.DateTime(2000,1,1,0,0,10)
    T=(D0,D1)
    F=FlowFields(u=uC/dx,v=vC/dx,period=T)

    np,nq=size(F.u0)
    x=np*(0.4 .+ 0.2*rand(100))
    y=nq*(0.4 .+ 0.2*rand(100))
    I=Individuals(F,x,y)
    solve!(I,T)

    fig=CairoMakie.plot( DriftersDataset( data=(I=I,ϕ=ϕ), options=(plot_type=:simple_plot1,) ) )
    @test isa(fig,CairoMakie.Figure)
end
