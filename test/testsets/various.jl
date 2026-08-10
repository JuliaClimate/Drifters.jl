@testset "various" begin
    u,v,w,pos=random_flow_field(format=:Array)
    F=FlowFields(u,u,v,v,[0,1.0])
    I=Individuals(F,pos...,(D=(problem_type=:default,),))
    ∫!(I)

    @suppress show(I)
    diff(I)
    size(I)
    J=similar(I)
    @test isa(J,Individuals)

    G=convert_to_FlowFields(u,v,10.0)
    tmp2=nearest_to_xy(G.u0,3.,3.,1.)
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
