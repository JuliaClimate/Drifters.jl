@testset "downloads" begin
    p0=Drifters.datadeps.getdata("global_ocean_circulation_inputs")
    Drifters.datadeps.getdata("flt_example")
    @test ispath(p0)
end
