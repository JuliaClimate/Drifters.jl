@testset "time axes" begin
    for time_unit in (:DateTime,:seconds)
        Drifters.start_times(3,time_unit=time_unit,period=:day)
        @suppress start_times=Drifters.start_times(3, time_unit=time_unit, direction=:forward)
        @suppress backward_times=Drifters.start_times(3, time_unit=time_unit, direction=:backward)
    end
    @test true

    mon=30*86400.00
    TA=Drifters.TimeAxis(0.0,mon,0.0,mon,true)
    D0=Drifters.Dates.DateTime(2000,1,1); D1=D0+Drifters.Year(1)
    TA=Drifters.TimeAxis(D0,D1,D0,D1,true)
    @test isa(TA,Drifters.TimeAxis)

    Drifters.ECCO.set_times(:DateTime,:forward)
    TA,T,times=Drifters.ECCO.set_times(:seconds,:backward)
    @test isa(TA,Drifters.TimeAxis)
end
