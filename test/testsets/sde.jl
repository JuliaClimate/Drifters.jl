@testset "SDE" begin
    import Drifters: ex_SDE
    SDE = Base.get_extension(Drifters, :DriftersStochasticDiffEqExt)
    MK = Base.get_extension(Drifters, :DriftersMakieExt)

    IC=ex_SDE.initial_conditions(100)
    SDE.main_loop(IC; ex_SDE.default_parameters()...)

    # Another run of the dispersion model, just for plotting trajectories
    tmp=SDE.demo_paths(IC; ex_SDE.default_parameters()...)
	fb=MK.plot_paths(za=tmp.za,zb=tmp.zb)

    # Eulerian model for comparison
    T,T0=ex_SDE.EulerianModel(10);
	f=MK.plot_EulerianModel(T,T0)

    # Compute population statistics
	st=ex_SDE.gridded_stats(IC)
	fs=MK.plot_stats(st,T=T)

    # Compute histogram
    z,sol=SDE.solve_paths(IC.u₀a)
    ρ,x_centers=ex_SDE.particle_density(sol)

    ex_SDE.g_erf(IC.u₀a,(0.5,0.1,0.001),0.0)
    ex_SDE.f_gauss(IC.u₀a,(0.5,0.1,0.001),0.0)

    @test isa(IC.u₀a,Vector)
end
