module DriftersStochasticDiffEqExt

import StochasticDiffEq, Drifters
import StochasticDiffEq: EnsembleProblem, SDEProblem, EM, solve, stack, remake
import Drifters: _SDEProblem, default_solver, ensemble_solver, ex_SDE

notes()=println("""
- For now : `step!(u₀)` is the sequence of `solve_paths` and `fold_tails`.
- For deeper integration into Drifters (WIP), we need `_SDEProblem`, `default_solver`, `ensemble_solver`.

```
function ∫!(I::Individuals,T::Tuple)
prob = ODEProblem(🚄,📌, T ,P)

function ensemble_solver(prob;solver=Tsit5(),reltol=1e-8,abstol=1e-8,safetycopy=false)
indiv_prob = ODEProblem(prob.f,u0[1],prob.tspan,prob.p)
ensemble_prob = EnsembleProblem(indiv_prob,prob_func=prob_func,safetycopy=safetycopy)
```
""")

"""
    main_loop(IC;p=0.5,nt=5)

```
IC is a NamedTuple providing u₀a,u₀b,ca,cb,np
p=0.5 # fraction of mass exchanged with neighbors every time step
nt=5 # number of time steps
```
"""
function main_loop(IC;p=0.5,nt=5,do_mix_neighbors=true)
    (; u₀a,u₀b,ca,cb,np) = IC
    for tt in 1:nt
        step!(u₀a)
        step!(u₀b)
        do_mix_neighbors ? ex_SDE.mix_neighbors(u₀a,ca,u₀b,cb,p) : false
    end
    return "done with model run"
end

"""
    step!(u₀)

Sequence of `solve_paths` and `fold_tails`.

```
using Drifters, StochasticDiffEq
SDE = Base.get_extension(Drifters, :DriftersStochasticDiffEqExt)
?SDE.step!
```
"""
function step!(u₀; do_fold_tails=true)
    za,sol_a=solve_paths(u₀)
    do_fold_tails ? ex_SDE.fold_tails(za) : nothing
    u₀[:]=za[:,end]
end

"""
From the SciML docs: Ito interpretation is the default.

```
# Itô interpretation (default)
prob = SDEProblem(f!, g!, u0, tspan, interpretation = :Ito)

# Stratonovich interpretation  
prob = SDEProblem(f!, g!, u0, tspan, interpretation = :Stratonovich)

# Or at solver level
sol = solve(prob, RKMil(interpretation = :Stratonovich))
```
"""

function solve_paths(u₀; P=ex_SDE.default_parameters())
    f(u,p,t) = 0.0
    g(u,p,t) = 0.1
    params=(P.mldepth,P.thickness,P.mlkappa,P.seafloor,P.depthscale)
    config=P.configuration
    if config==:basic
        prob = SDEProblem(f,g,u₀,P.tspan)
        sol=solve(prob,EM(),dt=P.dt)
    elseif config==:piecewise
        prob = SDEProblem(ex_SDE.f_piecewise,ex_SDE.g_piecewise,u₀,P.tspan,params)
        sol=solve(prob,EM(),dt=P.dt,callback = ex_SDE.surface_reflect(P.seafloor))
	elseif config==:ggl90
		rf = [0.0, -10.0, -20.0, -30.0, -40.0, -50.0, -60.0, -70.0, -80.01, -90.04, -100.15, 
			-110.47, -121.27, -133.03, -146.45, -162.48999999999998, -182.30999999999997, 
			-207.15999999999997, -238.25999999999996, -276.67999999999995, -323.17999999999995, 
			-378.17999999999995, -441.67999999999995, -513.26, -592.16, -677.31, -767.49, 
			-861.45, -958.0300000000001, -1056.2800000000002, -1155.5300000000002, -1255.5400000000002, 
			-1356.8700000000001, -1461.43, -1572.76, -1695.59, -1834.6799999999998, -1993.62, -2174.45, 
			-2378.0, -2604.5, -2854.0, -3126.5, -3422.0, -3740.5, -4082.0, -4446.5, 
			-4834.0, -5244.5, -5678.0, -6134.5]
		kappa = ex_SDE.DEFAULT_KAPPA
		kappa_itp = ex_SDE.build_interpolator(Γ.RF,kappa)
		f_ggl,g_ggl = ex_SDE.build_kappa_profile(kappa_itp)
		prob = SDEProblem(f_ggl,g_ggl,u₀,P.tspan,params)
        sol=solve(prob,EM(),dt=P.dt,callback = ex_SDE.ocean_reflect(P.seafloor))
	elseif config==:dimensional_idealized
		rf = [0.0, -10.0, -20.0, -30.0, -40.0, -50.0, -60.0, -70.0, -80.01, -90.04, -100.15, 
			-110.47, -121.27, -133.03, -146.45, -162.48999999999998, -182.30999999999997, 
			-207.15999999999997, -238.25999999999996, -276.67999999999995, -323.17999999999995, 
			-378.17999999999995, -441.67999999999995, -513.26, -592.16, -677.31, -767.49, 
			-861.45, -958.0300000000001, -1056.2800000000002, -1155.5300000000002, -1255.5400000000002, 
			-1356.8700000000001, -1461.43, -1572.76, -1695.59, -1834.6799999999998, -1993.62, -2174.45, 
			-2378.0, -2604.5, -2854.0, -3126.5, -3422.0, -3740.5, -4082.0, -4446.5, 
			-4834.0, -5244.5, -5678.0, -6134.5]
		kappa = zeros(50)
		iz_bottom = 10
		kappa[1:iz_bottom] .= 0.02
		println("Mixed layer depth:", -Γ.RF[iz_bottom+1], "m")
		println("Transition thickness:", Γ.DRF[iz_bottom], "m")
		kappa_itp = ex_SDE.build_interpolator(Γ.RF,kappa)
		f_di,g_di = ex_SDE.build_kappa_profile(kappa_itp)
		prob = SDEProblem(f_di,g_di,u₀,P.tspan,params)
        sol=solve(prob,EM(),dt=P.dt,callback = ex_SDE.ocean_reflect(P.seafloor))
    else
        error("unknown config")
    end
    stack(sol(0:0.01:1)),sol
end

function demo_paths(IC::NamedTuple; do_fold_tails=true)
    (; u₀a,u₀b,ca,cb,np) = IC
    za,sol_a=solve_paths(u₀a)
    do_fold_tails ? ex_SDE.fold_tails(za) : nothing
    zb,sol_b=solve_paths(u₀b)
    do_fold_tails ? ex_SDE.fold_tails(zb) : nothing
    (za=za,zb=zb)
end

"""
    _SDEProblem(f::Function,g::Function,u₀,tspan)
    
Calls `StochasticDiffEq.SDEProblem(f,g,u₀,tspan,p)``

```
using Drifters, StochasticDiffEq
?Drifters._SDEProblem
```
"""
_SDEProblem(f::Function,g::Function,u₀,tspan,p)=SDEProblem(f,g,u₀,tspan,p)

"""
Assuming that `∫!` calls `_SDEProblem` instead of `ODEProblem`.

```
function plot_one_traj(sol)
    tmp=[u[1] for u in sol.u]; lines(tmp)
    tmp=[u[2] for u in sol.u]; lines!(tmp)
    tmp=[u[3] for u in sol.u]; lines!(tmp)
    current_figure()
end

n=10; u=zeros(1,1,n); v=zeros(1,1,n); w=zeros(1,1,n+1);
F=FlowFields(u,u,v,v,w,w,[0,1.0])

p=100; x=zeros(100); y=zeros(100); z=rand(p); 
I=Individuals(F,x,y,z,(D=(problem_type=:SDE,),))
sol=∫!(I)

plot_one_traj(sol[1])
```
"""

P=ex_SDE.default_parameters()
cb=ex_SDE.surface_and_bottom_reflect(P.seafloor)

default_solver(prob::SDEProblem) = solve(prob,EM(),dt=P.dt,callback = cb)

function ensemble_solver(prob::SDEProblem;safetycopy=false)
	u0 = prob.u0
	prob_func(prob,i,repeat) = remake(prob,u0=u0[i])
	indiv_prob = SDEProblem(prob.f,prob.g,u0[1],prob.tspan,prob.p)

	ensemble_prob = EnsembleProblem(indiv_prob,prob_func=prob_func,safetycopy=safetycopy)
	solve(ensemble_prob,EM(),trajectories=length(u0),dt=P.dt,callback = cb)
end

end

