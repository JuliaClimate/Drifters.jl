module DriftersStochasticDiffEqExt

import StochasticDiffEq, Drifters
import StochasticDiffEq: EnsembleProblem, SDEProblem, EM, solve, stack
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
        sol=solve(prob,EM(),dt=P.dt,callback = ex_SDE.surface_and_bottom_reflect(P.seafloor))
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
    _SDEProblem(f::Function,g::Vector,u₀,tspan)
    
Calls `StochasticDiffEq.SDEProblem(f,g,u₀,tspan)``

```
using Drifters, StochasticDiffEq
?Drifters._SDEProblem
```
"""
_SDEProblem(f::Function,g::Vector,u₀,tspan)=SDEProblem(f,g,u₀,tspan)

"""
    solve_prob(prob,dt=dt)=solve(prob,EM(),dt=dt)
    
```
using Drifters, StochasticDiffEq
SDE = Base.get_extension(Drifters, :DriftersStochasticDiffEqExt)
?SDE.solve_prob
```
"""
default_solver(prob::SDEProblem) = error("yet to be implemented")
#solve(prob,EM(),dt=dt)

function ensemble_solver(prob::SDEProblem;solver=Tsit5(),reltol=1e-8,abstol=1e-8,safetycopy=false)
    error("yet to be implemented")
end

end

