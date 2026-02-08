module DriftersStochasticDiffEqExt

import StochasticDiffEq, Drifters
import StochasticDiffEq: EnsembleProblem, SDEProblem, EM, solve, stack
import Drifters: _SDEProblem, default_solver, ensemble_solver

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
    step!(u₀)

Sequence of `solve_paths` and `fold_tails`.

```
using Drifters, StochasticDiffEq
SDE = Base.get_extension(Drifters, :DriftersStochasticDiffEqExt)
?SDE.step!
```
"""
function step!(u₀)
    za=solve_paths(u₀)
    fold_tails(za)
    u₀[:]=za[:,end]
end

function solve_paths(u₀)
    dt=1e-3
    f(u,p,t) = 0.0
    g(u,p,t) = 0.1
    #dt = 1//2^(4)
    tspan = (0.0,1.0)
    prob = SDEProblem(f,g,u₀,(0.0,1.0))
    sol=solve(prob,EM(),dt=dt)
    stack(sol(0:0.01:1))
end

function fold_tails(z)
    while !isempty(findall( xor.(z.>1.0,z.<0.0) ))
        z[findall(z.<0.0)].=-z[findall(z.<0.0)]
        z[findall(z.>1.0)].=(2.0 .- z[findall(z.>1.0)])
    end
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

