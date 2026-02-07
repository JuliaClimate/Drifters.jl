module DriftersStochasticDiffEqExt

import StochasticDiffEq, Drifters
import StochasticDiffEq: EnsembleProblem, SDEProblem, EM, solve, stack

notes()=println("""
- For now : `step!(u₀)` is the sequence of `solve_paths` and `fold_tails`.
- For deeper integration into Drifters, these will need to overloaded.

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
    set_prob(f,g,u₀,tspan)=SDEProblem(f,g,u₀,tspan)

```
using Drifters, StochasticDiffEq
SDE = Base.get_extension(Drifters, :DriftersStochasticDiffEqExt)
?SDE.set_prob
```
"""
set_prob(f,g,u₀,tspan)=SDEProblem(f,g,u₀,tspan)

"""
    solve_prob(prob,dt=dt)=solve(prob,EM(),dt=dt)
    
```
using Drifters, StochasticDiffEq
SDE = Base.get_extension(Drifters, :DriftersStochasticDiffEqExt)
?SDE.solve_prob
```
"""
solve_prob(prob,dt=dt)=solve(prob,EM(),dt=dt)

end

