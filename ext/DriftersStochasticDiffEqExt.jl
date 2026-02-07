module DriftersStochasticDiffEqExt

import StochasticDiffEq, Drifters
import StochasticDiffEq: EnsembleProblem, SDEProblem, EM, solve

notes()=println("""
- For now : `set_prob` and `solve_prob` can be used via the extenstion.
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

