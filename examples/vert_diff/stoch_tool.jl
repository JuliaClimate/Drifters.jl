_at(x, t) = x(t)                 # for functions/functors
_at(x::Number, t) = x            # for constants

function hovmoller_density(sol; nbins=20, xmin=nothing, xmax=nothing, normalize=true,absolute = false)
    ts = sol.t
    if absolute
        us = [abs.(u) for u in sol.u]
    else
        us = sol.u
    end

    allmins = minimum(minimum, us)
    allmaxs = maximum(maximum, us)
    xmin = isnothing(xmin) ? allmins : xmin
    xmax = isnothing(xmax) ? allmaxs : xmax

    x_edges = collect(range(xmin, xmax; length=nbins+1))
    x_centers = (x_edges[1:end-1] .+ x_edges[2:end]) ./ 2
    dx = x_edges[2] - x_edges[1]

    nt = length(ts)
    ρ = zeros(nbins, nt)

    for (k, uk) in enumerate(us)
        h = fit(Histogram, uk, x_edges)
        counts = Float64.(h.weights)
        if normalize
            ρ[:, k] .= counts ./ (sum(counts) * dx)
        else
            ρ[:, k] .= counts
        end
    end

    return ρ,x_centers
end

function surface_reflect()
    condition(u,t,integrator) = true
    function affect!(integrator)
        integrator.u .= abs.(integrator.u)   # reflect all components (vector/array-safe)
    end
    cb = DiscreteCallback(condition, affect!)
    return cb
end

"""
    ocean_reflect(bottom) -> DiscreteCallback

Build a `DifferentialEquations.jl` callback that enforces no-flux (reflecting)
boundary conditions at the sea surface (`u = 0`) and the seafloor (`u = bottom`)
during SDE integration.

Rather than true reflection, this clamps the integrator state into the valid
depth range `[bottom, 0.0]` after every step, which is a cheap and adequate
approximation for a random-walk particle tracker as long as the step size is
small relative to the distance from the boundary.

# Arguments
- `bottom`: the (negative) depth of the seafloor (NOT MIXED LAYER BOTTOM!)

# Returns
A `DiscreteCallback` (condition always `true`) intended to be passed via
`callback = ocean_reflect(bottom)` to `solve`.

# Used by
[`create_sde_updator`](@ref), as the boundary condition for the Euler–Heun
SDE integration of individual particle depths.
"""
function ocean_reflect(bottom)
    condition(u,t,integrator) = true
    function affect!(integrator)
        integrator.u .= clamp.(integrator.u,bottom,0.0)   # reflect all components (vector/array-safe)
    end
    cb = DiscreteCallback(condition, affect!)
    return cb
end

"""
    build_interpolator(rf, kappa) -> Interpolations.Extrapolation

Construct a 1-D linear interpolant `κ(z)` for the vertical diffusivity profile,
suitable for evaluating `κ` at arbitrary (non-grid) depths during SDE
integration.

`kappa` is defined on the depth grid `rf` (grid faces, ordered from surface
`0` to seafloor, most negative). Two adjustments are made before building the
interpolant:
- The diffusivity is extended one cell above the surface by duplicating the
  first sub-surface value (`kappa[2]`), and padded with `0.0` at the deepest
  face, since the underlying `rf`/`kappa` are staggered.
- The corresponding above surface depth grid is set as `+1000` (very 
  far above surface), so that `κ` is treated as constant above
  the surface rather than extrapolated.

# Arguments
- `rf::AbstractVector`: depth grid faces (surface to seafloor, decreasing).
- `kappa::AbstractVector`: diffusivity values at `rf` (same length as `rf`).

# Returns
An interpolation object supporting `itp(z)` and `Interpolations.gradient`,
with flat (constant) extrapolation outside the domain.

"""
function build_interpolator(rf,kappa)
	kap = [kappa[2];kappa[2:end]; 0.0]
	rff = [1000;rf[2:end]]# it is all the same diffusivity above surface
	itp = linear_interpolation(reverse(rff), reverse(kap),
				  extrapolation_bc = Flat())
	return itp
end

"""
    build_kappa_profile(itp) -> (κ_func, f_func, g_func)

Given a diffusivity interpolant `itp` (from [`build_interpolator`](@ref)),
build the three functions needed to define a 1-D SDE for particle depth `z`
undergoing diffusion with spatially-varying diffusivity κ(z):

    dz = f(z) dt + g(z) dW

using the standard Itô-corrected random-walk formulation for
depth-dependent mixing (the "gradient" or Visser correction), so that the
particle distribution converges to the correct steady state rather than
spuriously accumulating in low-diffusivity regions.

# Returns a tuple of three functions, each with signature `(u, p, t)` as
expected by `StochasticDiffEq.jl`:
- `kappa_func`: κ(z), the diffusivity itself.
- `f_func`: the deterministic drift term, `∂κ/∂z` (the Itô correction term).
- `g_func`: the stochastic diffusion term, `sqrt(2κ(z))`, clamped to be
  non-negative to guard against small negative values from interpolation.

# Used by
[`create_sde_updator`](@ref), to construct the `SDEProblem` solved with
`EulerHeun()`.
"""
function build_kappa_profile(itp)
	kappa_func(u,p,t) = itp.(u)
	f_func(u,p,t) = Float64[Interpolations.gradient(itp, z)[1] for z in u]
	g_func(u,p,t) = sqrt.(max.(itp.(u), 0.0) .* 2.0)
	return kappa_func,f_func,g_func
end

"""
    create_sde_updator(kappa, rf) -> updator::Function

Build a direct SDE-based particle updator: given a diffusivity profile
`kappa` on depth grid `rf`, returns a function `updator(u_0, tf; dt=10)` that
integrates each particle's depth forward by time `tf` (seconds) using
Euler–Heun integration of the vertical diffusion SDE, reflecting particles
off the surface and seafloor via [`ocean_reflect`].

This is the fallback path used by
[`create_hybrid_updator`] when the fast eigen-decomposition/PDE method
(see [`create_always_pde_updator`] would be inaccurate  when a
particle's displacement is small relative to the grid spacing, so the
discrete PDE solution is poorly resolved. 

This method should not be used under large diffusivity or long time, 
even if you plan to use very short time step. 

# Arguments
- `kappa::AbstractVector`: diffusivity at each face of `rf`.
- `rf::AbstractVector`: depth grid faces (surface to seafloor).

# Returns
`updator(u_0, tf; dt=10) -> Vector`, where `u_0` are initial particle depths,
`tf` is the integration time (seconds), and `dt` is the internal SDE step
size. Returns the updated particle depths after time `tf`.

"""
function create_sde_updator(kappa,rf)
	itp = build_interpolator(rf,kappa)
	κ,𝒻,ℊ = build_kappa_profile(itp)
	function updator(u_0,tf;dt=10)
		tspan = (0.0,tf)
		p = nothing
		prob = SDEProblem(𝒻,ℊ,u_0,tspan,p)
		sol=solve(prob,EulerHeun(),dt=dt,callback = ocean_reflect(rf[end]));
		return sol.u[end]
	end
	return updator
end

"""
    build_diffusion_matrix(kappa_edges, drc, drf) -> A::Matrix{Float64}

Construct the finite-volume operator matrix `A` such that `dc/dt = A * c`
for the 1-D vertical diffusion equation

    ∂c/∂t = ∂/∂z (κ ∂c/∂z)

on a staggered grid, with implicit no-flux boundary conditions (achieved by
`kappa_edges` being `0` at the top and bottom faces, i.e. no flux b.c.)).

# Grid layout
- `kappa_edges[k]`: diffusivity at grid face `k` (must vanish at the first
  and last face for no-flux boundaries).
- Concentrations/probabilities live at cell centers; `drc[k]` is the
  center-to-center spacing and `drf[k]` is the cell thickness.

# Arguments
- `kappa_edges::AbstractVector`: diffusivity at each face, length `Nc+1`.
- `drc::AbstractVector`: center-to-center spacing, length `Nc`.
- `drf::AbstractVector`: cell thickness, length `Nc`.

# Returns
- `A::Matrix{Float64}` of size `(Nc, Nc)`, the discrete diffusion operator.

# Usage
    c(t) = exp(A * t) * c₀

"""
function build_diffusion_matrix(kappa_edges::AbstractVector, 
								drc::AbstractVector,
								drf::AbstractVector,
							   )
    Nc = length(drf)

    A = zeros(Nc, Nc)

    for k in 1:Nc
        if k < Nc
            coeff = kappa_edges[k+1] / (drc[k+1] * drf[k])
            A[k, k]   -= coeff
            A[k, k+1] += coeff 
        end

        if k > 1
            coeff = kappa_edges[k] / (drc[k] * drf[k])
            A[k, k]   -= coeff
            A[k, k-1] += coeff
        end
    end

    return A
end

"""
    apply_matrix_exp(eigs, c0, t) -> Vector

Propagate a concentration/probability vector `c0` forward by time `t` under
the linear diffusion operator whose eigen-decomposition is `eigs` (as
returned by `LinearAlgebra.eigen` applied to a matrix from
[`build_diffusion_matrix`](@ref)), i.e. compute `exp(A * t) * c0` without
forming the matrix exponential explicitly.

Using the eigen-decomposition once and reapplying it here is what makes the
PDE-based updators ([`create_always_pde_updator`](@ref),
[`create_hybrid_updator`](@ref)) fast: the expensive `eigen(A)` call happens
once per diffusivity profile, and every particle group's forward propagation
is then just a change of basis, a diagonal scaling, and a change back.

# Arguments
- `eigs`: an `Eigen` factorization object (`eigs.values`, `eigs.vectors`).
- `c0::AbstractVector`: initial probability/concentration vector.
- `t`: propagation time.

# Returns
The propagated vector `exp(A * t) * c0`.

"""
function apply_matrix_exp(eigs, c0, t)
    Vt_c0 = eigs.vectors \ c0             # V⁻¹ * c0 (not V' * c0)
    Vt_c0 .*= exp.(eigs.values .* t)
    return eigs.vectors * Vt_c0
end

"""
    pt_depth_from_prob(n, pdf, rf, drf) -> depths::Vector

Draw `n` random depths by inverse-transform sampling from a discretized
probability density `pdf` defined on grid cells with faces `rf` and
thicknesses `drf`.

The per-cell mass (`pdf .* drf`, clamped non-negative) is accumulated into a
CDF, degenerate (zero-mass) intervals are masked out so the interpolant is
well-defined, and `n` uniform random numbers are mapped through the inverse
CDF (linearly interpolated) to produce sampled depths.

# Arguments
- `n::Integer`: number of samples to draw.
- `pdf::AbstractVector`: probability density per cell (not necessarily
  normalized; the CDF is renormalized internally).
- `rf::AbstractVector`: grid face depths corresponding to `pdf`'s cells.
- `drf::AbstractVector`: cell thicknesses.

# Returns
A `Vector` of `n` sampled depths.

"""
function pt_depth_from_prob(n, pdf, rf, drf)
    mass = max.(pdf, 0.0) .* drf

    cdf = [0.0; cumsum(mass)]
    cdf ./= cdf[end]

    mask = [true; diff(cdf) .> 0]

    itp = linear_interpolation(cdf[mask], rf[mask])

    depths = itp.(rand(n))

    return depths
end

"""
    find_bin(z_arr, rf) -> bin_indices::Vector{Int}

Locate the grid cell index (1-based, into the `drf`/cell-center arrays) that
each depth in `z_arr` falls into, given face depths `rf` ordered from
surface (`0`) to seafloor (most negative). Essential step for creating onehot vector. 

Internally reverses `rf` to increasing order to use `searchsortedlast`, then
maps the result back to the original (surface-to-seafloor) cell indexing,
clamped to the valid range `[1, Nbins]`.

# Arguments
- `z_arr::AbstractVector`: particle depths to bin.
- `rf::AbstractVector`: grid face depths (surface to seafloor, decreasing).

# Returns
A `Vector{Int}` of cell indices, same length as `z_arr`.

"""
function find_bin(z_arr, rf)
    rf_inc = reverse(rf)
    Nbins = length(rf) - 1
    return Nbins .+ 1 .- clamp.(searchsortedlast.(Ref(rf_inc), z_arr), 1, Nbins)
end

"""
    onehot(k, n) -> Vector{Float64}

Return a length-`n` vector that is `1.0` at index `k` and `0.0` elsewhere —
i.e. a point-mass probability distribution concentrated in cell `k`.

Used as the initial condition `p0` when propagating a delta-function
distribution of particles (all starting in the same grid cell) forward
under diffusion via [`apply_matrix_exp`](@ref).

# Arguments
- `k::Integer`: index of the "hot" entry.
- `n::Integer`: length of the output vector.

"""
onehot(k, n) = [i == k ? 1.0 : 0.0 for i in 1:n]

"""
    create_always_pde_updator(kappa, coords) -> updator::Function

Build a particle-depth updator that *always* advances particles using the
exact finite-volume PDE / eigen-decomposition method (no SDE fallback),
given a diffusivity profile `kappa` and grid `coords` (as produced by
`make_coords`, with fields `RF`, `RC`, `DRF`, `DRC`).

The diffusion matrix is built and eigen-decomposed once, up front, since
`kappa` and `coords` are fixed for the lifetime of the returned closure.
Each call to the returned `updator(u_0, tf)` then:
1. Bins the initial particle depths `u_0` into grid cells via
   [`find_bin`](@ref).
2. For each occupied cell, propagates the corresponding one-hot probability
   distribution forward by time `tf` via [`apply_matrix_exp`](@ref).
3. Resamples new depths for all particles that started in that cell from the
   propagated distribution via [`pt_depth_from_prob`](@ref).

This is efficient because particles sharing a starting bin share a single
eigen-propagation step, but it can be inaccurate when a particle's
displacement is small relative to the grid spacing (see
[`create_hybrid_updator`](@ref) for an adaptive alternative).

# Arguments
- `kappa::AbstractVector`: diffusivity at each face of `coords.RF`.
- `coords::NamedTuple`: grid definition with fields `RF`, `RC`, `DRF`, `DRC`.

# Returns
`updator(u_0, tf) -> Vector`: given initial particle depths `u_0` and a
(possibly negative — sign is ignored) elapsed time `tf`, returns updated
depths after diffusive mixing.
"""
function create_always_pde_updator(kappa, coords)
	M_d = build_diffusion_matrix(kappa,coords.DRC,coords.DRF)
	eig_diff = eigen(M_d);
	# @assert all(eig_diff.values .≤ 1e-14)
	function updator(u_0,tf)
		tf = abs(tf)
		u_f = similar(u_0)
		lookup = find_bin(u_0,coords.RF)
		for iz in unique(lookup)
			# println(iz)
			p0 = onehot(iz, length(coords.DRF))
			pf = apply_matrix_exp(eig_diff, p0, tf)
			mask = lookup.==iz
			n_pt = sum(mask)
			u_f[mask] = pt_depth_from_prob(n_pt, pf, coords.RF, coords.DRF)
		end
		return u_f
	end
end

"""
    create_hybrid_updator(kappa, coords) -> updator::Function

Build a particle-depth updator that combines the fast PDE/eigen-decomposition
method with a direct SDE fallback, switching between them per-particle-group
based on whether the PDE solution can be trusted.

Like [`create_always_pde_updator`](@ref), the diffusion matrix for `kappa` on
`coords` is built and eigen-decomposed once. For each call to
`updator(u_0, tf)`, particles are grouped by starting depth bin
([`find_bin`](@ref)), and for each group the propagated probability
distribution is computed via [`apply_matrix_exp`](@ref). That distribution is
then checked for two conditions before it's trusted:

- **Well-resolved spread**: the distribution's standard deviation (first and
  second moments) must be at least twice the starting cell's thickness —
  otherwise the particles haven't diffused far enough for the discrete grid
  to resolve their new distribution accurately.
- **Well-conditioned diffusivity**: among cells with non-negligible
  probability mass, the ratio of max to min diffusivity must be below `2` —
  otherwise sharp local variation in κ makes the discrete solution
  unreliable.

Groups that fail either check are collected into a `revisit` list and
instead integrated directly with [`create_sde_updator`](@ref) (with
`dt = tf/100`), which is slower but accurate regardless of displacement
size. All other groups are resolved via the fast PDE/resampling path
([`pt_depth_from_prob`](@ref)), exactly as in
[`create_always_pde_updator`](@ref).

# Arguments
- `kappa::AbstractVector`: diffusivity at each face of `coords.RF`.
- `coords::NamedTuple`: grid definition with fields `RF`, `RC`, `DRF`, `DRC`.

# Returns
`updator(u_0, tf) -> Vector`: given initial particle depths `u_0` and a
(possibly negative — sign is ignored) elapsed time `tf`, returns updated
depths after diffusive mixing, using the SDE fallback only where needed for
accuracy.
"""
function create_hybrid_updator(kappa, coords)
    M_d = build_diffusion_matrix(kappa, coords.DRC, coords.DRF)
    eig_diff = LinearAlgebra.eigen(M_d)

    function updator(u_0, tf)
        tf = abs(tf)
        u_f = similar(u_0)
        lookup = find_bin(u_0, coords.RF)
        revisit = Int[]

        for iz in unique(lookup)
            p0 = onehot(iz, length(coords.DRF))/coords.DRF[iz]
            pf = apply_matrix_exp(eig_diff, p0, tf)

			final_kappa = kappa[pf.*coords.DRF.>0.005]
			cond_number = maximum(final_kappa)/minimum(final_kappa)
			kappa_well_conditioned = cond_number<2

			first_moment = sum(pf./coords.DRF.*coords.RC)
			second_moment = sum(pf./coords.DRF.*(coords.RC.^2))
			spread = sqrt(second_moment - first_moment^2)
            if spread<2*coords.DRF[iz] && kappa_well_conditioned
                push!(revisit, iz)
				continue
            end
            mask = lookup .== iz
            n_pt = sum(mask)
            u_f[mask] = pt_depth_from_prob(n_pt, pf, coords.RF, coords.DRF)
        end

        if length(revisit) > 0
            sde_updator = create_sde_updator(kappa, coords.RF)
            revisit_mask = [iz in revisit for iz in lookup]
            u_f[revisit_mask] = sde_updator(u_0[revisit_mask], tf; dt=tf/100)
        end

        return u_f
    end

    return updator
end

