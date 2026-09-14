# Frozen pre-optimization solver for interleaved benchmark comparisons.
module SolverReference
using RheologyCalculator, StaticArrays
import RheologyCalculator: AbstractCompositeModel, _direct_leaf_elastic_correction,
    second_invariant_value, correct_xnorm, mynorm, branch_strain_rate_mask,
    backsolve, max_feasible_step

function solve_reference(c::AbstractCompositeModel, x::SVector, vars0, others; xnorm0 = nothing, atol::Float64 = 1.0e-12, rtol::Float64 = 1.0e-12, itermax = 1.0e4, verbose::Bool = false)
    # Pre-correct ONLY the direct elastic leafs of the outer composite
    # (simple Maxwell backstress).  Tensor arithmetic is used here so that
    # second_invariant(ε + τ0/(2G·dt)) is evaluated correctly even for
    # non-coaxial ε/τ0 pairs.
    # ParallelModel branch corrections are handled implicitly inside
    # compute_residual via subtract_elastic_correction, so they must NOT be
    # included here to avoid double-counting.
    ε_corr = _direct_leaf_elastic_correction(c, vars0.ε, others)
    εII = second_invariant_value(vars0.ε .+ ε_corr)
    vars = merge(vars0, (; ε = εII))

    # vars = merge((; ε = εII), vars0)
    xnorm = correct_xnorm(x, xnorm0)
    r = compute_residual(c, x, vars, others)   # initial residual
    it = 0
    er0 = mynorm(r, xnorm)
    # `oftype` keeps the residual a single type across the loop, so that the
    # value stored in the returned `RCSolution` is inferrable.
    er = oftype(er0, Inf)

    nonneg = branch_strain_rate_mask(c)

    α = 1.0e0
    stagnant_iters = 0
    while er > atol && er > rtol * er0
        it += 1

        J = jacobian(c, x, vars, others)
        Δx = backsolve(J, r)
        α = max_feasible_step(x, Δx, nonneg)
        α = bt_line_search_reference(
            Δx, x, c, vars, others, xnorm, er;
            α = α, ρ = 0.5, lstol = 0.95, α_min = 0.1
        )
        x_next = x + α .* Δx

        # check convergence
        r = compute_residual(c, x_next, vars, others)
        er = mynorm(r, xnorm)

        # Once the update is below floating-point resolution, continuing the
        # Newton iteration cannot change either the iterate or its residual.
        # Also require the residual to be finite: non-finite residuals follow
        # the existing diagnostic path below. A decreasing but slow residual
        # must remain an iteration-limit failure, not a stagnation failure.
        if isfinite(er) && x_next == x
            stagnant_iters += 1
        else
            stagnant_iters = 0
        end
        x = x_next

        it > itermax && break

        if stagnant_iters ≥ 3
            throw(NonConvergenceError(it, er, x, :stagnation))
        end

        # ε_corr = effective_strain_rate_correction(c, vars0.ε, others.τ0, others)
        # ε_eff  = vars0.ε .+ ε_corr
        # εII    = second_invariant_value(ε_eff)
        # vars   = merge(vars0, (; ε = εII)) # this mames it type unstable; TODO

    end
    if verbose && it > 1
        println("Iterations: $it, Error: $er, α = $α")
    end
    # A NaN residual compares false against both tolerances and so exits the loop
    # by the same door as a converged one; `isfinite` is what separates them.
    isfinite(er) && (er ≤ atol || er ≤ rtol * er0) || throw(NonConvergenceError(it, er, x))
    return RCSolution(x, it, er)
end

function bt_line_search_reference(Δx, x, composite, vars, others, xnorm, rnorm; α = 1.0, ρ = 0.5, lstol = 0.9, α_min = 1.0e-8)

    α_initial = α
    best_α = α
    best_rnorm = Inf

    while α ≥ α_min
        # Apply scaled update
        perturbed_x = @. x + α * Δx

        # Get updated residual
        perturbed_r = compute_residual(composite, perturbed_x, vars, others)
        perturbed_rnorm = mynorm(perturbed_r, xnorm)

        if isfinite(perturbed_rnorm) && perturbed_rnorm < best_rnorm
            best_α, best_rnorm = α, perturbed_rnorm
        end

        # The first feasible trial may be limited by non-negativity, so do not
        # demand artificial decrease from a step that merely avoids growth.
        target = α == α_initial ? 1.0 : lstol
        if isfinite(perturbed_rnorm) && perturbed_rnorm ≤ target * rnorm
            return α
        end

        # Bisect step length
        α *= ρ
    end

    return best_α
end


end
