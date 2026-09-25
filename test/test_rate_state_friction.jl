using RheologyCalculator, Test
using RheologyCalculator.RheologyModels
using RheologyCalculator.RheologyModels: RateStateFriction, update_Ω, vars_2D, second_invariant_2D, zero_stress_tensor_2D
import RheologyCalculator: compute_stress

@testset "RateStateFriction in a parallel branch" begin
    rs = RateStateFriction(0, 0.5, 4.0e-9, 0.011, 0.015, 0.0047, 0, 500)
    G, η, P, dt, Ω_old = 20.0e9, 1.0e12, 50.0e6, 1.0e6, -5.0
    c = SeriesModel(IncompressibleElasticity(G), ParallelModel(rs, LinearViscosity(η)))
    @test x_keys(c) === (:τ, :ε)

    vars = vars_2D(1.0e-9, 1.0e-20)
    εII = second_invariant_2D(vars.ε)
    others = (; dt, τ0 = (zero_stress_tensor_2D(),), P0 = (0.0,), P, Ω_old)
    x0 = initial_guess_x(c, vars, (; τ = 0.0, P), others)
    x = solve(c, x0, vars, others; xnorm0 = normalisation_x(c, 1.0e6, εII))

    # Reference by bisection on the branch strain rate: the branch stress is the
    # sum of both elements' stresses, and the elastic strain rate makes up the rest.
    τ_branch(εb) = compute_stress(rs; ε = εb, Ω_old, P, dt) + 2η * εb
    f(εb) = τ_branch(εb) / (2G * dt) + εb - εII
    lo, hi = 0.0, εII
    for _ in 1:200
        m = (lo + hi) / 2
        f(m) > 0 ? (hi = m) : (lo = m)
    end
    @test 0 < x[2] < εII
    @test x[1] ≈ τ_branch(lo) rtol = 1.0e-12
    @test x[2] ≈ lo rtol = 1.0e-12
end

@testset "first Newton step is line-searched" begin
    # From this seed the first Newton step lands far past the root, where the
    # slip rate blows up (a residual of 1.5e36). If that step is accepted,
    # Newton needs 89 iterations to walk back instead of 10.
    rs = RateStateFriction(0.0, 0.2, 4.0e-9, 0.011, 0.017, 0.05, 0.0, 625.0)
    c = SeriesModel(IncompressibleElasticity(3.0e10), rs)
    vars = (; ε = 1.0e-9)
    others = (; dt = 1.0e6, τ0 = (0.0,), P0 = (0.0,), P = 5.0e7, Ω_old = -5.0)
    xnorm = normalisation_x(c, 1.0e6, 1.0e-9)
    x0 = 0.7 .* initial_guess_x(c, vars, (; τ = 0.0, P = 5.0e7), others)

    # itermax = 0 with zero tolerances stops after exactly one iteration
    first_step = try
        solve(c, x0, vars, others; xnorm0 = xnorm, itermax = 0, atol = 0.0, rtol = 0.0)
    catch e
        e
    end
    @test first_step isa NonConvergenceError
    @test first_step.iterations == 1
    @test first_step.residual < 1.0e3

    sol = solve(c, x0, vars, others; xnorm0 = xnorm)
    @test sol.iterations ≤ 15
    @test sol.x[1] ≈ 8.822062639363287e6 rtol = 1.0e-12
end

@testset "RateStateFriction state update stays finite" begin
    rs = RateStateFriction(0, 0.5, 4.0e-9, 0.011, 0.015, 0.0047, 0, 500)
    # An overflowed slip rate and a fully relaxed state at dt = 0 both drive
    # Θ = exp(Ω) to zero before the log.
    @test isfinite(update_Ω(rs; ε = Inf, Ω_old = 0.0, dt = 1.0))
    @test isfinite(update_Ω(rs; ε = 0.0, Ω_old = -800.0, dt = 0.0))
    # Away from the floor the update is unchanged.
    @test update_Ω(rs; ε = 1.0e-12, Ω_old = 2.0, dt = 1.0e3) ≈ 1.9999024186269614
    @test update_Ω(rs; ε = 1.0e-6, Ω_old = 2.0, dt = 1.0e3) ≈ -12.429216196844383
end
