# A zero characteristic strain rate must not make residual rows unmeasurable.
using Test, StaticArrays
using RheologyCalculator.RheologyModels
import RheologyCalculator: SeriesModel, initial_guess_x, normalisation_x, solve, mynorm,
    compute_residual

@testset "zero characteristic strain rate" begin
    viscous = LinearViscosity(1.0e23)
    elastic = Elasticity(1.0e10, 2.0e11)

    @testset "mynorm measures rows with a zero factor" begin
        r = SA[1.0e-3, 2.0]
        @test mynorm(r, SA[0.0, 1.0]) == 1.0e-3 + 2.0
        @test mynorm(r, SA[0.0, 0.0]) == 1.0e-3 + 2.0
        @test mynorm(r, SA[1.0e-3, 2.0]) == 2.0
    end

    @testset "normalisation_x floors a vanishing scale" begin
        c = SeriesModel(viscous, elastic)
        xn = normalisation_x(c, 1.0e6, 0.0)
        @test all(!iszero, xn)
        @test normalisation_x(c, 1.0e6, 3.0) == normalisation_x(c, 1.0e6, 3.0)
        @test any(==(3.0), normalisation_x(c, 1.0e6, 3.0))
    end

    @testset "undeformed point converges to a true root" begin
        c = SeriesModel(viscous, elastic)
        vars = vars_2D(0.0, 0.0)                       # zero strain rate
        others = (; dt = 1.0e5, τ0 = ((5.0e5, -5.0e5, 0.0),), P0 = (0.3e6,))
        char_ε = second_invariant_2D(vars.ε) + abs(vars.θ)
        @test iszero(char_ε)
        xnorm = normalisation_x(c, 1.0e6, char_ε)

        x0 = initial_guess_x(c, vars, (; τ = 0.0, P = 0.3e6), others)
        sol = solve(c, x0, vars, others; xnorm0 = xnorm, itermax = 50)
        @test sol.residual ≤ 1.0e-12

        r = compute_residual(c, sol.x, (; ε = 0.0, θ = vars.θ), others)
        @test mynorm(r, SA[1.0, 1.0]) ≤ 1.0e-9

        sol2 = solve(c, SA[1.0e9, 1.0e9], vars, others; xnorm0 = xnorm, itermax = 50)
        r2 = compute_residual(c, sol2.x, (; ε = 0.0, θ = vars.θ), others)
        @test mynorm(r2, SA[1.0, 1.0]) ≤ 1.0e-9
        @test sol2.x ≈ sol.x rtol = 1.0e-10
    end
end
