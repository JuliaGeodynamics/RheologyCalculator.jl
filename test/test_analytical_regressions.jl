using Test
using RheologyCalculator
using RheologyCalculator.RheologyModels

@testset "analytical regressions" begin
    ε = 1.0e-14

    @testset "linear viscosity" begin
        η = 3.0e21
        c = SeriesModel(LinearViscosity(η))
        vars = (; ε)
        others = (;)
        x0 = initial_guess_x(c, vars, (; τ = 0.0), others)
        sol = solve(c, x0, vars, others)
        @test sol.x[1] ≈ 2η * ε
        @test compute_residual(c, sol.x, vars, others)[1] ≈ 0.0 atol = 1.0e-28
    end

    @testset "elasticity" begin
        G = 2.0e10
        dt = 4.0e6
        c = SeriesModel(IncompressibleElasticity(G))
        vars = (; ε)
        others = (; dt, τ0 = (0.0,))
        x0 = initial_guess_x(c, vars, (; τ = 0.0), others)
        sol = solve(c, x0, vars, others)
        @test sol.x[1] ≈ 2G * dt * ε
    end
end
