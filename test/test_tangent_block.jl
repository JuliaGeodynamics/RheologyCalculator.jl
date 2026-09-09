using Test
using StaticArrays
using RheologyCalculator
using RheologyCalculator.RheologyModels

@testset "deviatoric/volumetric tangent block" begin
    G, K, dt = 2.0e10, 3.0e10, 4.0e6
    c = SeriesModel(Elasticity(G, K))
    vars = (; ε = 1.0e-14, θ = 2.0e-15)
    others = (; dt, τ0 = (0.0,), P0 = (0.0,))
    x = initial_guess_x(c, vars, (; τ = 0.0, P = 0.0), others)
    sol = solve(c, x, vars, others)
    D = tangent_block(c, sol, vars, others)

    @test D isa SMatrix{2, 2}
    @test D ≈ SMatrix{2, 2}(2G * dt, 0.0, 0.0, -K * dt)
    @test isbitstype(typeof(D))
end
