using Test
using LinearAlgebra
using ForwardDiff
using RheologyCalculator
using RheologyCalculator.RheologyModels

@testset "Jacobian backend" begin
    c = SeriesModel(LinearViscosity(1.0e22), IncompressibleElasticity(1.0e10))
    vars = (; ε = 1.0e-14, θ = 0.0)
    others = (; dt = 1.0e10, τ0 = (0.0,), P0 = (0.0,))
    x = initial_guess_x(c, vars, (; τ = 0.0, P = 0.0), others)

    J = jacobian(c, x, vars, others)
    reference = ForwardDiff.jacobian(y -> compute_residual(c, y, vars, others), x)
    @test J == reference
    @test size(J) == (length(x), length(x))
end
