using Test
using LinearAlgebra
using StaticArrays
using RheologyCalculator
using RheologyCalculator.RheologyModels
import RheologyCalculator.RheologyModels: tensor_strain_rate_2D, tensor_strain_rate_3D

@testset "isotropic deviatoric tensor tangent" begin
    η = 3.0e3
    c = SeriesModel(LinearViscosity(η))

    @testset "2D Voigt" begin
        vars = (; ε = tensor_strain_rate_2D(1.0e-3), θ = 0.0)
        x = initial_guess_x(c, vars, (; τ = 0.0), (;))
        sol = solve(c, x, vars, (;))
        D = tangent_tensor(c, sol, vars, (;))
        @test D isa SMatrix{3, 3}
        @test D ≈ 2η * I(3)
        @test isbitstype(typeof(D))
    end

    @testset "3D Voigt" begin
        vars = (; ε = tensor_strain_rate_3D(1.0e-3), θ = 0.0)
        x = initial_guess_x(c, vars, (; τ = 0.0), (;))
        sol = solve(c, x, vars, (;))
        D = tangent_tensor(c, sol, vars, (;))
        @test D isa SMatrix{6, 6}
        @test D ≈ 2η * I(6)
        @test isbitstype(typeof(D))
    end
end

# Outer stress tensor τ = τII E / EII of a converged solve at strain rate ε.
function solved_stress(c, ε, others, x0, xnorm0)
    vars = (; ε, θ = 0.0)
    x = solve(c, x0, vars, others; xnorm0, rtol = 1.0e-14, itermax = 200).x
    E = RheologyCalculator.assembled_strain_rate(c, generate_equations(c), x, ε, others)
    EII = RheologyCalculator.second_invariant_value(E)
    return SVector(map(e -> x[stress_index(c)] * e / EII, E))
end

@testset "tensor tangent with non-coaxial history matches finite differences" begin
    dt = 1.0e10
    η₁, η₂ = LinearViscosity(1.0e21), LinearViscosity(5.0e20)
    G = IncompressibleElasticity(1.0e10)
    PL = PowerLawViscosity(2.5e34, 3)
    ε = (1.0e-14, -1.0e-14, 3.0e-15)
    # History rotated away from ε, so E = ε + H is not coaxial with ε.
    others = (; dt, τ0 = ((2.0e6, 1.0e6, -3.0e6),))
    for c in (
            SeriesModel(η₁, ParallelModel(η₂, G)),
            SeriesModel(η₁, ParallelModel(PL, G)),
            SeriesModel(η₁, ParallelModel(η₂, SeriesModel(PL, G))),
        )
        vars = (; ε, θ = 0.0)
        xnorm0 = normalisation_x(c, 1.0e7, 1.0e-14)
        x = solve(c, initial_guess_x(c, vars, (; τ = 1.0e7), others), vars, others; xnorm0, rtol = 1.0e-14, itermax = 200).x
        D = tangent_tensor(c, x, vars, others)

        h = 1.0e-6 * 1.0e-14
        D_fd = hcat(
            ntuple(3) do j
                δ = ntuple(k -> k == j ? h : 0.0, 3)
                (solved_stress(c, ε .+ δ, others, x, xnorm0) - solved_stress(c, ε .- δ, others, x, xnorm0)) / (2h)
            end...
        )
        @test D ≈ D_fd rtol = 1.0e-5

        # `tangent` is dτII/ds along ε + s ε/εII.
        εII = RheologyCalculator.second_invariant_value(ε)
        τII(s) = RheologyCalculator.second_invariant_value(Tuple(solved_stress(c, ε .* (1 + s / εII), others, x, xnorm0)))
        @test tangent(c, x, vars, others) ≈ (τII(h) - τII(-h)) / (2h) rtol = 1.0e-5
    end
end
