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
