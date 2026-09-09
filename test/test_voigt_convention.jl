using Test
using RheologyCalculator.RheologyModels
import RheologyCalculator.RheologyModels: zero_stress_tensor_2D,
    zero_stress_tensor_3D, second_invariant_2D, second_invariant_3D

@testset "Voigt tensor convention" begin
    @test zero_stress_tensor_2D() == (0.0, 0.0, 0.0)
    @test zero_stress_tensor_3D() == (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    @test second_invariant_2D((1.0, -1.0, 0.0)) ≈ 1.0
    @test second_invariant_3D((1.0, -1.0, 0.0, 0.0, 0.0, 0.0)) ≈ 1.0
end
