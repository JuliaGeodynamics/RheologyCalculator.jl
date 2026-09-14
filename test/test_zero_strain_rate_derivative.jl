# The second invariant is a norm and so is non-smooth at the zero tensor. Its
# VALUE at zero was already correct (√0 == 0), but the ForwardDiff derivative
# came back NaN, because `d√s/ds = 1/(2√s)` is Inf at s = 0 and Inf·0 = NaN in
# the dual part. A run started from rest is exactly at that point, so every
# quadrature point produced NaN tangents.
using Test, ForwardDiff, StaticArrays
using RheologyCalculator.RheologyModels
import RheologyCalculator: second_invariant, second_invariant_value

@testset "derivatives at a vanishing strain rate" begin

    @testset "values are unchanged" begin
        @test second_invariant(0.0, 0.0, 0.0) == 0.0
        @test second_invariant(0.0, 0.0, 0.0, 0.0, 0.0, 0.0) == 0.0
        @test second_invariant(1.0, -1.0, 0.0) == 1.0
        @test second_invariant_value((0.0, 0.0, 0.0)) == 0.0
    end

    @testset "derivative at exactly zero is finite" begin
        g2 = ForwardDiff.gradient(v -> second_invariant(v[1], v[2], v[3]), zeros(3))
        @test all(isfinite, g2)
        @test all(iszero, g2)

        g3 = ForwardDiff.gradient(v -> second_invariant(v...), zeros(6))
        @test all(isfinite, g3)
        @test all(iszero, g3)

        d = ForwardDiff.derivative(x -> second_invariant_value((x, -x, zero(x))), 0.0)
        @test isfinite(d)
        @test iszero(d)
    end

    @testset "away from zero the derivative is untouched" begin
        # the 2D form reconstructs zz = -xx - yy, so at pure shear (1, -1, 0)
        # the invariant is 1 and its gradient is (0.5, -0.5, 0)
        g = ForwardDiff.gradient(v -> second_invariant(v[1], v[2], v[3]), [1.0, -1.0, 0.0])
        @test all(isfinite, g)
        @test g ≈ [0.5, -0.5, 0.0] atol = 1.0e-12

        d = ForwardDiff.derivative(x -> second_invariant_value((x, -x, zero(x))), 2.0)
        @test d ≈ 1.0 atol = 1.0e-12

        # and a tiny but nonzero argument is still the ordinary derivative
        @test all(isfinite, ForwardDiff.gradient(v -> second_invariant(v[1], v[2], v[3]), [1.0e-30, 0.0, 0.0]))
    end

    @testset "stress reconstruction at a zero effective strain rate" begin
        zero_eps = (0.0, 0.0, 0.0)
        @test stress_tensor_from_invariant_2D(1.0e6, zero_eps) == (0.0, 0.0, 0.0)
        dτ = ForwardDiff.derivative(t -> sum(stress_tensor_from_invariant_2D(t, zero_eps)), 1.0e6)
        @test isfinite(dτ)
        dε = ForwardDiff.derivative(
            e -> sum(stress_tensor_from_invariant_2D(1.0e6, (e, -e, zero(e)))), 0.0
        )
        @test isfinite(dε)
    end
end
