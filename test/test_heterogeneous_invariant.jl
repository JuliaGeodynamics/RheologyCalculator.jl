# `second_invariant_value` must accept tuples of mixed component type.
using Test, ForwardDiff, StaticArrays
using RheologyCalculator.RheologyModels
import RheologyCalculator: second_invariant_value, second_invariant

@testset "heterogeneous invariant tuples" begin

    @testset "mixed component types" begin
        @test second_invariant_value((1.0, -1, 0)) == 1.0
        @test second_invariant_value((1, -1, 0)) == 1.0
        @test second_invariant_value((1.0, -1.0, 0.0)) == 1.0
        @test second_invariant_value((1.0, -1, 0, 0, 0, 0.0)) ==
            second_invariant_value((1.0, -1.0, 0.0, 0.0, 0.0, 0.0))
    end

    @testset "homogeneous tuples are unchanged" begin
        for t in ((1.0, -1.0, 0.0), (0.5, 0.25, -0.75), (0.0, 0.0, 0.0))
            @test second_invariant_value(t) == second_invariant(t...)
        end
        @test second_invariant_value(2.0) == 2.0      # the scalar method still applies
    end

    @testset "partial-width ForwardDiff" begin
        d = ForwardDiff.derivative(x -> second_invariant_value((x, -1.0, 0.0)), 1.0)
        @test isfinite(d)
        full = ForwardDiff.gradient(v -> second_invariant_value((v[1], v[2], v[3])), [1.0, -1.0, 0.0])
        @test d ≈ full[1] atol = 1.0e-12

        g = ForwardDiff.gradient(v -> second_invariant_value((v[1], v[2], 0.0)), [1.0, -1.0])
        @test g ≈ full[1:2] atol = 1.0e-12

        d3 = ForwardDiff.derivative(x -> second_invariant_value((x, -1.0, 0.0, 0.0, 0.0, 0.0)), 1.0)
        @test isfinite(d3)
    end
end
