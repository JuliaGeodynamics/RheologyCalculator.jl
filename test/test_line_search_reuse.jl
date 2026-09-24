using Test, StaticArrays
using RheologyCalculator, RheologyCalculator.RheologyModels, ForwardDiff
import RheologyCalculator: compute_residual, _bt_line_search_result, bt_line_search

struct LineSearchProbe{F}
    f::F
    calls::Base.RefValue{Int}
end

@testset "Residual reuse preserves numeric and AD behavior" begin
    for T in (Float32, Float64)
        c = SeriesModel(LinearViscosity(T(5)))
        x = SVector(T(0.7))
        stress = rate -> solve(c, x, (; ε = rate), (;))[1]
        @test stress(T(0.1)) ≈ 10 * T(0.1)
        @test ForwardDiff.derivative(stress, T(0.1)) ≈ 10
        @test ForwardDiff.derivative(rate -> ForwardDiff.derivative(stress, rate), T(0.1)) ≈ 0 atol = 1.0e-12
    end
end

function compute_residual(probe::LineSearchProbe, x::SVector, vars, others)
    probe.calls[] += 1
    return SVector(probe.f(x[1]))
end

@testset "Line-search residual reuse" begin
    cases = (
        ("initial acceptance", identity, 1.0, 0.1, 1.0, 1),
        ("backtracked acceptance", x -> 2x, 1.0, 0.1, 0.25, 3),
        ("best trial is not last", x -> (x - 0.5)^2 + 2, 1.0, 0.1, 0.5, 4),
        ("nonfinite first trial", x -> x == 1 ? NaN : 0.5, 1.0, 0.1, 0.5, 2),
        ("no finite trial", x -> NaN, 1.0, 0.1, 1.0, 4),
        ("feasible step below minimum", identity, 0.05, 0.1, 0.05, 1),
    )
    for (name, f, initial, minimum_step, expected_step, calls) in cases
        @testset "$name" begin
            probe = LineSearchProbe(f, Ref(0))
            α, x, r, norm = _bt_line_search_result(
                SA[1.0], SA[0.0], probe, (;), (;), SA[1.0], 1.0;
                α = initial, α_min = minimum_step, lstol = 0.95
            )
            @test α == expected_step
            @test x == SA[expected_step]
            @test isequal(r, SA[f(expected_step)])
            @test isequal(norm, abs(f(expected_step)))
            @test probe.calls[] == calls
            probe.calls[] = 0
            @test bt_line_search(
                SA[1.0], SA[0.0], probe, (;), (;), SA[1.0], 1.0;
                α = initial, α_min = minimum_step, lstol = 0.95
            ) == expected_step
            @test probe.calls[] == (initial < minimum_step ? 0 : calls)
        end
    end
end
