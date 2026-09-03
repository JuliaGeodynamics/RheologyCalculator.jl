using Test
using RheologyCalculator
using RheologyCalculator.RheologyModels
using StaticArrays
using ForwardDiff
using SparseConnectivityTracer

import RheologyCalculator: initial_guess_x, normalisation_x, solve, SeriesModel,
    safe_inv, safe_inv_one
import RheologyCalculator.RheologyModels: DruckerPragerCap

# The tracer guards replaced three inline value-dependent expressions. Pin the
# Float64 semantics of each against the expression it replaced, including the
# zero cases that motivated the guards in the first place.
@testset "tracer guards preserve Float64 semantics" begin
    for v in (4.0, -2.5, 1.0e-30)
        @test safe_inv(v)     === (iszero(v) ? zero(v) : inv(v))
        @test safe_inv_one(v) === (iszero(v) ? one(v)  : inv(v))
    end
    @test safe_inv(0.0)     === 0.0    # guard: no Inf
    @test safe_inv_one(0.0) === 1.0    # guard: neutral element, no Inf
    @test safe_inv(4.0)     === 0.25
    @test safe_inv_one(4.0) === 0.25
end

@testset "SparseConnectivityTracer extension" begin
    @test Base.get_extension(
        RheologyCalculator, :RheologyCalculatorSparseConnectivityTracerExt,
    ) !== nothing

    viscous = LinearViscosity(1.0e23)
    elastic = Elasticity(1.0e10, 2.0e11)
    plastic = DruckerPragerCap(; C = 1.0e6, ϕ = 30.0, ψ = 10.0, η_vp = 1.0e18, Pt = -5.0e5)
    c = SeriesModel(viscous, elastic, plastic)

    others = (; dt = 1.0e5, τ0 = (zero_stress_tensor_2D(),), P0 = (0.3e6,))
    args = (; τ = 0.0e3, P = 0.3e6, λ = 0)
    v0 = [7.0e-14, 7.0e-15]

    function solve_outputs(input)
        strain_rate, volumetric_strain_rate = input
        vars = vars_2D(strain_rate, volumetric_strain_rate)
        x = initial_guess_x(c, vars, args, others)
        return solve(c, x, vars, others;
            xnorm0 = normalisation_x(c, plastic.C, 7.0e-14),
        )
    end

    @testset "global sparsity detection succeeds" begin
        # Without the extension this throws
        # "Function iszero requires primal value(s)".
        S = SparseConnectivityTracer.jacobian_sparsity(
            solve_outputs, v0, TracerSparsityDetector(),
        )
        @test size(S) == (length(normalisation_x(c, plastic.C, 7.0e-14)), length(v0))

        # The solve short-circuits under a tracer
        # each output must depend on every traced input.
        @test all(S)

        # The essential correctness property of any sparsity pattern: it must
        # be a superset of the true nonzeros, never miss one.
        J = ForwardDiff.jacobian(solve_outputs, v0)
        @test all(S[J .!= 0])
    end

    @testset "short-circuit returns one entry per unknown" begin
        n = length(normalisation_x(c, plastic.C, 7.0e-14))
        S = SparseConnectivityTracer.jacobian_sparsity(
            solve_outputs, v0, TracerSparsityDetector(),
        )
        @test size(S) == (n, length(v0))
        @test all(S)
    end

    @testset "traced auxiliary values survive local extraction" begin
        function extract_temperature(temperature_input)
            others = (; T = temperature_input[1], dt = 1.0,
                τ0 = ((0.0, 0.0, 0.0),))
            local_others = RheologyCalculator.extract_local_kwargs(others, (:τ0,), 1)
            return local_others.T
        end

        initial_temperature = [1.0]
        detector = TracerSparsityDetector()
        sparsity = SparseConnectivityTracer.jacobian_sparsity(
            extract_temperature, initial_temperature, detector,
        )

        # The local temperature remains connected to its single input even
        # though `others` also contains Float64 and history-tuple fields.
        expected_sparsity = trues(1, 1)
        @test sparsity == expected_sparsity
    end

    @testset "short-circuit includes auxiliary dependencies" begin
        n = length(normalisation_x(c, plastic.C, 7.0e-14))

        function solve_with_temperature(inputs)
            stress_seed = inputs[1]
            temperature = inputs[2]

            # `x` must be tracer-valued to dispatch to the tracer-only solve.
            # It depends only on `stress_seed`; temperature is supplied only
            # through `others`.
            x = SVector{n}(fill(stress_seed, n))
            vars = (; ε = 0.0)
            others = (; T = temperature)
            return solve(c, x, vars, others)
        end

        initial_inputs = [1.0, 1.0]
        detector = TracerSparsityDetector()
        sparsity = SparseConnectivityTracer.jacobian_sparsity(
            solve_with_temperature, initial_inputs, detector,
        )

        # Every local solution component conservatively depends on both the
        # traced solver seed and the temperature stored in `others`.
        expected_sparsity = trues(n, 2)
        @test sparsity == expected_sparsity
    end
end
