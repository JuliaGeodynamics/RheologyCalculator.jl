using RheologyCalculator, Test
using RheologyCalculator.RheologyModels
# core-internal helpers used across the tests
import RheologyCalculator: compute_stress_elastic, compute_pressure_elastic, compute_residual,
    mynorm, _direct_leaf_elastic_correction, second_invariant_value
# tensor helpers live in the RheologyModels submodule and are not exported
import RheologyCalculator.RheologyModels: second_invariant_2D, tensor_strain_rate_2D,
    vars_2D, zero_stress_tensor_2D, stress_tensor_from_invariant_2D, elastic_stress_history_2D

function runtests()
    files = readdir(@__DIR__)
    test_files = filter(f -> startswith(f, "test_") && endswith(f, ".jl"), files)

    # Explicit `test_*.jl` names in ARGS select individual files, so that
    # `Pkg.test(; test_args=["test_foo.jl"])` runs only those. Going through
    # `Pkg.test` rather than including a file directly is what materializes the
    # weak test dependencies (SparseConnectivityTracer), which are absent under
    # the bare package project.
    requested = filter(arg -> startswith(arg, "test_") && endswith(arg, ".jl"), ARGS)
    if !isempty(requested)
        unknown = setdiff(requested, test_files)
        isempty(unknown) || throw(ArgumentError("unknown test file(s): $(join(unknown, ", "))"))
        test_files = filter(in(requested), test_files)
    elseif "--allocations-only" in ARGS
        test_files = filter(==("test_allocations.jl"), test_files)
    elseif Base.JLOptions().code_coverage != 0
        filter!(!=("test_allocations.jl"), test_files)
    end

    for f in test_files
        !isdir(f) && include(f)
    end
    return
end

runtests()
