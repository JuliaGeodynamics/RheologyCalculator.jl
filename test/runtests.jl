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
    test_files = filter(startswith("test_"), files)

    allocations_only = "--allocations-only" in ARGS
    if allocations_only
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
