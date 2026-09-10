using Test
using RheologyCalculator
using RheologyCalculator.RheologyModels
import RheologyCalculator: compute_strain_rate, compute_stress

const _R_GK = 8.314

@testset "Goldsby–Kohlstedt mechanisms" begin
    r = GoldsbyKohlstedtCreep(1.7, 1.4, 1.0e-10, 49.0e3, _R_GK)
    ε = compute_strain_rate(r; τ = 1.0e5, T = 250.0, d = 1.0e-3)
    @test compute_stress(r; ε = ε, T = 250.0, d = 1.0e-3) ≈ 1.0e5

    d = GoldsbyKohlstedtDiffusion(_R_GK, 9.1e-4, 59.4e3, 1.0e-4, 59.4e3,
        1.97e-5, 1.0e-9)
    @test isfinite(compute_strain_rate(d; τ = 1.0e5, T = 250.0, d = 1.0e-3))

    c = SeriesModel(ParallelModel(d, DislocationCreep(4, 0, 1.0e-20, 64.0e3, 0.0, _R_GK), r))
    @test length(x_keys(c)) == 2
end
