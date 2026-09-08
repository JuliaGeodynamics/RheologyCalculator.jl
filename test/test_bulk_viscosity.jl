import RheologyCalculator: compute_pressure, compute_volumetric_strain_rate

@testset "Bulk viscosity sign convention" begin
    bulk = BulkViscosity(4.0)

    @test compute_volumetric_strain_rate(bulk; P = 12.0) == -3.0
    @test compute_pressure(bulk; θ = 3.0) == -12.0
end