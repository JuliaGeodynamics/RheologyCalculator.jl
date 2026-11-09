# Each fixture is (name, composite, vars, args, others). The topologies are the
# ones whose equation generation is unrolled per composite type: losing
# inference on any of them turns a static residual evaluation into a dynamic
# one without changing a single number.
function type_stability_fixtures()
    viscous = LinearViscosity(5.0e19)
    viscous2 = LinearViscosity(1.0e20)
    elastic = Elasticity(1.0e10, 1.0e12)
    elasticinc = IncompressibleElasticity(1.0e10)
    bulkvisc = BulkViscosity(1.0e18)

    return (
        (
            "series elastic-viscous",
            SeriesModel(elastic, viscous),
            (; ε = 1.0e-15, θ = 1.0e-20),
            (; τ = 1.0e2, P = 1.0e6),
            (; dt = 1.0e10, τ0 = (0.0,), P0 = (0.0,)),
        ),
        (
            "Kelvin-Voigt branch",
            SeriesModel(elasticinc, ParallelModel(elasticinc, viscous)),
            (; ε = 1.0e-15),
            (; τ = 1.0e2),
            (; dt = 1.0e10, τ0 = (0.0, 0.0), P0 = (0.0, 0.0)),
        ),
        (
            "generalized Maxwell branch",
            SeriesModel(viscous, ParallelModel(SeriesModel(viscous2, elasticinc), viscous2)),
            (; ε = 1.0e-15),
            (; τ = 1.0e2),
            (; dt = 1.0e10, τ0 = (0.0,), P0 = (0.0,)),
        ),
        (
            "plastic series",
            SeriesModel(LinearViscosity(1.0e22), elasticinc, DruckerPrager(1.0e6, 30.0, 0.0)),
            (; ε = 1.0e-14),
            (; τ = 1.0e6),
            (; dt = 1.0e10, τ0 = (0.0,), P = 1.0e6, P0 = (0.0,)),
        ),
        (
            "volumetric series",
            SeriesModel(elastic, bulkvisc),
            (; ε = 1.0e-15, θ = 1.0e-20),
            (; τ = 1.0e2, P = 1.0e6),
            (; dt = 1.0e10, τ0 = (0.0,), P0 = (0.0,)),
        ),
    )
end
