using Test

@testset "allocations" begin
    viscous1 = LinearViscosity(5.0e19)
    viscous2 = LinearViscosity(1.0e20)
    elastic = Elasticity(1.0e10, 1.0e12)
    elasticinc = IncompressibleElasticity(1.0e10)
    bulkvisc = BulkViscosity(1.0e18)

    # flat composite: no nested branches at all
    c_flat = SeriesModel(elastic, viscous1)
    vars_flat = (; ε = 1.0e-15, θ = 1.0e-20)
    others_flat = (; dt = 1.0e10, τ0 = (0.0,), P0 = (0.0,))
    x_flat = initial_guess_x(c_flat, vars_flat, (; τ = 1.0e2, P = 1.0e6), others_flat)

    # single branch nested one level deep (pre-existing regression coverage)
    s1 = SeriesModel(viscous1, viscous2)
    p = ParallelModel(s1, viscous2)
    c_single_branch = SeriesModel(viscous1, p)
    vars_sb = (; ε = 1.0e-15)
    others_sb = (;)
    x_sb = initial_guess_x(c_single_branch, vars_sb, (; τ = 1.0e2), others_sb)

    # sibling branches where a non-last branch nests further: the topology
    # that exposed the equation-indexing bug in generate_offsets_parallel
    p1 = ParallelModel(viscous1, viscous2)
    c_multi_branch = SeriesModel(viscous1, p, p1)
    vars_mb = (; ε = 1.0e-15)
    others_mb = (;)
    x_mb = initial_guess_x(c_multi_branch, vars_mb, (; τ = 1.0e2), others_mb)

    # Kelvin-Voigt branch: an elastic leaf in parallel with a viscous one, so
    # the branch carries elastic history and reaches _η_KV
    c_kv = SeriesModel(elasticinc, ParallelModel(elasticinc, viscous1))
    vars_kv = (; ε = 1.0e-15)
    others_kv = (; dt = 1.0e10, τ0 = (0.0, 0.0), P0 = (0.0, 0.0))
    x_kv = initial_guess_x(c_kv, vars_kv, (; τ = 1.0e2), others_kv)

    # Generalized Maxwell branch: a viscoelastic sub-series inside a parallel
    # branch, the only topology that exercises _η_eff_maxwell / _η_eff_elastic
    # and the η_star = η_eff_M/η_el backstress weighting
    c_maxwell = SeriesModel(viscous1, ParallelModel(SeriesModel(viscous2, elasticinc), viscous2))
    vars_gm = (; ε = 1.0e-15)
    others_gm = (; dt = 1.0e10, τ0 = (0.0,), P0 = (0.0,))
    x_gm = initial_guess_x(c_maxwell, vars_gm, (; τ = 1.0e2), others_gm)

    # Plastic composite, loaded past yield so compute_lambda returns a nonzero
    # multiplier and the lambda paths in add_child / subtract_parent are live
    c_plastic = SeriesModel(LinearViscosity(1.0e22), elasticinc, DruckerPrager(1.0e6, 30.0, 0.0))
    vars_pl = (; ε = 1.0e-14)
    others_pl = (; dt = 1.0e10, τ0 = (0.0,), P = 1.0e6, P0 = (0.0,))
    x_pl = initial_guess_x(c_plastic, vars_pl, (; τ = 1.0e6), others_pl)

    # Volumetric composite: exercises the compute_pressure /
    # compute_volumetric_strain_rate equation elimination path
    c_volumetric = SeriesModel(elastic, bulkvisc)
    vars_vol = (; ε = 1.0e-15, θ = 1.0e-20)
    others_vol = (; dt = 1.0e10, τ0 = (0.0,), P0 = (0.0,))
    x_vol = initial_guess_x(c_volumetric, vars_vol, (; τ = 1.0e2, P = 1.0e6), others_vol)

    cases = (
        (c_flat, vars_flat, others_flat, x_flat),
        (c_single_branch, vars_sb, others_sb, x_sb),
        (c_multi_branch, vars_mb, others_mb, x_mb),
        (c_kv, vars_kv, others_kv, x_kv),
        (c_maxwell, vars_gm, others_gm, x_gm),
        (c_plastic, vars_pl, others_pl, x_pl),
        (c_volumetric, vars_vol, others_vol, x_vol),
    )

    @testset "generate_equations" begin
        for (c, _, _, _) in cases
            generate_equations(c)
            @test (@allocated generate_equations(c)) == 0
        end
    end

    @testset "compute_residual" begin
        for (c, vars, others, x) in cases
            compute_residual(c, x, vars, others)
            @test (@allocated compute_residual(c, x, vars, others)) == 0
        end
    end

    @testset "solve" begin
        for (c, vars, others, x) in cases
            solve(c, x, vars, others)
            @test (@allocated solve(c, x, vars, others)) == 0
        end
    end

    @testset "tangent" begin
        for (c, vars, others, x) in cases
            xc = solve(c, x, vars, others)
            tangent(c, xc, vars, others)
            @test (@allocated tangent(c, xc, vars, others)) == 0
        end
    end
end
