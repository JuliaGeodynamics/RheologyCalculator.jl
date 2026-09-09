using RheologyCalculator, Test
using RheologyCalculator.RheologyModels

@testset "component partition reports series components" begin
    η, G, ε, dt = 1.0e19, 1.0e10, 1.0e-13, 1.0e10
    c = SeriesModel(
        LinearViscosity(η),
        IncompressibleElasticity(G),
        DruckerPrager(1.0e6, 30.0, 0.0),
    )
    vars, others = (; ε), (; dt, τ0 = (0.0,))
    x = initial_guess_x(c, vars, (; τ = 1.0e6), others)
    sol = solve(c, x, vars, others)
    p = component_partition(c, sol, vars, others)

    @test p isa ComponentPartition
    @test p.viscous_elements == (c.leafs[1],)
    @test p.elastic_elements == (c.leafs[2],)
    @test p.plastic_elements == (c.leafs[3],)
    @test p.viscous_indices == (1,)
    @test p.elastic_indices == (1,)
    @test p.plastic_indices == (1,)

    τ = sol[stress_index(c)]
    @test p.viscous_τ == (τ,)
    @test p.elastic_τ == (τ,)
    @test p.plastic_τ == (τ,)
    @test p.viscous_ε[1] + p.elastic_ε[1] + p.plastic_ε[1] ≈ ε
end

@testset "component partition reports parallel component stresses" begin
    η1, η2, ε, dt = 1.0e19, 1.0e20, 1.0e-13, 1.0e10
    c = SeriesModel(ParallelModel(LinearViscosity(η1), LinearViscosity(η2)))
    vars, others = (; ε), (; dt)
    x = initial_guess_x(c, vars, (; τ = 1.0e6), others)
    sol = solve(c, x, vars, others)
    p = component_partition(c, sol, vars, others)

    εb = sol[2]
    @test p.viscous_ε == (εb, εb)
    @test p.viscous_τ == (2η1 * εb, 2η2 * εb)
    @test p.viscous_elements == c.branches[1].leafs
    @test p.viscous_indices == (1, 2)
end
