# Dilatancy couples the plastic multiplier to the volumetric balance: a dilatant
# Drucker-Prager element expands as it shears, and the surrounding volumetric
# elements must absorb that expansion.

import RheologyCalculator: series_state_functions, parallel_state_functions,
    compute_volumetric_strain_rate, compute_volumetric_plastic_strain_rate, _isvolumetric

const C_dp, ϕ_dp, dt_dp = 2.0, 30.0, 0.1
const G_dp, K_dp, ε̇_dp = 1.0e3, 1.0e3, 1.0e-2

dilatancy_others() = (; dt = dt_dp, τ0 = ((0.0, 0.0, 0.0),), P0 = (0.0,))
dilatancy_vars() = (; ε = (ε̇_dp, -ε̇_dp, 0.0), θ = 0.0)

function solve_dilatant(c)
    others = dilatancy_others()
    x0 = initial_guess_x(c, (; ε = ε̇_dp, θ = 0.0), (; τ = 0.0), others)
    sol = solve(c, x0, dilatancy_vars(), others; xnorm0 = normalisation_x(c, 1.0e6, ε̇_dp), atol = 1.0e-12, rtol = 1.0e-12)
    return NamedTuple{x_keys(c)}(Tuple(sol.x))
end

series_dilatant(ψ) = SeriesModel(LinearViscosity(2.5e3), Elasticity(G_dp, K_dp), DruckerPrager(C_dp, ϕ_dp, ψ))
branch_dilatant(ψ) = SeriesModel(
    LinearViscosity(2.5e3), Elasticity(G_dp, K_dp),
    ParallelModel(LinearViscosity(25.0), DruckerPrager(C_dp, ϕ_dp, ψ))
)

@testset "dilatancy is a type-level property" begin
    plain, dilatant = DruckerPrager(C_dp, ϕ_dp, 0.0), DruckerPrager(C_dp, ϕ_dp, 5.0)
    @test !_isvolumetric(plain)
    @test _isvolumetric(dilatant)

    @test compute_volumetric_strain_rate ∉ series_state_functions(plain)
    @test compute_volumetric_strain_rate ∈ series_state_functions(dilatant)
    @test compute_volumetric_plastic_strain_rate ∉ parallel_state_functions(plain)
    @test compute_volumetric_plastic_strain_rate ∈ parallel_state_functions(dilatant)

    # A non-dilatant element must not pay for the unknowns it does not need: the
    # pressure here is the compressible elasticity's, present either way.
    @test x_keys(series_dilatant(0.0)) == (:τ, :λ, :P)
    @test x_keys(series_dilatant(30.0)) == (:τ, :λ, :P)
    @test x_keys(branch_dilatant(0.0)) == (:τ, :ε, :λ, :τ_pl, :P)
    @test x_keys(branch_dilatant(30.0)) == (:τ, :ε, :λ, :τ_pl, :P, :θ, :P_pl)
end

@testset "volumetric flow rule" begin
    dp = DruckerPrager(C_dp, ϕ_dp, 30.0)
    # θ_pl = 2 sinψ ε̇_II_pl and ε̇_II_pl = λ/2, so the multiplier enters undivided.
    @test compute_volumetric_strain_rate(dp; λ = 2.0) ≈ 2.0 * sind(30.0)
    @test compute_volumetric_plastic_strain_rate(dp; λ = 2.0) ≈ 2.0 * sind(30.0)
    # Dilation is positive, so a dilatant element expands.
    @test compute_volumetric_strain_rate(dp; λ = 2.0) > 0
    @test iszero(compute_volumetric_strain_rate(DruckerPrager(C_dp, ϕ_dp, 0.0); λ = 2.0))
end

@testset "dilatancy enters the solve" begin
    for build in (series_dilatant, branch_dilatant)
        s10 = solve_dilatant(build(10.0))
        s30 = solve_dilatant(build(30.0))
        @test s10.τ != s30.τ

        for (ψ, s) in ((10.0, s10), (30.0, s30))
            θ_pl = s.λ * sind(ψ)
            # Mass conservation: the elastic volumetric strain rate absorbs the
            # plastic dilation, since the prescribed total θ is zero.
            @test -s.P / (K_dp * dt_dp) + θ_pl ≈ 0 atol = 1.0e-14
            # Yield is satisfied at the pressure the volumetric balance produced.
            P_yield = haskey(s, :P_pl) ? s.P_pl : s.P
            τ_yield = haskey(s, :τ_pl) ? s.τ_pl : s.τ
            @test τ_yield - P_yield * sind(ϕ_dp) - C_dp * cosd(ϕ_dp) ≈ s.λ atol = 1.0e-12
        end
    end
end

@testset "dilatancy hardens" begin
    # Expansion against a compressible surrounding raises the pressure, which
    # raises the Drucker-Prager yield stress.
    τ = map(ψ -> solve_dilatant(series_dilatant(ψ)).τ, (0.0, 10.0, 20.0, 30.0))
    P = map(ψ -> get(solve_dilatant(series_dilatant(ψ)), :P, 0.0), (0.0, 10.0, 20.0, 30.0))
    @test issorted(τ)
    @test issorted(P)
    @test allunique(τ)
end

@testset "the branch carries the plastic dilation" begin
    s = solve_dilatant(branch_dilatant(30.0))
    # The branch volumetric strain rate is the plastic one: the parallel
    # viscosity has no volumetric response.
    @test s.θ ≈ s.λ * sind(30.0)
    # Nothing else in the branch carries volumetric stress, so the element
    # pressure is the branch pressure.
    @test s.P_pl ≈ s.P
end

@testset "zero dilatancy leaves the volumetric balance untouched" begin
    for build in (series_dilatant, branch_dilatant)
        s0 = solve_dilatant(build(0.0))
        @test s0.λ > 0
        # No dilation to absorb, so the prescribed θ = 0 costs no pressure.
        @test s0.P ≈ 0 atol = 1.0e-14
    end
end

@testset "the solved dilation is readable" begin
    others = dilatancy_others()
    for build in (series_dilatant, branch_dilatant), ψ in (0.0, 10.0, 30.0)
        c = build(ψ)
        x0 = initial_guess_x(c, (; ε = ε̇_dp, θ = 0.0), (; τ = 0.0), others)
        sol = solve(c, x0, dilatancy_vars(), others; xnorm0 = normalisation_x(c, 1.0e6, ε̇_dp), atol = 1.0e-12, rtol = 1.0e-12)
        s = NamedTuple{x_keys(c)}(Tuple(sol.x))
        @test volumetric_plastic_strain_rate(c, sol, others) ≈ s.λ * sind(ψ) atol = 1.0e-16
    end

    # A composite with no plastic element has no plastic dilation to report.
    c = SeriesModel(LinearViscosity(2.5e3), Elasticity(G_dp, K_dp))
    x0 = initial_guess_x(c, (; ε = ε̇_dp, θ = 0.0), (; τ = 0.0), others)
    sol = solve(c, x0, dilatancy_vars(), others)
    @test iszero(volumetric_plastic_strain_rate(c, sol, others))
end
