# The cap's consistency residual must be dimensionally homogeneous.
using Test, StaticArrays
using RheologyCalculator.RheologyModels
import RheologyCalculator.RheologyModels: DruckerPragerCap, compute_F, compute_lambda
import RheologyCalculator: SeriesModel, initial_guess_x, normalisation_x, solve

@testset "cap consistency residual" begin

    @testset "residual is dimensionally homogeneous" begin
        # Scaling λ and η_vp by a common factor (a change of time unit) must
        # leave the residual unchanged. The unregularized below-yield branch is
        # excluded: it keeps a `oneunit` factor purely to pin λ at zero.
        for ηvp in (1.0e18, 1.0e19, 1.0e20)
            pl = DruckerPragerCap(; C = 1.0e6, ϕ = 30.0, ψ = 10.0, η_vp = ηvp, Pt = -5.0e5)
            for (τ, P) in ((8.0e5, 1.0e5), (2.0e6, -1.0e5), (1.0e5, -4.0e5))
                λ = 1.0e-11
                r = compute_lambda(pl; τ = τ, λ = λ, P = P)
                for s in (1.0e3, 1.0e-3)
                    pls = DruckerPragerCap(; C = 1.0e6, ϕ = 30.0, ψ = 10.0, η_vp = ηvp * s, Pt = -5.0e5)
                    rs = compute_lambda(pls; τ = τ, λ = λ / s, P = P)
                    @test rs ≈ r rtol = 1.0e-12
                end
            end
        end
    end

    @testset "the yielding branch is exactly Duvaut-Lions" begin
        for ηvp in (0.0, 1.0e19)
            pl = DruckerPragerCap(; C = 1.0e6, ϕ = 30.0, ψ = 10.0, η_vp = ηvp, Pt = -5.0e5)
            τ, P = 2.0e6, -1.0e5
            F = compute_F(pl, τ, P)
            @test F > -1.0e-8                            # this state is yielding
            for λ in (0.0, 1.0e-11, 1.0e-3)
                @test compute_lambda(pl; τ = τ, λ = λ, P = P) ≈ -F + λ * ηvp rtol = 1.0e-12
            end
        end
    end

    @testset "F = 0 at convergence with η_vp = 0" begin
        viscous = LinearViscosity(1.0e23)
        elastic = Elasticity(1.0e10, 2.0e11)
        pl = DruckerPragerCap(; C = 1.0e6, ϕ = 30.0, ψ = 10.0, η_vp = 0.0, Pt = -5.0e5)
        c = SeriesModel(viscous, elastic, pl)
        others = (; dt = 1.0e8, τ0 = (zero_stress_tensor_2D(),), P0 = (0.0,))

        for rate in (1.0e-11, 1.0e-10)
            vars = vars_2D(rate, 0.0)
            x0 = initial_guess_x(c, vars, (; τ = 0.0, P = 0.0, λ = 0.0), others)
            xn = normalisation_x(c, pl.C, second_invariant_2D(vars.ε) + abs(vars.θ))
            sol = solve(c, x0, vars, others; xnorm0 = xn, itermax = 200)
            τ, λ, P = sol.x[1], sol.x[2], sol.x[3]
            @test λ > 0                                  # this state really is yielding
            # F vanishes on the scale of the cohesion, not of λ
            @test abs(compute_F(pl, τ, P)) / pl.C ≤ 1.0e-12
        end
    end

    @testset "below yield the multiplier is driven to zero" begin
        pl = DruckerPragerCap(; C = 1.0e6, ϕ = 30.0, ψ = 10.0, η_vp = 0.0, Pt = -5.0e5)
        τ, P = 1.0e5, 1.0e5
        @test compute_F(pl, τ, P) < -1.0e-8              # safely below yield
        @test compute_lambda(pl; τ = τ, λ = 0.0, P = P) == 0.0
        @test compute_lambda(pl; τ = τ, λ = 1.0e-11, P = P) > 0
        @test compute_lambda(pl; τ = τ, λ = -1.0e-11, P = P) < 0
    end
end
