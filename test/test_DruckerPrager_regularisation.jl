import RheologyCalculator: compute_lambda, compute_lambda_parallel

@testset "Drucker-Prager regularisation" begin
    default = DruckerPrager(10.0, 30.0, 0.0)
    regularised = DruckerPrager(10.0, 30.0, 0.0, 25.0)

    @test default.η_vp == 1.0
    @test regularised.η_vp == 25.0
    @test DruckerPrager(10.0f0, 30.0f0, 0.0f0) isa DruckerPrager{Float32}

    τ = 20.0
    P = 0.0
    λ = 0.25
    F = τ - regularised.C * regularised.cosϕ

    @test compute_lambda(regularised; τ, P, λ) == F - λ * regularised.η_vp
    @test compute_lambda_parallel(regularised; τ_pl = τ, P, λ) == F - λ * regularised.η_vp
    @test compute_lambda(DruckerPrager(10.0, 0.0, 0.0, 0.0); τ, P, λ) == 10.0

    mixed = DruckerPrager(10, 30.0f0, 0, 25.0)
    @test mixed isa DruckerPrager{Float64}
    @test @inferred(compute_lambda(regularised; τ, P, λ)) == F - λ * regularised.η_vp
end
