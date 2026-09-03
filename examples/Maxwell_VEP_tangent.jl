using RheologyCalculator
using RheologyCalculator.RheologyModels, StaticArrays, ForwardDiff, Chairmarks
using RheologyCalculator.RheologyModels: second_invariant_2D, tensor_strain_rate_2D, zero_stress_tensor_2D, elastic_stress_history_2D
import ForwardDiff: ForwardDiff



@inline function compute_stress_tensor(ε::SVector{3, T}, c) where T
    θ = 0e0
    # input variables (constant)
    εᵢⱼ    = Tuple(ε)
    τ0ᵢⱼ   = (zero_stress_tensor_2D(),)
    vars   = (; ε = εᵢⱼ, θ)
    # guess variables (we solve for these, differentiable)
    args   = (; τ = 0e0, λ = 0)
    # other non-differentiable variables needed to evaluate the state functions
    others = (; dt = 1.0e8, P = 1.0e6, τ0 = τ0ᵢⱼ, P0 = (0.0,))

    x       = initial_guess_x(c, vars, args, others)
    char_τ  = c.leafs[3].C
    char_ε  = second_invariant_2D(vars.ε)
    xnorm   = normalisation_x(c, char_τ, char_ε)

    τII = solve(c, x, vars, others, verbose=false, xnorm0=xnorm)[1]
    τᵢⱼ = elastic_stress_history_2D(c, τII, vars.ε, τ0ᵢⱼ, others)[1]
    return SVector{3}(τᵢⱼ)
end

@inline function compute_stress_tensor(ε::SVector{3, T}, c, index::Int) where T
    # input variables (constant)
    εᵢⱼ    = Tuple(ε)
    τ0ᵢⱼ   = (zero_stress_tensor_2D(),)
    vars   = (; ε = εᵢⱼ, θ = 0e0)
    # guess variables (we solve for these, differentiable)
    args   = (; τ = 0e0, λ = 0)
    # other non-differentiable variables needed to evaluate the state functions
    others = (; dt = 1.0e8, P = 1.0e6, τ0 = τ0ᵢⱼ, P0 = (0.0,))

    x       = initial_guess_x(c, vars, args, others)
    char_τ  = c.leafs[3].C
    char_ε  = second_invariant_2D(vars.ε)
    xnorm   = normalisation_x(c, char_τ, char_ε)

    τII = solve(c, x, vars, others, verbose=false, xnorm0=xnorm)[1]
    τᵢⱼ = elastic_stress_history_2D(c, τII, vars.ε, τ0ᵢⱼ, others)[1][index]
end

@inline function ∇σij(εxx, εyy, εxy, c, index)
    ∂σij∂εxx = ForwardDiff.derivative(εxx -> compute_stress_tensor(SA[εxx, εyy, εxy], c, index), εxx)
    ∂σij∂εyy = ForwardDiff.derivative(εyy -> compute_stress_tensor(SA[εxx, εyy, εxy], c, index), εyy)
    ∂σij∂εxy = ForwardDiff.derivative(εxy -> compute_stress_tensor(SA[εxx, εyy, εxy], c, index), εxy)
    return ∂σij∂εxx, ∂σij∂εyy, ∂σij∂εxy
end
 
@inline tangent_operator(ε, c) = tangent_operator(ε..., c)

@inline function tangent_operator(εxx, εyy, εxy, c)
    ∂σxx∂εxx, ∂σxx∂εyy, ∂σxx∂εxy = ∇σij(εxx, εyy, εxy, c, 1)
    ∂σyy∂εxx, ∂σyy∂εyy, ∂σyy∂εxy = ∇σij(εxx, εyy, εxy, c, 2)
    ∂σxy∂εxx, ∂σxy∂εyy, ∂σxy∂εxy = ∇σij(εxx, εyy, εxy, c, 3)
    return SA[
        ∂σxx∂εxx  ∂σxx∂εyy  ∂σxx∂εxy
        ∂σyy∂εxx  ∂σyy∂εyy  ∂σyy∂εxy
        ∂σxy∂εxx  ∂σxy∂εyy  ∂σxy∂εxy
    ]
end

@inline tangent_operator_diagonal(ε, c) = tangent_operator_diagonal(ε..., c)

@inline function tangent_operator_diagonal(εxx, εyy, εxy, c)

    ∂σxx∂εxx = ForwardDiff.derivative(εxx -> compute_stress_tensor(SA[εxx, εyy, εxy], c, 1), εxx)
    ∂σyy∂εyy = ForwardDiff.derivative(εyy -> compute_stress_tensor(SA[εxx, εyy, εxy], c, 2), εyy)
    ∂σxy∂εxy = ForwardDiff.derivative(εxy -> compute_stress_tensor(SA[εxx, εyy, εxy], c, 3), εxy)

    return SA[
        ∂σxx∂εxx  0e0       0e0
        0e0       ∂σyy∂εyy  0e0
        0e0       0e0       ∂σxy∂εxy
    ]
end

viscous = LinearViscosity(1e22)
elastic = IncompressibleElasticity(10e9)
plastic = DruckerPrager(1e6, 30, 0)

# Maxwell visco-elasto-plastic model
# elastic --- viscous --- plastic
c  = SeriesModel(viscous, elastic, plastic)
ε =  SA[tensor_strain_rate_2D(1.0e-14)...]

@assert tangent_operator(ε, c) == tangent_operator_diagonal(ε, c)

@b tangent_operator($ε, $c)

@b tangent_operator_diagonal($ε, $c)
