"""
    ComponentPartition

Decomposition of the deviatoric stress and strain rate of every viscous,
elastic, and plastic component from a converged local solution.

For each category, the `ε`, `τ`, `elements`, and `indices` tuples are index-aligned. In a
series equation, components share stress and have individual strain rates. In
parallel equation, they share strain rate and have individual stresses.
"""
struct ComponentPartition{EV, TV, RV, IV, EE, TE, RE, IE, EP, TP, RP, IP}
    viscous_ε::EV
    viscous_τ::TV
    viscous_elements::RV
    viscous_indices::IV
    elastic_ε::EE
    elastic_τ::TE
    elastic_elements::RE
    elastic_indices::IE
    plastic_ε::EP
    plastic_τ::TP
    plastic_elements::RP
    plastic_indices::IP
end

# `_show_round` is defined in dissipation_partition.jl, which is included first.

function Base.show(io::IO, ::MIME"text/plain", p::ComponentPartition)
    println(io, "ComponentPartition:")
    println(io, "  viscous  ε = ", _show_round(p.viscous_ε), "  τ = ", _show_round(p.viscous_τ), "  indices = ", p.viscous_indices)
    println(io, "  elastic  ε = ", _show_round(p.elastic_ε), "  τ = ", _show_round(p.elastic_τ), "  indices = ", p.elastic_indices)
    println(io, "  plastic  ε = ", _show_round(p.plastic_ε), "  τ = ", _show_round(p.plastic_τ), "  indices = ", p.plastic_indices)
    return nothing
end

"""
    component_partition(c::AbstractCompositeModel, sol::AbstractVector, vars, others)

Return the per-component deviatoric strain-rate and stress contributions of a
converged local solution. The `vars` argument is accepted to preserve the
package's conventional post-processing calling sequence; the component values are
reconstructed from `sol` and `others`.
"""
function component_partition(c::AbstractCompositeModel, sol::AbstractVector, _vars, others)
    eqs = generate_equations(c)
    n = length(eqs)
    x = SVector{n}(ntuple(i -> @inbounds(sol[i]), n))
    args = generate_args_template(eqs, x, others)

    viscous = _component_contributions(eqs, args, others, Val(:viscous))
    elastic = _component_contributions(eqs, args, others, Val(:elastic))
    plastic = _component_contributions(eqs, args, others, Val(:plastic))
    return ComponentPartition(viscous..., elastic..., plastic...)
end

@generated function _component_contributions(
        eqs::NTuple{N, CompositeEquation}, args, others, ::Val{category},
    ) where {N, category}
    N == 0 && return :(((), (), (), ()))
    return quote
        @inline
        ε_0 = (); τ_0 = (); elements_0 = (); indices_0 = ()
        Base.@nexprs $N i -> begin
            ε_i, τ_i, elements_i, indices_i = _component_equation(
                eqs[i].fn,
                eqs[i],
                args[i],
                others,
                Val(category),
                ε_{i - 1},
                τ_{i - 1},
                elements_{i - 1},
                indices_{i - 1},
            )
        end
        ($(Symbol(:ε_, N)), $(Symbol(:τ_, N)), $(Symbol(:elements_, N)), $(Symbol(:indices_, N)))
    end
end

@inline _component_equation(
        ::F, _eq, _args, _others, _category, ε, τ, elements, indices,
    ) where {F} = (ε, τ, elements, indices)

@inline function _component_equation(
        fn::Union{typeof(compute_strain_rate), typeof(compute_stress)},
        eq,
        args,
        others,
        category,
        ε,
        τ,
        elements,
        indices,
    )
    values = evaluate_state_function_perleaf(fn, eq.rheology, args, others, eq.el_number)
    return _component_leaves(
        eq.rheology, values, eq.el_number, fn, args, category, ε, τ, elements, indices,
    )
end

@generated function _component_leaves(
        rheology::NTuple{N, AbstractRheology}, values, element_numbers, fn, args, category,
        ε, τ, elements, indices,
    ) where {N}
    N == 0 && return :((ε, τ, elements, indices))
    return quote
        @inline
        ε_0 = ε; τ_0 = τ; elements_0 = elements; indices_0 = indices
        Base.@nexprs $N i -> begin
            ε_i, τ_i, elements_i, indices_i = _component_step(
                rheology_category(rheology[i]),
                category,
                rheology[i],
                element_numbers[i],
                fn,
                values[i],
                args,
                ε_{i - 1},
                τ_{i - 1},
                elements_{i - 1},
                indices_{i - 1},
            )
        end
        ($(Symbol(:ε_, N)), $(Symbol(:τ_, N)), $(Symbol(:elements_, N)), $(Symbol(:indices_, N)))
    end
end

@inline _component_pair(::typeof(compute_strain_rate), value, args) = (value, args.τ)
@inline _component_pair(::typeof(compute_stress), value, args) = (args.ε, value)

@inline function _component_step(
        ::Val{c}, ::Val{c}, rheology, index, fn, value, args, ε, τ, elements, indices,
    ) where {c}
    ε_component, τ_component = _component_pair(fn, value, args)
    return (
        (ε..., ε_component),
        (τ..., τ_component),
        (elements..., rheology),
        (indices..., index),
    )
end

@inline _component_step(
        ::Val, ::Val, _rheology, _index, _fn, _value, _args, ε, τ, elements, indices,
    ) = (ε, τ, elements, indices)
