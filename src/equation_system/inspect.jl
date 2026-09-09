"""
    ModelInspection

Description of the unknowns of a composite model, one entry per entry of the
solver vector `x`, as returned by [`inspect`](@ref).

It is an `AbstractVector` of `NamedTuple`s, so entries can be indexed, iterated,
and searched by position; it displays as a table.
"""
struct ModelInspection{T, N, E <: NTuple{N, Any}} <: AbstractVector{T}
    entries::E
end

ModelInspection(entries::Tuple) = ModelInspection{eltype(entries), length(entries), typeof(entries)}(entries)

Base.size(::ModelInspection{T, N}) where {T, N} = (N,)
Base.IndexStyle(::Type{<:ModelInspection}) = IndexLinear()
Base.@propagate_inbounds Base.getindex(insp::ModelInspection, i::Int) = insp.entries[i]

"""
    inspect(c::AbstractCompositeModel)

Describe the unknowns of composite model `c`, one entry per entry of the solver
vector `x` and hence of the [`RCSolution`](@ref) returned by [`solve`](@ref).

Each entry is a `NamedTuple` with fields

- `var`: the name of the unknown, as in [`x_keys`](@ref);
- `equation`: the state function whose residual the entry solves;
- `isglobal`: whether that equation belongs to the outermost series model
  rather than to a parallel branch;
- `elements`: the rheology elements the equation spans, each paired with its
  number in the per-type numbering that `display(c)` draws.
- `parent` and `children`: equation-graph dependencies.
- `inputs`: prescribed differentiable input fields used by the equation.
- `history`: auxiliary/history fields requested by its rheology elements.

A name repeats when several equations share the same physical unknown: a
composite with a parallel branch has one `:τ` per branch, of which the global
one is the stress of the model as a whole, and a Kelvin chain has one `:ε` per
branch. The `equation` and `elements` fields are what tell those apart.

# Example

```julia
julia> c = SeriesModel(
           LinearViscosity(1.0e18),
           ParallelModel(
               SeriesModel(LinearViscosity(1.0e19), IncompressibleElasticity(1.0e10)),
               LinearViscosity(1.0e18),
           ),
       );

julia> inspect(c)
3-element ModelInspection:
  index  var  equation             scope   elements
      1  τ    compute_strain_rate  global  LinearViscosity 1
      2  ε    compute_stress       branch  LinearViscosity 2
      3  τ    compute_strain_rate  branch  LinearViscosity 3, IncompressibleElasticity 1
```
"""
inspect(c::AbstractCompositeModel) = ModelInspection(inspection_entries(generate_equations(c)))

@generated function inspection_entries(eqs::NTuple{N, CompositeEquation}) where {N}
    return quote
        @inline
        e = Base.@ntuple $N i -> _inspection_entry(eqs[i])
        superflatten(e)
    end
end

# One entry per key of the equation, so the description keeps step with
# `x_keys` even for an equation that contributes more than one unknown.
@inline function _inspection_entry(eq::CompositeEquation)
    ks = keys(differentiable_kwargs(eq.fn))
    equation = nameof(eq.fn)
    isglob = isglobal(eq) === Val(true)
    elements = map(=>, map(r -> nameof(typeof(r)), eq.rheology), eq.el_number)
    inputs = keys(residual_kwargs(eq.fn))
    history = _history_keys(eq.rheology)
    return ntuple(Val(length(ks))) do i
        (;
            var = ks[i], equation, isglobal = isglob, elements,
            parent = eq.parent, children = eq.child, inputs, history,
        )
    end
end

@generated function _history_keys(rheology::NTuple{N, AbstractRheology}) where {N}
    keys = Symbol[]
    for T in rheology.parameters
        local_keys = T <: AbstractElasticity ? (:τ0, :P0) :
            T <: AbstractViscosity ? (:d,) : ()
        for key in local_keys
            key ∉ keys && push!(keys, key)
        end
    end
    return Expr(:tuple, (QuoteNode(key) for key in keys)...)
end

Base.summary(io::IO, insp::ModelInspection) = print(io, length(insp), "-element ", nameof(typeof(insp)))

const INSPECTION_HEADER = ("index", "var", "equation", "scope", "elements")

function Base.show(io::IO, ::MIME"text/plain", insp::ModelInspection)
    rows = [_inspection_row(i, entry) for (i, entry) in pairs(insp)]
    widths = ntuple(length(INSPECTION_HEADER)) do j
        maximum(textwidth(row[j]) for row in (INSPECTION_HEADER, rows...))
    end
    summary(io, insp)
    print(io, ":")
    for row in (INSPECTION_HEADER, rows...)
        println(io)
        _print_inspection_row(io, row, widths)
    end
    return nothing
end

function _inspection_row(i, entry)
    elements = join((string(name, " ", number) for (name, number) in entry.elements), ", ")
    return (
        string(i),
        string(entry.var),
        string(entry.equation),
        entry.isglobal ? "global" : "branch",
        elements,
    )
end

# The last column is left unpadded so that no row ends in whitespace.
function _print_inspection_row(io::IO, row, widths)
    print(io, "  ", lpad(row[1], widths[1]))
    for j in 2:lastindex(row)
        print(io, "  ", j == lastindex(row) ? row[j] : rpad(row[j], widths[j]))
    end
    return nothing
end
