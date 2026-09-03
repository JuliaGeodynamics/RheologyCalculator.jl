module RheologyCalculatorSparseConnectivityTracerExt

using RheologyCalculator
using RheologyCalculator: AbstractCompositeModel
using SparseConnectivityTracer: AbstractTracer
using StaticArrays

# Global tracers carry a sparsity pattern but no primal value, so the
# value-dependent guards in the solver must fall back to their unguarded form.
# The guard only ever avoided Inf/NaN in the *value*; the dependency pattern is
# identical either way.
@inline RheologyCalculator.safe_inv(v::AbstractTracer)     = inv(v)
@inline RheologyCalculator.safe_inv_one(v::AbstractTracer) = inv(v)

# Union the sparsity patterns of every numeric leaf of a nested structure
tracer_union(acc, v::Number) = acc + v
tracer_union(acc, v::NamedTuple) = tracer_union(acc, values(v))
tracer_union(acc, ::Any) = acc

function tracer_union(acc, v::Union{Tuple, AbstractArray})
    for value in v
        acc = tracer_union(acc, value)
    end
    return acc
end

"""
Return a conservative sparsity pattern for the local Newton
solve.

A global tracer records dependencies but has no numerical
values, so the Newton
iteration cannot run with one: it needs numerical Jacobians and
a linear solve.

Instead, every output is marked as depending on every traced
input found in `x`, `vars`, and `others`. This may include dependencies that are not
present numerically,
but it never omits one, which is the required property of a
sparsity pattern.
"""
function RheologyCalculator.solve(c::AbstractCompositeModel, x::SVector{N, T}, vars, others; kwargs...) where {N, T <: AbstractTracer}
    dep = zero(T)
    for value in (x, vars, others)
        dep = tracer_union(dep, value)
    end
    return SVector{N, T}(ntuple(_ -> dep, N))
end

end
