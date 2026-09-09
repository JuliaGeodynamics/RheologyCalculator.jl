# Tuple walkers used throughout equation generation and residual assembly.
#
# These replace hand-written `@generated` bodies. Written as recursion over
# `first`/`Base.tail`, they unroll and infer for tuples of mixed element types at
# any length; `ntuple(f, Val(N))` is explicitly unrolled only up to N = 10 and
# falls back to a path that does not infer for heterogeneous results, which the
# equation tuples are — a `CompositeEquation` carries its state function in a
# type parameter, so no two entries share a type.

"""
    maptuple(f, t::Tuple)

Apply `f` to each element of `t`, returning a tuple of the results.
"""
@inline maptuple(f::F, t::Tuple) where {F} = (f(first(t)), maptuple(f, Base.tail(t))...)
@inline maptuple(::F, ::Tuple{}) where {F} = ()

"""
    maptuple(f, t1::Tuple, t2::Tuple)

Apply `f` to corresponding elements of two tuples walked in lockstep, returning
a tuple of the results. The tuples must be the same length.
"""
@inline maptuple(f::F, t1::Tuple, t2::Tuple) where {F} =
    (f(first(t1), first(t2)), maptuple(f, Base.tail(t1), Base.tail(t2))...)
@inline maptuple(::F, ::Tuple{}, ::Tuple{}) where {F} = ()

"""
    foldtuple(op, init, f, t::Tuple)

Left-fold `op` over `f` applied to each element of `t`, starting from `init`.
"""
@inline foldtuple(op::O, init, f::F, t::Tuple) where {O, F} =
    foldtuple(op, op(init, f(first(t))), f, Base.tail(t))
@inline foldtuple(::O, acc, ::F, ::Tuple{}) where {O, F} = acc

"""
    flatmaptuple(f, t::Tuple)

Apply `f` to each element of `t`, where `f` returns a tuple, and splice the
results into a single flat tuple.
"""
@inline flatmaptuple(f::F, t::Tuple) where {F} = (f(first(t))..., flatmaptuple(f, Base.tail(t))...)
@inline flatmaptuple(::F, ::Tuple{}) where {F} = ()
