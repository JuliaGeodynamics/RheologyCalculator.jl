"""
    validate(r::AbstractRheology)
    validate(c::AbstractCompositeModel, vars, others)

Validate a rheology or a local solve input before entering the numerical hot
path. This is intentionally host-side; the validated model and inputs remain
immutable and suitable for device execution.

Validation is opt-in and does not run inside `solve`, residual, Jacobian, batch,
or tangent kernels. Viscosity history `d` is accepted when present but is not
required for viscosity laws that do not use it.
"""
function validate(r::AbstractRheology)
    isbitstype(typeof(r)) || throw(ArgumentError("$(typeof(r)) must be an isbits type"))
    for name in fieldnames(typeof(r))
        value = getfield(r, name)
        value isa Real || continue
        isfinite(value) || throw(ArgumentError("$(typeof(r)).$name must be finite"))
        name in (:η, :G, :K, :χ, :A, :ε0, :Q, :σb, :σr) && value > 0 ||
            name in (:η, :G, :K, :χ, :A, :ε0, :Q, :σb, :σr) &&
            throw(ArgumentError("$(typeof(r)).$name must be positive"))
    end
    return r
end

function validate(c::AbstractCompositeModel, vars, others)
    isbitstype(typeof(c)) || throw(ArgumentError("$(typeof(c)) must be an isbits type"))
    isbitstype(typeof(vars)) || throw(ArgumentError("vars must be an isbits type"))
    isbitstype(typeof(others)) || throw(ArgumentError("others must be an isbits type"))
    hasproperty(vars, :ε) || throw(ArgumentError("vars must contain `ε`"))
    all(isfinite, (vars.ε,)) || throw(ArgumentError("vars.ε must be finite"))
    _validate_elements(c.leafs, others)
    _validate_elements(c.branches, others)
    return c
end

function _validate_elements(elements::Tuple, others)
    for element in elements
        _validate_elements(element, others)
    end
    return nothing
end

function _validate_elements(r::AbstractRheology, others)
    validate(r)
    for key in history_kwargs(r)
        # Viscosity history `d` is model-specific; basic viscosity laws do not
        # need it. Elastic history is required whenever an elastic element is
        # present because those methods consume τ0/P0 directly.
        key === :d && !hasproperty(others, key) && continue
        hasproperty(others, key) || throw(ArgumentError("others must contain `$(key)` for $(typeof(r))"))
        history = getproperty(others, key)
        history isa Tuple || throw(ArgumentError("others.$key must be a tuple for $(typeof(r))"))
        isempty(history) && throw(ArgumentError("others.$key must not be empty for $(typeof(r))"))
        all(isfinite, history) || throw(ArgumentError("others.$key must contain finite values"))
    end
    return nothing
end

function _validate_elements(c::AbstractCompositeModel, others)
    _validate_elements(c.leafs, others)
    _validate_elements(c.branches, others)
    return nothing
end
