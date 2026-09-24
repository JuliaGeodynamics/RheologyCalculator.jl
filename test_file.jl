# MWE: generalised Maxwell — a composite where the solution vector contains two
# distinct `:τ` unknowns, and the second one is a REQUIRED OUTPUT rather than an
# internal detail.
#
#     SeriesModel
#       ├── LinearViscosity(η₁)                      bulk viscous response
#       └── ParallelModel
#             ├── SeriesModel(LinearViscosity(η₂),   ← Maxwell element:
#             │               IncompressibleElasticity(G))   viscous + elastic in series
#             └── LinearViscosity(η₃)                 parallel dashpot
#
# The nested `SeriesModel` is physically necessary here: it is the Maxwell
# element sitting inside a parallel branch. Both it and the composite as a whole
# carry a deviatoric stress, so `x_keys` contains `:τ` twice.
#
# Crucially, the nested `:τ` is the stress on the elastic spring — the value that
# must be carried to the next timestep as `τ0`. It is not something you would
# want the API to hide.

using RheologyCalculator, Test
using RheologyCalculator.RheologyModels
import RheologyCalculator: compute_stress_elastic, generate_equations, differentiable_kwargs

η₁, η₂, G, η₃ = 1.0e18, 1.0e19, 1.0e10, 1.0e18
ε, dt = 1.0e-14, 1.0e10

c      = SeriesModel(
    LinearViscosity(η₁),
    ParallelModel(SeriesModel(LinearViscosity(η₂), IncompressibleElasticity(G)),
                  LinearViscosity(η₃)),
)
vars   = (; ε)
others = (; dt, τ0 = (0.0,))          # one history entry per elastic source

x0  = initial_guess_x(c, vars, (; τ = 1.0), others)
sol = solve(c, x0, vars, others)

# MWE: `inspect` reports `branch` for equations that are not in any branch.
#
# `_inspection_row` (src/equation_system/inspect.jl:105) renders the scope column as
#
#     entry.isglobal ? "global" : "branch"
#
# which assumes "not global" == "lives in a parallel branch". That is false. The
# plastic consistency equations added by `add_local_equations` (equations.jl:80) are
# built with `Val(false)` while belonging to the *outer* SeriesModel, so a model with
# no ParallelModel at all still prints `branch`.
#
# This contradicts the `inspect` docstring ("whether that equation belongs to the
# outermost series model rather than to a parallel branch") and docs/src/composites.md.

using RheologyCalculator, Test
using RheologyCalculator.RheologyModels

has_parallel(c) = occursin("ParallelModel", string(typeof(c)))

# ---------------------------------------------------------------- the bug
# A flat viscoelastoplastic series. No branches anywhere.
vep = SeriesModel(LinearViscosity(1.0e23), Elasticity(1.0e10, 2.0e11), DruckerPrager(1.0e6, 30.0, 0.0))
@test !has_parallel(vep)

insp = inspect(vep)
println("VEP model (no ParallelModel):")
show(stdout, MIME"text/plain"(), insp); println("\n")

# the λ entry is reported as living in a branch that does not exist
λ_entry = only(filter(e -> e.var === :λ, collect(insp)))
@test λ_entry.equation === :compute_lambda
@test λ_entry.isglobal == false                    # -> rendered as "branch"
@test !has_parallel(vep)                           # ...but there is no branch

# ------------------------------------------------- why the column is unusable
# A model that DOES have a parallel branch produces an entry that is
# indistinguishable from the one above: same `isglobal`, same printed label.
kv = SeriesModel(LinearViscosity(1.0e19), ParallelModel(LinearViscosity(1.0e20), IncompressibleElasticity(1.0e10)))
@test has_parallel(kv)

branch_entry = only(filter(e -> !e.isglobal, collect(inspect(kv))))
println("Kelvin-Voigt model (has a real ParallelModel):")
show(stdout, MIME"text/plain"(), inspect(kv)); println("\n")

# identical scope flag for a real branch and a non-branch:
@test branch_entry.isglobal == λ_entry.isglobal == false

# and the rendered tables agree on the label, so a reader cannot tell them apart
scope_of(c, var) = begin
    line = only(filter(l -> occursin(" $var ", l), split(sprint(io -> show(io, MIME"text/plain"(), inspect(c))), "\n")))
    occursin("global", line) ? "global" : "branch"
end
@test scope_of(vep, "λ") == "branch"    # wrong: no branch in the model
@test scope_of(kv, "ε")  == "branch"    # right: genuinely in a branch

println("scope(VEP, :λ) = ", scope_of(vep, "λ"), "   <- model has ParallelModel? ", has_parallel(vep))
println("scope(KV,  :ε) = ", scope_of(kv, "ε"),  "   <- model has ParallelModel? ", has_parallel(kv))
println("\nBoth print `branch`; only one is in a branch.")


# ============================================================================
# The three-way scope, and what the two-way column currently prints instead.
#
# Balance equations belong to a network node:
#     compute_strain_rate, compute_volumetric_strain_rate  (SeriesModel node)
#     compute_stress, compute_pressure                     (ParallelModel node)
# Constraint equations close an element internally and correspond to NO node:
#     compute_lambda, compute_lambda_parallel,
#     compute_plastic_strain_rate, compute_volumetric_plastic_strain_rate
#
# proposed:  scope = is_constraint(fn) ? "local" : isglobal(eq) ? "global" : "branch"
# current :  scope =                                isglobal(eq) ? "global" : "branch"

const CONSTRAINT_FNS = (
    :compute_lambda, :compute_lambda_parallel,
    :compute_plastic_strain_rate, :compute_volumetric_plastic_strain_rate,
)
is_constraint(fn::Symbol) = fn in CONSTRAINT_FNS

current_scope(e) = e.isglobal ? "global" : "branch"
proposed_scope(e) = is_constraint(e.equation) ? "local" : e.isglobal ? "global" : "branch"

function compare_scopes(name, c)
    println("== ", name, "   (model contains ParallelModel: ", has_parallel(c), ")")
    println("   idx  var    equation                          current   proposed")
    for (i, e) in enumerate(collect(inspect(c)))
        cur, prop = current_scope(e), proposed_scope(e)
        flag = cur == prop ? "" : "   <-- changes"
        println("   ", rpad(i, 5), rpad(string(e.var), 7), rpad(string(e.equation), 34),
                rpad(cur, 10), rpad(prop, 9), flag)
    end
    println()
    return nothing
end

v(η) = LinearViscosity(η)
e_(G) = IncompressibleElasticity(G)
ec = Elasticity(1.0e10, 2.0e11)
dp = DruckerPrager(1.0e6, 30.0, 0.0)

compare_scopes("VEP (no ParallelModel)", SeriesModel(v(1.0e19), ec, dp))
compare_scopes("parallel + plastic", SeriesModel(v(1.0e19), ParallelModel(v(1.0e20), dp)))
compare_scopes("nested Maxwell (real branches)", SeriesModel(v(1.0e18), ParallelModel(SeriesModel(v(1.0e19), e_(1.0e10)), v(1.0e18))))
compare_scopes("Kelvin-Voigt", SeriesModel(v(1.0e19), ParallelModel(v(1.0e20), e_(1.0e10))))

# Entries whose label is wrong today are exactly the constraint equations.
for c in (SeriesModel(v(1.0e19), ec, dp), SeriesModel(v(1.0e19), ParallelModel(v(1.0e20), dp)))
    for e in collect(inspect(c))
        if is_constraint(e.equation)
            @test current_scope(e) == "branch"      # what it prints now
            @test proposed_scope(e) == "local"      # what it should print
        else
            @test current_scope(e) == proposed_scope(e)   # balance rows are unaffected
        end
    end
end
println("All balance rows keep their current label; only constraint rows change.")
