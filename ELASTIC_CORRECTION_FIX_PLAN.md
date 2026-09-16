# Plan: fixing the elastic strain-rate correction

Delivered as **one PR**.

References:

- Findings: `docs/derivations/tensor_reduction_audit.md` (S*, D* IDs)
- Math: `docs/derivations/tensor_reduction.typ` (up to `<series-eff>`), `docs/derivations/noncoaxial_reduction.typ` (`<eff>`, `<constitutive>`, `<scalar>`).
  The remaining formulas (recovery, history update, nonlinear closure, tangent) are stated in §1 so this plan stands on its own.

---

## 1. Target formulation

Composite: a `SeriesModel` with leafs *l* (any rheology) and `ParallelModel` blocks *b*. A block holds direct elements (dashpots, springs *s*) and `SeriesModel` sub-branches *j* made of leaf elements only. A block is **history-bearing** if it contains a spring.

### 1.1 Coefficients

Per history-bearing block *b* and sub-branch *j*, with every element represented by its secant viscosity η (spring: η = GΔt):

- η_M,j = (Σ_{elements of j} 1/η)⁻¹
- η*ᵢ = 1 for a direct spring; η*ᵢ = η_M,j/(GᵢΔt) for spring *i* of sub-branch *j* (its **own** G)
- η_KV,b = Σ_{direct elements} η + Σ_j η_M,j

### 1.2 History tensor and scalar equation

- **H** = Σ_{outer-series springs k} τᵒ_k/(2G_kΔt) + Σ_b (1/(2η_KV,b)) Σ_{i∈b} η*ᵢ τᵒᵢ
- **E** = ε̇ + **H** (assembled as a tensor); E_II = second_invariant(**E**); **n** = **E**/E_II
- The global equation is the invariant of the assembled tensor: (series-leaf strain rates at τ_II) + Σ_b (block unknown) = E_II. It is **never** ε̇_II + Σ(...)_II.

### 1.3 Shifted local unknowns

Elastic state functions ignore τ0. So inside a history-bearing block the graph unknowns are **shifted** values: x̃_b (block strain rate) and x̃_j (sub-branch stress) of the history-free problem with the same τ and the same coefficients. The shifted tensors x̃_b **n** and x̃_j **n** are all coaxial with τ, which is why the graph's scalar "children − parent" equations are exact.

The physical tensors are:

- ε̇ᵖ_b = x̃_b **n** − (1/(2η_KV,b)) Σ_{i∈b} η*ᵢ τᵒᵢ
- τ_j = x̃_j **n** − (η_M,j/η_KV,b) Σ_{i∈b} η*ᵢ τᵒᵢ + Σ_{k∈j} η*_jk τᵒ_jk

### 1.4 Nonlinear elements in history-bearing blocks (exact)

A state-dependent element must be evaluated at its **physical** state, but the coefficients in §1.3 depend on those same states. The loop is closed with physical invariants as extra unknowns:

- **e_b** ≈ (ε̇ᵖ_b)_II for a block with a state-dependent direct element
- **t_j** ≈ (τ_j)_II for a sub-branch with a state-dependent element

with closure residuals built from assembled tensors:

- r_b = e_b − ‖x̃_b **n** − Σ η*ᵢ τᵒᵢ/(2η_KV,b)‖_II
- r_j = t_j − ‖x̃_j **n** − (η_M,j/η_KV,b) Σ_{i∈b} η*ᵢ τᵒᵢ + Σ_{k∈j} η*_jk τᵒ_jk‖_II

Secant viscosities come from the existing state functions, evaluated at the physical invariants (so no `compute_viscosity` is needed and nothing is implicit):

- direct element of block *b*: η = compute_stress(r; ε = e_b)/(2 e_b)
- element of sub-branch *j*: η = t_j / (2 compute_strain_rate(r; τ = t_j))

The graph equations use the same secant viscosities on the shifted unknowns. The element's contribution is its state function at the physical argument, rescaled:

- direct element: compute_stress(r; ε = e_b) · x̃_b/e_b
- sub-branch element: compute_strain_rate(r; τ = t_j) · x̃_j/t_j

Because the graph, **H**, and the closures all use one set of coefficients evaluated at (e_b, t_j), the converged solution satisfies the tensor equations exactly. Linear elements give ratios that don't depend on their argument, so they reduce to today's evaluation. **Zero-argument guard:** the ratios x̃/e and x̃/t use the package's `safe_inv` pattern, which stays compatible with the sparsity tracer.

Blocks and sub-branches where every element is state-independent get **no** extra unknowns. Their coefficients are constant within the step.

### 1.5 Recovery and history update (after the solve)

- τ = τ_II **n**
- ε̇ᵖ_b and τ_j from §1.3, using the converged coefficients
- one tensor per spring, in `global_eltype_numbering` order:
  - outer-series spring: τ
  - direct spring of block *b*: τᵒ_s + 2G_sΔt ε̇ᵖ_b
  - spring of sub-branch *j*: τ_j
- check (in tests): τ = Σ_{direct} (stress at ε̇ᵖ_b) + Σ_j τ_j

### 1.6 Tangent

τ = τ_II(x(ε̇)) **n**(x(ε̇), ε̇). Differentiate the residual with respect to the ε̇ tensor components (ForwardDiff over an `SVector`), solve J dx/dε̇ = −∂R/∂ε̇, and chain through **n**. With coefficients that don't depend on the state this reduces to the symmetric form dτ = τ_II' **n**⟨**n**, dε̇⟩ + (τ_II/E_II)(dε̇ − **n**⟨**n**, dε̇⟩). With state-dependent coefficients it is non-symmetric in general.

### 1.7 Assumptions (documented, enforced where possible)

- isotropic laws; inelastic strain rates of outer-series leafs coaxial with τ (J2-type flow)
- sub-branches contain only leaf elements (error otherwise, S3)
- τ0 is supplied in the current frame (the caller rotates it)

---

## 2. Constraints

| Constraint | Test |
|---|---|
| zero allocations in `solve` / residual / Jacobian | `test_allocations.jl` |
| type stability | `test_type_stability.jl` |
| GPU invariants | `test_gpu_invariants.jl` |
| SparseConnectivityTracer compatibility | `test_sparsity_tracer.jl` |
| ForwardDiff Jacobians | `test_jacobians.jl`, `test_jacobian_backend.jl` |
| compile-time structure (`@generated`, literal τ0 indices, extra unknowns decided by type) | design |
| fail-fast on unsupported configurations | `test_simplification_plan.jl` |
| scalar-ε results unchanged, history signed | `test_VE*.jl`, `test_analytical_regressions.jl` |
| unchanged `x` layout for models without state-dependent elements in history-bearing blocks | `test_equation_graphs.jl`, partition tests |

---

## Step 1: References and failing tests (no source changes)

- [x] `test/reference/tensor_reference.jl` (test-only). A monolithic full-tensor solve per step.
  - Unknowns: the τ tensor, the block strain-rate tensors, the sub-branch stress tensors.
  - Residuals: the constitutive laws element by element.
  - Newton with a ForwardDiff Jacobian, which also covers nonlinear laws. No package solver code; only the package's scalar state functions for element laws.
- [x] Topologies:
  - linear: Maxwell `S(η,G)`; KV `S(η₁,P(η₂,G))`; mixed `S(η₁,P(η₂,S(η₃,G)))`; Burgers `S(η₁,G₁,P(η₂,G₂))`; two blocks `S(η₁,P(η₂,G₁),P(η₃,S(η₄,G₂)))`; direct spring plus sub-branch `S(η₁,P(G₁,S(η₃,G₂)))`; two springs in a sub-branch `S(η₁,P(η₂,S(η₃,G_a,G_b)))` with G_a ≠ G_b
  - nonlinear: power law as a direct KV element `S(η₁,P(PL,G))`; power law in a Maxwell sub-branch `S(η₁,P(η₂,S(PL,G)))`; power-law outer leaf with a linear block `S(PL,P(η₂,G))`
- [x] `test/test_elastic_history_tensor.jl`: time loop of `solve` + `elastic_stress_history_2D/3D` against the reference, comparing τ and every spring tensor at every step (`rtol ~1e-9`). Load cases:
  1. constant pure shear
  2. load reversal after n steps
  3. change of direction, pure shear → simple shear
  4. 3D, all six components
  5. scalar ε, constant and reversed
- [x] `@test_broken` for each confirmed failure; record confirmed / not reproduced per ID in the audit notes.

## Step 2: Local fixes and guards

- [x] **S2.** Per-spring η* using each spring's own G. Replace `_η_eff_elastic` (`findfirst`) with a version indexed by position.
- [x] **S3.** `_assert_kv_nesting_supported` rejects **any** composite nested inside a sub-branch of a history-bearing block. Add a `@test_throws "…"`.
- [x] **S13.** The "no spring found" fallback becomes an error.
- [x] **S15, S16, S19.** Remove dead code (`_direct_leaf_correction_scalar`, `strain_rate_correction.jl:97`, `solver.jl:175-178`).
- [x] **S22.** Use `zero` of the element type for fold seeds where the tracer allows it.
- [x] **Trait (exported):** `viscosity_depends_on_state(::AbstractRheology)`, documented for user rheologies.
  - default `true`: selects the exact nonlinear path (§1.4), which is also correct for linear laws, just more expensive;
  - `false` for audited linear rheologies (`LinearViscosity`, `Elasticity`, `IncompressibleElasticity`, `BulkElasticity`, …).
  - The trait only selects between two exact paths, so a wrong `true` costs performance, never correctness. A wrong `false` on a nonlinear law would be silently wrong; say so in the docstring.

## Step 3: Tensor strain rate reaches the residual; single tensor **H** (S1)

- [x] **Confirm the readers of `vars.ε`.** It is expected to be read only by `subtract_parent` for the global equation (`vars[eq.ind_input]`), since `generate_args_template` builds element arguments from `x` and `others`. List any other reader and adapt it.
- [x] `solve`, `tangent`, `tangent_block` and `tangent_tensor` stop reducing ε before the solve. `vars.ε` stays a tensor (or a `Number`) down to the residual.
- [x] **`_history_tensor(c, eqs, x, ε, others)`** (`@generated`): **H** from §1.2, coefficients from §1.1. Constant blocks use constant coefficients; blocks with extra unknowns use coefficients at (e_b, t_j) from `x`. It reuses `_accumulate_weighted_backstress(identity, …)`, `_checked_η_KV` and the literal τ0 offsets (`_branch_tau0_offsets`). For scalar ε it returns a signed scalar.
- [x] **Global residual:** `subtract_parent` for the global deviatoric equation subtracts `second_invariant_value(ε .+ H)`. This replaces both the pre-solve `_direct_leaf_elastic_correction` and `subtract_elastic_correction`.
- [x] **Remove** the scalar implicit path: `subtract_elastic_correction`, `_subtract_first`, `_implicit_elastic_correction`, `_kv_implicit_corrections_scalar`, `_kv_implicit_branch_correction_scalar`, `_nth_true_position`, `_weighted_backstress_scalar`, `_direct_leaf_elastic_correction`. Update the imports in `test/runtests.jl` and `test/test_solver_convergence.jl`.
- [x] **Zero E:** E_II has no derivative at **E** = 0, which a load reversal can pass through. Use a tracer-compatible guard (in the style of `safe_inv`) for the gradient, and add a test that passes through reversal and a test that starts from rest.
- [x] **S24.** For a signed scalar ε, block strain-rate unknowns must not be bounded to be non-negative (`branch_strain_rate_mask`): bound them only for tensor ε, where they are invariants.
- [x] **`effective_strain_rate_correction`** (exported): keep `(c, ε, τ0, others)` for models whose history-bearing blocks have constant coefficients. Add a `(c, x, ε, others)` method for the general case. The old method throws a clear error for models that need `x`.

Exit: Step 1 linear cases match the reference for τ.

## Step 4: Exact nonlinear elements in history-bearing blocks (S4, S5)

- [ ] **Equation graph.** For each history-bearing block with a state-dependent direct element, add one unknown e_b; for each sub-branch with a state-dependent element, add one unknown t_j. Decide at the type level (`@generated`) so other models keep their current graph.
  - New `CompositeEquation` state functions: `compute_block_strain_rate_invariant` (r_b) and `compute_subbranch_stress_invariant` (r_j), with no parent/children contributions, handled like `compute_lambda` in `add_children` / `subtract_parent`.
  - Their residuals (§1.4) need **E**, so they receive `vars.ε` and `others.τ0`.
- [ ] **Secant evaluation** (§1.4). In `evaluate_state_function_perleaf`, for equations of a history-bearing block with extra unknowns, evaluate each state-dependent element at the physical argument (e_b or t_j) and rescale by the shifted/physical ratio. Linear elements are evaluated as today.
- [ ] **One coefficient function** `_block_coefficients(branch, x, others)` used by the graph rescaling, `_history_tensor`, the closures, and the post-processing. This removes the three different evaluation points (S5).
- [ ] **Layout-dependent API** for models that gain unknowns: `x_keys` (new keys, e.g. `:εII_block`, `:τII_sub`), `initial_guess_x` (e_b from the shifted guess; t_j from the shifted stress), `normalisation_x`, `branch_strain_rate_mask` (e_b and t_j are non-negative), `inspect`.
- [ ] Allocation, type-stability, sparsity-tracer and GPU tests for at least one nonlinear topology from Step 1.
- [ ] **Plastic elements in history-bearing blocks.** Plastic elements carry their own λ unknowns and local equations (`compute_lambda`, `compute_plastic_strain_rate`), not a secant stress law, so §1.4 does not apply directly. Investigate `P(DruckerPrager, η_reg, G)` and `P(η, S(DruckerPrager, G))`:
  - derive the physical argument each plastic equation needs (yield on the physical τ_j or on the block stress share);
  - check whether the λ equations can read e_b / t_j and the physical tensors;
  - add a reference case to Step 1.

  If it can be made exact, implement it here. Otherwise this is the one remaining open decision (§Decisions).

Exit: Step 1 nonlinear cases match the reference.

## Step 5: Recovery, history update, tangents (S6, S7, S23)

- [ ] **`elastic_stress_history(c, x, ε, τ0, others)`** (dimension-agnostic, NTuple-based, allocation-free), following §1.5. `elastic_stress_history_2D/3D(c, x, ε, τ0, others)` delegate to it.
  - The coaxial overload `(c, τII, ε, τ0, others)` stays for models without history-bearing blocks and throws for the others.
- [ ] **`compute_stress_elastic`:**
  - new `compute_stress_elastic(c, x, ε, others)` returns the invariants of the §1.5 tensors (signed for scalar ε);
  - the existing `(c, x, others)` keeps working for models whose springs are all outer-series leafs, and throws with a pointer to the ε method otherwise;
  - remove `_branch_elastic_stress`, `_branch_elastic_info`, `_kv_corrections_elastic_stress`.
- [ ] **S23:** check whether `inspect`, `component_partition` and `dissipation_partition` read shifted block/sub-branch unknowns as physical values. Where they need physical values, use §1.3 (with ε) or the e_b/t_j unknowns. Add a reference comparison for the KV and nonlinear KV cases.
- [ ] **Tangents** (§1.6):
  - `tangent_tensor` differentiates the residual with respect to the ε̇ components and chains through **n**;
  - `tangent` / `tangent_block` return dτ_II/dε̇ along the applied direction, documented;
  - test all of them against finite differences of `solve` on non-coaxial linear and nonlinear KV cases.
- [ ] Update and run the examples that call these functions (`Burgers.jl`, `transient_creep.jl`, `Maxwell_KV_Maxwell.jl`, `elastic_corrections_playground*.jl`, and the `Maxwell_*` examples using `compute_stress_elastic`).

Exit: all Step 1 cases pass for τ and every spring tensor at every step; existing analytical regressions unchanged; tangents match finite differences.

## Step 6: Comments and docs

- [ ] Source comments and docstrings:
  - S8, S9, S14, S18 (describe **H** and the exact nonlinear closure);
  - S10 (sign), S11/S12 (sum / inverse sum), S17 (offset), S20 (parallel elements share strain rate), S21 (unsupported "proven");
  - remove the "Option 3" / "old approach" text with the removed code.
- [ ] `docs/src/strain_rate_correction.md`:
  - D1/D2 (tensor **E**, invariant of the assembled tensor), D3 (history update), D4 (per-spring η*), D5, D6, D8, D10;
  - supported topologies, the shifted unknowns, the e_b/t_j unknowns and their `x_keys`;
  - the `viscosity_depends_on_state` trait for user rheologies;
  - an updated implementation table;
  - an example with load reversal and a change of direction (linear Burgers and a nonlinear KV), checked against the reference.
- [ ] `docs/src/api.md`: the new and changed methods.
- [ ] Typst derivations: add §1.3–§1.5 (shifted unknowns, closure, recovery), since the docs cite them.
- [ ] Changelog: results change for tensor-ε models with blocks; new unknowns for nonlinear history-bearing blocks; the `compute_stress_elastic` / `effective_strain_rate_correction` method changes.
- [ ] Audit notes: final status per ID.

---

## Decisions

| ID | Decision |
|---|---|
| Q1 | Nonlinear elements in history-bearing blocks are **supported exactly** (Step 4) |
| Q2 | New `compute_stress_elastic(c, x, ε, others)`; the old method keeps working for outer-series springs only (pending confirmation) |
| Q3 | `viscosity_depends_on_state` exported, default `true` |
| Q4 | One PR |
| **Open** | Plastic elements in history-bearing blocks: exact support if the Step 4 investigation finds a formulation, otherwise to be decided |

## Risks

- **x-layout change** for nonlinear history-bearing blocks (new unknowns). This is breaking for users who index `x` by position in those models; the other models keep their layout.
- **Newton robustness** with the extra closure unknowns and near E = 0 (reversal). Watch iteration counts in the reversal tests.
- **Cost:** tensor ε in every residual evaluation, plus the extra unknowns where needed. Benchmark Burgers and a nonlinear KV case before and after.
- **Wrong `false` trait** on a user nonlinear law gives silently wrong results; the docstring must be explicit.
