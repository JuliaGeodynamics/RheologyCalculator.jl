# Tensor reduction / elastic strain-rate correction: audit notes

Sources audited:

| Tag | File |
|---|---|
| **TYP** | `docs/derivations/tensor_reduction.typ` |
| **SRC** | `src/post_processing/strain_rate_correction.jl`, plus call sites in `src/equation_system/{solver,equations,tangent}.jl`, `src/post_processing/post_calculations.jl`, `src/utils/tensor_helpers.jl` |
| **DOC** | `docs/src/strain_rate_correction.md` |
| **PDF** | `creep.pdf` ("Generalized visco-elastic element") |

Status legend:

- **✔ verified**: checked by hand algebra in this pass
- **? suspected**: follows from reading the code/derivation; needs a numerical test (Step 4)
- Severity: **[BUG]** wrong results possible · **[PHYS]** modeling difference / approximation with consequences · **[MATH]** error in a written equation · **[DRIFT]** docs/comments disagree with code · **[NIT]** typo, terminology, style

Notation used below: ε = strain rate (TYP/SRC/DOC write ε, PDF writes ε̇), τᵒ = τ0 = spring stress at the previous step,
η_M = 1/(1/η + 1/(GΔt)) (Maxwell effective viscosity), η_KV = Σ effective viscosities in a parallel block,
η* = η_M/(GΔt) = η/(η+GΔt).

---

## 1. Equations extracted from each source

### 1.1 TYP

| Case | Result |
|---|---|
| Maxwell `S(η, G)` | ε_eff = ε + τᵒ/(2GΔt) = τ/(2η) + τ/(2GΔt) |
| KV `S(η₁, P(η₂, G))` | εᵖ = (τ − τᵒ)/(2η_KV), η_KV = η₂ + GΔt; ε_eff = ε + τᵒ/(2η_KV) |
| Mixed `S(η₁, P(η₂, S(η₃, G)))` | τ^Max = 2η_M(εᵖ + τᵒ/(2GΔt)); εᵖ = (τ − 2η_M τᵒ/(2GΔt))/(2η₂ + 2η_M); ε = τ/(2η₁) + εᵖ |
| Generalized | τ/(2η_KV) = εᵖ + (1/(2η_KV))(**2** Σ η*_i τᵒ_i); η_KV = Σ η_eff,i; η*_i = η_i/(η_i + G_iΔt) (Maxwell) or 1 |
| Newton | r = ε̇ − f(τ) = 0, τⁿ⁺¹ = τⁿ − J⁻¹ r |

### 1.2 SRC

| Quantity | Code | Formula |
|---|---|---|
| Direct elastic leaf of outer series | `effective_strain_rate_correction(::Val{true}, c::AbstractRheology, …)` | τ0/(2·compute_viscosity) = τ0/(2GΔt), **tensor** |
| η_KV | `_η_KV` | Σ_direct leafs η_l (dashpot η, spring GΔt) + Σ_sub-branches η_M |
| η_M (sub-branch) | `_η_eff_maxwell` | 1/Σ_leafs(1/η_l) over *all* leafs of the sub-branch |
| η_el (sub-branch) | `_η_eff_elastic` | GΔt of the **first** spring in the sub-branch (`findfirst`) |
| ws | `_accumulate_weighted_backstress` | Σ_direct τ0 + Σ_sub (η_M/η_el) τ0 |
| Branch correction (tensor, public API) | `_kv_branch_correction` | ws/(2η_KV), viscosities evaluated at ε = ‖ε_total‖_II |
| Branch correction inside `solve` | `_kv_implicit_branch_correction_scalar` | Σ η* ‖τ0‖_II/(2η_KV), viscosities at ε = x[branch eq] |
| Global residual | `compute_residual` → `subtract_elastic_correction` | R₁ = f(x) − ‖ε + leafcor‖_II − Σ_branches ws_II/(2η_KV) |
| Post-solve spring stress (branch) | `_branch_elastic_stress` | εᵖ = (τ_II − ws_II)/(2η_KV); τ_spring = τ0_II + 2GΔt εᵖ (direct) or 2η_M(εᵖ + τ0_II/(2η_el)) (sub); viscosities at **ε = 0.0** |
| Tensor history | `_elastic_stress_history` | every spring tensor = coaxial with ε_eff (full tensor correction) |

How the pieces are wired:

- `solve`, `tangent`, `tangent_block` and `tangent_tensor` apply only the **direct-leaf** correction, in tensor form, before the solve.
- Branch corrections are applied **implicitly and in scalar form** inside the residual.

### 1.3 DOC

The Maxwell, KV and mixed cases match TYP. The generalized result is
Δε_eff = (1/(2η_KV)) Σ η*_i τᵒ_i, and DOC writes the total ε_eff as the leaf correction plus the branch corrections, all in tensor form.
DOC says the full correction is applied "before each Newton step".

### 1.4 PDF

The generalized element is a reference Maxwell element (s₀, G₀, η₀) in parallel with N Maxwell branches (sᵢ, Gᵢ, ηᵢ):

- (1) ṡ₀/(2G₀) + s₀/(2η₀) = ε̇ (reference); ṡᵢ/(2Gᵢ) + sᵢ/(2ηᵢ) = ε̇ (non-reference)
- (2) s₀ = τ − Σ sᵢ
- (3) α₀(τ − Σsᵢ) − β₀(τ* − Σsᵢ*) = ε̇, with α = β + 1/(2η) and β = 1/(2GΔt)
- (4) sᵢ = ε̇/αᵢ + (βᵢ/αᵢ) sᵢ*
- (5) ε̇ = (1 + Σ α₀/αᵢ)⁻¹ (α₀τ − β₀τ* + Σ(β₀ − βᵢα₀/αᵢ) sᵢ*)
- (6) τ = 2η* ε̇*, reduced to the scalar residual r(τ_II) = τ_II − 2η* ε̇*_II
- Nonlinear CASE 1: parameters depend on the top-level τ_II.
- Nonlinear CASE 2: parameters depend on the branch sᵢ, with extra unknowns (sᵢ)_II.

---

## 2. Errors and issues found

### 2.1 In TYP (`tensor_reduction.typ`)

**Status: T1–T11 all fixed in the rewritten `tensor_reduction.typ`.** The rewrite adds sections for notation and element laws,
Maxwell sub-branches with several elements, the generalized parallel block, the general series composite, and reduction to a scalar equation.
The table below describes the file as it was before the rewrite.
The non-coaxial case is derived in `noncoaxial_reduction.typ`, which covers tensor recovery, history update,
nonlinear viscosities with block-local unknowns, the consistent tangent, and history rotation. It is the
reference for S1, S5, S6 and S7 and for the `tangent_tensor` concern.

| # | Sev | Status | Location | Issue |
|---|---|---|---|---|
| T1 | MATH | ✔ | Generalized Maxwell, first equation | Spurious factor **2**: `(1/(2η_KV))(2 Σ η* τᵒ)` should be `(1/(2η_KV)) Σ η* τᵒ`. From εᵖ = (τ − Σ η*τᵒ)/(2η_KV) (the KV/mixed results) you get τ/(2η_KV) = εᵖ + Σ η*τᵒ/(2η_KV). SRC and DOC use the correct form. |
| T2 | MATH | ✔ | Generalized Maxwell, η*_M | η*_M = η_i/(η_i + G_iΔt) is only right for a sub-branch with **one** dashpot and **one** spring. With several leafs in series the weight of spring k is η_M/(G_kΔt), where η_M is the inverse sum over all leafs. The general form isn't stated. |
| T3 | MATH | ✔ | Generalized Maxwell | η_eff,i is never defined for each kind of element (dashpot → η, direct spring → GΔt, Maxwell sub-branch → η_M). The "1 else" case also covers non-elastic entries, which have no τᵒ. |
| T4 | NIT | ✔ | Generalized Maxwell | "where $k$ is the branch number": the formula uses index i, and k doesn't appear. |
| T5 | MATH | ✔ | Mixed KV–Maxwell, last line | Stops at ε = τ/(2η₁) + (τ − 2η_M τᵒ/(2GΔt))/(2η₂ + 2η_M) without writing ε_eff. It never names η_KV = η₂ + η_M or η* = η_M/(GΔt), so it doesn't connect to the generalized section. |
| T6 | NIT | ✔ | Headings | `= Generalized Maxwell Body` is a level-1 heading, but the other cases are level-2 (`==`). |
| T7 | NIT | ✔ | Throughout | Strain rate written as ε (no dot), which clashes with the PDF's ε̇. There's also a `bold(tau^o)` vs `bold(tau)^o` mix: the superscript ends up inside or outside the bold. |
| T8 | NIT | ✔ | End of file | Stray fragment `$cal(O)(4)$`. The Newton snippet doesn't define J (it should be J = ∂r/∂τ = −∂f/∂τ). |
| T9 | NIT | ✔ | Mixed section | "Lets get" → "Let's get". |
| T10 | MATH | ✔ | Whole file | Everything is written as tensors, but nothing says how the tensor result reduces to the scalar invariant the solver uses. That is exactly where SRC departs from it (see S1). |
| T11 | MATH | ✔ | Whole file | No outer-series Maxwell leaf appears alongside a branch (e.g. `S(G₁, P(η₂, S(η₃, G₂)))`), so the total ε_eff = ε + leaf + branch terms is never derived. DOC states it without derivation. |

### 2.2 In SRC

| # | Sev | Status | Location | Issue |
|---|---|---|---|---|
| S1 | BUG/PHYS | ✔ confirmed | `_kv_implicit_branch_correction_scalar`, `_weighted_backstress_scalar`; `subtract_elastic_correction` | Inside `solve` the branch correction is a **sum of invariants**: ‖ε + leafcor‖_II + Σ η*‖τ0_i‖_II/(2η_KV). The correct quantity is the **invariant of the tensor sum** ‖ε + leafcor + Σ η*τ0_i/(2η_KV)‖_II. (a) With non-coaxial τ0 and ε this overestimates ε_eff by the triangle inequality. (b) **On load reversal** (τ0 opposite to ε), ‖τ0‖_II ≥ 0 loses the sign: the correction *adds* where it should *subtract*. The KV spring then pushes in the wrong direction, so stress relaxation and recovery (e.g. Burgers creep recovery) come out wrong. The PDF's explicit requirement (p. 6 NOTE: compute the full effective strain-rate tensor, then reduce) is violated. The tests only use positive scalar τ0 and one-component ε, so they can't catch this. |
| S2 | BUG | ✔ confirmed | `_η_eff_elastic` (`findfirst`) used in `_accumulate_weighted_backstress` and `_branch_elastic_info` | A Maxwell sub-branch with **two or more springs** in series gives every spring the weight η_M/(G₁Δt) of the *first* spring. The correct weight is η_M/(G_kΔt) (derived: εᵖ = τ^M/(2η_M) − Σ_k τᵒ_k/(2G_kΔt)). Wrong whenever G_a ≠ G_b. It's also silently wrong in `_branch_elastic_stress`. |
| S3 | BUG | ? | `_η_KV` → `_η_eff_maxwell(sub.leafs, …)`; `_assert_kv_nesting_supported` | `_η_eff_maxwell` only looks at `sub.leafs` and ignores `sub.branches`. The nesting guard only raises when the deeper nesting contains a spring, so a **non-elastic** parallel block nested inside a Maxwell sub-branch (e.g. `P(η₂, S(η₃, G, P(ηa, ηb)))`) is silently left out of η_M. That gives wrong η_KV and η*, with no error. This breaks the fail-fast rule. |
| S4 | BUG/PHYS | ? not isolated (masked by S5) | `_kv_implicit_branch_correction_scalar` (args ε = `x[bpos]`) | For nonlinear viscosities, η_KV and η* are evaluated at ε = x[branch eq]. That value (a) is **shifted** from the physical branch strain rate: the branch equation uses τ0-free springs, so x[bpos] = εᵖ + ws/(2η_KV). (b) Even the physical εᵖ is **not** the strain rate of a dashpot inside a Maxwell sub-branch, which carries εᵖ minus the spring's share. The block comment ("correct for both linear and nonlinear viscosities", "true fully-implicit solution") is therefore unsupported. Compare PDF CASE 1 (parameters use top-level τ_II) and CASE 2 (extra unknowns (sᵢ)_II). |
| S5 | BUG/PHYS | ✔ confirmed | Nonlinear η in branches: `_kv_branch_correction` (ε = ‖ε_total‖), `_kv_implicit_branch_correction_scalar` (ε = x[bpos]), `_branch_elastic_stress` (**ε = 0.0**) | One time step evaluates the same η_KV/η* at **three different strain rates**. The public tensor correction (used by `elastic_stress_history_*`), the solve, and the post-solve spring-stress reconstruction disagree for any nonlinear branch viscosity. ε = 0.0 gives η = ∞ or 0 for a power law, so `_checked_η_KV` will throw or produce ∞. The comment at `post_calculations.jl:204-209` justifies ε = 0.0 by pointing at "the same level of rigor" as the pre-solve path, which the solve no longer uses. |
| S6 | PHYS | ✔ confirmed | `_elastic_stress_history` | Every spring's history tensor is rebuilt **coaxial with ε_eff**. For a nested (KV / sub-branch) spring its true stress is generally not coaxial with the total ε_eff, e.g. during rotation of loading or on reversal while the spring still carries the old orientation. Combined with S1, orientation and sign information about τ0 are both lost. |
| S7 | BUG | ✔ confirmed | `_branch_elastic_stress` | Same sum-of-invariants and sign issue as S1: εᵖ = (τ_II − Σ η*‖τ0‖_II)/(2η_KV), then τ_spring = ‖τ0‖_II + 2GΔt εᵖ. On reversal the reconstructed spring stress has the wrong magnitude. |
| S8 | DRIFT | ✔ | `strain_rate_correction.jl:19-28` header | "Called once per solve step (in solve()…)… ε_eff = ε + correction(τ0)… what the Newton solver actually sees". Only the direct-leaf part is applied pre-solve; branch parts are implicit and scalar. |
| S9 | DRIFT | ✔ | Docstring `effective_strain_rate_correction` (l. 43-48) | "Elastic elements contribute τ0/(2η)". Wrong for branch springs, which contribute η*τ0/(2η_KV). |
| S10 | MATH (comment) | ✔ | l. 110-113 | "…which must be **subtracted** from the total ε before solving": the code **adds** it (ε_eff = ε + τ0/(2GΔt)). The physical wording is also confused: τ0/(2GΔt) is a strain rate, not "strain … times dt". |
| S11 | NIT | ✔ | l. 171; `post_calculations.jl:190` | "arithmetic mean": η_KV is a **sum**, not a mean. |
| S12 | NIT | ✔ | `_η_eff_maxwell` docstring and comment (l. 386-395); DOC | "harmonic mean": 1/Σ(1/η) is the inverse sum (1/N of the harmonic mean), not the harmonic mean. |
| S13 | NIT | ✔ | l. 407-417 | `_η_eff_elastic` returns `0.0` when no spring is found ("safe fallback"). It's unreachable today, but if it were reached η_M/0 = Inf would silently propagate. That should be an error (fail-fast). |
| S14 | DRIFT | ✔ | l. 506-526 block comment | Refers to "Option 3" (a planning artifact) and "the old pre-correction approach" (history). It claims the correction includes "Σ direct-leaf τ0/(2η)", but `_implicit_elastic_correction` explicitly excludes it (l. 541-543). |
| S15 | NIT | ✔ | l. 546-564 `_direct_leaf_correction_scalar` | Dead code: never called anywhere in `src/` or `test/`. It also calls `compute_viscosity(leafs[pos], others)` without `ε`. |
| S16 | NIT | ✔ | l. 97 | Commented-out dead code. |
| S17 | NIT | ✔ | l. 315 vs `post_calculations.jl:196` | `el_idx_start` is described as the "1-based index" in one place and the "0-based index" in the other. It is a count, or offset, of entries already consumed. |
| S18 | DRIFT | ✔ | `solver.jl:92-94` docstring | "the solver applies `effective_strain_rate_correction`": it only applies `_direct_leaf_elastic_correction`. |
| S19 | NIT | ✔ | `solver.jl:175-178` | Commented-out dead code with a typo ("mames") and a TODO. |
| S20 | MATH (comment) | ✔ | `post_calculations.jl:130-134` | "For a ParallelModel branch, all sub-elements carry the same stress τ_shared": wrong. Sub-elements of a parallel block share **strain rate**, and the *block* carries τ_shared. Line 132 also calls εᵖ the "branch-local **inelastic** strain rate", but it's the total branch strain rate. l. 79 says "plastic-free strain rate". |
| S21 | DRIFT | ✔ | `post_calculations.jl:42-48` vs l. 49-51 | The comment says nested springs are "reconstructed separately"; that's true, but the comment near l. 102 says "This is **proven** to equal the true physical stress" with no reference to a proof. |
| S22 | NIT | ✔ | seeds `0.0` in `_η_KV`, `_η_eff_maxwell`, `_weighted_backstress_scalar`, `_kv_implicit_corrections_scalar` | Hard-coded `Float64` seeds promote `Float32` inputs to `Float64`. |
| S24 | BUG | ✔ confirmed | `branch_strain_rate_mask` / `max_feasible_step` in `solve` | Block strain-rate unknowns are bounded to be non-negative. With a signed scalar ε the physical block strain rates change sign after a load reversal, so the bounded Newton step cannot reach the solution: `S(η₁, P(η₂, G₁), P(η₃, S(η₄, G₂)))` throws `NonConvergenceError` under scalar reversal. |
| S25 | PHYS | ✔ confirmed | `PowerLawViscosity` scalar state functions | `τ^n` and `ε^(1/n)` are evaluated on the signed argument, so a power law in a block cannot follow a signed scalar ε once the block strain rate or stress turns negative: `S(η₁, P(PL, G₁))` throws `DomainError` under scalar reversal. Tensor ε is unaffected, since the block unknowns are invariants. |

### 2.3 In DOC (`strain_rate_correction.md`)

| # | Sev | Status | Location | Issue |
|---|---|---|---|---|
| D1 | DRIFT | ✔ | Intro and "Implementation" | "Before passing the prescribed strain rate to the Newton-Raphson solver… absorb those terms" and "called automatically by `solve` before each Newton step". Only direct Maxwell leafs are pre-corrected. Branch corrections enter the residual implicitly, in **scalar** form (S1). The table leaves out `_implicit_elastic_correction` and `_kv_implicit_*`. |
| D2 | DRIFT | ✔ | Generalized Maxwell, total formula | Presented as the tensor ε_eff the solver sees. The solver actually sees ‖ε + leafcor‖_II + Σ(scalar branch terms), which differs for non-coaxial or reversed τ0 (S1). |
| D3 | DRIFT/BUG | ? | Mixed body example, comments on `elastic_stress_history_2D` | "compute_stress_elastic must look at x[3], not x[1]" and "locate the spring stress at index 3". In the code the nested spring stress is **reconstructed** from τ_shared = x[1] by closed-form inversion (`_branch_elastic_stress`), and x[3] is explicitly skipped (`post_calculations.jl:49-51`). |
| D4 | MATH | ✔ | η*_i cases | η*_i = η_νᵢ/(η_νᵢ + GᵢΔt) assumes one dashpot and one spring per sub-branch (same as T2). SRC uses η_M/η_el over all leafs, and DOC doesn't document that. |
| D5 | NIT | ✔ | Generalized Maxwell | "For a purely viscous branch η* = 1 but τᵒ = 0": a purely viscous branch has **no** τ0 entry at all (it isn't counted in the τ0 tuple). The sentence suggests a zero entry is stored. |
| D6 | NIT | ✔ | Maxwell body | "the elastic backstress divided by the elastic stiffness": divided by 2GΔt, which is not the stiffness G. "harmonic mean" (S12). |
| D7 | NIT | ✔ | KV body | "arithmetic sum": fine, but it clashes with the "arithmetic mean" in SRC comments (S11). |
| D8 | NIT | ✔ | Example | τ₀ is used for the *initial stress* of the analytic solution while τ0/τᵒ means backstress everywhere else. |
| D9 | ✔ OK | ✔ | Example analytic solution | Re-derived: τ₀ = 2η₁η₂ε̇/(η₁+η₂), τ_∞ = 2η₁(η₂+η₃)ε̇/(η₁+η₂+η₃), t_relax = η₃(η₁+η₂)/(G(η₁+η₂+η₃)). **All correct.** |
| D10 | NIT | ✔ | Throughout | British spellings ("generalises", "specialisation", "catalogue") mixed with American ones ("normalization"). |
| D11 | ? | ? | Worked example | It should run (dispatch traced: 1-tuple ε → scalar invariant), but that needs confirming numerically. |

### 2.4 In PDF (`creep.pdf`)

| # | Sev | Status | Location | Issue |
|---|---|---|---|---|
| P1 | ✔ OK | ✔ | Eqs. (1)–(5) | Re-derived step by step; discretization, (3), (4), the substitution, the separation and (5) are all **correct**. |
| P2 | ✔ OK | ✔ | Use-case parameter settings | Spring (α₀=β₀, N=0), dashpot (β₀=0), Maxwell, KV (β₀=0, α₁=β₁), Jeffreys, SLS, generalized Maxwell and Burgers are all consistent with (1)–(5). |
| P3 | NIT | ✔ | p. 1–6 | The symbol `*` is overloaded: **history** (τ*, sᵢ*), **effective viscosity** (η*) and **effective strain rate** (ε̇*). Easy to misread. |
| P4 | NIT | ✔ | p. 4 | "quasi-harmonic mean": it's the inverse of the sum of inverses, not a mean. |
| P5 | NIT | ✔ | p. 6, 8 | Typos: "sightly" → slightly, "equaiton" → equation, "form step 1" → from step 1. |
| P6 | MATH | ✔ | p. 6, (6) | η* and ε̇* aren't given explicitly. From (5): **2η* = 1/α₀ + Σ 1/αᵢ** and **ε̇* = ε̇ + (β₀/α₀ s₀* + Σ βᵢ/αᵢ sᵢ*)/(2η*)**, with s₀* = τ* − Σsᵢ*. Stating these would make the link to the package explicit (see §3). |
| P7 | MATH | ✔ | p. 1, (3) | The reference history (τ* − Σsᵢ*) is correct only if the stored τ* and sᵢ* are mutually consistent. If τ* and sᵢ* are advected or rotated separately (e.g. Jaumann), s₀* is not stored directly and inconsistencies go straight into it. |
| P8 | NIT | ✔ | p. 1 | A rigid element (η, G → ∞) gives αᵢ = 0 in (4)–(5); the degenerate case isn't discussed. |
| P9 | NIT | ✔ | p. 8 | CASE 2 is announced as "to be defined in the final part" (p. 6) and only given as an algorithm, with no equivalent of (5) in terms of the extra unknowns. |

---

## 3. Cross-source consistency

### 3.1 TYP ↔ SRC ↔ DOC (linear, coaxial, same-sign τ0)

| Quantity | TYP | SRC | DOC | Match |
|---|---|---|---|---|
| Maxwell leaf correction | τᵒ/(2GΔt) | τ0/(2·GΔt) (tensor) | τᵒ/(2GΔt) | ✔ |
| KV correction | τᵒ/(2η_KV), η_KV = η₂ + GΔt | ws/(2η_KV), direct spring weight 1, spring adds GΔt to η_KV | same | ✔ |
| Mixed correction | η*τᵒ/(2η_KV), algebraically | η_M/η_el · τ0/(2η_KV) | same | ✔ |
| Generalized | extra factor 2 (**T1**) | Σ η*τ0/(2η_KV) | Σ η*τ0/(2η_KV) | ✘ TYP only |
| η* with >1 spring in sub-branch | undefined (T2) | first spring's G (**S2**) | undefined (D4) | ✘ |
| Tensor vs invariant | tensor | leaf: tensor; branch: sum of invariants (**S1**) | tensor | ✘ SRC |
| When applied | n/a | leaf pre-solve, branch in residual | all pre-solve (**D1**) | ✘ DOC |

### 3.2 PDF ↔ package (✔ verified algebraically)

Map the PDF element to `ParallelModel(SeriesModel(η₀,G₀), SeriesModel(η₁,G₁), …, SeriesModel(η_N,G_N))`. Using
1/αⱼ = 2η_M,j and βⱼ/αⱼ = ηⱼ/(ηⱼ + GⱼΔt) = η*ⱼ, dividing PDF (5) by α₀ and using s₀* = τ* − Σsᵢ* gives

  **2η_KV ε̇ = τ − Σ_{j=0..N} η*ⱼ sⱼ***,  η_KV = Σ_{j=0..N} η_M,j

This is **exactly** the package/TYP/DOC branch relation εᵖ = (τ − ws)/(2η_KV). The PDF's (6) is ε̇* = ε̇ + ws/(2η_KV)
with η* = η_KV. **For linear elements in the tensor form, the two formulations are equivalent.** A spring-only
branch (η → ∞) gives η* = 1 and η_M = GΔt, and a dashpot-only branch (G → ∞) gives η* = 0 and η_M = η. Both
match the SRC special cases.

### 3.3 Important differences (PDF vs package)

| # | Topic | PDF | Package | Consequence |
|---|---|---|---|---|
| X1 | Invariant reduction | Build the full ε̇* **tensor** (including all history), then take ε̇*_II; the PDF explicitly says to redo this every iteration in the nonlinear case | Direct leafs: tensor. Branches: **scalar** sum of ‖τ0‖_II | S1: wrong on non-coaxial loading and **wrong sign on reversal** |
| X2 | History storage | τ* (top-level) + sᵢ* (non-reference branches); s₀* implied | One tensor per spring, all rebuilt coaxial with ε_eff (S6) | Package loses orientation and sign of nested springs |
| X3 | Nonlinear parameters | CASE 1: functions of top-level τ_II (always available). CASE 2: extra unknowns (sᵢ)_II with its own residual | Functions of **strain rate**, evaluated at a shifted branch strain rate (S4), and inconsistently across solve, correction and post-processing (S5) | Neither CASE 1 nor CASE 2; results for nonlinear branch rheology unverified |
| X4 | Residual form | Stress form: r = τ_II − 2η*ε̇*_II | Strain-rate form with extra local unknowns per branch | Equivalent in the linear case; different Newton conditioning |
| X5 | Topology | Parallel block of Maxwell branches (any spring or dashpot can be switched off) in series with other such blocks; arbitrary sequences covered | Arbitrary Series/Parallel trees, but the correction is only valid to one alternation level; S3 silently drops a non-elastic nested parallel; S2 breaks for multi-spring sub-branches | Package allows more topologies than it correctly handles |
| X6 | Viscous-only reduction | Series → inverse sum, parallel → sum | `_η_eff_maxwell` / `_η_KV`, identical | ✔ |
| X7 | Plasticity | Mentions 3-invariant (quadratic τ²) case | Not part of this correction | Out of scope; noted |

---

## 4. Numerical checks

`test/test_elastic_history_tensor.jl` runs 12-step time loops of `solve` + `elastic_stress_history_2D/3D` and compares the outer τ_II and every spring tensor at every step (`rtol = 1e-9`) against `test/reference/tensor_reference.jl`, a monolithic Newton solve on the full tensors (outer τ, block strain rates, sub-branch stresses) with no reduction. The reference is itself checked against the closed-form backward-Euler Maxwell update for a non-coaxial τᵒ.

Load cases: constant pure shear; reversal after 6 steps; pure shear → simple shear; 3D with all six components and a change of direction; scalar ε, constant and reversed.

| Topology | Coaxial / scalar | Scalar reversal | Tensor reversal, rotation, 3D | IDs |
|---|---|---|---|---|
| Maxwell `S(η,G)` | pass | pass | pass | — |
| KV, mixed, Burgers, direct spring + sub-branch | pass | pass | fail (τ_II error 5–870 %) | S1, S6, S7 |
| two blocks | pass | `NonConvergenceError` | fail | S1, S6, S7, S24 |
| two springs in a sub-branch (G_a ≠ G_b) | fail (τ_II error 7 %) | fail | fail | S2 (+ S1) |
| power law as outer leaf, linear block | pass | pass | fail | S1, S6, S7 |
| power law as direct KV element / in a sub-branch | `ArgumentError` (NaN η_KV) | same | same | S5 |

- S1, S6, S7: a signed scalar ε keeps the sign of τ0, so these only show with tensor ε.
- S4: `solve` itself does not throw for the power-law block cases, but the post-solve spring reconstruction does (S5), so the solve result is not compared yet.
- Not covered by these tests: S3 (nested non-elastic parallel inside a sub-branch; the reference rejects that topology), D3, D11.
