# Changelog

## Unreleased

### Changed (results)

- Elastic history in `ParallelModel` blocks is carried as a tensor. The global
  deviatoric equation is posed on the second invariant of the assembled strain
  rate `ε + H`, where `H` is the history tensor of
  `effective_strain_rate_correction`. Previously the block terms entered as a
  sum of invariants, which was wrong for histories not coaxial with `ε` and
  had the wrong sign after a load reversal. Results change for tensor-`ε`
  models with springs in blocks; scalar-`ε` results and models whose springs
  are all outer-series leafs are unchanged.
- Each spring of a sub-branch is weighted by its own modulus, so sub-branches
  with several springs of different stiffness give correct results.
- Spring stresses inside blocks are recovered as full tensors, not rebuilt
  coaxial with `ε + H`. `elastic_stress_history_2D/3D(c, x, ε, τ0, others)`
  return them.
- Elements whose viscosity depends on the state (for example
  `PowerLawViscosity`) inside a block that carries elastic history are
  evaluated at their physical strain rate or stress. Such blocks add closure
  unknowns to the solver vector, listed by `x_keys` as `:εII_block` and
  `:τII_subbranch`. Models without such elements keep their layout.
- `component_partition` and `dissipation_partition` evaluate block elements at
  the physical block strain rate and sub-branch stresses.

### Added

- `elastic_stress_history(c, x, ε, τ0, others)`: dimension-agnostic history
  update, one tensor (or signed scalar) per spring.
- `compute_stress_elastic(c, x, ε, others)`: invariants of the recovered spring
  stresses.
- `effective_strain_rate_correction(c, x, ε, others)`: the history tensor at a
  solver vector.
- `viscosity_depends_on_state(r)`: trait that decides which elements of a block
  need closure unknowns. It defaults to `true`; a user rheology whose viscosity
  is independent of the stress and strain rate may return `false`.

### Breaking

- `compute_stress_elastic(c, x, others)` and `elastic_stress_history_2D/3D(c,
  τII, ε, τ0, others)` throw for a composite with springs inside blocks and a
  tensor history, since those spring stresses depend on `ε`; use the forms
  that take `x` and `ε`.
- `effective_strain_rate_correction(c, ε, τ0, others)` throws when a block that
  carries elastic history holds an element whose viscosity depends on the
  state; use the form that takes the solver vector `x`.
- Composites whose blocks carry elastic history and nest a composite inside a
  sub-branch throw an error.
- The solver vector of a model with a state-dependent element in a block that
  carries elastic history has extra entries.
