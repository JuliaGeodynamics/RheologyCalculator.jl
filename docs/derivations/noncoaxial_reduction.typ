#set math.equation(numbering: "(1)")
#set heading(numbering: "1.")

#let ed = $dot(bold(epsilon))$
#let tt = $bold(tau)$
#let to = $bold(tau)^o$
#let II = $"II"$
#let ip(a, b) = $lr(⟨ #a, #b ⟩)$

#align(center, text(16pt, weight: "bold")[Non-coaxial tensor reduction])

This note complements `tensor_reduction.typ`. It covers what the scalar
reduction needs when the history tensors $#to _i$ are not coaxial with the
prescribed strain rate $#ed$: recovering tensors from the scalar solution,
updating the history, nonlinear viscosities, the consistent tangent, and
rotating the stored history.

// = Assumptions

// + *Isotropy.* Every element parameter ($eta$, $G$) is a positive scalar that
//   may depend on the state only through tensor invariants. Anisotropic viscosity
//   or elasticity invalidates everything below.
// + *Time discretization.* Springs use the backward-Euler law
//   $#ed = (#tt - #to) \/ (2 G Delta t)$, with $#to$ expressed in the current
//   frame (see @sec-rotation).
// + *Topology.* A `SeriesModel` of dashpot leafs $l$, spring leafs $k$ and
//   `ParallelModel` blocks $b$. Each block holds dashpots, springs $s$ and
//   `SeriesModel` sub-branches $j$ made of leaf elements only.

= Equations carried over from `tensor_reduction.typ`

For block $b$, with $eta^("eff")_(M,j)$ the series effective viscosity of sub-branch $j$:

$
eta_("KV",b) = sum_l eta_l + sum_s G_s Delta t + sum_j eta^("eff")_(M,j),
quad
eta^star_i = cases(
  1 & "direct spring of the block",
  eta^("eff")_(M,j) \/ (G_i Delta t) & "spring of sub-branch" j,
)
$ <weights>

$
#ed^p_b = (#tt - sum_(i in b) eta^star_i #to _i) / (2 eta_("KV",b)),
$ <block>

$
#ed^("eff") = #ed + bold(H), quad
bold(H) = sum_k #to _k / (2 G_k Delta t)
+ sum_b 1 / (2 eta_("KV",b)) sum_(i in b) eta^star_i #to _i,
$ <eff>

$
#tt = 2 eta^("eff") #ed^("eff"), quad
1 / eta^("eff") = sum_l 1 / eta_l + sum_k 1 / (G_k Delta t) + sum_b 1 / eta_("KV",b) .
$ <constitutive>

Under assumption 1, $eta^("eff") > 0$ is a scalar, so @constitutive makes $#tt$
coaxial with $#ed^("eff")$, *not* with $#ed$. The history tensor $bold(H)$ is a
weighted sum of arbitrarily oriented tensors.

= Scalar equation

Taking the invariant of @constitutive:

$
tau_#II = 2 eta^("eff") (#ed^("eff"))_#II, quad
(#ed^("eff"))_#II = (#ed + bold(H))_#II .
$ <scalar>

// = Tensor recovery

// Given the scalar solution $tau_#II$, the stress tensor follows from @constitutive:

// $
// #tt = tau_#II / (#ed^("eff"))_#II #ed^("eff") = tau_#II bold(n),
// quad
// bold(n) = #ed^("eff") / (#ed^("eff"))_#II,
// $ <recovery>

// for $(#ed^("eff"))_#II > 0$; if $#ed^("eff") = bold(0)$ then $#tt = bold(0)$.
// The direction $bold(n)$ is that of the *effective* strain rate. Using the
// direction of $#ed$ instead is correct only when $bold(H)$ is coaxial with $#ed$.

// The block strain rates then follow from @block:

// $
// #ed^p_b = (tau_#II bold(n) - sum_(i in b) eta^star_i #to _i) / (2 eta_("KV",b)),
// $ <block-recovery>

// which is in general coaxial with neither $#tt$, $#ed$ nor $#ed^("eff")$.

// = History update

// The stress of every spring at the end of the step is the next step's history.
// All springs need their own full tensor.

// $
// #tt^("new")_k = #tt quad "spring leaf of the outer series",
// $ <hist-series>

// $
// #tt^("new")_s = #to _s + 2 G_s Delta t #ed^p_b quad "direct spring of block" b,
// $ <hist-direct>

// $
// #tt^("new")_(j k) = #tt _j = 2 eta^("eff")_(M,j) #ed^p_b + sum_(k') eta^star_(j k') #to _(j k')
// quad "spring" k "of sub-branch" j "of block" b .
// $ <hist-sub>

// All springs of one sub-branch share the stress $#tt _j$, and all spring leafs of
// the outer series share $#tt$. The stresses of a block satisfy the parallel
// balance

// $
// #tt = sum_l 2 eta_l #ed^p_b + sum_s #tt^("new")_s + sum_j #tt _j,
// $

// which follows from @weights and @block and can be used as a check.

// Reconstructing $#tt^("new")_s$ or $#tt^("new")_(j k)$ as a tensor coaxial with
// $#ed^("eff")$ (or with $#ed$), with only its invariant taken from the solve, discards the
// orientation carried by $#to$ and is wrong whenever $#to$ is not coaxial with $#ed^("eff")$.

// == Example: Kelvin-Voigt body under a change of loading direction

// Model `SeriesModel(η₁, ParallelModel(η₂, G))`, 2D Voigt order $(x x, y y, x y)$.
// The spring carries pure-shear history $#to = (s, -s, 0)$, and the new loading is
// simple shear $#ed = (0, 0, e)$. With $eta_("KV") = eta_2 + G Delta t$,
// $1 \/ eta^("eff") = 1 \/ eta_1 + 1 \/ eta_("KV")$ and $a = s \/ (2 eta_("KV"))$:

// $
// #ed^("eff") = (a, -a, e), quad
// (#ed^("eff"))_#II = sqrt(a^2 + e^2) < e + a quad (a, e > 0),
// $

// $
// #tt = 2 eta^("eff") (a, -a, e), quad
// #ed^p = (#tt - #to) / (2 eta_("KV"))
// = ( a (eta^("eff") - eta_("KV")) / eta_("KV"),
//   -a (eta^("eff") - eta_("KV")) / eta_("KV"),
//   eta^("eff") e / eta_("KV") ) .
// $

// Since $eta^("eff") < eta_("KV")$, the normal components of $#ed^p$ have the
// opposite sign to those of $#tt$: the spring unloads its pure-shear stress while
// the block shears. The new spring stress
// $#tt^("new") = #to + 2 G Delta t #ed^p$ has, in general, normal and shear parts in a ratio
// different from that of $#ed^("eff")$ and of $#ed$, so neither direction reproduces it.

// = Nonlinear viscosities

// When element viscosities depend on the state, $eta_("KV",b)$, $eta^star_i$ and
// $eta^("eff")$ are functions of the unknowns, and so is $bold(H)$. Consequently
// both the magnitude *and the direction* of $#ed^("eff")$ change during the
// iteration, and $bold(H)$ must be reassembled as a tensor at every residual
// evaluation. Only the scalar coefficients change; the tensors $#ed$ and
// $#to _i$ are fixed within a step, so each evaluation is a weighted sum of
// fixed tensors followed by one invariant.

// == Viscosities that depend on top-level invariants

// If every viscosity depends only on $tau_#II$, the single unknown is $tau_#II$ and

// $
// r_0(tau_#II) = (#ed + bold(H)(tau_#II))_#II - tau_#II / (2 eta^("eff")(tau_#II)) = 0 .
// $ <r0>

// == Viscosities that depend on the state inside a block

// A viscosity inside block $b$ may depend on a local invariant: the block
// strain-rate invariant $(#ed^p_b)_#II$ for a dashpot that is a direct leaf of the block,
// or the sub-branch stress invariant $(#tt _j)_#II$ for an element of sub-branch $j$
// (a dashpot $l$ in sub-branch $j$ has strain-rate invariant $(#tt _j)_#II \/ (2 eta_l)$).
// These invariants are not functions of $tau_#II$ alone in closed form, so they become
// additional unknowns. With the unknown vector

// $
// bold(x) = (tau_#II, {e_b}, {t_j}), quad e_b tilde (#ed^p_b)_#II, quad t_j tilde (#tt _j)_#II,
// $

// all viscosities are evaluated at $bold(x)$, the tensors are rebuilt from $bold(x)$ via
// @eff, @recovery, @block-recovery and @hist-sub, and the residuals are

// $
// r_0(bold(x)) = (#ed + bold(H)(bold(x)))_#II - tau_#II / (2 eta^("eff")(bold(x))),
// $ <rx0>

// $
// r_b (bold(x)) = e_b - (tau_#II bold(n)(bold(x)) - sum_(i in b) eta^star_i (bold(x)) #to _i)_#II / (2 eta_("KV",b)(bold(x))),
// $ <rxb>

// $
// r_j (bold(x)) = t_j - (2 eta^("eff")_(M,j)(bold(x)) #ed^p_b (bold(x)) + sum_k eta^star_(j k)(bold(x)) #to _(j k))_#II .
// $ <rxj>

// Each invariant is taken of an assembled tensor. The scalar shortcuts

// $
// e_b = (tau_#II - sum_i eta^star_i (#to _i)_#II) / (2 eta_("KV",b)), quad
// t_j = 2 eta^("eff")_(M,j) e_b + sum_k eta^star_(j k) (#to _(j k))_#II
// $

// are exact only for positively coaxial histories and, like the sum-of-invariants
// form of @scalar, give the wrong sign under load reversal.

// The derivative of an invariant of an assembled tensor, needed for the Jacobian, is

// $
// (partial bold(A)_#II) / (partial x) = #ip($bold(A)$, $(partial bold(A)) / (partial x)$) / bold(A)_#II,
// quad bold(A)_#II > 0 .
// $ <dinv>

// = Consistent tangent

// Consider the case of @r0. Differentiating @constitutive with respect to
// $#ed$, with $bold(H)$ and $eta^("eff")$ depending on $#ed$ only through $tau_#II$:

// $
// d #tt = 2 eta^("eff") d #ed + (2 eta^("eff") bold(H)' + 2 eta^("eff")' #ed^("eff")) d tau_#II,
// $

// where $(dot)' = d(dot) \/ d tau_#II$. Differentiating @r0 with @dinv and
// $g = d(tau_#II \/ (2 eta^("eff"))) \/ d tau_#II$:

// $
// d tau_#II = #ip($bold(n)$, $d #ed$) / (g - #ip($bold(n)$, $bold(H)'$)) .
// $

// Hence

// $
// d #tt = 2 eta^("eff") d #ed + bold(v) #ip($bold(n)$, $d #ed$), quad
// bold(v) = (2 eta^("eff") bold(H)' + 2 eta^("eff")' #ed^("eff")) / (g - #ip($bold(n)$, $bold(H)'$)) .
// $ <tangent>

// The tangent is the isotropic part $2 eta^("eff") bb(I)$ plus a rank-one update
// $bold(v) ⊗ bold(n)$ acting along the *effective* strain-rate direction:

// - linear elements ($bold(H)' = bold(0)$, $eta^("eff")' = 0$): $d #tt = 2 eta^("eff") d #ed$;
// - nonlinear, no history in the blocks ($bold(H)' = bold(0)$): $bold(v) parallel bold(n)$, a symmetric rank-one update along $bold(n)$;
// - nonlinear with block history ($bold(H)' != bold(0)$): $bold(v)$ has a component
//   along the $#to _i$, so the tangent is *not symmetric* in general.

// A tangent built from the direction of $#ed$, or one assuming $#tt parallel #ed$,
// is therefore incorrect whenever $bold(H)$ is not coaxial with $#ed$ and the viscosities are nonlinear.

// = Rotation of the stored history <sec-rotation>

// The backward-Euler spring law compares tensors from two instants, so $#to$
// must be expressed in the current frame before it enters @eff. For a material
// spin $bold(W)$ over the step, a first-order (Jaumann) update is

// $
// #to = #tt^n + Delta t (bold(W) #tt^n - #tt^n bold(W)),
// $

// or, with the incremental rotation $bold(R)$, $#to = bold(R) #tt^n bold(R)^T$.
// The same rotation must be applied to *every* stored spring history, because
// mixing rotated and unrotated histories in @eff corrupts both $bold(H)$ and the
// balance checks above. Under pure rotation the invariants $(#to _i)_#II$ are
// unchanged but their orientation relative to $#ed$ is not, so this effect is
// invisible to any purely scalar treatment of the history.
