#set math.equation(numbering: "(1)")

#let ed = $dot(bold(epsilon))$
#let tt = $bold(tau)$
#let to = $bold(tau)^o$
#let GDt = $G Delta t$

= Tensor reduction operations

== Notation and element laws

$#ed$ is the deviatoric strain-rate tensor, $#tt$ the deviatoric stress tensor, and
$#to$ the stress carried by a spring at the previous time step. All tensors are
symmetric and deviatoric; $(dot)_("II")$ denotes the second invariant.

A linear dashpot of viscosity $eta$ obeys $#ed = #tt \/ (2 eta)$. A spring of
shear modulus $G$, discretized with backward Euler, obeys

$
#ed = (#tt - #to) / (2 #GDt),
$

so a spring behaves as a dashpot of viscosity $#GDt$ plus a history term. Its
*effective viscosity* is $eta = #GDt$. Effective viscosities combine as

$
eta_("series") = (sum_l 1 / eta_l)^(-1), quad
eta_("parallel") = sum_l eta_l .
$

In what follows, the *effective strain rate* $#ed^("eff")$ is the prescribed
strain rate with every history term moved to the left-hand side, so that the
right-hand side depends on $#tt$ only.

== Maxwell body

Model: `SeriesModel(η, G)`.

Strain rates add in series:

$
#ed = 1/(2 eta) #tt + (#tt - #to)/(2 #GDt).
$

Moving the history term to the left-hand side:

$
#ed^("eff") = #ed + #to / (2 #GDt) = #tt (1/(2 eta) + 1/(2 #GDt)) = #tt / (2 eta^("eff")_M),
$ <maxwell-eff>

with the Maxwell effective viscosity

$
eta^("eff")_M = (1/eta + 1/(#GDt))^(-1).
$

== Kelvin-Voigt body

Model: `SeriesModel(η₁, ParallelModel(η₂, G))`.

The outer series splits the strain rate into the dashpot $eta_1$ and the
parallel block, whose strain rate is $#ed^p$:

$
#ed = 1/(2 eta_1) #tt + #ed^p .
$

Stresses add in the parallel block, and both elements carry the strain rate $#ed^p$:

$
#tt = 2 eta_2 #ed^p + (2 #GDt #ed^p + #to) = 2 (eta_2 + #GDt) #ed^p + #to .
$

With the Kelvin-Voigt effective viscosity $eta_("KV") = eta_2 + #GDt$:

$
#ed^p = (#tt - #to) / (2 eta_("KV")) .
$

Substituting into the series equation and moving the history term to the left-hand side:

$
#ed^("eff") = #ed + #to / (2 eta_("KV")) = #tt (1/(2 eta_1) + 1/(2 eta_("KV"))) .
$

== Mixed Kelvin-Voigt and Maxwell body

Model: `SeriesModel(η₁, ParallelModel(η₂, SeriesModel(η₃, G)))`.

$
#ed = 1/(2 eta_1) #tt + #ed^p, quad
#tt = 2 eta_2 #ed^p + #tt^("M"),
$

where $#tt^("M")$ is the stress in the Maxwell sub-branch, whose strain rate is also $#ed^p$:

$
#ed^p = #tt^("M") / (2 eta_3) + (#tt^("M") - #to) / (2 #GDt)
= #tt^("M") / (2 eta^("eff")_M) - #to / (2 #GDt),
quad
eta^("eff")_M = (1/eta_3 + 1/(#GDt))^(-1).
$

Solving for the sub-branch stress:

$
#tt^("M") = 2 eta^("eff")_M (#ed^p + #to / (2 #GDt)) .
$

Substituting into the parallel stress balance:

$
#tt = 2 (eta_2 + eta^("eff")_M) #ed^p + eta^("eff")_M / (#GDt) #to
= 2 eta_("KV") #ed^p + eta^star #to,
$

with

$
eta_("KV") = eta_2 + eta^("eff")_M, quad
eta^star = eta^("eff")_M / (#GDt) = eta_3 / (eta_3 + #GDt) .
$

Hence

$
#ed^p = (#tt - eta^star #to) / (2 eta_("KV")),
$

and substituting into the series equation:

$
#ed^("eff") = #ed + (eta^star #to) / (2 eta_("KV")) = #tt (1/(2 eta_1) + 1/(2 eta_("KV"))) .
$

The weight $eta^star in (0, 1)$ tends to $1$ for a soft spring ($G -> 0$) and to
$0$ for a stiff spring ($G -> infinity$), in which case the sub-branch acts as the
dashpot $eta_3$ alone.

== Maxwell sub-branch with several elements

Consider a `SeriesModel` sub-branch $j$ of a parallel block, made of dashpots
and springs $k$ with moduli $G_(j k)$ and histories $#to _(j k)$. All its elements carry
the sub-branch stress $#tt _j$, and their strain rates add up to $#ed^p$:

$
#ed^p = sum_("dashpots") #tt _j / (2 eta_l)
+ sum_k (#tt _j - #to _(j k)) / (2 G_(j k) Delta t)
= #tt _j / (2 eta^("eff")_(M,j)) - sum_k #to _(j k) / (2 G_(j k) Delta t),
$

where $eta^("eff")_(M,j)$ is the series effective viscosity of *all* elements in the sub-branch,
springs included through $eta = G_(j k) Delta t$. Solving for the sub-branch stress:

$
#tt _j = 2 eta^("eff")_(M,j) #ed^p + sum_k eta^star_(j k) #to _(j k),
quad
eta^star_(j k) = eta^("eff")_(M,j) / (G_(j k) Delta t) .
$ <subbranch>

Each spring is weighted by its *own* modulus. For a single dashpot $eta$ and a
single spring $G$ this reduces to $eta^star = eta \/ (eta + #GDt)$.

== Generalized parallel block

A `ParallelModel` block whose elements are dashpots $l$, springs $s$ (direct
leafs) and `SeriesModel` sub-branches $j$ all share the strain rate $#ed^p$,
and their stresses add. A direct spring carries $2 G_s Delta t #ed^p + #to _s$;
a sub-branch carries @subbranch. Therefore

$
#tt = 2 eta_("KV") #ed^p + sum_i eta^star_i #to _i,
$

where the sum runs over every spring $i$ in the block, and

$
eta_("KV") = sum_l eta_l + sum_s G_s Delta t + sum_j eta^("eff")_(M,j),
$

$
eta^star_i = cases(
  1 & "direct leaf of the block",
  eta^("eff")_(M,j) \/ (G_i Delta t) & "spring of sub-branch" j .
)
$

Dashpots contribute to $eta_("KV")$ but carry no history. Inverting:

$
#ed^p = (#tt - sum_i eta^star_i #to _i) / (2 eta_("KV")) .
$ <block>

The derivation assumes that each sub-branch contains only leaf elements (no
further nested parallel blocks).

== General series composite

Consider a `SeriesModel` with dashpot leafs $l$, spring leafs $k$ and parallel
blocks $b$. All carry the stress $#tt$, and strain rates add:

$
#ed = sum_l #tt / (2 eta_l) + sum_k (#tt - #to _k) / (2 G_k Delta t) + sum_b #ed^p_b .
$

Substituting @block for each block and moving every history term to the left-hand side:

$
#ed^("eff")
= #ed
+ sum_k #to _k / (2 G_k Delta t)
+ sum_b 1 / (2 eta_("KV",b)) sum_(i in b) eta^star_i #to _i
= #tt / (2 eta^("eff")),
$ <series-eff>

with

$
1 / eta^("eff") = sum_l 1 / eta_l + sum_k 1 / (G_k Delta t) + sum_b 1 / eta_("KV",b) .
$

The Maxwell, Kelvin-Voigt, and mixed bodies above are special cases of @series-eff.

// == Reduction to a scalar equation

// @series-eff is a tensor equation with a scalar coefficient, so $#tt$ is coaxial
// with $#ed^("eff")$, and taking the second invariant of both sides gives the
// equivalent scalar equation

// $
// tau_("II") = 2 eta^("eff") dot.op (#ed^("eff"))_("II") .
// $

// The invariant must be taken of the *full* effective strain-rate tensor, i.e. of the sum in @series-eff:

// $
// (#ed^("eff"))_("II") != (#ed)_("II") + sum_k (#to _k)_("II") / (2 G_k Delta t) + dots
// $

// in general. The two sides agree only when every history tensor $#to _i$ is
// coaxial with $#ed$ *and* has the same sign. Non-coaxial loading makes the
// right-hand side too large, and load reversal ($#to$ opposite to $#ed$) gives the
// history terms the wrong sign.

// When the viscosities are nonlinear, $eta^("eff")$, $eta_("KV",b)$ and $eta^star_i$
// depend on the state, and @series-eff keeps its form only if they are evaluated at
// quantities available at the top level (e.g. $tau_("II")$ or $(#ed^("eff"))_("II")$).
// Viscosities that depend on the stress or strain rate *inside* a block require the
// block's local invariants as additional unknowns.

// == Newton iteration

// With $f(tau_("II")) = tau_("II") \/ (2 eta^("eff"))$, the scalar residual and its Newton update are

// $
// r(tau_("II")) = (#ed^("eff"))_("II") - f(tau_("II")) = 0,
// quad
// tau_("II")^(n+1) = tau_("II")^n - r / J,
// quad
// J = (partial r) / (partial tau_("II")) = - (partial f) / (partial tau_("II")) .
// $
