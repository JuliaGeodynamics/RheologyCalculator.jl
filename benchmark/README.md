# Performance experiments

Run from the repository root with an isolated, lightweight environment:

```sh
julia --project=benchmark -e 'using Pkg; Pkg.instantiate()'
julia --project=benchmark --threads=1 benchmark/jacobian.jl baseline
julia --project=benchmark benchmark/report.jl benchmark/results/BASELINE.toml
```

`jacobian.jl` writes a new timestamped TOML result; it never overwrites a run.
Keep baseline and experiment results under version control. The environment's
Manifest is ignored; results record its hash and key dependency versions.
`RC_BENCH_SECONDS` (default 0.15) and `RC_BENCH_ROUNDS` (default 3) control duration.
Use the same values for comparable runs. Keep other CPU-heavy work idle.
On Windows, optionally pin only the benchmark process to a logical CPU (mask 4
selects logical CPU 2) before collecting a baseline and all its candidates:

```powershell
$env:RC_BENCH_AFFINITY = '4'
julia --project=benchmark --threads=1 benchmark/jacobian.jl baseline_pinned
```

The mask is recorded and compared. It does not change machine-wide settings.
On Linux use an external affinity tool such as `taskset` consistently instead.

Before each optimization, capture a baseline on the same machine/environment.
After the change, run the unchanged harness with a descriptive experiment label:

```sh
julia --project=benchmark --threads=1 benchmark/jacobian.jl experiment_name
julia --project=benchmark benchmark/report.jl benchmark/results/BASELINE.toml benchmark/results/CANDIDATE.toml
```

Save that Markdown comparison in `results/` and add a dated entry to
`EXPERIMENTS.md`: hypothesis, exact change, baseline/candidate filenames,
validation, improvements, regressions, and keep/reject decision. Repeat promising
or noisy results in a fresh process. Record rejected approaches too. A new
fixture or harness change requires a new baseline; the comparison tool rejects
different environments, harnesses, and fixture/kernel coverage.

The initial runs still show absolute timing drift. When a proposed gain is small
relative to unchanged-code repeat variability, add the old and new algorithms to
the same interleaved comparison (including full solves for solver changes) before
claiming an improvement. A cross-process percentage alone is not sufficient.

The harness measures residuals, Jacobians, their separate combined cost, static
backsolves, cold-start and warm-start solves. Three randomized rounds help expose
ordering/drift effects. The five Jacobian variants use Chairmarks' interleaved,
randomized comparative sampling within each round, with five times the per-kernel
time budget. Each result contains within-round quartiles, minima,
sample counts, and allocations; reports use the median of round medians.
Inputs are loaded through a Ref inside the timed function to resist constant
folding. Compilation and fixture setup are excluded from steady-state timings.
First preparation time is diagnostic only: it includes compilation plus setup,
and is affected by compilation shared with preceding fixtures.

Full solves now also use interleaved sampling: original/current, cold/warm.
`reference_solver.jl` freezes the original solve and line-search implementation;
the harness verifies exact solution, residual, and iteration-count agreement.
It shares lower-level residual/Jacobian functions with the package, so a future
experiment changing those functions needs an appropriate additional reference.
The reference file is included in the harness fingerprint. Single-run reports
show full-solve improvements and regressions against this frozen reference.

Current preprocessing experiments are confined to the benchmark harness:
cached callable, cached full-width config, cached chunk-1 config, and precomputed
equation metadata. No cache is introduced into production `solve`. Config and
equation construction have separate timings. Prepared residual code mirrors the
production pipeline and is checked at initial, converged, and perturbed points.
Every fixture must converge; numerical/shape checks fail the run on a mismatch.
Cold-start guesses are `0.7 * initial_guess_x(...)`, since the generated guess
already solves some simple models exactly. Full-width and chunk-1 config paths
also have an untimed assertion that each evaluates the residual once.

This initial suite covers 10 scalar, volumetric, nested, and plastic fixtures.
It does not yet cover cap models, a GPU runtime, or a full time-stepping workload.
Kernel speedups are hypotheses until integration improves complete solves.
Separate tolerances and convergence settings must be recorded if introduced.
