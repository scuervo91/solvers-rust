# solvers-rust — pensum (indefinite)

Goal: become **prolific in Rust** while learning **numerical methods** (Chapra) and shipping **petroleum / process engineering CLIs**.

Rhythm: **3 sessions / week**, ~60–90 min each. Never skip the Rust day for “just theory.”

| Day | Name | Focus |
| --- | --- | --- |
| **A** | Chapra lab | Read 1 section → implement or harden a solver in this crate |
| **B** | Rust craft | Same code: ownership, traits, errors, tests, benches, API shape |
| **C** | Domain CLI | Apply solvers to a real workflow (material balance, flash, decline, …) |

Mark progress in `PROGRESS.md` (one line per session).

---

## How each session works (always)

1. **Pick** the next unchecked item in the current Phase below.
2. **Write** failing tests first (Chapra textbook numbers / handmade balances).
3. **Implement** until `cargo test` is green.
4. **Log** in `PROGRESS.md`: date · day A/B/C · what shipped · Chapra § · Rust concept.
5. **Stop** when the checkbox is done — don’t binge three phases in one night.

---

## Phase 0 — Foundations (finish if incomplete) ✅ mostly done

- [x] Root-finding library skeleton (`RootResult`, errors, convergence)
- [x] Bracketing: bisection, false position, modified false position
- [x] Open: Newton–Raphson, secant, modified secant, fixed point
- [x] Polynomials helper
- [ ] Day B: unify public API (`solve_root` trait or enum of methods)
- [ ] Day B: property tests + Chapra example table as `tests/fixtures/`
- [ ] Day C: tiny CLI `solvers roots --method newton --expr ...` (or stdin JSON)

---

## Phase 1 — Linear systems (Chapra: linear algebraic equations)

**A — Chapra:** Gauss elimination, pivoting, LU, Jacobi / Gauss–Seidel, norms & condition.

**B — Rust:** `ndarray` layouts, `Result` for singular matrices, `thiserror` variants, generics over `f64` first then trait `Float` later.

**C — Domain:** **steady material balance** as `A x = b` (species × units).

Deliverable crate path:

```text
src/linear/          # elimination, lu, iterative
bins/mbcli or src/bin/mb.rs
examples/material_balance_mixer.rs
```

CLI sketch (build toward this):

```bash
cargo run --bin mb -- solve --file cases/mixer_two_stream.json
# → component flows out, closure error, method used
```

Case JSON: feeds, specs, unknowns tagged; build `A,b` in Rust; solve; report.

**Done when:** LU + one iterative method + 2 textbook tests + 1 mixer/separator MB case via CLI.

---

## Phase 2 — Nonlinear systems & recycling (Chapra: nonlinear eqs)

**A:** Newton for systems, Jacobian (analytic + finite difference), successive substitution / Wegstein (process engineering classic).

**B:** Traits `NonlinearSystem { fn residual(&self, x) -> Array1; fn jacobian(...) }`, lifetimes on borrowed specs, avoid clones in hot loops.

**C:** **Recycle material balance** (tear stream + Wegstein/Newton). CLI: `mb recycle --file cases/flash_recycle.json`.

**Done when:** tear-stream example converges; CLI prints iteration history (optional `--trace`).

---

## Phase 3 — Curve fit & data (Chapra: least squares / interpolation)

**A:** Linear regression, polynomial fit, spline or Lagrange (pick one), linearization for exponential/power models.

**B:** Polars for CSV in → `ndarray` for math → tidy report struct; feature-gate heavy deps.

**C:** **Decline curve / production fit** CLI: `fit decline --csv wells.csv --model arps` (exp or hyperbolic simplified).

**Done when:** fit + R² / residuals exported; one real-looking synthetic CSV in `cases/`.

---

## Phase 4 — Integration & ODEs (Chapra: numerical integration, IVPs)

**A:** Trapezoid, Simpson, Romberg; Euler, Heun, RK4; stiff taste (qualitative only).

**B:** Iterator-style integrators, `step()` API, criterion benches for RK4 vs Euler.

**C:** **Tank / accumulator material balance dynamic**: `d(VC)/dt = in − out`. CLI: `mb tank --file cases/tank.yaml --t-end 10 --dt 0.1` → CSV time series.

**Done when:** RK4 tank matches analytic exponential washout case within tol.

---

## Phase 5 — Optimization & VLE lite (Chapra: optimization; domain stretch)

**A:** Golden section, Newton opt, gradient descent / simple Nelder–Mead.

**B:** Feature `opt`, clear separation `solvers_rust::opt` vs domain crates.

**C:** Mini **flash / bubble-T** (Raoult + ideal) as nonlinear solve — or **loss minimization** on MB reconciliation (gross error soft start).

**Done when:** one optimization + one thermo-or-recon CLI subcommand.

---

## Phase 6 — Workspace & “prolific” Rust (ongoing every 4–6 weeks)

Restructure when Phase 2–3 feel cramped:

```text
solvers-rust/           # library: roots, linear, nlin, integrate, ode, fit, opt
mb-cli/                 # material balance workflows
fit-cli/                # decline / regression
py-solvers/             # optional later: PyO3 export for dune/Python
```

Rust themes to rotate on **Day B** forever:

1. Ownership & borrowing in numeric loops  
2. Traits & dyn vs static dispatch  
3. Error design (`thiserror` / `anyhow` only in bins)  
4. Testing: unit, integration, proptest, fixtures  
5. `cargo bench` + criterion  
6. CLI: `clap` derive, exit codes, `--format json`  
7. Serde for cases  
8. Logging / tracing  
9. Feature flags & workspace crates  
10. Unsafe-free SIMD later / `rayon` for embarrassingly parallel cases  
11. Docs.rs-quality module docs  
12. PyO3 bridge (link to your Python stack)

---

## Phase ∞ — Repeat with harder Chapra + harder O&G

After Phase 5, loop:

| Cycle | Chapra depth | Domain CLI |
| --- | --- | --- |
| ∞.1 | PDE intro / finite differences | 1D diffusion / pressure diffusivity toy |
| ∞.2 | Stiff ODEs | CSTR / simple kinetics |
| ∞.3 | Sparse linear | Big MB networks |
| ∞.4 | Uncertainty | Monte Carlo on MB inputs |
| ∞.5 | Parallel cases | batch CLI over `cases/*.json` |

Each cycle = ~4–6 weeks at 3 days/week.

---

## Weekly template (copy into calendar)

- **Day A (e.g. Mon):** Chapra §X · implement · tests from book  
- **Day B (e.g. Wed):** Refactor API · docs · benches · clippy/`cargo fmt`  
- **Day C (e.g. Fri):** Domain case + CLI UX · write `cases/…` · short note in `PROGRESS.md`

Miss a day → do a **30‑min Day B** (tests/clippy) so the streak doesn’t die.

---

## North-star projects (creativity targets)

1. **`mb` CLI** — material balance (steady → recycle → dynamic tank)  
2. **`fit` CLI** — decline / lab correlation from CSV  
3. **`flash` CLI** — ideal flash (optional thermo)  
4. **Notebook parity** — same case in Python vs Rust timings (blog for yourself)  
5. **PyO3** — call `solvers_rust` from dune/Python later  

---

## Definition of “prolific”

Not “finished Chapra.” You are prolific when you can, without fear:

- open a new binary crate in the workspace in <15 min  
- express a balance as residuals + pick a solver  
- ship CLI + fixtures + tests the same week  
- read unfamiliar Rust numeric code and change it safely  

---

## References

- Chapra & Canale — *Numerical Methods for Engineers* (your edition’s chapter map)  
- The Rust Book + *Rust by Example* (Day B ammo)  
- `ndarray` book / docs  
- Optional: Himmelblau / Felder for MB problem ideas; Ahmed or Craft for decline forms  

