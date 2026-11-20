# Schatzman–Paoli time-stepping scheme (`SchatzmanPaoliOSI`)

This note describes the Schatzman–Paoli scheme of the kernel, how it maps onto the code,
the bugs fixed in November 2026 and how the scheme was validated.

References: M. Schatzman, *A class of nonlinear differential equations of second order
in time*, Nonlinear Analysis 2 (1978); L. Paoli, M. Schatzman, *A numerical scheme for
impact problems I & II*, SIAM J. Numer. Anal. 40 (2002); V. Acary, B. Brogliato,
*Numerical Methods for Nonsmooth Dynamical Systems*, Springer (2008), §10.

---

## 1. Problem and scheme

Linear time-invariant Lagrangian system (`LagrangianLinearTIDS` only) with unilateral
constraints and a Newton impact law:

```
M q'' + C q' + K q = F_ext(t) + Gᵀ λ,      g(q) = H q + b ≥ 0,   0 ≤ g(q) ⊥ λ ≥ 0
```

The scheme is a two-step (central difference) scheme on the positions. The impact law
is not written on the velocities but on the positions: the constraint is imposed on a
combination of the new and the old gaps. With `Δ²q = q_{k+1} − 2 q_k + q_{k−1}`:

```
M Δ²q + h² C (θ (q_{k+1} − q_{k−1}) / (2h) + (1 − θ) v_k)
      + h² K (θ² q_{k+1} + (1 − θ²) q_k)
      − h² (θ F_ext(t_{k+1}) + (1 − θ) F_ext(t_k))  =  Gᵀ r_{k+1}

0 ≤ g(q_{k+1}) + e g(q_{k−1})  ⊥  r_{k+1} ≥ 0

v_{k+1} = (q_{k+1} − q_{k−1}) / (2h)
```

For a ball arriving on the ground (`g(q_{k−1}) > 0`, `g(q_{k+1}) = −e g(q_{k−1}) < 0`),
the restitution is obtained on the positions over two steps: this is the discrete
counterpart of `g'(t⁺) = −e g'(t⁻)`. The scheme lets the gaps become slightly negative,
by O(h) at an impact.

## 2. Mapping onto the code

| step | formula | code (`SchatzmanPaoliOSI.cpp`) |
|---|---|---|
| iteration matrix | `W = M + (hθ/2) C + h² θ² K` | `initializeIterationMatrix` |
| free residual | left-hand side above at `q_{k+1} = q_k` | `computeResidu` (`RESIDU_FREE`) |
| free position | `q_free = q_k − W⁻¹ RESIDU_FREE` | `computeFreeState` (`FREE`) |
| LCP | `w = H W⁻¹ Hᵀ r + g(q_free) + e g(q_{k−1})` | `computeFreeOutput`, `_NSLEffectOnFreeOutput` |
| update | `q_{k+1} = q_free + W⁻¹ Hᵀ r`, `v_{k+1}` centred | `updateState` |

Notes:

- The LCP is written on index set 0 (all the constraints are always in the LCP).
- `g(q_{k−1})` is read in the memory of the output: `Interaction::yMemory(0)`, entry 1
  (entry 0 is `g(q_k)`). The memories have depth `_steps = 2`.
- r is stored in `p[0]`, `λ = r` at level 0.
- `θ` weighs the external forces, the damping and the stiffness. `θ = 1/2` gives a
  second-order scheme.

### Start of the two-step scheme

At `t_0`, the position memory must hold `(q_0, q_{−1})`. `initializeWorkVectorsForDS`
sets

```
q_{−1} = q_0 − h v_0 + h²/2 M⁻¹ (F_ext(t_0) − K q_0 − C v_0)
```

so that the first free step gives the second-order Taylor expansion of `q_1` (exact for
constant forces). The output memory is filled with `g(q_0)` twice. This is only an
approximation of `g(q_{−1})`, and it matters only if a contact is active at `t_0` with
`e > 0`.

### Memory of the output

`TimeStepping::initializeNewtonSolve` calls `OneStepIntegrator::updateAndSwapAllOutput`
at the beginning of every step, after the event of the previous step has pushed the
outputs in memory (`TimeDiscretisationEvent::process` → `pushInteractionsInMemory`). For
a one-step integrator (memory of depth 1), pushing the output again is harmless. For the
two-step integrator, it would replace `g(q_{k−1})` by `g(q_k)`. Once the memory holds
`_steps` values, `updateAndSwapAllOutput` therefore only refreshes the latest entry
(`Interaction::updateLatestInMemory`).

## 3. Bugs fixed (November 2026)

1. **One-step shift at the start.** The memories were started with
   `(q_{−1}, q_0) := (q_0, q_0 + h v_0)`: the whole trajectory was one step ahead (the
   value computed at `t_k` was the one of `t_{k+1}`).
2. **Impact law on `g(q_k)` instead of `g(q_{k−1})`.** The output was pushed three times
   per step (once by the event, twice by `updateAndSwapAllOutput`), so both entries of
   the memory held `g(q_k)`. The scheme solved `0 ≤ g(q_{k+1}) + e g(q_k) ⊥ r ≥ 0`, which
   is not consistent with the impact law: the restitution was wrong.
3. **Stiffness replaced by the identity.** `computeResidu` used
   `h² hasStiffnessMatrix() q_k` (a bool) instead of `h² K q_k` (regression of the
   refactoring; master used K).

Bugs 1 and 2 were found by comparing with hots-python. They were present on master too
(the reference of `BouncingBallTS-SchatzmanPaoliOSI` reproduced them).

## 4. Validation

Python script: `kernel/src/simulation/test/validate_with_hots.py` (run the CppUnit test
`TimeSteppingSchemesTest` with `SICONOS_TEST_DUMP_TRAJECTORIES=1` first).

| case | reference | max \|q − q_ref\| |
|---|---|---|
| `beadsSchatzmanPaoli` (3 beads, 3 contacts, e = 0.8) | `hots.paoli_schatzman.paoli_schatzman_fixed` | 5.6e-12 |
| `springChainSchatzmanPaoli` (5 masses, K, C, 2 contacts) | `paoli_with_kc` (the formulas of §1 in Python) | 3.2e-15 |
| `BouncingBallTS-SchatzmanPaoliOSI` (tutorials), t ≤ 6 s | `paoli_schatzman_fixed` | 1.4e-13 |
| same, 6 s < t ≤ 10 s (ball at rest, \|q\| ≈ 1e-11) | same | 5.8e-7 (round-off in the micro-bounces) |

The hots scheme (a port of the Scilab `PaoliLinearOnestep`) starts with one Moreau step.
That start is the same as the Taylor start for constant forces. Its stiffness and damping
terms differ from §1: `h² K (θ q_{k+1} + (1 − θ) q_k)` and a damping term in
`(q_{k+1} + q_k)`. On `springChainSchatzmanPaoli` before the first impact (t < 0.0525), the
error against the exact solution (matrix exponential) is:

- Siconos (§1): 7.3e-5, which is O(h²) for h = 2.5e-3;
- hots: 3.5e-3.

So hots is not used as the reference for systems with K or C.

Kernel refs: `beadsSchatzmanPaoli.ref` (regenerated), `springChainSchatzmanPaoli.ref`
(new). Tutorials: `BouncingBallTS-SchatzmanPaoliOSI.ref` (regenerated).
