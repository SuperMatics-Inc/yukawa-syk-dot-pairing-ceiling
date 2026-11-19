# Replay: A sharp bound on pair formation in the Yukawa–SYK dot

This bundle accompanies the paper by Davis Rens. The main theorem concerns the centered causal normal saddle at equal electronic and phonon flavor counts, with the large-N limit taken first. The temperature is the highest linear pairing onset over eligible saddles; it is not an exclusion of first-order transitions. Esterlis and Schmalian derived the limiting kernel and the factor 3/8; this proof excludes any overshoot at finite coupling. The flavor-interval extension is specified below.

Use Python 3.13 and an environment containing the packages below. From this directory:

```sh
uv run --no-project --with python-flint==0.9.0 --with numpy==2.5.5 --with scipy==1.18.1 python -B run.py
```

The driver pins the covered sources by SHA-256, uses the exact rational rectangle endpoints and checks the leaf counts. It checks 17 expected refusals: one endpoint normalization, and five wide-guard and three finite-cover mutations for each of the two covers. A full replay takes several minutes on one CPU. Do not use Python's optimization mode. There are no private input paths. Dependencies need to be installed only if not already present. The floating trial in `two_tight.py` is not used as proof data until it has passed interval checks.

| File | What it checks |
|---|---|
| `cover.py` | Whole finite-width strip, by Arb interval enclosures, analytic normal tails, and positive head Schur witnesses |
| `replay.py` | Finite inequalities in the analytic wide-width estimate, with an infinite-tail gamma-two residual |
| `two_tight.py`, `gamma-two.md` | Separate tight enclosure for gamma two, depending on the stated analytic tail estimates, self-adjoint spectral theory and Arb |
| `FiniteCore.lean` | Optional Lean check of finite rational folding, mass, and weak-width identities |

For Lean, from a mathlib project configured with `leanprover/lean4:v4.35.0`, run `lake env lean /path/to/replay/FiniteCore.lean`. The Lean check is not a formalization of the infinite-dimensional arguments or the Arb program. To build the manuscript locally, enter `paper/` and run `tectonic main.tex`.

The finite-cover SHA-256 is `7deea3c26b7eff520888693176c0ae57b8b9b3bd47a0024fab3a8a2989d57151`. In the replay the canonical five argument sets (column order exactly as shown) are:

| `--sigma` | `--p` | `--n` | `--k` | `--passes` | `--depth` | Mass / weak / empty / pairing |
|---|---|---:|---:|---:|---:|---|
| `15/32 4` | `0 33` | 12 | 6 | 3 | 28 | 35 / 3 / 0 / 35 |
| `4 16` | `0 130` | 24 | 10 | 6 | 28 | 98 / 0 / 0 / 163 |
| `16 40` | `0 324` | 32 | 12 | 6 | 28 | 235 / 0 / 0 / 381 |
| `40 100` | `0 810` | 40 | 16 | 6 | 28 | 647 / 0 / 0 / 1155 |
| `100 180` | `0 1457` | 48 | 20 | 6 | 31 | 1050 / 0 / 0 / 1592 |

Run the small focus with `python -B cover.py --sigma 99 100 --p 3 16/5 --n 40 --k 16 --passes 6 --depth 16`: it must pass at 8 / 0 / 0 / 14. Add `--mutant crossing-half` or `--mutant inflate-tail` to produce an uncovered cell, or `--mutant double-polarization` to fail the independent signed-shift normalization anchor. A mutant rejection checks that a specific change is detected, not that every possible coding error is detectable.

## Closed flavor interval

`flavor/RIGOR.md` gives the changes for the full saddle-coefficient interval \(99/50\le r=b/a\le101/50\). Its references to the earlier wide-width and finite-cover arguments correspond to Appendices A and B of `../paper/main.tex`, respectively; the gamma-two tail input is `gamma-two.md` in this directory. The sharp limiting coefficient is \(C_\infty\sqrt{3/(4r)}\) at each fixed ratio, and the bound covers all finite positive couplings under the same centered, causal large-\(N\) hypotheses. The closed interval is in \(r\), not the number ratio before including the spin trace.

| `--sigma` | `--p` | `--n` | `--k` | `--passes` | `--depth` | Mass / weak / empty / pairing |
|---|---|---:|---:|---:|---:|---|
| `15/32 4` | `0 33` | 16 | 8 | 4 | 28 | 32 / 3 / 0 / 36 |
| `4 16` | `0 131` | 28 | 12 | 6 | 29 | 99 / 0 / 0 / 160 |
| `16 40` | `0 327` | 36 | 14 | 7 | 30 | 230 / 0 / 0 / 394 |
| `40 100` | `0 818` | 56 | 20 | 8 | 32 | 608 / 0 / 0 / 1096 |
| `100 180` | `0 1472` | 64 | 24 | 8 | 34 | 1078 / 0 / 0 / 1619 |

The focus command is `python -B flavor/cover.py --sigma 99 100 --p 3 16/5 --n 56 --k 20 --passes 8 --depth 16` (12 / 0 / 0 / 14). The same three cover mutations are refused. The driver also runs `flavor/replay.py` and challenges all five wide-width mutations. For Lean, under the toolchain above, run `lake env lean /path/to/replay/flavor/FiniteCore.lean`. This checks the extra rational flavor guards and the normalization with symbolic \(r\).
