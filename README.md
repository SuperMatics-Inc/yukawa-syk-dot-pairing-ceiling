# Replay supplement for A sharp bound on pair formation in the Yukawa-SYK dot

Davis Rens, SuperMatics Inc. ([ORCID](https://orcid.org/0000-0001-5839-5197)).

This repository accompanies the preprint [A sharp bound on pair formation in the Yukawa-SYK dot](https://doi.org/10.5281/zenodo.24371778). The replay supplement is [10.5281/zenodo.24371782](https://doi.org/10.5281/zenodo.24371782). The source repository is [yukawa-syk-dot-pairing-ceiling](https://github.com/SuperMatics-Inc/yukawa-syk-dot-pairing-ceiling). The paper's interval enclosure for the gamma-two constant relies on the separate [gamma-model replay certificate](https://doi.org/10.5281/zenodo.24366341). It cites the related [continuum Yukawa-SYK bound](https://doi.org/10.5281/zenodo.24366469), which treats a different saddle.

## Scope and evidence

The paper proves a strict upper bound on the strongest linear pairing instability at every positive finite coupling in the centered, causal, equal-flavor large-N Yukawa-SYK dot saddle with real couplings and full phonon polarization. The large-N limit is taken first. The conclusion is not a bound on a first-order transition, a noncentered saddle, a finite-N transition, or a material.

The finite rational identities and interval enclosures are certified by the replay and its rejection checks. The gamma-two enclosure is conditional on the analytic tail estimates stated in `replay/gamma-two.md`. The Lean files check finite rational lemmas only; they do not formalize the infinite-dimensional proof or the interval calculations.

| Result | Label | Evidence |
|---|---|---|
| Strict finite-coupling upper bound for the stated centered, causal large-N saddle | PROVED | `paper/main.tex` and its stated hypotheses |
| Finite rational identities, interval covers, and expected refusal paths | CERTIFIED | `replay/run.py`, the pinned replay sources, and the CI workflow |
| Finite algebra used in the rational folding and flavor guards | CERTIFIED | `replay/FiniteCore.lean` and `replay/flavor/FiniteCore.lean` |
| Gamma-two constant enclosure | CERTIFIED, conditional on the stated analytic tail estimates | `replay/two_tight.py` and `replay/gamma-two.md` |

## Contents and replay

- `paper/` contains the manuscript source and bibliography.
- `replay/` contains the exact-rational identities, interval certificates, mutation checks, and finite Lean sources. Its README records the rational strip counts and the separate focus checks.
- `lean/` pins Lean 4.35.0 and the Mathlib dependency used by the Lean sources.
- `CITATION.cff` gives citation metadata. `MANIFEST.sha256` pins the released source tree.

With Python 3.13, install the pinned replay dependencies and run the driver from the repository root:

```sh
python -m pip install python-flint==0.9.0 numpy==2.5.5 scipy==1.18.1
python -B replay/run.py
sha256sum -c MANIFEST.sha256
```

The driver checks the equal-flavor and closed flavor-interval certificates, the finite-cover counts, and the expected mutation refusals. It takes several minutes on one CPU. Do not use Python's optimization mode.

For the finite Lean checks, from `lean/` with elan installed, run:

```sh
lake exe cache get
lake env lean ../replay/FiniteCore.lean
lake env lean ../replay/flavor/FiniteCore.lean
```

The manuscript can be built from `paper/` with Tectonic:

```sh
tectonic --keep-intermediates main.tex
```

The GitHub Actions workflow runs the replay and its mutation checks, verifies the source manifest, compiles both Lean files, and builds the manuscript. It uploads the resulting PDF as a workflow artifact.

## Licenses and provenance

The replay code and Lean sources are licensed under Apache-2.0; the preprint is distributed under CC-BY-4.0 through its Zenodo record. The related gamma-model certificate retains its own record and license.

AI agents contributed proofs, certificates, code, Lean sources, and manuscript drafts or edits under the author's direction and review.