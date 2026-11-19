"""Run the dot ceiling certificate and its rejection checks, without private paths."""
import hashlib
import os
from pathlib import Path
import subprocess
import sys

from flint import arb, ctx
import two_tight

ROOT = Path(__file__).resolve().parent
PINS = {
    "cover.py": "7deea3c26b7eff520888693176c0ae57b8b9b3bd47a0024fab3a8a2989d57151",
    "replay.py": "717230167af33ef942704babdcd6320ae00aed469580b346ee5b4876510534c4",
    "FiniteCore.lean": "7dcadd7a55d10cacba01c6e8e66feb0905da5c3f22a528554dfd2bfa291c2e1e",
    "two_tight.py": "78634b231b0519bceb73827eb157e57933d55cafe81b6f93691fc0cd646d8411",
    "flavor/cover.py": "c0b7c15099a8edb14a07ae7679e589c47115f105318aacf1d47d62212d825fcd",
    "flavor/replay.py": "c5261548f6cd43e8930c2d2e2bc3c170b1e1793bca00c576fc7907da6a925531",
    "flavor/FiniteCore.lean": "980587b3ab4f0b060bf2eab18b478c9c2f98a95994a22e8df354f86a0d30de1d",
}
STRIPS = (
    ("15/32 4 0 33 12 6 3 28", (35, 3, 0, 35)),
    ("4 16 0 130 24 10 6 28", (98, 0, 0, 163)),
    ("16 40 0 324 32 12 6 28", (235, 0, 0, 381)),
    ("40 100 0 810 40 16 6 28", (647, 0, 0, 1155)),
    ("100 180 0 1457 48 20 6 31", (1050, 0, 0, 1592)),
)
FOCUS = "99 100 3 16/5 40 16 6 16"
FLAVOR_STRIPS = (
    ("15/32 4 0 33 16 8 4 28", (32, 3, 0, 36)),
    ("4 16 0 131 28 12 6 29", (99, 0, 0, 160)),
    ("16 40 0 327 36 14 7 30", (230, 0, 0, 394)),
    ("40 100 0 818 56 20 8 32", (608, 0, 0, 1096)),
    ("100 180 0 1472 64 24 8 34", (1078, 0, 0, 1619)),
)
FLAVOR_FOCUS = "99 100 3 16/5 56 20 8 16"
WIDE_MUTANTS = ("crossing-half", "mass-sign", "weight-sign", "compression-sign", "skip-edge")
FINITE_MUTANTS = {
    "crossing-half": "uncovered cell",
    "inflate-tail": "uncovered cell",
    "double-polarization": "folded polarization normalization failed",
}


def require(condition, message):
    if not condition:
        raise ArithmeticError(message)


def invoke(*args, refusal=None):
    env = os.environ.copy()
    env.pop("PYTHONOPTIMIZE", None)
    print("+ python -B", *args, flush=True)
    out = subprocess.run([sys.executable, "-B", *args], cwd=ROOT, env=env,
                         text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    print(out.stdout, end="", flush=True)
    if refusal is None:
        require(out.returncode == 0, f"certificate failed: exit {out.returncode}")
    else:
        require(out.returncode == 1 and refusal in out.stdout,
                f"mutant was not refused for the expected reason: exit {out.returncode}")
        print("expected refusal: exit 1", flush=True)
    return out.stdout


def cover_args(row, script):
    sl, su, pl, pu, n, k, passes, depth = row.split()
    return (script, "--sigma", sl, su, "--p", pl, pu, "--n", n,
            "--k", k, "--passes", passes, "--depth", depth)


def cover(row, counts, script):
    result = invoke(*cover_args(row, script))
    expected = dict(zip(("mass", "weak", "empty", "pair"), counts))
    require(f"counts={expected}" in result, "cover passed with a different receipt")
    if script == "flavor/cover.py":
        require("r=[99/50,101/50]" in result, "unexpected flavor interval")


def plateau(lo, hi, normalization=1):
    s = (arb(3) / 8).sqrt() * normalization
    require(arb("0.1118965174138946") < s * lo < s * hi < arb("0.1118965174138948"),
            "plateau endpoint or normalization failed")


def main():
    require(__debug__, "replay requires Python assertions")
    for name, expected in PINS.items():
        require(hashlib.sha256((ROOT / name).read_bytes()).hexdigest() == expected,
                f"SHA-256 mismatch: {name}")
    ctx.prec = 160
    actual_lo, actual_hi, _ = two_tight.certify()
    lo = arb("0.1827262477723295078215179")
    hi = arb("0.1827262477723295176168720")
    require(lo < actual_lo < actual_hi < hi, "gamma-two interval was not reproduced")
    require(4 / (3 * arb.pi() ** 2 * lo ** 2) < arb("4.046104720"),
            "finite-strip mass endpoint failed")
    plateau(lo, hi)
    try:
        plateau(lo, hi, normalization=2)
    except ArithmeticError as exc:
        print(f"expected doubled-temperature refusal: {exc}", flush=True)
    else:
        raise ArithmeticError("doubled-temperature mutant survived")
    print("PASS: pinned gamma-two certificate and one expected refusal")
    cover_lane("replay.py", "cover.py", STRIPS, FOCUS, (8, 0, 0, 14))
    print("PASS: pinned equal-flavor certificate and eight expected refusals")
    cover_lane("flavor/replay.py", "flavor/cover.py",
               FLAVOR_STRIPS, FLAVOR_FOCUS, (12, 0, 0, 14))
    print("PASS: pinned flavor-interval certificate and eight expected refusals")


def cover_lane(wide_script, cover_script, strips, focus, focus_counts):
    invoke(wide_script)
    for mutant in WIDE_MUTANTS:
        invoke(wide_script, "--mutant", mutant, refusal="AssertionError")
    for row, counts in strips:
        cover(row, counts, cover_script)
    cover(focus, focus_counts, cover_script)
    for mutant, reason in FINITE_MUTANTS.items():
        invoke(*cover_args(focus, cover_script), "--mutant", mutant, refusal=reason)


if __name__ == "__main__":
    main()
