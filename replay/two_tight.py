"""Gamma=2 KAY-operator certificate, using a finite trial and an infinite tail.

Run from this directory:
    uv run --no-project --with numpy --with scipy --with python-flint python -u two_tight.py

The floating-point eigensolver and Cholesky factor only propose data.  All
certificate arithmetic, including the full finite-trial residual and the
factorization check, uses 160-bit Arb balls.  No input numerical eigenvalue is
assumed to be an upper bound.  The gamma=2 operator and the tail and Temple
lemmas are the ones stated in RIGOR.md; they are not formalized here.
"""

import os

# Set before NumPy/SciPy import, even when the inherited environment differs.
for _name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
              "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[_name] = "1"

import numpy as np
from scipy.linalg import matmul_toeplitz
from scipy.signal import fftconvolve
from scipy.sparse.linalg import LinearOperator, eigsh
from flint import arb, arb_mat, arb_poly, ctx

ctx.prec = 160
TRIAL_SIZE = 8192
GAP_SIZE = 192


def floating_trial(n):
    """Propose a positive trial, without using floats as certificate bounds."""
    indices = np.arange(n, dtype=float)
    roots = 1 / np.sqrt(2 * indices + 1)
    powers = 1 / np.arange(1, 2 * n + 1, dtype=float) ** 2
    diagonal = -2 * np.r_[0.0, np.cumsum(powers[:n - 1])] / (2 * indices + 1)
    column = np.r_[0.0, powers[:n - 1]]

    def apply(v):
        x = roots * v
        toeplitz = matmul_toeplitz((column, column), x, check_finite=False)
        hankel = fftconvolve(powers, x[::-1])[n - 1:2 * n - 1]
        return roots * (toeplitz + hankel) + diagonal * v

    op = LinearOperator((n, n), matvec=apply, dtype=np.float64)
    _, vectors = eigsh(op, k=1, which="LA", tol=2e-15,
                       v0=np.ones(n), maxiter=1000)
    v = vectors[:, 0]
    if v[0] < 0:
        v = -v
    if not np.all(np.isfinite(v)) or not np.all(v > 0):
        raise ArithmeticError("floating trial is not finite and positive")
    return v


def _gap_bound(v, rho):
    """Upper bound lambda_2 of the infinite operator by checked deflation."""
    n = GAP_SIZE
    roots = [1 / arb(2 * i + 1).sqrt() for i in range(n)]
    pp = [arb(0)] + [arb(1) / (k * k) for k in range(1, 2 * n + 1)]
    matrix = arb_mat(n, n)
    harmonic = arb(0)
    for i in range(n):
        if i:
            harmonic += pp[i]
        matrix[i, i] = (pp[2 * i + 1] - 2 * harmonic) * roots[i] ** 2
        for j in range(i + 1, n):
            x = (pp[j - i] + pp[i + j + 1]) * roots[i] * roots[j]
            matrix[i, j] = matrix[j, i] = x

    y = arb_mat([[arb(float(x))] for x in v[:n]])
    y = y / (y.transpose() * y)[0, 0].sqrt()
    shift = 4 * float(rho.mid())
    if not arb(shift) > 0:
        raise ArithmeticError("deflation shift must be positive")
    deflated = matrix - arb(shift) * (y * y.transpose())

    # Build only a floating-point proposal to a factor of cutoff*I-P.
    ii = np.arange(n)
    distance = np.abs(ii[:, None] - ii[None, :])
    hankel = ii[:, None] + ii[None, :] + 1
    weights = np.sqrt(2 * ii + 1)
    proposal = (np.where(distance == 0, 0.0, 1.0 / np.maximum(distance, 1) ** 2)
                + 1.0 / hankel ** 2) / np.outer(weights, weights)
    proposal[np.diag_indices(n)] -= (
        2 * np.r_[0.0, np.cumsum(1.0 / np.arange(1, n, dtype=float) ** 2)]
        / (2 * ii + 1))
    yf = v[:n] / np.linalg.norm(v[:n])
    proposal -= shift * np.outer(yf, yf)
    cutoff = float(np.linalg.eigvalsh(proposal)[-1] + 1e-6)
    margin = arb("0.0000005")
    chol = arb_mat(np.linalg.cholesky(
        (cutoff - 0.0000005) * np.eye(n) - proposal).tolist())
    error = -deflated - chol * chol.transpose()
    for i in range(n):
        error[i, i] += arb(cutoff) - margin
    for i in range(n):
        if not sum((abs(error[i, j]) for j in range(n)), arb(0)) < margin:
            raise ArithmeticError("Arb deflation factorization was not certified")

    # The nonnegative tail-majorant Schur bound c_n and the cross-block
    # Hilbert-Schmidt bound beta^2 are valid for the entire infinite tail.
    t = arb(n)
    c = (4 / t ** 2 * (1 + 1 / (2 * arb(1).exp()))
         + 2 * arb(2).zeta() / t + 2 / t ** 2 + 1 / (2 * t ** 2))
    beta2 = sum((2 * (1 / arb(n - i) ** 5 + 1 / (4 * arb(n - i) ** 4))
                 / (2 * i + 1) for i in range(n)), arb(0))
    s = arb(cutoff)
    return (s + c) / 2 + (((s - c) / 2) ** 2 + beta2).sqrt()


def gap_bound(v, rho):
    # Reserve a strict, outward cushion for the separate Temple recheck.
    return _gap_bound(v, rho).upper() + arb("1e-20")


class _SecondBoundError(ArithmeticError):
    def __init__(self, actual, validated):
        super().__init__("Temple second-eigenvalue bound was not independently revalidated")
        self.actual = actual
        self.validated = validated


def _require_strict_bound(actual, validated):
    # Overlap and equality are not sufficient: Arb must prove the ordering.
    if not actual > validated:
        raise _SecondBoundError(actual, validated)


def _temple_upper(rho, residual2, second, proposed):
    # Rebuild and check the all-infinite form bound from the operator, not
    # from the value supplied to Temple or from a published enclosure.
    validated = _gap_bound(proposed, rho).upper()
    _require_strict_bound(second, validated)
    return rho + residual2 / (rho - second)


def certify(n=TRIAL_SIZE):
    if n < GAP_SIZE:
        raise ValueError("trial size must be at least GAP_SIZE")
    proposed = floating_trial(n)
    v = [arb(float(x)) for x in proposed]
    if not all(x > 0 for x in v):
        raise ArithmeticError("trial positivity could not be certified")
    norm2 = sum((x * x for x in v), arb(0))
    roots = [1 / arb(2 * i + 1).sqrt() for i in range(n)]
    x = [v[i] * roots[i] for i in range(n)]

    # Finite exact-kernel action by Arb polynomial convolution.  With
    # d_0=0, d_k=k^-2, h_k=(k+1)^-2 and x^R reversed,
    # (Gv)_i = roots_i * [(d*x)_i + (d*x^R)_{n-1-i}
    #                       + (h*x^R)_{n-1+i} - 2 H_i*x_i].
    d = arb_poly([arb(0)] + [arb(1) / (k * k) for k in range(1, n)])
    h = arb_poly([arb(1) / (k * k) for k in range(1, 2 * n)])
    xp = arb_poly(x)
    xr = arb_poly(x[::-1])
    forward = d * xp
    backward = d * xr
    hankel = h * xr
    harmonic = arb(0)
    action = []
    for i in range(n):
        if i:
            harmonic += arb(1) / (i * i)
        action.append(roots[i] * (forward[i] + backward[n - 1 - i]
                                  + hankel[n - 1 + i] - 2 * harmonic * x[i]))
    rho = sum((v[i] * action[i] for i in range(n)), arb(0)) / norm2
    second = gap_bound(proposed, rho)
    if not (rho > 0 and second < rho):
        raise ArithmeticError("positive Rayleigh quotient and spectral gap not certified")

    finite_r2 = sum(((action[i] - rho * v[i]) ** 2 for i in range(n)), arb(0)) / norm2
    # For each i<n, sum_{m>=n} G_mi^2 <= R_i^2.  The triangle
    # inequality bounds the full infinite tail of Gv, not a finite truncation.
    tail = sum((v[i] * (2 * (1 / arb(n - i) ** 5
                              + 1 / (4 * arb(n - i) ** 4))
                        / (2 * i + 1)).sqrt() for i in range(n)), arb(0))
    tail_r2 = tail ** 2 / norm2
    upper_rho = _temple_upper(rho, finite_r2 + tail_r2, second, proposed)
    lower_t = rho.sqrt() / (2 * arb.pi())
    upper_t = upper_rho.sqrt() / (2 * arb.pi())
    width = upper_t - lower_t
    if not (lower_t > 0 and upper_t > lower_t and width < arb("1e-15")):
        raise ArithmeticError("temperature interval width was not certified")
    print("gamma=2; trial_size=", n, "; gap_size=", GAP_SIZE, sep="")
    print("rho Rayleigh ball:", rho.str(25))
    print("second-eigenvalue upper bound ball:", second.str(25))
    print("finite residual squared ball:", finite_r2.str(12))
    print("infinite tail residual squared upper-bound ball:", tail_r2.str(12))
    print("Tc/g lower endpoint ball:", lower_t.str(25))
    print("Tc/g upper endpoint ball:", upper_t.str(25))
    print("Tc/g enclosure uses lower endpoint's lower bound and upper endpoint's upper bound")
    print("Tc/g width bound ball:", width.str(12))
    return lower_t, upper_t, second


if __name__ == "__main__":
    certify()
