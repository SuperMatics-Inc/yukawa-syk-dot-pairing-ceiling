# The gamma-two eigenvalue enclosure

**PROVED, subject to the indicated operator and interval analysis.** The gamma-two matrix is Eq. (2) of the manuscript, with indices `n,m >= 0`, weights `w_n=2n+1`, and zero on the *Toeplitz* diagonal only. The diagonal Hankel term and the finite negative normal term remain. This appendix gives the tail and eigenvalue ingredients behind `two_tight.py`, which replays the gamma-two part of [the separate certificate](https://doi.org/10.5281/zenodo.24366341). The computation depends on the infinite estimates, standard self-adjoint spectral theory and Arb outward arithmetic. A finite floating-point solve alone is not an enclosure.

For the block with indices `n,m>=N` drop the negative diagonal and use positive Schur weights `w_n^(-1/2)`. The weighted row sum splits at `m <= n/2`, `n/2 < m <= 2n`, and `m > 2n`. In the first range both transfers are at most `4/n²` and `2 sum_(N<=m<=n/2) 1/w_m <= 1+log(n/N)`. Taking the supremum and bounding `u^-2 log u <= 1/(2e)` gives `4/N² (1+1/(2e))`. In the middle range `1/w_m<=1/n`, the Toeplitz part sums to at most `2 zeta(2)`, and the Hankel part is at most `2/n` before the weighted denominator, since it has at most `2n` terms bounded by `1/n²`. In the last range, both transfers are at most `4/m²`, so the combined contribution is at most `sum_(m>2n) 4/m³ <= 1/(2n²)`. The resulting upper bound for the spectral top of the tail is

```
c_N = 4/N² (1+1/(2e)) + 2 zeta(2)/N + 2/N² + 1/(2N²).
```

These estimates also bound the norm of the nonnegative tail majorant. The discarded negative diagonal tends uniformly to zero on the tail. The head-to-tail operator is finite rank at every fixed split, which proves compactness of the self-adjoint limit. For each `i<N`, `k=m-i>=N-i` gives two transfers at most `1/k²` each, and `2m+1>=2k`. Thus the full coupling row norm squared is at most

```
R_i² = 2/(2i+1) * [ (N-i)^(-5) + (N-i)^(-4)/4 ].
```

Its derivation uses the first summand plus the `k^-5` integral. In the gap bound use `N=192` and `beta²=sum_(i<N) R_i²` as an upper bound on the squared cross norm. In the trial residual use `N=8192` and `||B* v|| <= sum_i |v_i| R_i`, not a finite tail sample. This bounds the part of the residual outside the trial support.

To bound the second eigenvalue, deflate the head by a positive rank-one term with head unit vector `y`. The original operator agrees with the deflated operator on the quadratic forms of vectors orthogonal to `y` extended by zero, so min--max bounds the original second eigenvalue by the deflated spectral top. For cutoff `s` and `delta=5e-7`, the code proposes a floating Cholesky factor `L`, then encloses `E=(s-delta)I-A_deflated-LL^T` with Arb and requires each absolute row sum strictly below `delta`. Hence `sI-A_deflated` is positive definite. Completing the two block bounds gives

```
u_2 = (s+c_N)/2 + sqrt(((s-c_N)/2)^2+beta²).
```

The value passed to Temple is checked strictly above a freshly evaluated all-tail bound. If `q` is the exact trial Rayleigh quotient and `v` the exact finite-supported trial of unit norm, then `q <= rho <= q+||Av-qv||²/(q-u_2)` provided `q>u_2`. To see the upper direction, the spectral measure has `(lambda-u_2)(lambda-rho)>=0` everywhere on the spectrum; integration gives `||Av-qv||² >= (rho-q)(q-u_2)`. The finite action in `q` and the residual is a pair of Toeplitz convolutions and one Hankel convolution, checked with Arb. The only role of SciPy is to propose a finite vector; its binary64 entries are then treated as exact dyadic data. `run.py` verifies the gamma enclosure before using it in the finite-strip mass cap and the interval for `C_infty sqrt(3/8)`. The arithmetic mutant doubles the latter normalization and must be rejected.
