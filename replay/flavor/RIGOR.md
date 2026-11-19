# Uniform extension around equal flavor

**PROVED, with the finite strip dependent on the Arb replay:** Let the exact centered causal normal and real-channel singlet equations be those in the dot manuscript (`../../paper/main.tex`), with the large-N limit taken first, bare frequency \(\Omega_0>0\), and \(r=b/a\) in the closed rational interval \([99/50,101/50]\). For all finite positive \(a\),

\[
 k_BT_{\rm lin}<\tau_*(r):=C_\infty\sqrt{3/(4r)}\,\Omega_0,
 \qquad \lim_{a\to\infty}k_BT_{\rm lin}=\tau_*(r).
\]

Here \(T_{\rm lin}\) is the supremum of linear-instability temperatures over eligible centered causal normal saddles. No uniqueness, free-energy ordering, finite-N claim or exclusion of first-order transitions is assumed. The interval is for the ratio in the saddle, including the spin trace; the equal-flavor point is \(r=2\). This is a uniform *interval* proof, not a scan at finitely many ratios. It does not decide the intermediate-coupling question at arbitrary positive \(r\).

## Normalization and weak width

Use the dimensionless variables of the manuscript: \(B=a/(\pi^4\tau^3)\), \(\sigma^2=B d_0\), \(p=B/\sigma\), \(\epsilon=1/\sigma\), \(q_n=\sigma h_n\), \(K_j=p d_j\), \(\alpha=p/4\), \(m=\Omega_0^2/[r(\pi\tau)^2]\). At \(\tau\ge\tau_*(r)\), \(m\le16\kappa/3\) with \(\kappa=1/\rho\) the inverse gamma-two eigenvalue. In particular the bound on \(m\) is independent of \(r\). The exact equations, with \(Q\) the signed full-line reflection of \(q\), give

\[
 d_j^{-1}=4j^2+p\epsilon+(rp\epsilon/2)\|Q-S_jQ\|^2,
 \qquad m=p\epsilon/r+2p\epsilon\sum_{n\ge0}q_n^2.
\]

Use \(r_- = 99/50\), \(r_+ = 101/50\). The weak Hankel bound in the manuscript does not depend on \(r\), because its only transfer input is \(0\le d_j\le (4j^2)^{-1}\). It excludes pairing if \(B\le25/13\). For \(\sigma\le15/32\) we have \(B\le r_+ m\sigma^2 < 25/13\); the last rational inequality is checked in `cover.py` and `FiniteCore.lean`.

## Finite width

Use the completely folded norm of Appendix B in `../../paper/main.tex` on the \(H=Q/\sigma\) scale. Given \(L,U\) as outward lower/upper bounds of that norm square and \(x_{\rm lo},x_{\rm hi}\) bounds for \(p/\sigma\), replace the refined transfer interval by

\[
 [4j^2+x_{\rm hi}+(r_+/2) B_{\rm hi}U]^{-1}\le d_j\le
 [4j^2+x_{\rm lo}+(r_-/2) B_{\rm lo} L]^{-1}.
\]

This follows by positivity of all factors. For the initial all-transfer lower bound, \(\|H-S_j H\|^2\le4\|H\|^2\le4m/B\); replacing the resulting \(2r_+m\) by \(4r_+ m_{\rm cap}\) is conservative. The static-mass lower bound uses \(x_{\rm lo}/r_+\) rather than \(x_{\rm lo}/2\). The normal root upper/lower bounds, Stieltjes difference and pairing block bounds do not otherwise change. In particular the bound \(q_n^2\le1+(p/\sigma)/4\) uses the static *residual* \(p/\sigma\) in \(d_1^{-1}\), which does not acquire an extra \(r\) factor. Caps on \(p\) cover to \(r_+m_{\rm cap}\sigma_{\rm hi}\) in the rectangles. The cap \(m_{\rm cap}=4.046104720\) is checked from the certified lower endpoint for \(C_\infty\); `../gamma-two.md` states the tail assumptions and normalization.

## Wide widths

Apply the large-width estimate in Appendix A of `../../paper/main.tex` with \(\epsilon\le1/180\). Here are all flavor-dependent changes. For the initial bootstrap use \(q_n^2\le1+r_+m/4<4\). The midpoint and normal-sum arguments after this bound are independent of \(r\): \(p<16/5\), the \(q_n\) decrease, \(0\le1-q_n\le\epsilon b_n\), \(b_n=w_n/2+(4/5)H_n\), and \(2\epsilon\sum q_n^2\ge4/3-(5p/9)\epsilon-10\epsilon^2\). In the improved static estimate put

\[
 c=\frac34\left(\frac59\frac{16}5-\frac1{r_+}\right),\quad
 m\ge\frac43p[1-c\epsilon-(15/2)\epsilon^2].
\]

The first term in the exact mass equation is positive; the lower bound uses \(1/r_+\). The guard in `replay.py` proves \(\alpha/\kappa\le[1-c\epsilon-(15/2)\epsilon^2]^{-1}\le1+(101/100)\epsilon\) at the worst endpoint \(\epsilon=1/180\), by increasing the correction in both numerator and denominator.

In the transfer bound replace \(D_j\) by

\[
 D_j=(4/5)^2[j^{-4}+2r_+j^{-3}+(13r_+/500)j^{-2}].
\]

The crossing part of the norm is at most \(4j\); the same-sign part, including both halves, is at most \(2j\cdot j(13/5)\epsilon\), with \(\epsilon\le1/100\) in this inequality. After multiplication by \(r/2\) these give the displayed coefficients. The rational trial in `replay.py` uses these constants in the finite projection and all-infinite residual tail, giving \(\|\Delta z\|<7\epsilon\) and \(\|\Delta\|<19\epsilon\) for the exact perturbation \(\Delta=S-\alpha A\) and the infinite gamma-two unit top vector \(z\).

For the first transfer, the crossed term is at least \(4q_0^2\) and \(q_0\ge1-\epsilon/2\). Let \(e=1/180\) and \(A_*=4/5\). Monotonicity of \(x/(1+\alpha\epsilon x)\) gives the correction at least \(\epsilon\alpha^2\) times

\[
 \frac{1+2r_-(1-e/2)^2}{1+(1+2r_+)A_*e}>19/4.
\]

The denominator is a valid upper bound **after** replacing \(x\) by its smaller crossing-only value in the increasing function. The simple bound \(\|Q-S_jQ\|^2\le8j\) gives \(K_1\ge\alpha/[1+(1+4r_+)A_*e]\) and \(K_2\ge\alpha/[4+(1+8r_+)A_*e]\). These retain the strict \((93/100)\alpha,(6/25)\alpha\) thresholds. The signed normal weight contractions and the critical-vector bounds are unchanged. Thus the same Rayleigh correction \(F(\alpha)=(98/100)\alpha-(793/200)\alpha^2\) and the same three spectral-gap regimes in Appendix A of `../../paper/main.tex` give \(\sup\operatorname{spec}S<1\) on all \(\sigma\ge180\).

The argument in the manuscript excluding approach to the ceiling at fixed finite \(a\) is unchanged: use \(m=B/(r\sigma^2)+2B\sum h_n^2\), \(r\) fixed positive, to rule out both zero and infinite limiting width. Dominated convergence then carries all necessary inequalities to the limit in a common scaled space. Existence at positive \(\tau\) and convergence to the gamma-two top as \(a\to\infty\) hold at every positive finite \(r\) by the existence appendix and strong-coupling section of the manuscript, both written with \(b=ra\) arbitrary. At \(\tau<\tau_*(r)\) the limiting eliminated pairing top exceeds one, completing sharpness without extrapolating in \(r\) or in the Matsubara index.
