"""Branch-free Arb replay for a closed rectangle in (static width, B/width).

A failed cell is uncovered, never accepted. Coordinates are dyadic and all
claims depend on interval predicates, not the subdivision heuristic.
"""
import argparse
from fractions import Fraction as F

from flint import arb, ctx

ctx.prec = 80
MASS = arb('4.046104720')


def Q(x):
    return arb(x.numerator)/x.denominator


def span(lo, hi):
    return lo.union(hi)


def bounds(b, lo, hi):
    a, z = max(b.lower(), lo.lower()), min(b.upper(), hi.upper())
    return (a, z) if not a > z else None


def root(u, sigma):
    return 2/((u*u+4*sigma*sigma).sqrt()+u)


def transfer(j, a, B, norm2, factor):
    return 1/(4*j*j+a+factor*B*norm2)


def guard_polarization(factor):
    # h_0=1, h_{i>0}=0, reflected with a minus sign. For j=1 the
    # independent full-line bubble difference is 2-(-1)=3.
    # At r=2 the shift term is 2B*3; the folded squared norm is 4+2=6.
    if not abs(transfer(1, arb(2), arb(3), arb(6), factor)-1/arb(24)) < 1/arb(10)**20:
        raise ArithmeticError('folded polarization normalization failed')


def cell(sl, su, pl, pu, n, k, passes, mutant):
    factor = 2 if mutant == 'double-polarization' else 1
    if pu == 0:
        return 'weak'
    si0, si1, p0, p1 = map(Q, (sl,su,pl,pu))
    B0, B1 = si0*p0, si1*p1
    if B1 <= Q(F(25,13)):
        return 'weak'
    w = [2*i+1 for i in range(n+1)]
    a0, a1 = p0/si1, p1/si0
    # Upper mass at the critical plateau, rounded outward via GAMMA/RESULT.md.
    mstar = MASS
    if a0/2 > mstar:
        return 'mass'
    maxq = (1+a1/4).sqrt()
    monotone = p1*maxq/4 < 1
    if monotone:
        maxq = arb(1)
    lip = (1+p1/2)/si0**2
    R = 5*p1*maxq/6
    top = [min(1/arb(wi),maxq/si0).upper() for wi in w]
    tail_top = top[k]  # unrefined monotone envelope, valid at every index >= k
    low = [root(arb(wi)+R,si1).lower() for wi in w]
    # Midpoint lower bound for the unrepresented positive mass tail.
    # The primitive of t(x)^2 from x to infinity is t(x)+t(x)^3/3.
    qtail = si0*root(arb(2*n)+R,si0)
    mass_tail = max(arb(0),(qtail+qtail**3/3-1/(2*si0**2))/(2*si1))
    if a0/2+2*B0*(sum(x*x for x in low[:n])+mass_tail) > mstar:
        return 'mass'
    d0 = [arb(0)]+[1/(4*arb(j*j)+a1+8*mstar) for j in range(1,2*n)]
    d1 = [arb(0)]+[1/(4*arb(j*j)+a0) for j in range(1,2*n)]
    pU = [[(arb(1)/abs(i-j)**2-arb(1)/(i+j+1)**2)/4 if i!=j else arb(0)
           for j in range(n)] for i in range(n)]
    for _ in range(passes):
        # Folded polarization: j crossing terms and twice the same-sign sum.
        for j in range(1,n+1):
            L = U = arb(0)
            for i in range(j):
                l, u = low[i]+low[j-1-i], top[i]+top[j-1-i]
                L += l*l
                U += u*u
            for i in range(n-j):
                l = max(arb(0),low[i]-top[i+j],low[i+j]-top[i])
                if not l > 0:
                    l = arb(0)
                u = (top[i]-low[i+j]).union(top[i+j]-low[i]).upper()
                if monotone:
                    u = min(u,j*lip)
                L += 2*l*l
                U += 2*u*u
            # h_i, h_{i+j} both lie in [0,1/w_i].
            N = w[n-j]
            tail_shift = arb(1)/N**2+arb(1)/(2*N)
            if monotone:
                tail_shift = min(tail_shift,j*top[n-j]*min(top[n-j],j*lip))
            U += 2*tail_shift
            if mutant == 'crossing-half':
                L /= 2  # weaker, can only cause unresolved cells
            d0[j] = max(d0[j],transfer(j,a1,B1,U,factor)).lower()
            d1[j] = min(d1[j],transfer(j,a0,B0,L,factor)).upper()
        plow, pup = [], []
        for i in range(n):
            lm = hm = arb(0)
            for j in range(n):
                if i == j:
                    continue
                l, h = abs(i-j), i+j+1
                lm += max(arb(0), (d0[l]-d1[h]).lower())*low[j]
                hm += min(pU[i][j],(d1[l]-d0[h]).upper())*top[j]
            v = n-i
            J0 = (B0*lm-B1*d1[w[i]]*top[i]).lower()
            tail_normal = top[n]*min(arb(1)/v**2+arb(1)/v, arb(w[i])/v**2)/4
            J1 = (B1*(hm+tail_normal)-B0*d0[w[i]]*low[i]).upper()
            found = bounds(span(root(arb(w[i])+J1,si1).lower(),root(arb(w[i])+J0,si0).upper()),low[i],top[i])
            if found is None:
                return 'empty'
            l, u = found
            plow.append(l)
            pup.append(u)
        low[:n], top[:n] = plow, pup
        if a0/2+2*B0*(sum(x*x for x in low[:n])+mass_tail) > mstar:
            return 'mass'
    # Head majorant; no omitted zero-transfer term. Drop only negative normal
    # contributions beyond n, and use Schur with a positive head test vector.
    head = []
    for i in range(k):
        row = [B1*(d1[abs(i-j)]+d1[i+j+1])*(top[i]*top[j]/(w[i]*w[j])).sqrt()
               for j in range(k)]
        normal = sum(max(arb(0),(d0[abs(i-j)]-d1[i+j+1]).lower())*low[j]
                     for j in range(n) if i!=j)
        row[i] = (2*B1*d1[w[i]]*top[i]-B0*normal)/w[i]
        head.append(row)
    # Tail Schur with sqrt(h_i/w_i). Head-tail coupling is bounded in HS norm.
    tail = (5*B1/6+B1/(4*arb(w[k])**2))*tail_top/w[k]
    if not tail < 1:
        return None
    cross = [(B1**2*tail_top*top[i]/(4*w[k]*w[i])*
              (arb(1)/(k-i)**4+arb(1)/(3*(k-i)**3))).sqrt() for i in range(k)]
    if mutant == 'inflate-tail':
        cross = [32*c for c in cross]
    for i in range(k):
        for j in range(k):
            head[i][j] += cross[i]*cross[j]/(1-tail)
    # The shift makes the head iteration positive. The vector is merely a
    # witness; only the subsequent verified weighted inequalities matter.
    vec = [arb(1)/(i+1)**3 for i in range(k)]
    for _ in range(40):
        hv = [sum(a*b for a,b in zip(row,vec)) for row in head]
        if all(v < z for v,z in zip(hv,vec)):
            return 'pair'
        if not hv[0]+vec[0] > 0:
            return None
        ref = hv[0]+vec[0]
        vec = [max(arb('1e-12'),((v+z)/ref).mid()) for v,z in zip(hv,vec)]
    return None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--sigma',nargs=2,type=F,required=True)
    parser.add_argument('--p',nargs=2,type=F,required=True)
    parser.add_argument('--n',type=int,default=12)
    parser.add_argument('--k',type=int,default=6)
    parser.add_argument('--passes',type=int,default=3)
    parser.add_argument('--depth',type=int,default=18)
    parser.add_argument('--mutant',choices=('crossing-half','inflate-tail','double-polarization'))
    args = parser.parse_args()
    sl,su = args.sigma
    pl,pu = args.p
    if not 0 < sl < su or not 0 <= pl < pu or not 2 <= args.k < args.n or args.passes < 1:
        parser.error('invalid rectangle or truncation')
    guard_polarization(2 if args.mutant == 'double-polarization' else 1)
    assert 4/(3*arb.pi()**2*arb('0.1827262477723295078215179')**2) < MASS
    if pl == 0 and not Q(pu)/(2*Q(su)) > MASS:
        parser.error('upper p must include the positive-mass bound')
    todo = [(sl,su,pl,pu,0)]
    counts = {key:0 for key in ('mass','weak','empty','pair')}
    while todo:
        a,b,c,d,dep = todo.pop()
        if result := cell(a,b,c,d,args.n,args.k,args.passes,args.mutant):
            counts[result] += 1
        elif dep < args.depth:
            if (b-a)*max(d,F(1)) > (d-c)*b:
                mid = (a+b)/2
                todo.extend([(a,mid,c,d,dep+1),(mid,b,c,d,dep+1)])
            else:
                mid = (c+d)/2
                todo.extend([(a,b,c,mid,dep+1),(a,b,mid,d,dep+1)])
        else:
            raise ArithmeticError(f'uncovered cell sigma=[{a},{b}] p=[{c},{d}], {counts=}')
    print(f'PASS: sigma=[{sl},{su}] p=[{pl},{pu}], {counts=}, n={args.n}, k={args.k}, normal passes={args.passes}')


if __name__ == '__main__':
    main()
