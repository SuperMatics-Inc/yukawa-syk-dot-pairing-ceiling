"""Arb guards for all widths >=180 uniformly on the stated flavor interval."""
import argparse
from fractions import Fraction as F

from flint import arb, ctx

ctx.prec = 96
MUTANTS = ('crossing-half', 'mass-sign', 'weight-sign', 'compression-sign', 'skip-edge')
RLO, RHI = arb(99)/50, arb(101)/50


def Q(a, b=1):
    return arb(a) / b


def gamma_guards():
    n = 48
    w = [2*i+1 for i in range(2*n)]
    k = [F(0)]+[F(1, j*j) for j in range(1,3*n)]
    H = [sum(k[1:i+1]) for i in range(n)]
    mat = [[k[abs(i-j)]+k[i+j+1]-2*(H[i] if i==j else 0)
            for j in range(n)] for i in range(2*n)]
    shift = [[mat[i][j]/w[i]+2*(i==j) for j in range(n)] for i in range(n)]
    v = [F(1)]+[F(0)]*(n-1)
    d = 10**12
    for _ in range(70):
        update = [sum(a*b for a,b in zip(row,v)) for row in shift]
        v = [F(round(d*a/update[0]),d) for a in update]
    assert min(v) >= 0
    y = [Q(a.numerator, a.denominator) for a in v]
    def row_value(row):
        return sum(Q(a.numerator,a.denominator)*b for a,b in zip(row,y))
    mv = list(map(row_value,mat))
    s = sum(a*a*wi for a,wi in zip(y,w))
    rho = sum(a*b for a,b in zip(y,mv))/s
    tail = 25*sum(y)**2/(8*(2*n-1)**4)
    e2 = (sum((b-rho*wi*a)**2/wi for a,b,wi in zip(y,mv,w))
          +sum(b*b/wi for b,wi in zip(mv[n:],w[n:]))+tail)/s
    if not (rho-Q(1,10) > Q(1,5) and e2 < Q(1,1500)**2):
        raise ArithmeticError(f'failure isolating eigenvector: {rho=} {e2=}')
    assert rho > Q(25,19) and rho+Q(1,1500) < Q(500,379)
    assert Q(3,2)*Q(1,1500)/(rho-Q(1,10)) < Q(1,200)
    assert y[0]**2/s-Q(1,100) > Q(4,5)
    c1 = (y[0]**2+2*sum(a*b for a,b in zip(y,y[1:]))-2*sum(a*a for a in y[1:]))/s
    assert c1-6*Q(1,200) > Q(11,10)
    assert (y[0]*y[1]/s)-Q(1,100)>Q(9,50)
    projected_guard(y,s)
    return rho, e2.sqrt()


def projected_guard(y,s):
    # Uniform bounds for |Delta v|/epsilon, where v_i=sqrt(w_i)y_i/sqrt(s).
    e,A,z2 = Q(1,100),Q(4,5),Q(5,3)
    n=len(y)
    def D_at(j):
        return A*A*(1+2*RHI*j+Q(13,500)*RHI*j*j)/j**4
    D = [arb(0)]+[D_at(j) for j in range(1,3*n)]
    b=[Q(2*i+1,2)+A*sum((Q(1,j*j) for j in range(1,i+1)),arb(0)) for i in range(2*n)]
    q=[max(arb(0),1-e*x) for x in b]
    H=arb(0)
    norm2=arb(0)
    for i in range(2*n):
        w=2*i+1
        if i:
            H+=D[i]
        diaglo=-(D[w] if w<len(D) else D_at(w))-2*A*b[i]/w**2
        diaghi=2*H+A*(z2+Q(1,2*w*w))*w+(A*A*Q(11,10) if not i else 2*A*A*z2*z2)
        lo = hi = arb(0)
        for j in range(n):
            if j==i:
                lo+=diaglo*y[j]
                hi+=diaghi*y[j]
            else:
                l,h=abs(i-j),i+j+1
                lo-=y[j]*(D[l]+D[h]+A*(Q(1,l*l)+Q(1,h*h))*(b[i]+b[j])/(1+(q[i]*q[j]).sqrt()))
        norm2 += (lo.union(hi).abs_upper())**2/w
    # Beyond 2n, i-j>=i/2 and w_i<=3*i; no trial diagonal there.
    cap=5*(D[1]+2*A*Q(23,30))
    norm2 += 3*cap*cap*sum(y)**2/(2*(2*n-1)**2)
    assert (norm2/s).sqrt()+Q(19,200)<Q(7)


def wide_guards(mutant):
    e, eb = Q(1,180), Q(1,100)
    A,z2,z3,z4=Q(4,5),Q(5,3),Q(121,100),Q(11,10)
    m = Q(304,75)
    # Compression of the gamma-two limit to n>=1.
    subtraction = -1 if mutant=='compression-sign' else 1
    assert (2*z2-subtraction*2-1-Q(1,4))/3 < Q(1,10)
    assert (2*z2-subtraction*Q(5,2)-Q(1,4)-Q(1,9))/5 < Q(1,10)
    assert Q(2,3*7) < Q(1,10)
    assert sum(Q(1,j*j) for j in range(1,6))+Q(1,5) < z2
    assert sum(Q(1,j**4) for j in range(1,4))+Q(1,81) < z4
    assert sum(Q(1,j**3) for j in range(1,5))+Q(1,32)<z3
    # Global bootstrap and midpoint mass, uniformly in the flavor ratio.
    assert arb("4.046104720") < m
    assert 1+RHI*m/4 < 4
    assert m*(1+eb)/(1-Q(5,3)*m*eb)<Q(9,2)
    assert (1/(1-Q(9,2)*eb/4)).sqrt() < Q(101,100)
    assert Q(5,6)*Q(9,2)*Q(101,100)<4
    assert m/(Q(4,3)-4*eb-eb**2/2) < 4*A
    assert A*Q(101,100)<1
    assert A*(1+2/(1-eb))+Q(5,6)*4*A*Q(4,3) < 6
    assert Q(5,6)*4*A+6*(1+eb)+Q(1,2) < 10
    # Transfer, diagonal and off-diagonal contractions at the uniform endpoints.
    assert z2/7 < Q(1,3)
    transfer = 2*A*A*(z4+2*RHI*z3+Q(13,500)*RHI*z2)
    weight = 2*A*Q(23,30)*(Q(7,4)*3+z2)+A*Q(23,30)*(2*z2+1)
    assert transfer+weight<19
    crossing = RLO if mutant=='crossing-half' else 2*RLO
    assert (1+crossing*(1-e/2)**2)/(1+(1+2*RHI)*A*e)>Q(19,4)
    assert 1/(1+(1+4*RHI)*A*e)>Q(93,100)
    assert Q(1,2)-(Q(26,15)+Q(1,8))*eb>Q(12,25)
    assert 2*Q(93,100)*Q(12,25)>Q(89,100)
    decay = -Q(89,100) if mutant=='weight-sign' else Q(89,100)
    assert Q(4,5)*(z2+Q(1,2)-decay)+Q(1,5)*(z2+Q(1,18))<Q(137,100)
    assert Q(4,5)*z4+Q(1,5)*2*z2*z2/3<Q(63,50)
    assert (3-A/9)/2-9*e/8>Q(7,5)
    assert 1/(4+(1+8*RHI)*A*e)>Q(6,25)
    assert Q(9,50)*(Q(12,25)+Q(7,5))*(Q(93,100)+Q(6,25))>Q(39,100)
    loss = 0 if mutant == 'skip-edge' else Q(39,100)
    coef = Q(137,100)-loss
    f = lambda a: coef*a-Q(793,200)*a*a
    assert coef-Q(793,200)*Q(11,20)<0
    assert coef-2*Q(793,200)*Q(37,50)<0
    signed_residual = -1/RHI if mutant=='mass-sign' else 1/RHI
    static = Q(3,4)*(Q(5,9)*4*A-signed_residual)+Q(15,2)*e
    assert static/(1-e*static)<Q(101,100)
    assert Q(11,20)/Q(379,500)+19*e<1
    gap = lambda a:a*(Q(25,19)-Q(1,10))-26*e
    assert gap(Q(11,20))>0
    assert Q(37,50)/Q(379,500)+49*e*e/gap(Q(11,20))<1
    assert Q(101,100)+f(Q(37,50))+49*e/gap(Q(37,50))<0


def main():
    cli = argparse.ArgumentParser()
    cli.add_argument('--mutant',choices=MUTANTS)
    args = cli.parse_args()
    rho, res = gamma_guards()
    if not __debug__:
        raise ArithmeticError("replay requires Python assertions")
    wide_guards(args.mutant)
    print(f'PASS: sigma >= 180, r=[99/50,101/50], infinite-tail Rayleigh={rho}, residual bound={res}')


if __name__=='__main__':
    main()
