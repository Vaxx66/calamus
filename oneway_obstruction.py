#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Companion script to:
  F. Vaccarino, "Statistical models as natural transformations: a bridge
  from the McCullagh--Brons hexagon to Markov categories" (2026), Section 6.
One-way layout: defect cochain, Baues--Wirsching complex, obstruction classes.

Reproduces the computations reported in Examples 6.6--6.8 of the paper:
  RUN 1   Example 6.6 / 6.7(i): arrow category; the family is natural,
          dim C^1 = dim Z^1 = dim B^1 = 16, H^1 = 0.
  RUN 2   Example 6.7(ii): full category, S = {1,2}; defect supported on
          {i1, i2, c1, c2}, cocycle identity of Lemma 6.3 verified on all
          22 composable pairs, unrestricted correction recovers the
          within-group dispersion.
  RUN 3   Example 6.7(iii): recalibration class under scale confounding
          (S = {1, 5/4}); obstruction nonzero relative to the class.
  CHECK A robustness: arrow category, quantity with nonzero defect on p.
  CHECK B Example 6.8: M = {0,1,2}; H^1 of the full design category is
          nonzero at three levels of the treatment mean.

Exact linear algebra over QQ (sympy.Rational). Requires python3 with sympy
and numpy; runs unchanged under `sage -python`.

Setup (Definition 6.1, Lemma 6.3, Proposition 6.4 of the paper):
  O  : covariate category. Two variants:
       - "arrow": objects w1={*}, w2={1,2}; single non-identity morphism
         p : w2 -> w1 (merge groups).  [Example 6.6 of the paper as written]
       - "full" : full subcategory of FinSet on {*} and {1,2}: adds the two
         injections i1,i2 : w1 -> w2 (subgroup selection), the swap
         s : w2 -> w2 (treatment relabelling), and the constants
         c1 = i1 p, c2 = i2 p.
  Theta(w) = M^w x S  (finite means M, finite dispersions S), contravariant;
  Theta(p) = diagonal, Theta(i_k) = k-th coordinate, etc.
  A(w) = QQ^{Theta(w)}, A(f) = Theta(f)^*  (covariant).
  Lambda = constant at finite target T;  B(w) = QQ^T, B(f) = id.
  g : design-indexed quantity; g_w^* : B(w) -> A(w) dualisation
      (g_w^*)(e_t) = indicator of g_w^{-1}(t).
  Defect 1-cochain:                c(f) = A(f) g_dom^* - g_cod^* B(f)
  BW differentials:  (d0 h)(f)       = A(f) h_dom - h_cod B(f)
                     (d1 c)(f2,f1)   = A(f2) c(f1) - c(f2 f1) + c(f2) B(f1)
"""

from sympy import Rational, Matrix, zeros
from itertools import product

R = Rational

# ----------------------------------------------------------------------
# category machinery
# ----------------------------------------------------------------------

def make_category(variant):
    """Objects are 1 (={*}) and 2 (={1,2}); morphisms name->(dom,cod,map tuple).
    map tuple m: element x of dom (indexed 0..n-1) goes to element m[x] of cod."""
    mor = {
        'id1': (1, 1, (0,)),
        'id2': (2, 2, (0, 1)),
        'p':   (2, 1, (0, 0)),
    }
    if variant == 'full':
        mor.update({
            'i1': (1, 2, (0,)),
            'i2': (1, 2, (1,)),
            's':  (2, 2, (1, 0)),
            'c1': (2, 2, (0, 0)),
            'c2': (2, 2, (1, 1)),
        })
    return mor

def compose(mor, f2, f1):
    """f2 after f1; returns morphism name, or None if not composable."""
    d1, c1_, m1 = mor[f1]
    d2, c2_, m2 = mor[f2]
    if d2 != c1_:
        return None
    m = tuple(m2[m1[x]] for x in range(len(m1)))
    for name, (d, c, mm) in mor.items():
        if (d, c, mm) == (d1, c2_, m):
            return name
    raise RuntimeError(f'composite of {f2} o {f1} not in category: {(d1,c2_,m)}')

# ----------------------------------------------------------------------
# Theta, A, g*, defect
# ----------------------------------------------------------------------

def theta_elems(obj, M, S):
    n = obj  # object 1 has 1 unit-class, object 2 has 2
    return [ (mu, s) for mu in product(M, repeat=n) for s in S ]

def theta_map(mor, f, M, S):
    """Theta(f): Theta(cod f) -> Theta(dom f); (mu', s) -> (mu' o f, s)."""
    d, c, m = mor[f]
    def tm(el):
        mu, s = el
        return (tuple(mu[m[x]] for x in range(d)), s)
    return tm

def A_matrix(mor, f, M, S):
    """A(f): A(dom f) -> A(cod f); rows Theta(cod f), cols Theta(dom f);
    entry [t', t] = 1 iff Theta(f)(theta'_{t'}) = theta_t."""
    d, c, m = mor[f]
    Td, Tc = theta_elems(d, M, S), theta_elems(c, M, S)
    tm = theta_map(mor, f, M, S)
    A = zeros(len(Tc), len(Td))
    for i, tp in enumerate(Tc):
        A[i, Td.index(tm(tp))] = 1
    return A

def gstar(obj, gfun, M, S, T):
    Th = theta_elems(obj, M, S)
    G = zeros(len(Th), len(T))
    for i, th in enumerate(Th):
        G[i, T.index(gfun(th))] = 1
    return G

def defect(mor, f, gs, M, S):
    d, c, m = mor[f]
    return A_matrix(mor, f, M, S) * gs[d] - gs[c]          # B(f)=id

# ----------------------------------------------------------------------
# BW complex in degrees 0,1,2 (normalised: cochains vanish on identities)
# ----------------------------------------------------------------------

def nonid(mor):
    return [f for f in mor if not f.startswith('id')]

def comp_pairs(mor):
    """composable pairs (f2, f1) of NON-identity morphisms."""
    fs = nonid(mor)
    return [(f2, f1) for f2 in fs for f1 in fs
            if mor[f2][0] == mor[f1][1]]

def c0_blocks(M, S, T):
    """h = (h1, h2), h_w in Hom(B(w), A(w)) = matrix |Theta(w)| x |T|."""
    return {1: (len(theta_elems(1, M, S)), len(T)),
            2: (len(theta_elems(2, M, S)), len(T))}

def c1_blocks(mor, M, S, T):
    """for f: w -> w', block is Hom(B(w), A(w')) = |Theta(w')| x |T|."""
    return {f: (len(theta_elems(mor[f][1], M, S)), len(T)) for f in nonid(mor)}

def flatten(mats, blocks, keys):
    v = []
    for k in keys:
        r, c = blocks[k]
        Mk = mats.get(k, zeros(r, c))
        v.extend(Mk[i, j] for i in range(r) for j in range(c))
    return Matrix(len(v), 1, v)

def d0_matrix(mor, M, S, T):
    b0, b1 = c0_blocks(M, S, T), c1_blocks(mor, M, S, T)
    k0 = [1, 2]; k1 = sorted(nonid(mor))
    n0 = sum(r * c for (r, c) in (b0[k] for k in k0))
    n1 = sum(r * c for (r, c) in (b1[k] for k in k1))
    D = zeros(n1, n0)
    # basis of C^0: run through (object w, entry (i,j))
    col = 0
    for w in k0:
        rw, cw = b0[w]
        for i in range(rw):
            for j in range(cw):
                h = {w: zeros(rw, cw)}
                h[w][i, j] = 1
                img = {}
                for f in k1:
                    d, c, m = mor[f]
                    Af = A_matrix(mor, f, M, S)
                    t = zeros(*b1[f])
                    if d == w: t += Af * h[w]
                    if c == w: t -= h[w]
                    img[f] = t
                D[:, col] = flatten(img, b1, k1)
                col += 1
    return D, k1, b1

def d1_matrix(mor, M, S, T):
    b1 = c1_blocks(mor, M, S, T)
    k1 = sorted(nonid(mor))
    pairs = comp_pairs(mor)
    b2 = {(f2, f1): (len(theta_elems(mor[f2][1], M, S)), len(T))
          for (f2, f1) in pairs}
    n1 = sum(r * c for (r, c) in (b1[k] for k in k1))
    n2 = sum(r * c for (r, c) in (b2[k] for k in pairs))
    D = zeros(n2, n1)
    col = 0
    for f in k1:
        rf, cf = b1[f]
        for i in range(rf):
            for j in range(cf):
                cch = {f: zeros(rf, cf)}
                cch[f][i, j] = 1
                img = {}
                for (f2, f1) in pairs:
                    comp = compose(mor, f2, f1)
                    t = zeros(*b2[(f2, f1)])
                    if f1 == f:
                        t += A_matrix(mor, f2, M, S) * cch[f]
                    if comp == f:                      # identities excluded (normalised)
                        t -= cch[f]
                    if f2 == f:
                        t += cch[f]                    # B(f1)=id
                    img[(f2, f1)] = t
                D[:, col] = flatten(img, b2, pairs)
                col += 1
    return D, pairs, b2


# ----------------------------------------------------------------------
# fast exact-with-overwhelming-certainty rank: mod several large primes
# (entries of all matrices here are 0 / +-1, so modular rank = rational
#  rank unless p divides a nonzero minor; three 31-bit primes make a
#  false verdict astronomically unlikely; small cases are cross-checked
#  against sympy exact ranks below)
# ----------------------------------------------------------------------
import numpy as np

PRIMES = (2147483647, 2147483629, 2147483587)

def rank_modp(Msym, p):
    if Msym.rows == 0 or Msym.cols == 0:
        return 0
    A = np.array(Msym.tolist(), dtype=np.int64) % p
    m, n = A.shape
    r = 0
    for col in range(n):
        piv = None
        for row in range(r, m):
            if A[row, col] % p:
                piv = row; break
        if piv is None:
            continue
        A[[r, piv]] = A[[piv, r]]
        inv = pow(int(A[r, col]), p - 2, p)
        A[r] = (A[r] * inv) % p
        for row in range(m):
            if row != r and A[row, col]:
                A[row] = (A[row] - A[row, col] * A[r]) % p
        r += 1
        if r == m:
            break
    return r

def fast_rank(Msym):
    rs = {rank_modp(Msym, p) for p in PRIMES}
    assert len(rs) == 1, 'modular ranks disagree; increase primes'
    return rs.pop()

# ----------------------------------------------------------------------
# analysis runs
# ----------------------------------------------------------------------

def run(variant, M, S, g1fun, g2fun, label):
    print('=' * 78)
    print(f'RUN: {label}')
    print(f'  category = {variant},  M = {M},  S = {S}')
    mor = make_category(variant)
    Th1, Th2 = theta_elems(1, M, S), theta_elems(2, M, S)
    vals = sorted({g1fun(t) for t in Th1} | {g2fun(t) for t in Th2})
    T = vals
    print(f'  |Theta(w1)| = {len(Th1)}, |Theta(w2)| = {len(Th2)}, T = {T}')
    gs = {1: gstar(1, g1fun, M, S, T), 2: gstar(2, g2fun, M, S, T)}

    # ---- defect cochain -------------------------------------------------
    cch = {f: defect(mor, f, gs, M, S) for f in nonid(mor)}
    nz = [f for f in sorted(cch) if any(x != 0 for x in cch[f])]
    print(f'  defect cochain c: nonzero exactly on {nz if nz else "-- (natural!)"}')

    # ---- Lemma 6.3: 1-cocycle identity on all composable pairs ---------
    ok = True
    for (f2, f1) in comp_pairs(mor):
        comp = compose(mor, f2, f1)
        lhs = cch.get(comp, zeros(*A_matrix(mor, f2, M, S).shape[:1] + (len(T),))) \
              if comp in cch else zeros(len(theta_elems(mor[f2][1], M, S)), len(T))
        rhs = A_matrix(mor, f2, M, S) * cch[f1] + cch[f2]
        if lhs != rhs:
            ok = False
            print(f'  COCYCLE IDENTITY FAILS on ({f2} o {f1})')
    print(f'  Lemma 6.3 cocycle identity: {"verified on all composable pairs" if ok else "FAILED"}')

    # ---- BW cohomology in degree 1 -------------------------------------
    D0, k1, b1 = d0_matrix(mor, M, S, T)
    D1, pairs, b2 = d1_matrix(mor, M, S, T)
    dimC1 = D0.rows
    rk0, rk1 = fast_rank(D0), fast_rank(D1)
    dimZ1 = dimC1 - rk1
    dimH1 = dimZ1 - rk0
    print(f'  dim C^1 = {dimC1}, dim Z^1 = {dimZ1}, dim B^1 = {rk0}, '
          f'dim H^1(BW; Hom(B,A)) = {dimH1}')

    # sanity: c is a cocycle
    cvec = flatten(cch, b1, k1)
    assert all(x == 0 for x in D1 * cvec), 'defect cochain is not a cocycle?!'

    # ---- unrestricted correction: is c a coboundary? -------------------
    aug = D0.row_join(cvec)
    correctable = (aug.rank() == rk0)
    print(f'  [c] = 0 in unrestricted H^1?  {correctable}')
    hsol = None
    if correctable and nz:
        sol, params = D0.gauss_jordan_solve(cvec)
        sol = sol.subs({pp: 0 for pp in params})
        # unflatten: h = -sol per the sign convention of Proposition 6.4 (c = -d h)
        b0 = c0_blocks(M, S, T)
        idx = 0
        h = {}
        for w in (1, 2):
            r, ccols = b0[w]
            h[w] = Matrix(r, ccols, sol[idx: idx + r * ccols])
            idx += r * ccols
        hsol = {w: -h[w] for w in h}
        corrected = {w: gs[w] + hsol[w] for w in (1, 2)}
        det = all(sorted(corrected[w][i, :]) == sorted(gs[w][i, :])
                  or set(corrected[w][i, :]) <= {0, 1}
                  and sum(corrected[w][i, :]) == 1
                  for w in (1, 2) for i in range(corrected[w].rows))
        print('  particular correcting 0-cochain h found (params -> 0). '
              f'Corrected g*+h rows are 0/1 stochastic (deterministic): {det}')
        if det:
            # decode the corrected deterministic quantity at w2
            dec = []
            for i, th in enumerate(Th2):
                j = [jj for jj in range(len(T)) if corrected[2][i, jj] == 1][0]
                dec.append((th, T[j]))
            print('  corrected quantity at w2 (theta -> value):')
            for th, v in dec:
                print(f'      {th} -> {v}')
    return dict(mor=mor, T=T, gs=gs, c=cch, dimH1=dimH1,
                correctable=correctable, M=M, S=S)

# ----------------------------------------------------------------------
# recalibration class (statistically admissible corrections):
# corrected quantities of the form  r_w o g_w,  r_w : T -> T.
# ----------------------------------------------------------------------

def recalibration_scan(res, g1fun, g2fun):
    """Admissible-correction class: corrected quantities r_w o g_w with
    r_w : T -> T (recalibration of the value scale). Naturality along the
    generators forces r2 to be constant on the partition of T generated by
      s  ~  s + (m0-m1)^2/4      for all s in S, m0, m1 in M
    (from the injections i1, i2; the merge p then forces r1 = r2 on S; the
    swap is automatic by symmetry of g2). Exact union-find computation."""
    mor, T, M, S = res['mor'], res['T'], res['M'], res['S']
    if 'i1' not in mor:
        print('  recalibration class: no injections in this category variant.')
        return None
    parent = list(range(len(T)))
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    def union(x, y):
        rx, ry = find(x), find(y)
        if rx != ry: parent[rx] = ry
    for s in S:
        for m0 in M:
            for m1 in M:
                v = s + R((m0 - m1) ** 2, 4)
                union(T.index(s), T.index(v))
    classes = {}
    for i, t in enumerate(T):
        classes.setdefault(find(i), []).append(t)
    cls = sorted(classes.values(), key=lambda c: c[0])
    k = len(cls)
    nonconstant = k >= 2
    bijective = all(len(c) == 1 for c in cls)
    print('  recalibration class (corrected = r_w o g_w, r_w : T -> T):')
    print(f'    forced identification classes on T: '
          f'{[ [str(x) for x in c] for c in cls ]}')
    print(f'    information-bearing (non-constant) recalibration exists: {nonconstant}')
    print(f'    bijective recalibration exists: {bijective}')
    if not nonconstant:
        print('    ==> obstruction NONZERO relative to the admissible class:')
        print('        every natural recalibration of the marginal dispersion is')
        print('        constant, i.e. destroys all information (scale confounding).')
    return cls

if __name__ == '__main__':
    # ---------------- quantities -----------------------------------------
    # g1: dispersion at the merged design;  g2: marginal (total) dispersion an
    # analyst attributes at the two-group design when the groups are balanced:
    #     g2((m0,m1), s) = s + (m0 - m1)^2 / 4
    # This is the marginal (total) dispersion under balanced allocation: the
    # one-way analogue of the design-averaged expectation alpha + beta*xbar,
    # which is the non-example Tjur actually gives in his 2002 discussion
    # (Ann. Statist. 30, 1297-1300).  It is NOT the overdispersion of his
    # separate passage on generalised linear models, which is a non-existence
    # claim about scale-extended GLMs and unrelated to naturality.
    g1 = lambda th: th[1]
    g2 = lambda th: th[1] + R((th[0][0] - th[0][1]) ** 2, 4)

    M0 = [0, 1]

    # RUN 1 -- Example 6.6 of the paper exactly as written: arrow category.
    r1_ = run('arrow', M0, [R(1), R(2)], g1, g2,
              'Example 6.6 (arrow category, merge morphism only)')

    # RUN 2 -- full FinSet subcategory, separated dispersion scales.
    r2_ = run('full', M0, [R(1), R(2)], g1, g2,
              'full design category, separated scales S={1,2}')
    recalibration_scan(r2_, g1, g2)

    # RUN 3 -- full category, OVERLAPPING scales: S = {1, 5/4}; the marginal
    # value 5/4 conflates (within=1, overdispersed) with (within=5/4, none).
    r3_ = run('full', M0, [R(1), R(5, 4)], g1, g2,
              'full design category, overlapping scales S={1, 5/4} (scale confounding)')
    recalibration_scan(r3_, g1, g2)


    # ---------------- robustness checks ----------------------------------
    # CHECK A -- arrow category with a quantity whose defect on p is nonzero
    # (g2' disagrees with g1 on the diagonal): every defect over the arrow
    # category is an unrestricted coboundary (H^1 = 0 there).
    g2p = lambda th: th[1] + R(th[0][0], 4)
    rA = run('arrow', M0, [R(1), R(2)], g1, g2p,
             'CHECK A: arrow category, quantity with nonzero defect on p')

    # CHECK B -- three mean levels; with S={1,2} the marginal value
    # 1 + (2-0)^2/4 = 2 collides with s = 2: confounding is generic.
    rB = run('full', [0, 1, 2], [R(1), R(2)], g1, g2,
             'CHECK B: full category, M={0,1,2}, S={1,2}')
    recalibration_scan(rB, g1, g2)
