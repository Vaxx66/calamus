#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Companion script to:
  F. Vaccarino, "Statistical models as natural transformations:
  meaningfulness, coherence and priors as states in Markov categories"
  (2026), Section 8 (Propositions 8.4 and 8.5).
One-way layout: coherent priors as natural families of states.

Companion to oneway_obstruction.py; same conventions (Example 6.7 of the
paper), exact arithmetic over QQ.

  Theta(w_n) = M^n x S,  contravariant.
  Design-reduction (insertion) morphisms i_1, i_2 : w1 -> w2 induce
      Theta(i_k) : Theta(w2) -> Theta(w1),  ((mu_1,mu_2), s) |-> (mu_k, s).
  Merge p : w2 -> w1 induces
      Theta(p) : Theta(w1) -> Theta(w2),   (mu, s) |-> ((mu,mu), s).

  A PRIOR is a family of states pi_w : I ~> Theta(w), i.e. a probability
  vector on Theta(w).
  Naturality of pi : Delta_I => delta Theta as functors Omega^op -> Stoch
  reads, for psi : w -> w' in Omega,
      pi_w = Theta(psi)_* pi_{w'} .

  Four strata of the two-object scheme (identifying w1 = {*} with {1}):
      Omega_red  = <i_1>             inclusions only (Tjur's class, Section 5
                                     of the paper)
      Omega_sel  = <i_1, i_2>        both selections into the single one-group
                                     object: identifies {1} and {2} with w1
      Omega_inj  = <i_1, i_2, sigma> all injections (McCullagh's class as
                                     stated); sigma = swap of the two groups
      Omega_full = all eight morphisms, merge p included (degenerate)
  COHERENT over a stratum = natural on it.  Closed forms (Proposition 8.4):
      red  : |M|^2|S| - 1
      sel  : |M|^2|S| - 1 - |S|(|M|-1)
      inj  : |M|(|M|+1)|S|/2 - 1
      full : |M||S| - 1
  sel and inj coincide iff |M| <= 2; their difference is |S| binom(|M|-1, 2).

  Jeffreys analogue (Proposition 8.5): pi_{w_n} = unif(M^n) tensor rho_n with
  rho_n(s) proportional to s^{-(n+2)/2}.  Its defect along i_1, i_2 is the
  total variation between rho_1 and rho_2; its defect along the merge p is
      (1/2) [ (1 - 1/|M|) + sum_s | rho_1(s) - rho_2(s)/|M| | ]  >=  1 - 1/|M|,
  with equality iff  Z_1/Z_2 <= |M| sqrt(s_min),  Z_n = sum_s s^{-(n+2)/2}.

Runs under plain python3 or `sage -python`.
"""

from sympy import Rational, Matrix, zeros, nsimplify
from itertools import product

R = Rational


# ----------------------------------------------------------------------
# Theta and pushforward matrices
# ----------------------------------------------------------------------

def theta_elems(n, M, S):
    """Theta(w_n) = M^n x S, as an ordered list."""
    return [(mu, s) for mu in product(M, repeat=n) for s in S]


def pushforward_matrix(f, M, S, dom_n, cod_n):
    """Matrix of pushforward along the measurable map
       f : Theta(w_{dom_n}) -> Theta(w_{cod_n}).
    Rows indexed by Theta(w_{cod_n}), cols by Theta(w_{dom_n});
    entry [t', t] = 1 iff f(t) = t'.  So (f_* pi) = P pi."""
    Td, Tc = theta_elems(dom_n, M, S), theta_elems(cod_n, M, S)
    P = zeros(len(Tc), len(Td))
    for j, t in enumerate(Td):
        P[Tc.index(f(t)), j] = 1
    return P


def theta_i(k):
    """Theta(i_k) : Theta(w2) -> Theta(w1)."""
    return lambda t: ((t[0][k],), t[1])


def theta_sigma(t):
    """Theta(sigma) : Theta(w2) -> Theta(w2), swap of the two mean coordinates."""
    (m1, m2), s = t
    return ((m2, m1), s)


def theta_p(t):
    """Theta(p) : Theta(w1) -> Theta(w2)."""
    return ((t[0][0], t[0][0]), t[1])


# ----------------------------------------------------------------------
# the coherence system
# ----------------------------------------------------------------------

STRATA = ('red', 'sel', 'inj', 'full')

PREDICTED = {
    'red':  lambda m, s: m * m * s - 1,
    'sel':  lambda m, s: m * m * s - 1 - s * (m - 1),
    'inj':  lambda m, s: m * (m + 1) * s // 2 - 1,
    'full': lambda m, s: m * s - 1,
}


def coherence_system(M, S, stratum):
    """Linear system in the concatenated unknown (pi_1 | pi_2).
    Returns (A, b, n1, n2) with A x = b the coherence identities
        pi_w = Theta(psi)_* pi_{w'}     (psi : w -> w' a generator of the stratum)
    together with the two normalisations (positivity is handled separately:
    the uniform prior is an interior feasible point of every stratum).
    Generators: red <i_1>, sel <i_1,i_2>, inj <i_1,i_2,sigma>,
    full <i_1,i_2,sigma,p>; the constants c_k = i_k p follow by pasting."""
    assert stratum in STRATA
    n1, n2 = len(theta_elems(1, M, S)), len(theta_elems(2, M, S))
    rows, rhs = [], []

    # insertions / selections  i_k : w1 -> w2,   pi_1 = Theta(i_k)_* pi_2
    for k in ((0,) if stratum == 'red' else (0, 1)):
        P = pushforward_matrix(theta_i(k), M, S, 2, 1)      # n1 x n2
        for r in range(n1):
            row = [0] * (n1 + n2)
            row[r] = 1
            for c in range(n2):
                row[n1 + c] -= P[r, c]
            rows.append(row)
            rhs.append(0)

    # relabelling  sigma : w2 -> w2,   pi_2 = Theta(sigma)_* pi_2
    if stratum in ('inj', 'full'):
        P = pushforward_matrix(theta_sigma, M, S, 2, 2)      # n2 x n2
        for r in range(n2):
            row = [0] * (n1 + n2)
            row[n1 + r] += 1
            for c in range(n2):
                row[n1 + c] -= P[r, c]
            rows.append(row)
            rhs.append(0)

    # merge  p : w2 -> w1,   pi_2 = Theta(p)_* pi_1
    if stratum == 'full':
        P = pushforward_matrix(theta_p, M, S, 1, 2)          # n2 x n1
        for r in range(n2):
            row = [0] * (n1 + n2)
            row[n1 + r] = 1
            for c in range(n1):
                row[c] -= P[r, c]
            rows.append(row)
            rhs.append(0)

    rows.append([1] * n1 + [0] * n2); rhs.append(1)
    rows.append([0] * n1 + [1] * n2); rhs.append(1)

    return Matrix(rows), Matrix(rhs), n1, n2


def solution_dimension(A, b):
    """Affine dimension of {x : A x = b}; -1 if infeasible."""
    aug = A.row_join(b)
    if aug.rank() != A.rank():
        return -1
    return A.cols - A.rank()


# ----------------------------------------------------------------------
# named prior rules
# ----------------------------------------------------------------------

def uniform_prior(n, M, S):
    Th = theta_elems(n, M, S)
    return Matrix([R(1, len(Th))] * len(Th))


def jeffreys_analogue(n, M, S):
    """Finite analogue of Jeffreys' rule for the Gaussian one-way layout.
    For N(X beta, sigma^2 I) with p mean parameters, Jeffreys' rule gives
        p(beta, sigma^2)  ~  (sigma^2)^{-(p+2)/2},  flat in beta.
    Here p = n = |w|, and s plays the role of sigma^2.  Proper version:
    uniform in mu, weight s^{-(n+2)/2} in s, then normalised."""
    Th = theta_elems(n, M, S)
    w = [nsimplify(s) ** R(-(n + 2), 2) for (_, s) in Th]
    Z = sum(w)
    return Matrix([wi / Z for wi in w])


def defect_report(name, pri, M, S):
    """Check coherence at the selections i_1, i_2 (stratum Omega_sel) of a
    prior RULE n |-> pi_{w_n}.  For the product rules tested here this is the
    same as coherence at i_1 alone (Omega_red) and at sigma as well
    (Omega_inj), since a product prior with uniform means is symmetric."""
    p1, p2 = pri(1, M, S), pri(2, M, S)
    out = {}
    for k in (0, 1):
        P = pushforward_matrix(theta_i(k), M, S, 2, 1)
        d = P * p2 - p1
        out[f'i{k+1}'] = d
    tv = {k: sum(abs(x) for x in d) / 2 for k, d in out.items()}
    coherent = all(all(x == 0 for x in d) for d in out.values())
    print(f'    {name}: coherent at i1,i2 = {coherent}'
          + ('' if coherent else f',  TV defect = '
             + ', '.join(f'{k}: {v}' for k, v in tv.items())))
    return coherent, out, tv


# ----------------------------------------------------------------------
# runs
# ----------------------------------------------------------------------

def run(M, S, label):
    print('=' * 78)
    print(f'RUN: {label}')
    print(f'  M = {M},  S = {[str(s) for s in S]}')
    n1, n2 = len(theta_elems(1, M, S)), len(theta_elems(2, M, S))
    print(f'  |Theta(w1)| = {n1}, |Theta(w2)| = {n2}; '
          f'ambient dim of (pi_1, pi_2) = {n1 + n2}')

    dims = {}
    for st in STRATA:
        A, b, _, _ = coherence_system(M, S, st)
        dims[st] = solution_dimension(A, b)
        pred = PREDICTED[st](len(M), len(S))
        head = f'Omega_{st}-coherent priors:'
        print(f'  {head:<28} affine dimension = {dims[st]:>3}'
              f'   (Proposition 8.4 predicts {pred})')
    print(f'  (unconstrained probability pairs: {n1 - 1 + n2 - 1})')

    # what does the merge condition force?
    Th2 = theta_elems(2, M, S)
    diag = [i for i, ((m1, m2), s) in enumerate(Th2) if m1 == m2]
    print(f'  (diagonal subset of Theta(w2) has {len(diag)} of {n2} points; '
          f'off-diagonal mass forced to 0 by the merge)')

    print('  named prior rules:')
    defect_report('uniform             ', uniform_prior, M, S)
    defect_report('Jeffreys analogue   ', jeffreys_analogue, M, S)

    # exchangeability forced by i1, i2
    P1 = pushforward_matrix(theta_i(0), M, S, 2, 1)
    P2 = pushforward_matrix(theta_i(1), M, S, 2, 1)
    Ex = P1 - P2
    print(f'  marginals-agree constraint (i1 vs i2) has rank {Ex.rank()} '
          f'on Theta(w2)-space of dim {n2}')
    return dims


def product_prior(n, M, S, rho):
    """(uniform on M^n) tensor rho[n], with rho[n] a probability vector on S."""
    Th = theta_elems(n, M, S)
    return Matrix([R(1, len(M)) ** n * rho[n][S.index(s)] for (_, s) in Th])


def closed_form_check():
    print('=' * 78)
    print('CLOSED-FORM CHECK (Proposition 8.4), |M| = 2,3,4 and |S| = 1,2,3')
    print('  predicted: red  |M|^2|S| - 1')
    print('             sel  |M|^2|S| - 1 - |S|(|M|-1)')
    print('             inj  |M|(|M|+1)|S|/2 - 1')
    print('             full |M||S| - 1')
    ok = True
    for m in (2, 3, 4):
        for s in (1, 2, 3):
            M, S = list(range(m)), [R(k + 1) for k in range(s)]
            got = {st: solution_dimension(*coherence_system(M, S, st)[:2])
                   for st in STRATA}
            pred = {st: PREDICTED[st](m, s) for st in STRATA}
            gap = got['sel'] - got['inj']
            good = (got == pred and gap == s * (m - 1) * (m - 2) // 2)
            ok &= good
            print(f'    |M|={m} |S|={s}:  '
                  + '  '.join(f'{st} {got[st]} (pred {pred[st]})' for st in STRATA)
                  + f'   sel-inj = {gap}   {"OK" if good else "MISMATCH"}')
    print(f'  closed forms: {"confirmed" if ok else "FAILED"}')


def product_prior_check():
    print('=' * 78)
    print('PRODUCT-PRIOR CRITERION  (pi_w = unif(M^w) tensor rho_w)')
    M, S = [0, 1], [R(1), R(4)]
    for name, rho in [('rho design-independent',
                       {1: [R(3, 5), R(2, 5)], 2: [R(3, 5), R(2, 5)]}),
                      ('rho design-dependent  ',
                       {1: [R(3, 5), R(2, 5)], 2: [R(1, 2), R(1, 2)]})]:
        p1, p2 = product_prior(1, M, S, rho), product_prior(2, M, S, rho)
        tv = []
        for k in (0, 1):
            d = pushforward_matrix(theta_i(k), M, S, 2, 1) * p2 - p1
            tv.append(sum(abs(x) for x in d) / 2)
        print(f'    {name}: coherent = {all(t == 0 for t in tv)},  '
              f'TV = {[str(t) for t in tv]}')


def _tv(p, q):
    keys = set(p) | set(q)
    return sum(abs(p.get(k, 0) - q.get(k, 0)) for k in keys) / 2


def jeffreys_dict(n, M, S):
    """The Jeffreys analogue as a dict t -> mass (S consists of squares of
    rationals, given by their square roots, so that everything stays in QQ)."""
    w = {s: R(1) / (s ** (n + 2)) for s in S}           # s = sqrt(sigma^2)
    Z = sum(w.values())
    Th = [(mu, s) for mu in product(M, repeat=n) for s in S]
    return {(mu, s): R(1, len(M) ** n) * w[s] / Z for (mu, s) in Th}, \
           {s: w[s] / Z for s in S}, Z


def closed_form_tv_check():
    """Proposition 8.5(2).  The defect of the Jeffreys analogue along i_1, i_2
    lives in the dispersion marginal: it equals TV(rho_1, rho_2) and does not
    depend on M.  The defect along the merge p is
        (1/2) [ (1 - 1/|M|) + sum_s | rho_1(s) - rho_2(s)/|M| | ]  >= 1 - 1/|M|,
    with equality iff Z_1/Z_2 <= |M| sqrt(s_min).  The last case below is one
    where the inequality is strict."""
    print('=' * 78)
    print('CLOSED-FORM DEFECT CHECK (Jeffreys analogue, Proposition 8.5(2))')
    print('  S is listed through the square roots r of its elements s = r^2')
    roots_strict = [R(1, 2)] + [R(k, 100) for k in range(105, 230)]
    cases = [([0, 1], [R(1), R(2)], 'S={1,4}'),
             ([0, 1], [R(1), R(2), R(3)], 'S={1,4,9}'),
             ([0, 1, 2], [R(1), R(2)], 'S={1,4}'),
             ([0, 1], [R(1), R(3, 2)], 'S={1,9/4}'),
             ([0, 1], [R(1), R(2), R(3), R(4)], 'S={1,4,9,16}'),
             ([0, 1], roots_strict, 'S={1/4} u {(k/100)^2 : 105<=k<=229}')]
    ok = True
    for M, roots, name in cases:
        m = len(M)
        p1, rho1, Z1 = jeffreys_dict(1, M, roots)
        p2, rho2, Z2 = jeffreys_dict(2, M, roots)
        tv_i = []
        for k in (0, 1):
            push = {}
            for ((mu, s), x) in p2.items():
                key = ((mu[k],), s)
                push[key] = push.get(key, 0) + x
            tv_i.append(_tv(push, p1))
        push_p = {((mu[0], mu[0]), s): x for ((mu, s), x) in p1.items()}
        tv_p = _tv(push_p, p2)
        pred_i = sum(abs(rho1[s] - rho2[s]) for s in roots) / 2
        pred_p = (1 - R(1, m) + sum(abs(rho1[s] - rho2[s] / m) for s in roots)) / 2
        bound = 1 - R(1, m)
        attained = (Z1 / Z2 <= m * min(roots))
        good = (tv_i[0] == tv_i[1] == pred_i and tv_p == pred_p
                and tv_p >= bound and ((tv_p == bound) == attained))
        ok &= good
        show = (lambda x: str(x) if len(str(x)) < 40 else f'{float(x):.6f}...')
        print(f'    |M|={m}, {name}: TV(i1)=TV(i2)={show(tv_i[0])} '
              f'(pred {show(pred_i)}), TV(p)={show(tv_p)} (pred {show(pred_p)}), '
              f'bound 1-1/|M|={bound} '
              f'{"attained" if tv_p == bound else "STRICT"}'
              f'   {"OK" if good else "MISMATCH"}')
    print(f'  defect formulas: {"confirmed" if ok else "FAILED"}')


if __name__ == '__main__':
    # Dispersion scales are taken to be squares of rationals so that the
    # Jeffreys analogue s^{-(n+2)/2} stays inside QQ.  The computation on the
    # priors does not depend on the scale confounding of Example 6.10(iii),
    # which concerns the quantity and not the prior.
    run([0, 1], [R(1), R(4)], 'M={0,1}, S={1,4}          (base case)')
    run([0, 1], [R(1), R(4), R(9)], 'M={0,1}, S={1,4,9}')
    run([0, 1, 2], [R(1), R(4)], 'M={0,1,2}, S={1,4}')
    run([0, 1], [R(1)], 'M={0,1}, S={1}            (degenerate control)')
    closed_form_check()
    product_prior_check()
    closed_form_tv_check()
