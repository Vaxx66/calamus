"""Exact (over Q) Baues-Wirsching H^1 of the two-object design scheme O_full / O_arr
of F. Vaccarino, "Statistical models as natural transformations: a bridge from the
McCullagh--Brons hexagon to Markov categories" (2026), Examples 6.6 and 6.8, with Lambda constant
at T = image of g and B(psi) = id.  Because B(psi)=id and B(omega)=k^T, the complex with
coefficients Hom(B(dom),A(cod)) is the |T|-fold direct sum of the complex with
coefficients psi |-> A(cod psi); we compute the latter and multiply by |T|."""
from fractions import Fraction
from itertools import product
import sympy as sp

def theta(n, M, S):
    return [(mu, s) for mu in product(M, repeat=n) for s in S]

def build(M, S, full=True):
    T1, T2 = theta(1, M, S), theta(2, M, S)
    idx1 = {t: i for i, t in enumerate(T1)}
    idx2 = {t: i for i, t in enumerate(T2)}
    # Theta(psi): Theta(cod) -> Theta(dom)   (contravariant)
    Th = {
        'p':  ('w2', 'w1', lambda t: ((t[0][0], t[0][0]), t[1])),          # p: w2->w1, Theta(p): Theta(w1)->Theta(w2)
        'i1': ('w1', 'w2', lambda t: ((t[0][0],), t[1])),                  # Theta(i1): Theta(w2)->Theta(w1)
        'i2': ('w1', 'w2', lambda t: ((t[0][1],), t[1])),
        's':  ('w2', 'w2', lambda t: ((t[0][1], t[0][0]), t[1])),
        'c1': ('w2', 'w2', lambda t: ((t[0][0], t[0][0]), t[1])),
        'c2': ('w2', 'w2', lambda t: ((t[0][1], t[0][1]), t[1])),
    }
    if not full:
        Th = {'p': Th['p']}
    obj = {'w1': (T1, idx1), 'w2': (T2, idx2)}
    # A(psi): A(dom) -> A(cod), (A(psi) f)(theta) = f(Theta(psi) theta); matrix [theta, theta0]
    def Amat(name):
        dom, cod, f = Th[name]
        Td, idxd = obj[dom]; Tc, idxc = obj[cod]
        Mx = sp.zeros(len(Tc), len(Td))
        for th in Tc:
            Mx[idxc[th], idxd[f(th)]] = 1
        return Mx
    A = {n: Amat(n) for n in Th}
    dims = {'w1': len(T1), 'w2': len(T2)}
    # composition table (psi2 after psi1): returns name or 'id'
    comp = {}
    for a in Th:
        for b in Th:
            if Th[b][1] != Th[a][0]:  # cod(b) must equal dom(a)  for a∘b
                continue
            # compute composite by its Theta action: Theta(a∘b) = Theta(b)∘Theta(a)
            fa, fb = Th[a][2], Th[b][2]
            dom, cod = Th[b][0], Th[a][1]
            Tc = obj[cod][0]
            act = tuple(fb(fa(t)) for t in Tc)
            if dom == cod and act == tuple(Tc):
                comp[(a, b)] = 'id'
            else:
                found = None
                for n in Th:
                    if Th[n][0] == dom and Th[n][1] == cod and tuple(Th[n][2](t) for t in Tc) == act:
                        found = n
                assert found is not None, (a, b)
                comp[(a, b)] = found
    names = list(Th)
    # C^0 = A(w1) (+) A(w2); C^1 = (+)_psi A(cod psi); C^2 = (+)_{(a,b)} A(cod a)
    off1, pos = {}, 0
    for n in names:
        off1[n] = pos; pos += dims[Th[n][1]]
    dimC1 = pos
    off0 = {'w1': 0, 'w2': dims['w1']}; dimC0 = dims['w1'] + dims['w2']
    d0 = sp.zeros(dimC1, dimC0)
    for n in names:
        dom, cod, _ = Th[n]
        r, k = off1[n], dims[cod]
        d0[r:r+k, off0[dom]:off0[dom]+dims[dom]] += A[n]            # A(psi) h_dom
        d0[r:r+k, off0[cod]:off0[cod]+dims[cod]] -= sp.eye(k)       # - h_cod
    pairs = list(comp)
    off2, pos = {}, 0
    for pr in pairs:
        off2[pr] = pos; pos += dims[Th[pr[0]][1]]
    dimC2 = pos
    d1 = sp.zeros(dimC2, dimC1)
    for (a, b) in pairs:
        cod = Th[a][1]; k = dims[cod]; r = off2[(a, b)]
        # A(a) c(b)
        d1[r:r+k, off1[b]:off1[b]+dims[Th[b][1]]] += A[a]
        # - c(a∘b)
        ab = comp[(a, b)]
        if ab != 'id':
            d1[r:r+k, off1[ab]:off1[ab]+k] -= sp.eye(k)
        # + c(a) B(b) = c(a)
        d1[r:r+k, off1[a]:off1[a]+k] += sp.eye(k)
    assert (d1 * d0).is_zero_matrix, "d1 d0 != 0"
    rk0, rk1 = d0.rank(), d1.rank()
    z1 = dimC1 - rk1
    h1 = z1 - rk0
    return dict(dimC0=dimC0, dimC1=dimC1, dimC2=dimC2, rank_d0=rk0, rank_d1=rk1, dimZ1=z1, dimH1_reduced=h1, pairs=len(pairs))

def Tset(M, S):
    return sorted({Fraction(s) + Fraction((m1 - m2) ** 2, 4) for m1 in M for m2 in M for s in S})

for M, S in [((0, 1), (1, 2)), ((0, 1, 2), (1, 2)), ((0, 1), (1, Fraction(5, 4))), ((0,1,2,3),(1,2))]:
    T = Tset(M, S)
    r = build(M, S, full=True)
    ra = build(M, S, full=False)
    print(f"M={M} S={tuple(str(x) for x in S)} |T|={len(T)} T={[str(t) for t in T]}")
    print("  O_full:", r, " => dim H^1 (all t) =", len(T) * r['dimH1_reduced'])
    print("  O_arr :", ra, " => dim H^1 (all t) =", len(T) * ra['dimH1_reduced'], " dim C^1(all t) =", len(T)*ra['dimC1'])
