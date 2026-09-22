# calamus

Computational supplement and manuscript for:

> F. Vaccarino, *Statistical models as natural transformations: meaningfulness, coherence and priors as states in Markov categories*. Politecnico di Torino, 2026. arXiv: **[identifier pending]**.

**Status (release `v34`, September 2026).** The manuscript (38 pages) is deposited in `manuscript/`. The arXiv identifier will be added here as soon as it is announced. This release is the one cited in the manuscript, and every number quoted there is read from the outputs committed here. The repository is the stable reference for the computational supplement of the CALAMUS research programme (*Categorical Languages for Models and Uncertainty in Statistics*).

## Abstract

We show that a statistical model in the sense of McCullagh, in the form given by Brøns, is a natural transformation between two functors from the category of designs to the Kleisli category Stoch of the Giry monad, provided that its components are measurable in the parameter. The condition is empty for finite models. A design-indexed quantity is a family of morphisms of Stoch defined on the parameter objects, called meaningful if it is natural. We prove that Tjur's criterion, imposed on parameter functions indexed by finite samples with multiplicities, forces the indexing by the support and then coincides with naturality over the insertions. For finite designs we show that a quantity can be corrected to a natural one within a given class of corrections if and only if a class vanishes in the first cohomology group of a Baues–Wirsching complex relative to that class, while its image in the absolute group is always zero. In the one-way layout the marginal dispersion is not meaningful, and the within-group dispersion is its unique correction which leaves the merged design unchanged. A prior is a family of states on the parameter objects, called coherent over a class of design morphisms if it is natural over that class. We show that coherence at a merge confines the prior to the image of the corresponding parameter map, that coherence over the insertions is Kolmogorov consistency, and that coherence over the injections adds the exchangeability assumed by the categorical de Finetti theorem. In the finite one-way scheme the coherent priors form polytopes of known dimension. The analogue of Jeffreys' general rule is not coherent, while the analogue of his rule for location–scale families is. Finally we show that ridge regression is the Bayesian inversion of the Gaussian linear model with respect to a Gaussian prior, which is coherent over the insertions and never over the injections.

## Contents

### Manuscript (`manuscript/`)

- `calamus_arkiv_draft_light_v02.pdf` — the manuscript (38 pages).
- `calamus_arkiv_draft_light_v02.tex`, `calamus_arkiv_draft_light_v02.bbl` — the source pair, identical to the one posted on arXiv; it builds with `pdflatex` alone (three passes, no BibTeX).
- `calamus_arkiv_draft_light_v02.bib` — the bibliography, exactly the 28 entries cited. The full chain `pdflatex → bibtex → pdflatex ×2` reproduces the committed `.bbl` byte for byte.

### Code for Section 6 (cohomological obstruction, one-way layout)

- `oneway_obstruction.py` — defect cochain, cocycle identity, Baues–Wirsching complex, corrections and recalibrations for the two-object design scheme of Example 6.7. Exact arithmetic over ℚ; about 10 s.
- `oneway_obstruction_output.txt` — committed output.
- `h1_check.py` — independent computation of H¹ through the reduced coefficient system ψ ↦ A(cod ψ) of Example 6.12, for both categories and |M| = 2, 3, 4; about 3 s.
- `h1_check_output.txt` — committed output.

### Code for Section 8 (priors as states)

- `oneway_prior_coherence.py` — the four strata of Proposition 8.4 and the defects of the Jeffreys analogue of Proposition 8.5. Exact arithmetic over ℚ; about 2 s.
- `oneway_prior_coherence_output.txt` — committed output.

## What the scripts compute

The notation is the one of the manuscript.

### `oneway_obstruction.py`

Input: a category (`O_arr` or `O_full`), two finite sets `M` (means) and `S` (dispersions), and the pair `(g_{ω1}, g_{ω2})`. For each run the script prints, in this order:

- **(a)** the defect cochain `c(ψ) = A(ψ) g*_ω − g*_ω' B(ψ)` of eq. (13) and the set of the morphisms on which it is not zero;
- **(b)** a check, entry by entry, of the cocycle identity of Lemma 6.3 on every composable pair of morphisms different from the identities (22 pairs in `O_full`, none in `O_arr`). This is a test of the implementation, not of the lemma;
- **(c)** `dim C¹`, `dim Z¹`, `dim B¹` and `dim H¹`. The two ranks are computed modulo three primes of 31 bits (all the entries of the matrices are 0 or ±1), with the same result, and they are computed again exactly over ℚ;
- **(d)** the verdict `[c] = 0` in `H¹`, which holds by Proposition 6.4(2), together with a particular solution `h` of `c = −d⁰h`, obtained by setting to zero the free parameters of the row-reduced system, a test that `g* + h` is again deterministic (rows with entries 0 and 1), and the table of the corrected quantity at `ω2`;
- **(e)** for the admissible class of Example 6.10(iii), the partition of `T` forced by naturality (Proposition 6.9), computed by union–find, and whether there is a natural recalibration which is not constant, or which is bijective.

The five runs carry in the output the labels used in Section 6.2 of the manuscript.

| label | category | M | S | quantity | used in |
|---|---|---|---|---|---|
| `RUN 1` | `O_arr` | {0,1} | {1,2} | marginal dispersion | Example 6.10(i) |
| `RUN 2` | `O_full` | {0,1} | {1,2} | marginal dispersion | Example 6.10(ii), (iii) |
| `RUN 3` | `O_full` | {0,1} | {1, 5/4} | marginal dispersion | Example 6.10(iii) |
| `CHECK A` | `O_arr` | {0,1} | {1,2} | `g'_{ω2}((μ1,μ2),s) = s + μ1/4`, which does not agree with `g_{ω1}` on the diagonal | Example 6.10(i) |
| `CHECK B` | `O_full` | {0,1,2} | {1,2} | marginal dispersion | Examples 6.10(iii), 6.12 |

### `h1_check.py`

Since `B(ψ) = id` and `B(ω) = k^T`, the complex with coefficients `Hom(B(dom), A(cod))` is the direct sum of `|T|` copies of the complex with coefficients `ψ ↦ A(cod ψ)`. The script computes the latter exactly over ℚ and multiplies by `|T|`. Its results must agree with output (c) of `oneway_obstruction.py`: `dim H¹ = 0` for `O_arr` and for `O_full` with `|M| = 2`, `dim H¹ = 2·5 = 10` for `|M| = 3` and `dim H¹ = 6·7 = 42` for `|M| = 4` (Example 6.12).

### `oneway_prior_coherence.py`

For each of the four strata `Ω_red = ⟨i₁⟩ ⊊ Ω_sel = ⟨i₁,i₂⟩ ⊊ Ω_inj = ⟨i₁,i₂,σ⟩ ⊊ Ω_full` of Proposition 8.4 the script builds the linear system in the unknown `(π_{ω1}, π_{ω2})` given by the identities `π_ω = Θ(ψ)_* π_{ω'}` along the generators of the stratum and by the two normalisations, and it computes the affine dimension of the set of the solutions. It then

- compares the four dimensions, and the difference `sel − inj = |S|·binom(|M|−1, 2)`, with the closed forms of Proposition 8.4 for `|M| = 2, 3, 4` and `|S| = 1, 2, 3`;
- checks the product-prior criterion of Proposition 8.5(1);
- computes exactly the total variation defects of the Jeffreys analogue along `i₁`, `i₂` and the merge `p`, and compares them with the closed forms of Proposition 8.5(2) in six cases. In the first five the bound `1 − 1/|M|` for the merge is attained. In the sixth, `|M| = 2` and `S = {1/4} ∪ {(k/100)² : 105 ≤ k ≤ 229}`, it is strict (`0.5078… > 1/2`). The base case `M = {0,1}`, `S = {1,4}` gives the marginals `(8/9, 1/9)`, `(16/17, 1/17)` and the defect `8/153` quoted in the manuscript.

## Running

Python ≥ 3.10 with `sympy` and `numpy`:

```
python3 oneway_obstruction.py      > oneway_obstruction_output.txt
python3 h1_check.py                > h1_check_output.txt
python3 oneway_prior_coherence.py  > oneway_prior_coherence_output.txt
```

The three scripts are self-contained and deterministic, and a fresh run reproduces the committed outputs byte for byte. Last reproduction check: 19 September 2026, with Python 3.12.3, sympy 1.14.0 and numpy 2.4.4.

## Changes with respect to the previous deposit (v33k)

- New manuscript: the Bayesian part (priors as states, Jeffreys' rules, ridge regression) is now in the same paper, Tjur's criterion is treated for sample-indexed quantities, and the obstruction is a class in a relative complex.
- New script `oneway_prior_coherence.py` with its committed output.
- `oneway_obstruction.py`: the headers of the five runs are now `RUN 1`, `RUN 2`, `RUN 3`, `CHECK A`, `CHECK B`, as in the manuscript, and the docstring follows the new numbering. No computation has changed: the committed output differs from the previous one in those five header lines only.
- `h1_check.py`: docstring only; the output is unchanged.

## License

MIT for the code in this repository (see `LICENSE`). The manuscript in `manuscript/` is © 2026 Francesco Vaccarino and is **not** covered by the MIT license; it is distributed under the terms selected for the arXiv posting.

## Citing

```bibtex
@misc{vaccarino_statmodels_2026,
  author        = {Vaccarino, Francesco},
  title         = {Statistical models as natural transformations: meaningfulness, coherence and priors as states in {Markov} categories},
  year          = {2026},
  eprint        = {2609.22929},
  archivePrefix = {arXiv},
  primaryClass  = {math.ST},
  note          = {Computational supplement: \url{https://github.com/Vaxx66/calamus}, release v34}
}
```

See also `CITATION.cff`.
