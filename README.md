# calamus

Computational supplement and manuscript for:

> F. Vaccarino, *Statistical models as natural transformations: a bridge from the McCullagh–Brøns hexagon to Markov categories*. Politecnico di Torino, 2026. arXiv: **[identifier pending]**.

**Status (10 September 2026).** The manuscript (version v33k, 36 pages) is deposited in `manuscript/` and has been submitted to arXiv; the identifier will be added here as soon as it is announced. This repository is the stable reference for the computational supplement cited in the manuscript and in the CALAMUS research programme (*Categorical Languages for Models and Uncertainty in Statistics*).

## Abstract

McCullagh's 2002 categorical definition of a statistical model, sharpened by Brøns into a hexagonal diagram of categories and functors, identified naturality as the mathematical content of the informal requirement that a parameter's meaning be "retained" across experimental designs. The programme supplied a language but not a calculus: no ambient categorical probability theory existed in which the naturality condition could be computed with, composed, or mechanically verified. We show that the Markov-categorical framework of Fritz, Cho–Jacobs and Perrone supplies exactly this calculus. Our main result recasts a McCullagh–Brøns model as a natural family of morphisms in the Kleisli category of the Giry monad, deterministic on the parameter side; under this translation McCullagh's criterion for a meaningful subparameter becomes a naturality condition, hence a type-level property under the functorial encoding, and Tjur's weaker elementary criterion, invariance under design reduction, is recovered as naturality restricted to the insertion morphisms. We prove the statement in full for finite (discrete) models, state the measurable case with its precise additional hypotheses, and formulate the coherence problem for finite designs in the cohomology of a diagram of algebras in the sense of Gerstenhaber–Schack and Baues–Wirsching, connecting it to recent work of Caputi and the author. We show that the obstruction so obtained is necessarily *class-relative*, the unrestricted class vanishes for every quantity, by construction, and discharge the resulting statement by an exact computation for the one-way layout, which certifies the non-meaningfulness of the marginal dispersion relative to the class of information-bearing recalibrations for scale-confounded designs, recovers the within-group dispersion as a correction by pure linear algebra, and exhibits H¹ of the design scheme as a computable invariant of independent interest.

## Contents

### Manuscript (`manuscript/`)

- `calamus_arxiv_draft_v33k.pdf` — the manuscript as submitted to arXiv (36 pages).
- `calamus_arxiv_draft_v33k.tex`, `calamus_arxiv_draft_v33k.bbl` — the arXiv source pair; builds with `pdflatex` (two passes, no BibTeX).
- `calamus_bridge_draft_clean_v33k.bib` — the bibliography (38 entries, exactly those cited).

### Code (Section 6 of the manuscript)

- `oneway_obstruction.py` — the Section 6 pipeline (Examples 6.6–6.8): defect cochain, cocycle identity, unrestricted-class vanishing, class-relative obstruction under scale confounding, and the within-group dispersion recovered as the correction. Deterministic, exact arithmetic; ≈10 s.
- `oneway_obstruction_output.txt` — committed output; a fresh run reproduces it byte for byte.
- `h1_check.py` — independent exact-over-ℚ certificate: ranks of the Baues–Wirsching cochain complexes for the design schemes of Examples 6.6 and 6.8 across treatment/scale configurations, confirming dim H¹ = 0 in the base cases and dim H¹ = 10 at three treatment levels (the invariant of Example 6.8). ≈3 s.
- `h1_check_output.txt` — committed output of the run above.

## Running

Python ≥ 3.10 with `sympy` and `numpy`:

```
python3 oneway_obstruction.py
python3 h1_check.py
```

Both scripts are self-contained and deterministic. Last reproduction check: 10 September 2026, both outputs identical to the committed files.

## License

MIT for the code in this repository (see `LICENSE`). The manuscript in `manuscript/` is © 2026 Francesco Vaccarino and is **not** covered by the MIT license; it is distributed under the terms selected for the arXiv posting.

## Citing

```bibtex
@misc{vaccarino_bridge_2026,
  author        = {Vaccarino, Francesco},
  title         = {Statistical models as natural transformations: a bridge from the {McCullagh}--{Br{\o}ns} hexagon to {Markov} categories},
  year          = {2026},
  eprint        = {XXXX.XXXXX},
  archivePrefix = {arXiv},
  primaryClass  = {math.ST},
  note          = {Computational supplement: \url{https://github.com/Vaxx66/calamus}}
}
```

Replace `XXXX.XXXXX` with the arXiv identifier once announced; see also `CITATION.cff`.
