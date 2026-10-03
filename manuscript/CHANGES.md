# Changes to the manuscript

## v3 (source `calamus_arkiv_draft_light_v04`, release `v36`) with respect to v2 (source `calamus_arkiv_draft_light_v03`, release `v35`)

Attribution of the diagram Markov category to Fritz. No numbered statement is added or removed; the numbering of theorems, propositions, definitions and equations is unchanged, the bibliography is unchanged (the `.bbl` file is identical), and the manuscript is still 39 pages.

1. **Section 1, related work.** The sentence on Fritz's sufficiency theorems is replaced. It now records that Fritz (Section 7, Proposition 7.1) shows that the functors from a small category D to the deterministic morphisms of a Markov category C, with all natural transformations as morphisms and the pointwise monoidal structure, form a Markov category Fun(D, C), applied there to stochastic processes, and that whether Fun(D, C) has conditionals when C does is his Problem 11.9; that Theorem 1.1 identifies a statistical model with a morphism of Fun(V × D^op, Stoch); and that what is specific to statistics is the indexing by the designs with the parameter functor as domain, Fritz's Remark 7.3 noting that processes require structure beyond the diagram category.
2. **Remark 4.6.** A paragraph is added: L and R are objects of Fun(V × D^op, Stoch) because L(m) and R(m) are deterministic, so Theorem 4.4 says that a measurable model is a morphism of that category, and a Q-coherent family of Definition 5.1 between functors with deterministic images is a morphism of Fun(Q^op, Stoch); the reading is descriptive, and the existence of conditionals in Fun(D, C) is Fritz's open Problem 11.9.
3. **Remark 8.6.** A sentence is added after the conjecture on the inversions: already for proper priors, whether the inversions of Definition 8.1 form a natural transformation is whether the diagram category of Remark 4.6 has Bayesian inversions, a form of Problem 11.9.
4. **Footnote of Section 6.2.** The release containing the scripts and this version of the manuscript is `v36`.
5. **Date.** October 2026.

## v2 (source `calamus_arkiv_draft_light_v03`, release `v35`) with respect to v1 (source `calamus_arkiv_draft_light_v02`, release `v34`)

All changes are in Sections 7–9 and in one footnote of Section 6.2. The numbering of the theorems, propositions, definitions and equations is unchanged; the bibliography is unchanged (the `.bbl` file is identical). The manuscript is now 39 pages instead of 38.

1. **Section 7.1.** The standard Borel spaces are described as the measurable spaces isomorphic to a Polish space with its Borel σ-algebra, equivalently to a Borel subset of ℝ, which is the wording of the cited sources.
2. **Definition 8.1.** The notion of a prior *coherent at* a single morphism ψ is now defined; it was used in Theorem 1.7 and in Section 9 without a definition.
3. **Proposition 8.2(3).** The naturality of the prior predictives is stated for the family indexed by the designs over the fixed unit object, as in Section 9, since the predictive depends on the design map and not only on the covariate space.
4. **Proposition 8.2(4).** The injectivity hypothesis may be imposed on any class of states which contains the prior and is stable under the parameter maps; the previous statement is the case of all the states. The clause on `FinStoch` is now proved: the kernel of the mixing map lies in the hyperplane of the vectors with zero sum, which is spanned by the differences of probability vectors.
5. **Remark 8.3.** The invariance under the automorphisms of the covariate space is the finite counterpart of the exchangeability of Fritz, Gonda and Perrone (Definition 4.1), which is stated there for a state on a countable Kolmogorov power.
6. **Proposition 8.4, proof.** The uniform prior on M²×S is not concentrated on the diagonal, so it is not a relative interior point of the polytope over `O_full`; that case is now handled through the affine embedding of the simplex on M×S. The statement and the four dimensions are unchanged.
7. **Remark 8.6.** "Not of the form π(dζ)f(z|ζ)" reads "not proportional to π(dζ)f(z|ζ)", as in Dawid, Stone and Zidek, Section 1.
8. **Definition 9.1.** The constant τ² does not depend on ω. The existence of compatible ridge data is asserted on the full subcategory of the subspaces of the coordinate spaces ℝ^p, which is the category of covariate spaces used in Section 9.
9. **Proposition 9.2(3).** Hoerl and Kennard print the prior covariance as (δ²/k)I with δ² undefined; the text now says so and reads δ² as the error variance σ².
10. **Proposition 9.4(4).** The identifiability hypothesis of Proposition 8.2(4) is imposed on the states with variance fixed at σ₀², which contain the ridge priors. On all the states of ω*×ℝ_{>0} it fails when n = p, as shown by an explicit pair of states with the same prior predictive; the proof for the fixed-variance states is given through the characteristic function of the location mixture.
11. **Proof of Proposition 9.4(3).** The annihilator of ker ψ is written (ker ψ)⁰ instead of (ker ψ)^⊥, since inner products are in play.
12. **Footnote of Section 6.2.** The release containing the scripts and this version of the manuscript is `v35`.
13. **Definition 8.1.** "Its second marginal" and "given its second output" read "its marginal on R(ω)" and "given R(ω)", the vocabulary of Section 7.2, which names marginals and conditionals by factor.
