# Counterexamples to Zoroa et al.'s strong-dominance conjecture for every number of draws above four

For four labels and five independent draws, the strictly positive vectors
`p = (48,48,2,2)/100` and `q = (50,46,3,1)/100` satisfy `p ≺ q`, but the
distinct-label mass differences have signs `-, +, -, +` at counts 1 through 4.
The positive difference at 2 followed by the negative difference at 3 rules
out every real single-crossing threshold. Ordinary stochastic dominance is
retained.
The same obstruction holds for every real `0 < t <= 1/100` for
`p_t=(1/2-2t,1/2-2t,2t,2t)` and `q_t=(1/2,1/2-4t,3t,t)`.

For every integer `n >= 5`, the same two mass signs hold for every real
`0 < t <= 1/(100*n*2^n)`. Adding sufficiently rare labels gives strictly
positive rational counterexamples for every `N >= 4` and `n >= 5`.
The proof gives an explicit choice of the additional probability mass.

## Read and build

- `note.tex` is the complete LaTeX source, with its bibliography included.
- Build with **Tectonic 0.17.0**: `tectonic note.tex`.
- Alternatively run `pdflatex note.tex` twice; the bibliography needs no
  separate BibTeX database or invocation.
- Python 3.9 or later suffices for `python3 verify.py`; there are no Python
  dependencies. It writes `verification.json` alongside the script.

The checker compares weighted enumeration of all 1024 draw words for each
vector, subset inclusion-exclusion, and the 56 multinomial occupancy vectors.
The concrete-witness calculations use integers; there is no floating-point rounding.
For the family it uses exact rational polynomials and interval bounds;
the resulting certificate covers the whole real interval, not a sampled grid.

## Scientific coverage

The written proof gives the inclusion-exclusion calculation and the
single-crossing obstruction and the one-parameter-family corollary.
The small-support lemma gives the three-point-support argument,
using the ordinary stochastic-order theorem of Wong and Yue (1973),
*A majorization theorem for the number of distinct outcomes in n independent
trials*, pp. 391–398, doi:10.1016/0012-365X(73)90070-8, explicitly credited
by Zoroa et al., Remark 6.1, pp. 25–26 of arXiv:1608.00155v1.
The all-draw-count theorem has a written Taylor-remainder proof, and the
all-dimension corollary has a written majorization and coupling proof.
The checker verifies the exact coefficient weights and derivative constants;
its endpoint checks for `n=5,...,20` are finite checks, not the all-n proof.
The Lean declaration `D5.S3.Analytic.Probability.StrongDominanceRefutation.result`
in the pinned source below proves `not claim` from the Theorem 1 vectors
at `N=4`, `n=5` only, and formalizes their decreasing-prefix inequalities
and the four mass values `a_p(2)`, `a_q(2)`, `a_p(3)`, `a_q(3)`.
The masses at counts 1 and 4, the cumulative differences, Theorem 3 and
Corollary 4 have written proofs and script checks only. The script checks
the theorem's coefficients and remainder constants and finite enlarged-vector
examples; the universal statements follow from the written proofs.

https://github.com/the-omega-institute/trureturing/blob/ad59550e53d607f9e9dbd8110a79414c46e97f7e/D5/S3/Analytic/Probability/StrongDominanceRefutation.lean

Posing source: Zoroa, Lesigne, Fernández-Sáez, Zoroa and Casas,
*The coupon collector urn model with unequal probabilities in ecology and
evolution*, arXiv:1608.00155v1, Section 8, Conjecture, page 30;
published in *Journal of the Royal Society Interface* 14 (2017), 20160643,
doi:10.1098/rsif.2016.0643.

## Authors and license

Haobo Ma — ChronoAI Pte Ltd, Singapore.

Wenlin Zhang — National University of Singapore; ChronoAI Pte Ltd, Singapore.

Copyright 2026 Haobo Ma and Wenlin Zhang. Creative Commons Attribution 4.0
International, **CC BY 4.0** (`cc-by-4.0`).
