# Nonlinear Dynamics of Pulsatile Blood Flow in Viscoelastic Vessels: A Dispersive Wave Approach

This repository contains the LaTeX source and supplementary files for the manuscript:

**Nonlinear Dynamics of Pulsatile Blood Flow in Viscoelastic Vessels: A Dispersive Wave Approach**

## Overview

The manuscript develops a systematic framework for modeling pulsatile blood flow in viscoelastic vessels, using asymptotic modeling techniques and a hierarchy of mathematical models:
- Full axisymmetric Euler equations
- Serre-Green-Naghdi (SGN) equations
- Cylindrical Boussinesq equations
- Korteweg–de Vries (KdV) and Benjamin–Bona–Mahony (BBM) equations
- Lumped (Windkessel) models (see Appendix)

The approach provides insight into the interplay between nonlinearity, dispersion, and fluid-structure interactions in large arteries.

## Authors

- Rim El Cheikh (Univ. Grenoble Alpes, Univ. Savoie Mont Blanc, CNRS, LAMA, Chambéry, France)
- Alexey Cheviakov (Department of Mathematics and Statistics, University of Saskatchewan, Saskatoon, Canada)
- Denys Dutykh (Mathematics Department, Khalifa University of Science and Technology, Abu Dhabi, United Arab Emirates)
- Dimitrios Mitsotakis (Victoria University of Wellington, School of Mathematics and Statistics, Wellington, New Zealand)

## How to Compile

The simplest way is the provided `Makefile`:

```
make           # build SerreCyl.pdf, then remove the intermediate files
make build     # build but keep the intermediate files (faster reruns)
make clean     # remove the intermediate files (keeps the PDF)
make distclean # remove the intermediate files and the PDF
```

Equivalently, the manual sequence (pdflatex + bibtex + makeglossaries, then two
further passes to resolve all references) is:

```
pdflatex SerreCyl.tex
bibtex SerreCyl
makeglossaries SerreCyl
pdflatex SerreCyl.tex
pdflatex SerreCyl.tex
```

`latexmk -pdf SerreCyl` also works, with one `makeglossaries SerreCyl` run for the acronym list.

## File Structure
- `SerreCyl.tex` — main LaTeX manuscript
- `biblio02.bib` — bibliography database
- `Makefile` — build automation (see above)
- `wiley-article.cls`, `WileyNJD-*.bst` — journal class and bibliography style
- `figs/` — figures used in the manuscript
- `scripts/` — Python (SymPy / NumPy / Matplotlib) scripts that symbolically verify
  the analytical results and generate the modulational-instability figure
- `maple/` — Maple worksheets for the conservation-law and symmetry computations
- `.gitignore` — ignores LaTeX build artifacts except the final PDF

## Reproducibility and symbolic verification

The analytical results are accompanied by Python scripts in `scripts/` that verify the
underlying algebra symbolically (SymPy) and generate the modulational-instability figure:

- `verify_consistency.py` — Laplace residual `O(δ⁶)`, the crux of the SGN ⇄ Euler consistency theorem
- `verify_dispersion.py` — asymptotic match of the Bessel and KdV dispersion relations; linear-stability signs
- `verify_energy.py` — energy / `L²` balance identities for the KdV and BBM reductions
- `verify_soliton_decay.py` — adiabatic amplitude-decay law of the solitary wave
- `verify_lie_algebra.py` — commutators of the point symmetries (the Galilei algebra)
- `verify_reductions.py` — travelling-wave and self-similar symmetry reductions
- `verify_poisson.py` — Hamiltonian / Poisson structure, Casimir, variational characterisation of solitary waves
- `verify_hodograph.py` — strict hyperbolicity of the dispersionless limit (hodograph / Euler–Poisson–Darboux)
- `mi_growth.py` — modulational-stability diagnostics and the KdV-vs-BBM figure

Run, for example, `python3 scripts/verify_dispersion.py`. Requirements: `sympy`, `numpy`, `matplotlib`.

## Code availability

All Maple (Maplesoft™) worksheets and Python symbolic-verification scripts developed and used in this work, as well as the complete LaTeX source files of this article, are openly available in this repository:

[https://github.com/dutykh/SerreCyl/](https://github.com/dutykh/SerreCyl/)

## License

This repository is licensed under the GNU Lesser General Public License v3.0 (LGPL-3.0).

You are free to use, modify, and distribute this work under the terms of the LGPL-3.0. For more details, see the [LICENSE.md](LICENSE.md) file in this repository or visit the [official license page](https://www.gnu.org/licenses/lgpl-3.0.html).

If you use this work or code, please cite the manuscript and acknowledge the authors.

## Contact

Corresponding author: Denys Dutykh (denys.dutykh@ku.ac.ae)

---

If you use this work or code, please cite the manuscript and acknowledge the authors.
