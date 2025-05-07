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

To compile the manuscript, run:

```
pdflatex SerreCyl.tex
bibtex SerreCyl
makeglossaries SerreCyl
pdflatex SerreCyl.tex
pdflatex SerreCyl.tex
```

You may also use `latexmk` for automated compilation:

```
latexmk -pdf SerreCyl.tex
```

## File Structure
- `SerreCyl.tex` — Main LaTeX manuscript
- `biblio.bib` — Bibliography file
- `figs/` — Figures used in the manuscript
- `*.sty` — Style files (if any)
- `.gitignore` — Ignores LaTeX build artifacts except the final PDF

## Code availability

All Maple (Maplesoft™) scripts developed and used in this work, as well as the complete LaTeX source files of this article, are openly available in this repository:

[https://github.com/dutykh/SerreCyl/](https://github.com/dutykh/SerreCyl/)

## License

This repository is intended for academic and research purposes. Please contact the corresponding author for reuse or collaboration inquiries.

## Contact

Corresponding author: Denys Dutykh (denys.dutykh@ku.ac.ae)

---

If you use this work or code, please cite the manuscript and acknowledge the authors.
