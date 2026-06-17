# Makefile — compile the SerreCyl manuscript.
#
# Requires a TeX distribution providing: pdflatex, bibtex, makeglossaries.
# (The optional `latexmk` target additionally needs latexmk.)
#
# Usage:
#   make            # build SerreCyl.pdf, then remove intermediate files
#   make build      # build SerreCyl.pdf but KEEP intermediates (fast reruns)
#   make latexmk    # alternative build driven by latexmk
#   make view       # open the resulting PDF
#   make clean      # remove intermediate files (keeps the PDF)
#   make distclean  # remove intermediate files AND the PDF

NAME       = SerreCyl
LATEX      = pdflatex
LATEXFLAGS = -interaction=nonstopmode -halt-on-error -file-line-error
BIBTEX     = bibtex
GLOSSARIES = makeglossaries

# Intermediate build artefacts (all are .gitignored). Listed once and
# reused by the default target, `clean`, and `distclean`.
AUXEXTS  = aux log out toc bbl blg acn acr alg glo gls glg ist \
           fls fdb_latexmk synctex.gz run.xml bcf lof lot
AUXFILES = $(addprefix $(NAME).,$(AUXEXTS))

.PHONY: all build latexmk view clean distclean

# Default target: build the PDF and then tidy up the compilation files, so
# that the repository is left containing only source files and the final
# PDF. The cleanup runs as a recipe step (after the `build` prerequisite),
# so it is never reordered before the build, even under `make -j`.
all: build
	$(RM) $(AUXFILES)

build: $(NAME).pdf

# Canonical full build: an initial pass to emit the .aux/.acn files, then
# the bibliography (bibtex) and the acronym glossary (makeglossaries),
# followed by TWO further LaTeX passes that resolve all cross-references
# and citations (the second pass settles forward references whose numbers
# shift after the bibliography and glossary are inserted).
$(NAME).pdf: $(NAME).tex biblio02.bib
	$(LATEX) $(LATEXFLAGS) $(NAME)
	$(BIBTEX) $(NAME)
	$(GLOSSARIES) $(NAME)
	$(LATEX) $(LATEXFLAGS) $(NAME)
	$(LATEX) $(LATEXFLAGS) $(NAME)

# Convenience target: let latexmk drive the LaTeX/BibTeX passes, with a
# single makeglossaries run in between to build the acronym list.
latexmk:
	latexmk -pdf -bibtex $(NAME)
	$(GLOSSARIES) $(NAME)
	latexmk -pdf -bibtex $(NAME)

view: $(NAME).pdf
	@( xdg-open $(NAME).pdf || open $(NAME).pdf ) >/dev/null 2>&1 &

# Remove intermediate build files; the tracked PDF is preserved.
clean:
	$(RM) $(AUXFILES)

# Remove everything, including the generated PDF.
distclean: clean
	$(RM) $(NAME).pdf
