# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line, and also
# from the environment for the first two.
SPHINXOPTS    ?=
SPHINXBUILD  ?= python -m sphinx
SOURCEDIR    = source
BUILDDIR     = build

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile

# Standard Sphinx targets
html:
	@$(SPHINXBUILD) -b html "$(SOURCEDIR)" "$(BUILDDIR)/html" $(SPHINXOPTS) $(O)
	@echo
	@echo "Build finished. The HTML pages are in $(BUILDDIR)/html."

dirhtml:
	@$(SPHINXBUILD) -b dirhtml "$(SOURCEDIR)" "$(BUILDDIR)/dirhtml" $(SPHINXOPTS) $(O)
	@echo
	@echo "Build finished. The HTML pages are in $(BUILDDIR)/dirhtml."

singlehtml:
	@$(SPHINXBUILD) -b singlehtml "$(SOURCEDIR)" "$(BUILDDIR)/singlehtml" $(SPHINXOPTS) $(O)
	@echo
	@echo "Build finished. The HTML page is in $(BUILDDIR)/singlehtml."

pickle:
	@$(SPHINXBUILD) -b pickle "$(SOURCEDIR)" "$(BUILDDIR)/pickle" $(SPHINXOPTS) $(O)
	@echo
	@echo "Build finished; now you can process the pickle files."

json:
	@$(SPHINXBUILD) -b json "$(SOURCEDIR)" "$(BUILDDIR)/json" $(SPHINXOPTS) $(O)
	@echo
	@echo "Build finished; now you can process the JSON files."

htmlhelp:
	@$(SPHINXBUILD) -b htmlhelp "$(SOURCEDIR)" "$(BUILDDIR)/htmlhelp" $(SPHINXOPTS) $(O)
	@echo
	@echo "Build finished; now you can run HTML Help Workshop with the" \
	      ".hhp project file in $(BUILDDIR)/htmlhelp."

qthelp:
	@$(SPHINXBUILD) -b qthelp "$(SOURCEDIR)" "$(BUILDDIR)/qthelp" $(SPHINXOPTS) $(O)
	@echo
	@echo "Build finished; now you can run \"qcollectiongenerator\" with the" \
	      ".qhcp project file in $(BUILDDIR)/qthelp, like this:"
	@echo "# qcollectiongenerator $(BUILDDIR)/qthelp/DeepChemServer.qhcp"
	@echo "To view the help file:"
	@echo "# assistant $(BUILDDIR)/qthelp/DeepChemServer.qhc"

devhelp:
	@$(SPHINXBUILD) -b devhelp "$(SOURCEDIR)" "$(BUILDDIR)/devhelp" $(SPHINXOPTS) $(O)
	@echo
	@echo "Build finished."
	@echo "To view the help file:"
	@echo "# mkdir -p $$HOME/.local/share/devhelp/DeepChemServer"
	@echo "# ln -s $(BUILDDIR)/devhelp $$HOME/.local/share/devhelp/DeepChemServer"
	@echo "# devhelp"

epub:
	@$(SPHINXBUILD) -b epub "$(SOURCEDIR)" "$(BUILDDIR)/epub" $(SPHINXOPTS) $(O)
	@echo
	@echo "Build finished. The epub file is in $(BUILDDIR)/epub."

epub3:
	@$(SPHINXBUILD) -b epub3 "$(SOURCEDIR)" "$(BUILDDIR)/epub3" $(SPHINXOPTS) $(O)
	@echo
	@echo "Build finished. The epub3 file is in $(BUILDDIR)/epub3."

latex:
	@$(SPHINXBUILD) -b latex "$(SOURCEDIR)" "$(BUILDDIR)/latex" $(SPHINXOPTS) $(O)
	@echo
	@echo "Build finished; the LaTeX files are in $(BUILDDIR)/latex."
	@echo "Run \`make' in that directory to run these through (pdf)latex" \
	      "(use \`make latexpdf' here to do that automatically)."

latexpdf:
	@$(SPHINXBUILD) -b latex "$(SOURCEDIR)" "$(BUILDDIR)/latex" $(SPHINXOPTS) $(O)
	@echo "Running LaTeX files through pdflatex..."
	@$(MAKE) -C $(BUILDDIR)/latex all-pdf
	@echo "pdflatex finished; the PDF files are in $(BUILDDIR)/latex."

latexpdfja:
	@$(SPHINXBUILD) -b latex "$(SOURCEDIR)" "$(BUILDDIR)/latex" $(SPHINXOPTS) $(O)
	@echo "Running LaTeX files through platex and dvipdfmx..."
	@$(MAKE) -C $(BUILDDIR)/latex all-pdf-ja
	@echo "pdflatex finished; the PDF files are in $(BUILDDIR)/latex."

text:
	@$(SPHINXBUILD) -b text "$(SOURCEDIR)" "$(BUILDDIR)/text" $(SPHINXOPTS) $(O)
	@echo
	@echo "Build finished. The text files are in $(BUILDDIR)/text."

man:
	@$(SPHINXBUILD) -b man "$(SOURCEDIR)" "$(BUILDDIR)/man" $(SPHINXOPTS) $(O)
	@echo
	@echo "Build finished. The manual pages are in $(BUILDDIR)/man."

texinfo:
	@$(SPHINXBUILD) -b texinfo "$(SOURCEDIR)" "$(BUILDDIR)/texinfo" $(SPHINXOPTS) $(O)
	@echo
	@echo "Build finished. The Texinfo files are in $(BUILDDIR)/texinfo."
	@echo "Run \`make' in that directory to run these through makeinfo" \
	      "(use \`make info' here to do that automatically)."

info:
	@$(SPHINXBUILD) -b texinfo "$(SOURCEDIR)" "$(BUILDDIR)/texinfo" $(SPHINXOPTS) $(O)
	@echo "Running Texinfo files through makeinfo..."
	@make -C $(BUILDDIR)/texinfo info
	@echo "makeinfo finished; the Info files are in $(BUILDDIR)/texinfo."

gettext:
	@$(SPHINXBUILD) -b gettext "$(SOURCEDIR)" "$(BUILDDIR)/locale" $(SPHINXOPTS) $(O)
	@echo
	@echo "Build finished. The message catalogs are in $(BUILDDIR)/locale."

changes:
	@$(SPHINXBUILD) -b changes "$(SOURCEDIR)" "$(BUILDDIR)/changes" $(SPHINXOPTS) $(O)
	@echo
	@echo "The overview file is in $(BUILDDIR)/changes."

linkcheck:
	@$(SPHINXBUILD) -b linkcheck "$(SOURCEDIR)" "$(BUILDDIR)/linkcheck" $(SPHINXOPTS) $(O)
	@echo
	@echo "Link check complete; look for any errors in the above output " \
	      "or in $(BUILDDIR)/linkcheck/output.txt."

doctest:
	@$(SPHINXBUILD) -b doctest "$(SOURCEDIR)" "$(BUILDDIR)/doctest" $(SPHINXOPTS) $(O)
	@echo "Testing of doctests in the sources finished, look at the " \
	      "results in $(BUILDDIR)/doctest/output.txt."

coverage:
	@$(SPHINXBUILD) -b coverage "$(SOURCEDIR)" "$(BUILDDIR)/coverage" $(SPHINXOPTS) $(O)
	@echo "Testing of coverage in the sources finished, look at the " \
	      "results in $(BUILDDIR)/coverage/python.txt."

xml:
	@$(SPHINXBUILD) -b xml "$(SOURCEDIR)" "$(BUILDDIR)/xml" $(SPHINXOPTS) $(O)
	@echo
	@echo "Build finished. The XML files are in $(BUILDDIR)/xml."

pseudoxml:
	@$(SPHINXBUILD) -b pseudoxml "$(SOURCEDIR)" "$(BUILDDIR)/pseudoxml" $(SPHINXOPTS) $(O)
	@echo
	@echo "Build finished. The pseudo-XML files are in $(BUILDDIR)/pseudoxml."

# Custom targets for DeepChem Server documentation
clean:
	rm -rf $(BUILDDIR)/*
	@echo "Removed everything under '$(BUILDDIR)'."

clean-cache:
	rm -rf $(BUILDDIR)/.doctrees
	@echo "Removed Sphinx cache files."

# Build with warnings as errors (strict mode)
strict:
	@$(SPHINXBUILD) -W -b html "$(SOURCEDIR)" "$(BUILDDIR)/html" $(SPHINXOPTS) $(O)
	@echo
	@echo "Build finished with warnings as errors. The HTML pages are in $(BUILDDIR)/html."

# Fast build (skip expensive operations)
fast:
	@$(SPHINXBUILD) -E -b html "$(SOURCEDIR)" "$(BUILDDIR)/html" $(SPHINXOPTS) $(O)
	@echo
	@echo "Fast build finished. The HTML pages are in $(BUILDDIR)/html."

# Auto-build and serve locally for development
serve: html
	@echo "Starting local server at http://localhost:8000"
	@echo "Press Ctrl+C to stop the server"
	@cd $(BUILDDIR)/html && python -m http.server 8000

# Watch for changes and rebuild automatically (requires sphinx-autobuild)
watch:
	@which sphinx-autobuild > /dev/null || (echo "sphinx-autobuild not found. Install with: pip install sphinx-autobuild" && exit 1)
	@sphinx-autobuild "$(SOURCEDIR)" "$(BUILDDIR)/html" $(SPHINXOPTS) $(O) \
		--host 0.0.0.0 \
		--port 8000 \
		--open-browser \
		--delay 1 \
		--ignore "*.tmp" \
		--ignore "*.swp" \
		--ignore "*~"

# Generate API documentation from source code
apidoc:
	@echo "Generating API documentation from source code..."
	@sphinx-apidoc -o $(SOURCEDIR)/api_reference/generated ../deepchem_server --force --module-first
	@echo "API documentation generated in $(SOURCEDIR)/api_reference/generated/"

# Full rebuild (clean + html)
rebuild: clean html

# Check for broken internal links
check-internal:
	@$(SPHINXBUILD) -b linkcheck "$(SOURCEDIR)" "$(BUILDDIR)/linkcheck" $(SPHINXOPTS) $(O) -D linkcheck_ignore='http.*'

# Generate spelling check report (requires sphinxcontrib-spelling)
spelling:
	@which sphinxcontrib-spelling > /dev/null || (echo "sphinxcontrib-spelling not found. Install with: pip install sphinxcontrib-spelling" && exit 1)
	@$(SPHINXBUILD) -b spelling "$(SOURCEDIR)" "$(BUILDDIR)/spelling" $(SPHINXOPTS) $(O)
	@echo "Spelling check complete. Results in $(BUILDDIR)/spelling/"

# Show documentation statistics
stats:
	@echo "Documentation Statistics:"
	@echo "========================="
	@find $(SOURCEDIR) -name "*.rst" | wc -l | sed 's/^/RST files: /'
	@find $(SOURCEDIR) -name "*.py" | wc -l | sed 's/^/Python files: /'
	@find $(SOURCEDIR) -name "*.rst" -exec wc -l {} + | tail -1 | sed 's/total/Total lines in RST files:/'
	@if [ -d $(BUILDDIR)/html ]; then \
		find $(BUILDDIR)/html -name "*.html" | wc -l | sed 's/^/Generated HTML files: /'; \
	else \
		echo "Generated HTML files: 0 (run 'make html' first)"; \
	fi

# Development setup help
dev-setup:
	@echo "Development Setup for DeepChem Server Documentation"
	@echo "=================================================="
	@echo ""
	@echo "1. Install required packages:"
	@echo "   pip install sphinx sphinx-rtd-theme sphinx-autobuild"
	@echo ""
	@echo "2. Optional packages for enhanced features:"
	@echo "   pip install sphinxcontrib-spelling"
	@echo ""
	@echo "3. Useful commands:"
	@echo "   make html        - Build HTML documentation"
	@echo "   make watch       - Auto-rebuild on changes (requires sphinx-autobuild)"
	@echo "   make serve       - Build and serve on localhost:8000"
	@echo "   make linkcheck   - Check for broken links"
	@echo "   make clean       - Clean build directory"
	@echo "   make strict      - Build with warnings as errors"
	@echo ""

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O) 