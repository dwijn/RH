# Build unit-tests.pdf: make -C docs/unit-tests
unit-tests.pdf: unit-tests.tex
	pdflatex -interaction=nonstopmode -halt-on-error unit-tests.tex >/dev/null
	pdflatex -interaction=nonstopmode -halt-on-error unit-tests.tex >/dev/null

clean:
	rm -f unit-tests.aux unit-tests.log unit-tests.out unit-tests.toc unit-tests.pdf
