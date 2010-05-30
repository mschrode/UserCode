# Makefile for note on resolution fit

SOURCES = ./figures/*.pdf \
	./*.tex \
	./*.cls
MAINFILE = Note_ResFit_Main.tex


# Default: Compile texfiles
Note_ResolutionFit.pdf: $(SOURCES)

	pdflatex $(MAINFILE)	

	@while ( grep "Rerun to get cross-references" $(MAINFILE:.tex=.log) > /dev/null ); do \
		pdflatex $(MAINFILE); \
	done

	@if ( grep "LaTeX Warning: Reference" $(MAINFILE:.tex=.log) > /dev/null ); then \
		echo ;\
		echo ;\
		echo '>>>>>> Undefined references:'; \
		grep "LaTeX Warning: Reference" $(MAINFILE:.tex=.log); \
		echo ;\
		echo ;\
	fi

	@if ( grep "LaTeX Warning: Citation" $(MAINFILE:.tex=.log) > /dev/null ); then \
		echo ;\
		echo ;\
		echo '>>>>>> Undefined citations'; \
		grep "LaTeX Warning: Citation" $(MAINFILE:.tex=.log); \
		echo ;\
		echo ;\
	fi

	mv $(MAINFILE:.tex=.pdf) Note_ResolutionFit.pdf

	rm -f ./*.aux ./*.toc ./*.log ./*.out ./*.nav ./*.snm


# Clean up
clean:
	rm -f ./*.aux ./*.toc ./*.log ./*.out ./*.nav ./*.snm ./*~
