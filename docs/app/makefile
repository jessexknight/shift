all:
	make tex && make bib && make tex && make tex
tex:
	pdflatex -interaction batchmode main
bib:
	biber --onlylog main
cls:
	cp ~/dev/tex/pkg/ms/ms/ms.cls .
clean:
	texclean .
