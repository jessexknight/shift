cat $1.tex > which.aux
pdflatex main ;\
cp main.pdf ../../out/fig/tikz/$1.pdf
convert -density 600 main.pdf ../../out/fig/tikz/png/$1.png
