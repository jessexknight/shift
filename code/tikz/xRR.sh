# compiles all 9 *RR figures (3 types x 3 shapes)
for xRR in tRR dRR nRR; do
  for shape in step ramp exp; do
    echo "\let\uRRshape\uRR$shape" | cat - ${xRR}.tex > which.aux
    pdflatex main ;
    cp main.pdf ../../out/fig/tikz/${xRR}.${shape}.pdf
  done
done
