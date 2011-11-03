set title "50x50 atom Graphene sheet with X wrapping enabled. Line plot to see how fast Eigenvalues increase."
set xlabel "n"
set ylabel "E"
plot "50x50wrapxpositivesort.txt" using 1:2
set terminal png
set output '50x50wrapxlineplot.png'
replot
