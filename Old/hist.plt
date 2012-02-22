binwidth = 0.037
bin(x,width)=width*floor(x/width) + binwidth/2.
set boxwidth binwidth
unset key
set xlabel "E" 0.0,0.7
set ylabel "Frequency" 2.5,0.0
set title "Density of states of 100x100 atom Graphene sheet with X wrapping enabled. Bin width of 0.037" 0.0,-1.0
set xrange [-3.1:3.1]
set yrange [0:0.5]
plot '100x100wrapx.txt' using (bin($1,binwidth)):(1/370.0) smooth freq with boxes lt -1,'dos-0.dat' u 1:($2) w l 1, 'dos-0.dat' u (-1*$1):($2) w l 1

set terminal png
set output '100x100wrapxhist3.png'
replot
