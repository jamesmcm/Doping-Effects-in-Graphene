binwidth = 0.07
bin(x,width)=width*floor(x/width) + binwidth/2.
set boxwidth binwidth
unset key
set xlabel "E bins"
set ylabel "Frequency"
set title "50x50 atom Graphene sheet with X wrapping disabled. Bin width of 0.07"
plot '50x50nowrap.txt' using (bin($1,binwidth)):(1.0) smooth freq with boxes lt -1
set terminal png
set output '50x50nowraphist.png'
replot
