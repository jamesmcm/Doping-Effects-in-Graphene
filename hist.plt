binwidth = 0.11
bin(x,width)=width*floor(x/width) + binwidth/2.
set boxwidth binwidth

unset key
set title "50x50 atom Graphene sheet with X wrapping enabled. Bin width of 0.11"
plot '50x50wrapx.txt' using (bin($1,binwidth)):(1.0) smooth freq with boxes lt -1
set terminal png
set output '50x50wrapxhist.png'
replot
