set yr [0:1]
set xr[-5:5]
set pointsize -1
set key off
plot "anan10test.dat" using 1:2 lt 2 with lines, "anan10test.dat" using 1:3 lt 3 with lines
#plot "chriscompare.dat" using 1:2 lt 4 with dots, "chriscompare.dat" using 1:3 lt 5 with dots