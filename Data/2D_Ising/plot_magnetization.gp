set terminal pngcairo dashed enhanced size 640,640
set output 'output.png'

set xlabel "T"
set ylabel "m"

set title "Magnetization of 2D Ising model"
set arrow from 2.2692,graph(0,0) to 2.2692,graph(1,1) nohead dashtype 2

plot "chi4_magnetization.dat" using 1:2 title "χ = 4" ps 2 pt 4, \
  "chi8_magnetization.dat" using 1:2 title "χ = 8" ps 2 pt 6, \
  "chi12_magnetization.dat" using 1:2 title "χ = 12" ps 2 pt 8, \
  "chi16_magnetization.dat" using 1:2 title "χ = 16" ps 2 pt 10
