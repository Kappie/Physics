set terminal pngcairo dashed enhanced size 640,640
set output 'output.png'

set xlabel "T"
set ylabel "m"

set title "Magnetization of 2D Ising model"
set arrow from 2.2692,graph(0,0) to 2.2692,graph(1,1) nohead dashtype 2

plot "chi4_adjustedreverse.dat" using 1:2
