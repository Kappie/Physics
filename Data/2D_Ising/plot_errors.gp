set terminal pngcairo dashed enhanced size 640,640
set output '2d_ising_environment_initialization.png'

#titles = "chi8tolerance1e-8random chi8tolerance1e-8adjusted chi8tolerance1e-8adjustedreverse"

set logscale y
set format y "%.2t*10^{%+03T}"
set xlabel "β"
set ylabel "|m - m_{exact}|"

set title "Error of m with three methods of initializing environment tensors for χ = 8"
set style arrow 1 lc 2
set arrow from 0.44069,graph(0,0) to 0.44069,graph(1,1) nohead dashtype 2

plot "chi8tolerance1e-8random.dat" using 1:2 title "random" ps 2 pt 4, \
  "chi8tolerance1e-8adjusted.dat" using 1:2 title "adjusted" ps 2 pt 6, \
  "chi8tolerance1e-8adjustedreverse.dat" using 1:2 title "adjusted reverse" ps 2 pt 8
