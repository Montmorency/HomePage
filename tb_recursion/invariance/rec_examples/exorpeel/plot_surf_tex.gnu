set terminal epslatex size 4.5 in, 4 in color colortext
set output 'exorpeel.tex'

set multiplot layout 2,1 rowsfirst 

set key right top
set format x '%.1f'

set ylabel 'Structural Energy (eV)'
set xlabel 'Energy (eV)'
set xr[-1.5:0.0]
set yr[-0.2:0.05]

pl 'll_7.dat' u ($1*13.605):($2*13.605) w lp pt 7 t '$E_{1}$',\
   'll_7.dat' u ($1*13.605):($3*13.605) w lp pt 8 t '$E_{2}$'

set ylabel '$[E_{2}-E_{1}]$ (meV)'
set xlabel '$E_{F}$'
set xr[-1.5:0.0]
set yr[0.0:15]
pl 'll_7.dat' u ($1*13.605):(1000*13.605*($3-$2)) w lp pt 9 t 'W'
