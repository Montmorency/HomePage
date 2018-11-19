set terminal epslatex size 4.5 in, 4 in color colortext
set output 'recal_dos.tex'

set multiplot layout 2,1 rowsfirst 
set ylabel 'N(E)'

set key left top
set format x '%.1f'

pl 'dos.dat' u ($1*13.605):2 w lp pt 7 t 'TermDOS',\
   'dos.dat' u ($1*13.605):3 w lp pt 7 t 'QuadDOS',\

set ylabel 'Structural Energy (eV)'
set xlabel 'Energy (eV)'

set key left
set yr [-2:6]
pl 'dos.dat' u ($1*13.605):($5*13.605) w lp pt 7 t 'Structural Energy',\
   'dos.dat' u ($1*13.605):4 w lp pt 7 t 'Integrated DOS'


