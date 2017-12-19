set multiplot layout 2,1 rowsfirst 
set ylabel 'N(E)'

pl 'dos.gnu' u ($1*13.605):2 w lp pt 7 t 'TermDOS',\
   'dos.gnu' u ($1*13.605):3 w lp pt 7 t 'QuadDOS',\

set ylabel 'Structural Energy (eV)'
set xlabel 'Energy (eV)'
pl 'dos.gnu' u ($1*13.605):($5*13.605) w lp pt 7 t 'Structural Energy',\
   'dos.gnu' u ($1*13.605):4 w lp pt 7 t 'Integrated Dos'


