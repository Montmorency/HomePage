set multiplot layout 3,1 rowsfirst 
set xr [-0.6:0.6]

set ylabel 'N(E)'
pl './exbccHrecal/dos.dat' u  1:3 w lp pt 7 t 'Fe-54 H Term',\
   './exbccrecal/53/dos.dat' u 1:3 w lp pt 7 t 'Fe-54 Term.'

set ylabel 'Energy (Ry.)'
pl './exbccHrecal/dos.dat' u  1:5 w lp pt 7 t 'Fe-54 H Term',\
   './exbccrecal/53/dos.dat' u 1:5 w lp pt 7 t 'Fe-54 Term.'

set arrow from -0.6, 3 to 0.6, 3
set arrow from -0.6, 2 to 0.6, 2
set xlabel 'N'
pl './exbccHrecal/dos.dat' u  1:4 w lp pt 7 t 'Fe-54 H Term',\
   './exbccrecal/53/dos.dat' u 1:4 w lp pt 7 t 'Fe-54 Term.'

