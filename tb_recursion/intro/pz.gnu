set terminal epslatex size 4.0 in, 2.75 in color colortext
set output 'pz.tex'

set multiplot layout 2,1 

set ylabel 'Energy (Ry.)'

set xr [0:1]
set yr [-0.1:0.0]
pl 'pz_vxc.dat' u 1:2 w lp pt 7 lc rgb 'red' t 'PZ Correlation Energy'

set xlabel '$\rho$ (e/a.u.$^{3}$)'
set yr [-1.0:0.0]
pl 'pz_vxc.dat' u 1:3 w lp pt 7 lc rgb 'blue' t 'PZ Exchange Energy'


