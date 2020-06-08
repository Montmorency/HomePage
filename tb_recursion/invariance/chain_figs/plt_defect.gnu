set terminal epslatex size 4.0 in, 4.0 in color colortext
set output 'adsorbate_plot.tex'
set multiplot layout 3,1 \
margins 0.2, 0.98, 0.05, 0.95 \
spacing 0.1, 0.1 

set xtics 
set ytics 2

unset xlabel
unset ylabel

set xr [0:2.5]
set yr [-2:2.5]
pl 'ge.dat' u (-$3):1 w l t 'Im[$G_{0}(E)$]'
set xr [0:0.6]
set yr [-4:2.5]
#set yr 
unset xtics
set ylabel 'E' rotate by 0
pl 'g_defect_1.dat' u (-$3):1 w l t 'b/10'
unset ylabel
set xtics
pl 'g_defect_3.dat' u (-$3):1 w l t '2b'


