set terminal epslatex size 4.0 in, 5.0 in color colortext
set output 'inband_plot.tex'
set multiplot layout 4,1 title 'Im $[G(E)]$'\
margins 0.2, 0.98, 0.05, 0.95 \
spacing 0.08, 0.025 

set xr [0:4]
set yr [-1.5:1.5]
unset xtics 
set key top right
set ytics 1

pl 'ge.dat' u (-1*$3):1 w l t 'Chain Im[G(E)]'

unset title

pl 'g_inband_1.dat' u (-1*$3):1 w l t 'b/10'

set ylabel 'E'
pl 'g_inband_2.dat' u (-1*$3):1 w l t 'b'
unset ylabel

set key center right
set xtics 1
pl 'g_inband_3.dat' u (-1*$3):1 w l t '2b'
