set terminal epslatex size 4.0 in, 3.5 in color colortext
set output 'adsorbate_plot.tex'

unset xtics 
unset ytics

set xlabel 'E'
set ylabel 'Im[G(E)]'

set yr [-3:1.5]

pl 'ge.dat' u 1:3 w l t '$G_{0}(E)$',\
'g_defect_1.dat' u 1:3 w l t 'b/10',\
'g_defect_3.dat' u 1:3 w l t '2b'
