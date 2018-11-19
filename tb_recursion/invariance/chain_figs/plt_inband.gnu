set terminal epslatex size 4.0 in, 3.5 in color colortext
set output 'inband_plot.tex'

set xlabel 'E'
set ylabel 'Im G(E)'

set yr [-5:2]

pl 'ge.dat' u 1:3 w l t 'Chain Im[G(E)]',\
  'g_inband_1.dat' u 1:3 w l t 'b/10',\
  'g_inband_2.dat' u 1:3 w l t 'b',\
  'g_inband_3.dat' u 1:3 w l t '2b'
