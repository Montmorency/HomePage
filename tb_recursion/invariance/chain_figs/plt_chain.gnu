set terminal epslatex size 3.25 in, 3.25 in color colortext
set output 'chain_plot.tex'

unset xtics 
unset ytics

set key at 0.1, 2.2
set xlabel 'E'
set ylabel 'G(E)'
set xr [-2:2]
set yr [-2.5:2.5]

#for b=0.4
set  label "b" at 0.8, -2.6
set  label "-b" at -0.8, -2.6

set label "1/b"  at -2.5, 2.4
set label "-1/b" at -2.5, -2.4

pl 'ge.dat' u 1:2 w l t 'Re[G(E)]', 'ge.dat' u 1:3 w l t 'Im[G(E)]'
