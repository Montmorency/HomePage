set terminal epslatex size 3.25 in, 3.25 in color colortext
set output 'chain_plot.tex'
set multiplot layout 2,1 

unset xtics 
unset ytics

unset key
set xr [-2.5:2.5]
set yr [-1.5:1.5]

#for b=0.4
#set label "1/b"  at 0.8, 2.4
#set label "-1/b" at 0.8, -2.4

set label 1 "$b$" at  -0.15, 0.9
set label 2 "$-b$" at -0.15, -0.80

set ylabel 'Im[$G(E)$]'
pl 'ge.dat' u (-$3):1 w l lt -1 lc rgb "blue"

unset label 1
unset label 2

set yr [-2.0:2.0]
set xlabel 'E'
set ylabel 'Re[$G(E)$]' 
pl 'ge.dat' u 2:1 w l lt -1 lc rgb "red"
