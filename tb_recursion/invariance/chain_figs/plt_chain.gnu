set terminal epslatex size 5.0 in, 3.25 in color colortext
set output 'chain_plot.tex'
set multiplot layout 2,1 \
margins 0.2, 0.98, 0.05, 0.95 \
spacing 0.08, 0.025 

unset xtics 
unset ytics

set xr [-2.5:2.5]
set yr [-1.5:1.5]

#for b=0.4
#set label "1/b"  at 0.8, 2.4
#set label "-1/b" at 0.8, -2.4

set ylabel 'E' offset -2,1 rotate by 0

set label 1 "$b$" at  -0.15, 0.9
set label 2 "$-b$" at -0.3, -0.80

pl 'ge.dat' u (-$3):1 w l lt -1 lc rgb "blue" t 'Im[$G(E)$]'

unset label 1
unset label 2

set yr [-2.0:2.0]

pl 'ge.dat' u 2:1 w l lt -1 lc rgb "red" t 'Re[$G(E)$]' 
