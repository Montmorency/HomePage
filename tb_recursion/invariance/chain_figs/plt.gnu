unset xtics 
unset ytics

set xlabel 'E'
set ylabel 'G(E)'
set xr [-2:2]
set yr [-2.5:2.5]

#for b=0.4
set  label "b" at 0.8, 0.0
set  label "-b" at -0.8, 0.0

set  label "1/b" at 0.0, 2.2
set  label "-1/b" at 0.0, -2.2

pl 'ge.dat' u 1:2 w l t 'Re[G(E)]', 'ge.dat' u 1:3 w l t 'Im[G(E)]'
