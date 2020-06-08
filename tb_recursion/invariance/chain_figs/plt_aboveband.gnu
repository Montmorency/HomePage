unset xtics 
unset ytics

set xlabel 'E'
set ylabel 'G(E)'

set yr [-5:0]

pl 'ge.dat' u 1:3 w l t 'Chain Im[G(E)]'
repl 'g_aboveband_1.dat' u 1:3 w l t 'b/10'
repl 'g_aboveband_2.dat' u 1:3 w l t 'b'
repl 'g_aboveband_3.dat' u 1:3 w l t '2b'
