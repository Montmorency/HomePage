unset xtics 
unset ytics

set xlabel 'E'
set ylabel 'G(E)'

set yr [-5:5]

pl 'ge.dat' u 1:2 w l t 'Chain Im[G(E)]'
repl 'g_inband_1.dat' u 1:2 w l t 'b/10'
repl 'g_inband_2.dat' u 1:2 w l t 'b'
repl 'g_inband_3.dat' u 1:2 w l t '2b'
