set multiplot layout 3,1 rowsfirst 
set xr [-0.6:0.6]
set ylabel 'N(E)'
set xlabel 'E (eV)'


pl './exbccrecal/53/dos.dat' u 1:2 w lp pt 7 t 'BCC 54 Term.',\
   './exbccrecal/53/dos.dat' u 1:3 w lp pt 7 t 'BCC 54 Quad.',\
   './exbccrecal/H_DFT/H2_dos.dat' u (($1-3.1)/13.605):(0.5*($2+$3)) w lp pt 7  t 'BCC Fe-H DFT'

pl './exbccHrecal/dos.dat' u 1:3 w lp pt 7 t 'Fe-54 H Quad',\
   './exbccrecal/53/dos.dat' u 1:3 w lp pt 7 t 'Fe-54 Quad.'

pl './exbccHrecal/dos.dat' u  1:2 w lp pt 7 t 'Fe-54 H Term',\
   './exbccrecal/53/dos.dat' u 1:2 w lp pt 7 t 'Fe-54 Term.'

#pl './exbccHrecal/dos.dat' u  1:4 w lp pt 7 t 'Fe-54 H Term',\
#   './exbccrecal/53/dos.dat' u 1:4 w lp pt 7 t 'Fe-54 Term.'
