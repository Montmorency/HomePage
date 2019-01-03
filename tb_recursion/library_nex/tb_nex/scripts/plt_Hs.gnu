set multiplot layout 3,1 rowsfirst 
set xr [-0.6:0.6]
set ylabel 'N(E)'
set xlabel 'E (eV)'


pl './exbccHrecal/dos.dat' u 1:3 w lp pt 7 t 'Fe-54 H Quad',\
   './exbccrecal/53/dos.dat' u 1:3 w lp pt 7 t 'Fe-54 Quad.',\
   './exbccHSrecal/dos.dat' u 1:3 w lp pt 7 t 'Fe-54 HS Quad.'

pl './exbccHSrecal/dos.dat' u  1:4 w lp pt 7 t 'Fe-54 HS',\
   './exbccHrecal/dos.dat' u  1:4 w lp pt 7 t 'Fe-54 H',\
   './exbccrecal/dos.dat' u 1:4 w lp pt 7 t 'Fe-54'

pl './exbccHSrecal/dos.dat' u  1:5 w lp pt 7 t 'Fe-54 HS',\
   './exbccHrecal/dos.dat' u  1:5 w lp pt 7 t 'Fe-54 H',\
   './exbccrecal/dos.dat' u 1:5 w lp pt 7 t 'Fe-54'

#pl './exbccHrecal/dos.dat' u  1:2 w lp pt 7 t 'Fe-54 H Term',\
#   './exbccrecal/53/dos.dat' u 1:2 w lp pt 7 t 'Fe-54 Term.',\
#   './exbccHSrecal/dos.dat' u 1:2 w lp pt 7 t 'Fe-54 Term'

