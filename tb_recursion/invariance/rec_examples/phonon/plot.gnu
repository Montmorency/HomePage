set terminal epslatex size 4.5 in, 4 in color colortext
set output 'phonon_dos.tex'

set key left
set xr [0:21]
set yr [0:13]

set xlabel 'E (meV)'
set ylabel 'N(E)'

pl 'fort.10' u($1*13.605):2 w lp pt 7 t 'Phonon LDOS',\
   'fort.10' u ($1*13.605):3 w lp pt 7 t 'Integrated LDOS'
