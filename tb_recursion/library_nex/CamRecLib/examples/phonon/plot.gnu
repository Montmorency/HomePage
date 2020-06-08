set xlabel 'E (eV)'
set ylabel 'N(E)'

pl 'fort.10' u($1*13.605):2 w lp pt 7 t 'Phonon LDOS',\
   'fort.10' u ($1*13.605):3 w lp pt 7 t 'INT LDOS'
