FC = gfortran
FCFLAGS = -g -fbounds-check
#FCFLAGS = -O2
#FCFLAGS +=
#LDFLAGS = -l_i_need

%.o: %.f
	$(FC) $(FCFLAGS) -c $<

exrecal.x: exrecal.o slkode.o equiv.o recal.o fccbnd.o out.o selfd.o setup.o skdd.o mmcal.o nncal.o fcclat.o
	$(FC) -o $@ $^ 

clean:
	rm -f *.o
