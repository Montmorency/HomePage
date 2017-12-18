FC = gfortran
FCFLAGS = -g -fbounds-check
#FCFLAGS = -O2
#FCFLAGS +=
#LDFLAGS = -l_i_need

%.o: %.f
	$(FC) $(FCFLAGS) -c $<

all : exrecal.x exproc.x 

exproc.x: exproc.o dencrs.o denqd.o fenval.o recsum.o termgn.o bndest.o bndref.o\
          bndwt.o bndcrd.o cfgen.o cfgpgn.o plyval.o recqd.o recrts.o recwt.o numc.o numd.o\
          taban.o wtmin.o
	$(FC) -o $@ $^ 
	cp ./exproc.x ../examples/exrecal/

exrecal.x: exrecal.o slkode.o equiv.o recal.o fccbnd.o out.o selfd.o setup.o skdd.o mmcal.o nncal.o fcclat.o
	$(FC) -o $@ $^ 
	cp ./exrecal.x ../examples/exrecal/


clean:
	rm -f *.o
