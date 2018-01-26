FC = gfortran
FCFLAGS = -g -fbounds-check -fdefault-real-8
#FCFLAGS = -O2
#FCFLAGS +=
#LDFLAGS = -l_i_need

%.o: %.f
	$(FC) $(FCFLAGS) -c $<

all: exphnds.x exphnplt.x exrecal.x exproc.x exorpeel.x expeeled.x

exorpeel.x: exorpeel.o slkode.o equiv.o recal.o fccbnd.o out.o selfd.o\
            setup.o skpd.o skpp.o selfp.o mmcal.o nncal.o fcclat.o addat.o peel.o orpeel.o slkode.o skdd.o
	$(FC) -o $@ $^
	cp ./exorpeel.x ../examples/exorpeel/

expeeled.x: expeeled.o ediff.o recrts.o numc.o numd.o
	$(FC) -o $@ $^
	cp ./expeeled.x ../examples/exorpeel/

exphnplt.x: exphnplt.o recsum.o cfgen.o denqd.o recqd.o recrts.o recwt.o numc.o numd.o
	$(FC) -o $@ $^ 
	cp ./exphnplt.x ../examples/phonon/

exphnds.x: exphnds.o recal.o scan1.o scan.o
	$(FC) -o $@ $^ 
	cp ./exphnds.x ../examples/phonon/

exproc.x: exproc.o dencrs.o denqd.o fenval.o recsum.o termgn.o bndest.o bndref.o\
          bndwt.o bndcrd.o cfgen.o cfgpgn.o plyval.o recqd.o recrts.o recwt.o numc.o numd.o\
          taban.o wtmin.o
	$(FC) -o $@ $^ 
	cp ./exproc.x ../examples/exrecal/

exrecal.x: exrecal.o slkode.o equiv.o recal.o fccbnd.o out.o selfd.o setup.o skdd.o mmcal.o nncal.o fcclat.o
	$(FC) -o $@ $^ 
	cp ./exrecal.x ../examples/exrecal/


clean:
	rm -f *.o *.x
