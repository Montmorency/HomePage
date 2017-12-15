FC = gfortran
FCFLAGS = -g -fbounds-check
#FCFLAGS = -O2
#FCFLAGS +=
#LDFLAGS = -l_i_need

programs = exrecal

all: $(programs) 

exrecal: slkode.o, fccbnd.o, equiv.o, fcclat.o, out.o, recal.o, setup.o

%: %.o
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.f
	$(FC) $(FCFLAGS) -c $<

clean:
	rm -f *.o
