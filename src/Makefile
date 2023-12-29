# using ifort and MKL
FC = ifort -r8 -O
LIBS = -qmkl=sequential

# using gfortran and MKL
#FC = gfortran -fdefault-real-8 -O
#LIBS = -L${MKLROOT}/lib/intel64 -lmkl_gf_lp64 -lmkl_sequential -lmkl_core

%.o: %.f90
	$(FC) -c $<

OBJFIT := basis.o gradient.o fit.o

OBJGET := basis.o gradient.o pes_shell.o test.o

fit.x : $(OBJFIT) 
	$(FC) -o $@ $^ $(LIBS)

test.x : $(OBJGET)
	$(FC) -o $@ $^ $(LIBS)

clean:
	rm -f *.o *.mod *.x	
