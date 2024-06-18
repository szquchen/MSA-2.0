# using ifort and MKL
FC = ifort
######## sequential ########
FLAGS = -r8 -O -qmkl=sequential ## for older versions of ifort compiler, you may need to use -mkl=sequential
######## parallel with openmp ########
#FLAGS = -r8 -O -qopenmp -qmkl=parallel ## for older versions of ifort compiler, you may need to use -openmp -mkl=parallel

# using gfortran and MKL. The example below uses MKL2024; for other MKL versions, please consult https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl-link-line-advisor.html
#FC = gfortran
######## sequential ########
#FLAGS = -fdefault-real-8 -fdefault-double-8 -m64 -I"${MKLROOT}/include" -L${MKLROOT}/lib -Wl,--no-as-needed -lmkl_gf_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl
######## parallel with openmp ########
#FLAGS = -fdefault-real-8 -fdefault-double-8 -m64 -I"${MKLROOT}/include" -L${MKLROOT}/lib -Wl,--no-as-needed -lmkl_gf_lp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl

%.o: %.f90
	$(FC) $(FLAGS) -c $<

OBJFIT := basis.o gradient.o fit.o

OBJGET := basis.o gradient.o pes_shell.o test.o

fit.x : $(OBJFIT) 
	$(FC) $(FLAGS) -o $@ $^

test.x : $(OBJGET)
	$(FC) $(FLAGS) -o $@ $^

clean:
	rm -f *.o *.mod *.x	
