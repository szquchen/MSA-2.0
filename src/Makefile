# using ifort and MKL
FC = ifort
FLAGS = -r8 -O
######## sequential ########
LIBS = -qmkl=sequential ## for older versions of ifort compiler, you may need to use -mkl=sequential
######## parallel with openmp ########
LIBS = -qopenmp -qmkl=parallel ## for older versions of ifort compiler, you may need to use -openmp -mkl=parallel

# using gfortran and MKL. The example below uses MKL2024; for other MKL versions, please consult https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl-link-line-advisor.html
#FC = gfortran
#FLAGS = -fdefault-real-8 -fdefault-double-8 -m64  -I"${MKLROOT}/include"
######## sequential ########
#LIBS = -L${MKLROOT}/lib -Wl,--no-as-needed -lmkl_gf_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl
######## parallel with openmp ########
#LIBS = -L${MKLROOT}/lib -Wl,--no-as-needed -lmkl_gf_lp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl

%.o: %.f90
	$(FC) $(FLAGS) -c $<

OBJFIT := basis.o gradient.o fit.o

OBJGET := basis.o gradient.o pes_shell.o test.o

fit.x : $(OBJFIT) 
	$(FC) $(FLAGS) -o $@ $^ $(LIBS)

test.x : $(OBJGET)
	$(FC) $(FLAGS) -o $@ $^ $(LIBS)

clean:
	rm -f *.o *.mod *.x	
