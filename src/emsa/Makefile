CC := g++
CFLAG := -O2 -Wall
HEADER_PATH = ./header
INCLUDE = -I ${HEADER_PATH}

LOCAL_LIB_PATH = ./lib
LLIB = -L${LOCAL_LIB_PATH} -lcppmol -lcppmat
GLIB = -lm -llapack

all:msa 
msa: msa.o monomial.o polynomial.o
	${CC} ${CFLAG} -o $@ $^ ${INCLUDE}; 

%:%.o
	${CC} ${CFLAG} -o $@ $^ ${INCLUDE} ${LLIB} ${GLIB}; 
-include .depend
depend: 
	g++ -MM *.cpp ${INCLUDE} > .depend


%.o:%.cpp
	${CC} ${CFLAG} -c $< ${INCLUDE};

# ========================================
.PHONEY:clean depend
clean:
	rm -rf *~  *.o  msa 

