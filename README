This is the MSA software that fits the potential energy surface using
fitting bases that are invariant with respect to permutations of like
atoms. The detailed theory is described in the paper
Xie, Z.; Bowman, J. M. J. Chem. Theory Comput. 2010, 6, 26-34
Nandi, A.; Qu, C.; Bowman, J. M. J. Chem. Theory Comput. 

Requirements:
To use this software, you need a
(A) a Fortran 90 and a C++ compiler;
(B) "dgelss" subroutine from LAPACK
(C) Perl and Python

We've tested on Linux operating system, using "ifort" and "g++" as the
Fortran and C++ compilers, and Intel MKL library for the "dgelss".
Other compilers should also work, but you may need to use different
flags. LAPACK is available at www.netlib.org, MKL is also freely
available at https://software.intel.com/en-us/mkl.

After downloading and unzipping the "MSA", the following folders and
files should be in the package:

(A) A folder called "src", which contains the source codes. Inside it
    there are
    (a) A folder called "emsa", which contains the C++ codes that
        generate the monomials and polynomials.
    (b) a test program "test.f90" that uses the fitted potential to
        calculate the potential of a given geometry.
    (c) A Makefile that is used to compile the Fortran codes.
    (d) A test data file "test.xyz" that contains an arbitrary geometry
        of CH4. This test data is used by the "getpot" program.
    (e) Two Perl scripts.

(B) A example data file "geom.inp" that contains 1000 Ch4 geometries
    and corresponding ab initio energies and gradients. This is
    an example dataset to fit the potential.

(C) A Python script "msa.py" to fit the potential

(D) A script called "reset" to clean up everything and reset the whole
    MSA directory to its original condition. (Be careful because all
    the output and results will be removed)

(E) An additional input file "param.inp" that specifies the values of 
    two parameters. This is optional and users could still use the
    default values.

Here are the procedure and a few explanations for using the software.
The procedure is also described in our tutorial video. You can watch
the video on our group page at:
https://scholarblogs.emory.edu/bowman/msa/
We recommend that you watch the video first and get familiar with this
software.

(0) First take a look at the "geom.inp" file.
    It contains geometries and corresponding energies and gradients.
    The first line is the number of atoms, and second line is the
    energy in hartree. Each of the following lines contains the x,y,z
    coordinates of an atom (in Angstrom) and then the gradient in x,y,z
    directions (in hartree/bohr). The order of atoms must be consistent
    with the permutation group. In this example, the permutation group
    is A4B for CH4, so the C should follow the four H's.

(1) Modify the two Makefile's, based on your compiler and LAPACK.

(2) Run the Python script:
    ./msa.py
    (or python msa.py)

    (A) Input the polynomial order you would like to use for the fitting.
        In our example we use 5. The number of coefficients increases
        rapidly when the polynomial order becomes larger. So you may want
        to start with a low polynomial order.
    (B) Input the permutation group. Our example is CH4, so we use "4 1"
        here.
    (C) Input the name of the data file, which is "geom.inp" in this
        example.

(3) The program tells you the number of coefficients and the number of
    points in the data file. And it asks you if you would like to continue.
    If the number of coefficients is too small (which leads to large
    fitting error), or too large (which may cause over-fitting), you can
    type "n" to terminate and then pick another polynomial order. If you
    would like continue, just enter y, and the program will continue.

(4) The program asks you if you want to include gradients in the fitting.
    Enter "y" or "n". Note that if you choose "y", the data file must
    contain the gradients. In this example just enter "y"

(5) The program asks you if you would like to specify the values of two
    additional parameters. You can use the default by entering "n", or
    enter the file name. Here an example "param.inp" is provided, and
    these two parameters are explained in this file.

(6) The program fits the potential energy surface and when it finishes,
    the root-mean-square fitting error is printed on the screen.
    The coefficients of the fit is written in "coeff.dat", and three
    Fortran code files "pes_shell.f90", "basis.f90", and "gradient.f90"
    are also generated.

(7) The test program "test.x" is compiled, and if you would like to run
    the test, use the command
    ./test.x test.xyz
    The ab initio and predicted values will be printed in this file, and
    if you use polynomial order of 5 with the default parameters, the
    differece between ab initio and prediction is small.

(8) If you would like to use the fit in your own program, pes_shell.f90,
    basis.f90, gradient.f90, and coeff.dat are necessary. Copy these four
    files to the folder that contains your own program, and in your own
    Fortran code, insert "use pes_shell", and "call pes_init()", (as we
    do in the "test.f90" example), and you can calculate the potential
    of any configuration using the "f" function, and the gradient using the
    "g" function.

If you have any questions, please contact Kee Wang or Chen Qu:
Kee : kee.wang@emory.edu
Chen: cqu3@emory.edu
