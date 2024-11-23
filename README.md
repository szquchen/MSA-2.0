# About MSA 2.0 and MSA 2.0.1 (latest version)
MSA 2.0.1 is a software that generates permutationally invariant polynomials (PIPs) and gradients and optionally does a fit to electronic energies and gradients using PIPs.

The detailed theory is discussed in the references listed at the end. More information can be found in the Tutorial.txt file and at https://scholarblogs.emory.edu/bowman/msa.

# Prerequisites
For PIPs
1. C++  compiler;
2. Python3 and Perl;
   
For fits

4. Fortran 90   
5. LAPACK library.

# Credits
Original MSA code: Zhen Xie

Using gradient in the fit: Chen Qu

Python Wrapper: Qingfeng (Kee) Wang


# References
1. Xie, Z., Bowman, J.M. Permutationally Invariant Polynomial Basis for Molecular Energy Surface Fitting via Monomial Symmetrization. J. Chem. Theory Comput. 2010, 6, 26-34.
2. Nandi, A. Qu, Chen, Bowman, J.M. Using Gradients in Permutationally Invariant Polynomial Potential Fitting: A Demonstration for CH4 Using as Few as 100 Configurations, J. Chem. Theory Comput. 2019, 15, 2826-2835
