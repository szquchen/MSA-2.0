# About MSA 2.0
MSA is a software that does a fit to electronic energies and gradients using a fitting basis that is invariant with respect to permutations of like atoms.

An older version can be downloaded at https://github.com/Kee-Wang/PES-Fitting-MSA. This new 2.0 version adds the function to include ab initio gradients in the fitting.

The detailed theory is discussed in the references listed at the end. More information can be found in the Tutorial.txt file and at https://scholarblogs.emory.edu/bowman/msa.

# Prerequisites
1. C++ and Fortran 90 compiler;
2. Python3 and Perl;
3. LAPACK library.

# Credits
MSA code: Zhen Xie

Gradient extension: Chen Qu

Python script: Qingfeng (Kee) Wang


# References
1. Xie, Z.; Bowman, J.M. Permutationally Invariant Polynomial Basis for Molecular Energy Surface Fitting via Monomial Symmetrization. J. Chem. Theory Comput. 2010, 6, 26-34.
2. Nandi, A.; Qu, C.; Bowman, J.M. An Assessment of the Effectiveness of Incorporating Gradients in Potential Energy Fitting Using Permutationally Invariant Polynomials. J. Chem. Theory Comput. (in preparation)
