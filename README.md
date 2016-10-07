# easyPCG

Hi there ! 

The project easyPCG provides a simple stand-alone C-library to solve a linear system of equations, say A.x=b, with help of 
the Preconditioned Conjugated Gradient (PCG) Method. 

I do not claim that this implementation is state-of-art numerics. Instead, the goal is to provide a solver that will allow 
the user to get it running in just 5 minutes without caring for any dependencies (BLAS,LAPACK, etc. ...). Just copy the 
*.c and *.h file in your project folder, include it, and have fun. 

Note that the matrix A must be symmetric and positive-definite. Otherwise, the algorithm may not converge. The matrix must be given in the commonly used Compressed Sparse Row format (CSR). Proconditioning is not necessary but will speed-up convergence. Some basic preconditioners are implemented, which can be selected by the user. If favored, you can also pass your own 
preconditioner and matrix-multiplication routine to the solver. 

The relevant files are "easyPCG.c" and ".h". The additional files "linalg_stuff.c" and "testPCG.c" are given to provide 
a little test/example program. To get it running, put all files *.c and *.h files in a folder with the makefile and 
type "make". You can then run "./testPCG arg1 arg2 ..."  that will solve the Poisson equation with a point source on a 2D square 
domain. 

The project files are documented in doxygen style. If you have doxygen installed, you can build the documentation yourself. Otherwise there is a "reference.pdf" in the repository as a manual. 
