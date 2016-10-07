/* lib: easyPCG  */

#ifndef ePCG_H
#define ePCG_H

/////////////////////////////////////////* includes *//////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <math.h>
#include <stdlib.h> 
#include <string.h>
// #include <omp.h> 

///////////////////////////////////////////* public types *//////////////////////////////////////////////////////////////////

/**
 * ID for preconditioner used in the PCG algorithm
 */
typedef enum PREC_ID{
	PCG_NONE,							/**< Preconditioner mode ID: no preconditioner used */ 
	PCG_USER_DEFINED,					/**< Preconditioner mode ID: user-defined preconditioner used */ 
	PCG_JACOBI,							/**< Preconditioner mode ID: Jacobi preconditioner used */ 
	PCG_ICHOL							/**< Preconditioner mode ID: incomplete Cholesky preconditioner used */ 
}preconditionerID;

//////////////////////////////////////* public functions */////////////////////////////////////////////////////////////

int PCG_solve(double* X,double* B,int n,double tol,int max_iter);
void PCG_set_preconditioner_mode(void* System_matrix,preconditionerID mode);
void PCG_set_default_matrix(void* A);
void PCG_set_precon(void (*Precon)(double* R,double* X,int n));
void PCG_set_mult(void (*Mult)(double* R,double* X,int n));
void PCG_clean();
//void CSR_incompleteCholeski(CSR_matrix* A,CSR_matrix* L);

#endif
